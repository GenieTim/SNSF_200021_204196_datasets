#!/usr/bin/env python
import argparse
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pyrender
import trimesh
from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import Universe
from warning_utils import format_Warning

warnings.formatwarning = format_Warning


def move_atoms_into_box(atom_coords, box_coords):
    """Ensure atoms are within the periodic box."""
    for i in range(len(atom_coords)):
        for j in range(3):  # x, y, z coordinates
            box_length = box_coords[j][1] - box_coords[j][0]
            atom_coords[i][j] = (
                atom_coords[i][j] - box_coords[j][0]
            ) % box_length + box_coords[j][0]
    return atom_coords


def create_glassy_materials(atom_colors, interactive=False):
    """Create glass-like materials with transparency."""
    alpha = 0.6 if interactive else 0.8
    return {
        atom_type: pyrender.MetallicRoughnessMaterial(
            metallicFactor=0.025,
            roughnessFactor=0.05,
            baseColorFactor=color + [alpha],
            alphaMode="BLEND",  # Enable transparency
        )
        for atom_type, color in atom_colors.items()
    }


def create_metallic_materials(atom_colors, interactive=False):
    """Create metallic materials with high reflectivity."""
    return {
        atom_type: pyrender.MetallicRoughnessMaterial(
            metallicFactor=0.8,
            roughnessFactor=0.2,
            baseColorFactor=color + [1.0],
            alphaMode="OPAQUE",
        )
        for atom_type, color in atom_colors.items()
    }


def create_matte_materials(atom_colors, interactive=False):
    """Create matte materials with no reflectivity."""
    return {
        atom_type: pyrender.MetallicRoughnessMaterial(
            metallicFactor=0.0,
            roughnessFactor=1.0,  # maximum roughness for diffuse look
            baseColorFactor=color + [1.0],
            alphaMode="OPAQUE",
        )
        for atom_type, color in atom_colors.items()
    }


def srgb_to_linear(rgb):
    """Convert sRGB [0,1] to linear RGB [0,1] for accurate color in 3D rendering."""
    rgb = np.array(rgb)
    linear = np.where(
        rgb <= 0.04045,
        rgb / 12.92,
        ((rgb + 0.055) / 1.055) ** 2.4,
    )
    return linear.tolist()


def visualize_universe_3d(
    universe: Universe,
    filename: str,
    camera_position: tuple,
    interactive: bool = False,
    style: str = "glassy",
):
    """
    Visualize LAMMPS universe in 3D using pyrender.

    In interactive mode, uses lower polygon counts for better performance during navigation.
    Style options: 'glassy', 'metallic', 'matte'
    """
    # Extract atomic data
    all_atoms = universe.get_atoms()
    atom_types = [a.get_type() for a in all_atoms]
    assert all(
        [isinstance(t, int) for t in atom_types]
    ), "All atom types must be integers."
    atom_coords = np.array([[a.get_x(), a.get_y(), a.get_z()] for a in all_atoms])

    # Get periodic box info
    box = universe.get_box()
    box_coords = np.array(
        [
            [box.get_low_x(), box.get_high_x()],
            [box.get_low_y(), box.get_high_y()],
            [box.get_low_z(), box.get_high_z()],
        ]
    )

    # Calculate box center and dimensions for camera positioning
    box_center = np.mean(box_coords, axis=1)
    box_dimensions = box_coords[:, 1] - box_coords[:, 0]
    max_dimension = np.max(box_dimensions)

    # If camera position is default (0,0,10), calculate optimal position
    if camera_position == (0.0, 0.0, 10.0):
        # Position camera to view the entire box with some margin
        camera_distance = max_dimension * 2.5  # Adjust multiplier as needed
        camera_position = (
            box_center[0] + camera_distance * 0.4,
            box_center[1] + camera_distance * 0.5,
            box_center[2] + camera_distance * 0.3,
        )

    atom_coords = move_atoms_into_box(atom_coords, box_coords)

    # Define atom colors using matplotlib's tab10 colormap
    tab10 = plt.get_cmap("tab10")
    atom_colors = {i + 1: list(tab10(i)[:3]) for i in range(10)}

    # Convert atom colors to linear RGB for accurate 3D rendering
    atom_colors = {k: srgb_to_linear(v) for k, v in atom_colors.items()}

    scene = pyrender.Scene()

    # Create shared mesh and materials for instanced atoms
    # Use lower subdivision for better interactive performance
    subdivisions = 2 if interactive else 3
    base_sphere = trimesh.creation.icosphere(subdivisions=subdivisions, radius=0.2)

    # Create materials based on specified style
    if style == "metallic":
        atom_materials = create_metallic_materials(atom_colors, interactive)
    elif style == "matte":
        atom_materials = create_matte_materials(atom_colors, interactive)
    else:  # default to glassy
        atom_materials = create_glassy_materials(atom_colors, interactive)

    # Group atoms by type and create instanced meshes
    unique_types = list(set(atom_types))
    for atom_type in unique_types:
        # Get indices of atoms of this type
        type_indices = [i for i, t in enumerate(atom_types) if t == atom_type]
        type_coords = atom_coords[type_indices]

        # Create transformation matrices for atoms of this type
        type_poses = np.tile(np.eye(4), (len(type_coords), 1, 1))
        type_poses[:, :3, 3] = type_coords

        # Use appropriate material for this atom type, fallback to type 1 if not found
        material = atom_materials.get(atom_type, atom_materials[1])

        # Add instanced atoms to the scene
        atom_mesh = pyrender.Mesh.from_trimesh(
            base_sphere, material=material, poses=type_poses
        )
        scene.add(atom_mesh)

    # Add periodic box frame
    box_min = box_coords[:, 0]
    box_max = box_coords[:, 1]

    # Define the 8 corners of the box
    corners = np.array(
        [
            [box_min[0], box_min[1], box_min[2]],  # 0: min corner
            [box_max[0], box_min[1], box_min[2]],  # 1
            [box_max[0], box_max[1], box_min[2]],  # 2
            [box_min[0], box_max[1], box_min[2]],  # 3
            [box_min[0], box_min[1], box_max[2]],  # 4
            [box_max[0], box_min[1], box_max[2]],  # 5
            [box_max[0], box_max[1], box_max[2]],  # 6: max corner
            [box_min[0], box_max[1], box_max[2]],  # 7
        ]
    )

    # Define the 12 edges of the box (connecting corner indices)
    box_edges = np.array(
        [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0],  # bottom face
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 4],  # top face
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],  # vertical edges
        ]
    )

    # Create line segments for the box frame
    box_vertices = []
    for edge in box_edges:
        box_vertices.extend([corners[edge[0]], corners[edge[1]]])
    box_vertices = np.array(box_vertices)

    # Create thicker box frame using cylinders instead of lines
    box_material = pyrender.MetallicRoughnessMaterial(
        metallicFactor=0.0,
        roughnessFactor=0.0,
        baseColorFactor=[0.1, 0.1, 0.1, 1.0],  # Dark gray frame
        alphaMode="OPAQUE",
    )

    # Create cylinder for each box edge with increased thickness
    box_cylinder_radius = 0.075
    box_cylinder = trimesh.creation.cylinder(
        radius=box_cylinder_radius, height=1.0, sections=8
    )

    box_poses = []
    for edge in box_edges:
        start = corners[edge[0]]
        end = corners[edge[1]]
        mid_point = (start + end) / 2
        edge_vector = end - start
        height = np.linalg.norm(edge_vector)

        if height > 0:
            # Create transformation matrix for the box edge
            transform = np.eye(4)
            transform[:3, 3] = mid_point

            # Orient cylinder along edge direction
            z_axis = edge_vector / height
            # Create orthonormal basis
            if abs(z_axis[2]) < 0.9:
                x_axis = np.cross([0, 0, 1], z_axis)
            else:
                x_axis = np.cross([1, 0, 0], z_axis)
            x_axis = x_axis / np.linalg.norm(x_axis)
            y_axis = np.cross(z_axis, x_axis)

            # Set rotation matrix and scale z-axis by edge length
            transform[:3, 0] = x_axis
            transform[:3, 1] = y_axis
            transform[:3, 2] = z_axis * height

            box_poses.append(transform)

    # Add instanced box frame to the scene
    box_mesh = pyrender.Mesh.from_trimesh(
        box_cylinder, material=box_material, poses=np.array(box_poses)
    )
    scene.add(box_mesh)

    # Handle periodic bonds using instancing
    edges = universe.get_edges()
    bond_from = edges["edge_from"]
    bond_to = edges["edge_to"]

    # Use lower polygon count for cylinders in interactive mode for better performance
    sections = 8 if interactive else 16
    base_cylinder = trimesh.creation.cylinder(
        radius=0.075, height=1.0, sections=sections
    )

    # Create bond material based on style
    if style == "metallic":
        bond_material = pyrender.MetallicRoughnessMaterial(
            metallicFactor=0.7,
            roughnessFactor=0.3,
            baseColorFactor=[0.6, 0.8, 1.0, 1.0],
            alphaMode="OPAQUE",
        )
    elif style == "matte":
        bond_material = pyrender.MetallicRoughnessMaterial(
            metallicFactor=0.0,
            roughnessFactor=0.8,
            baseColorFactor=[0.6, 0.8, 1.0, 1.0],
            alphaMode="OPAQUE",
        )
    else:  # glassy
        alpha = 0.6 if interactive else 0.8
        bond_material = pyrender.MetallicRoughnessMaterial(
            metallicFactor=0.025,
            roughnessFactor=0.05,
            baseColorFactor=[0.6, 0.8, 1.0, alpha],
            alphaMode="BLEND" if not interactive else "OPAQUE",
        )

    bond_poses = []
    omitted_bonds = 0
    omitted_half_bonds = 0
    primary_loop_bonds = 0
    min_bond_length = 1e-12  # Minimum bond length to consider valid

    for i, j in zip(bond_from, bond_to):
        if i == j:
            # Skip self-bonds
            primary_loop_bonds += 1
            continue
        start = atom_coords[i].copy()
        end = atom_coords[j].copy()
        box_lengths = box_coords[:, 1] - box_coords[:, 0]
        box_mins = box_coords[:, 0]
        box_maxs = box_coords[:, 1]

        # Calculate the minimum image bond vector
        bond_vector = end - start
        is_periodic = False

        # Check if bond crosses periodic boundary in any direction
        for axis in range(3):
            while abs(bond_vector[axis]) > box_lengths[axis] / 2:
                # Apply minimum image convention to find shortest bond path
                if bond_vector[axis] > 0:  # end is "ahead" of start
                    bond_vector[axis] -= box_lengths[axis]  # wrap backward
                else:  # end is "behind" start
                    bond_vector[axis] += box_lengths[axis]  # wrap forward
                is_periodic = True

        if not is_periodic:
            # No boundary crossing - draw normal bond
            mid_point = (start + end) / 2
            height = np.linalg.norm(end - start)

            if height > min_bond_length:
                transform = np.eye(4)
                transform[:3, 3] = mid_point

                z_axis = (end - start) / height
                if abs(z_axis[2]) < 0.9:
                    x_axis = np.cross([0, 0, 1], z_axis)
                else:
                    x_axis = np.cross([1, 0, 0], z_axis)
                x_axis = x_axis / np.linalg.norm(x_axis)
                y_axis = np.cross(z_axis, x_axis)

                transform[:3, 0] = x_axis
                transform[:3, 1] = y_axis
                transform[:3, 2] = z_axis * height

                bond_poses.append(transform)
            else:
                omitted_bonds += 1
        else:
            # Periodic bond - draw two segments showing actual bond direction
            # Calculate the corrected end position using minimum image convention
            corrected_end = start + bond_vector

            # Segment 1: From start atom to box boundary
            # Find which boundary the bond crosses first
            seg1_end = start.copy()
            min_t = float("inf")

            for axis in range(3):
                if abs(bond_vector[axis]) > 1e-10:  # Avoid division by zero
                    # Calculate intersection with both boundaries for this axis
                    if bond_vector[axis] > 0:  # Moving toward high boundary
                        t = (box_maxs[axis] - start[axis]) / bond_vector[axis]
                    else:  # Moving toward low boundary
                        t = (box_mins[axis] - start[axis]) / bond_vector[axis]

                    if 0 < t < min_t:
                        min_t = t
                        # Calculate intersection point
                        intersection = start + t * bond_vector
                        # Clamp to box boundaries to avoid floating point errors
                        intersection = np.clip(intersection, box_mins, box_maxs)
                        seg1_end = intersection

            # Create first segment (start to boundary) if it has reasonable length
            seg1_vector = seg1_end - start
            seg1_length = np.linalg.norm(seg1_vector)

            if seg1_length > min_bond_length:
                seg1_mid = (start + seg1_end) / 2
                transform1 = np.eye(4)
                transform1[:3, 3] = seg1_mid

                z_axis = seg1_vector / seg1_length
                if abs(z_axis[2]) < 0.9:
                    x_axis = np.cross([0, 0, 1], z_axis)
                else:
                    x_axis = np.cross([1, 0, 0], z_axis)
                x_axis = x_axis / np.linalg.norm(x_axis)
                y_axis = np.cross(z_axis, x_axis)

                transform1[:3, 0] = x_axis
                transform1[:3, 1] = y_axis
                transform1[:3, 2] = z_axis * seg1_length

                bond_poses.append(transform1)
            else:
                omitted_half_bonds += 1

            # Segment 2: From opposite boundary to end atom
            # Find where the bond would enter the box from the other side
            seg2_start = end.copy()
            min_t = float("inf")

            # Work backwards from end atom
            reverse_bond_vector = -bond_vector

            for axis in range(3):
                if abs(reverse_bond_vector[axis]) > 1e-10:  # Avoid division by zero
                    # Calculate intersection with boundaries
                    if reverse_bond_vector[axis] > 0:  # Moving toward high boundary
                        t = (box_maxs[axis] - end[axis]) / reverse_bond_vector[axis]
                    else:  # Moving toward low boundary
                        t = (box_mins[axis] - end[axis]) / reverse_bond_vector[axis]

                    if 0 < t < min_t:
                        min_t = t
                        # Calculate intersection point
                        intersection = end + t * reverse_bond_vector
                        # Clamp to box boundaries to avoid floating point errors
                        intersection = np.clip(intersection, box_mins, box_maxs)
                        seg2_start = intersection

            # Create second segment (boundary to end) if it has reasonable length
            seg2_vector = end - seg2_start
            seg2_length = np.linalg.norm(seg2_vector)

            if seg2_length > min_bond_length:
                seg2_mid = (seg2_start + end) / 2
                transform2 = np.eye(4)
                transform2[:3, 3] = seg2_mid

                z_axis = seg2_vector / seg2_length
                if abs(z_axis[2]) < 0.9:
                    x_axis = np.cross([0, 0, 1], z_axis)
                else:
                    x_axis = np.cross([1, 0, 0], z_axis)
                x_axis = x_axis / np.linalg.norm(x_axis)
                y_axis = np.cross(z_axis, x_axis)

                transform2[:3, 0] = x_axis
                transform2[:3, 1] = y_axis
                transform2[:3, 2] = z_axis * seg2_length

                bond_poses.append(transform2)
            else:
                omitted_half_bonds += 1

    # Print statistics about omitted bonds
    if omitted_bonds > 0 or omitted_half_bonds > 0:
        warnings.warn(
            f"{omitted_bonds} + {omitted_half_bonds} + {primary_loop_bonds} bonds were omitted "
            + f"due to very short lengths or periodicity issues when rendering {filename}."
        )

    # Add instanced bonds to the scene
    bond_mesh = pyrender.Mesh.from_trimesh(
        base_cylinder, material=bond_material, poses=np.array(bond_poses)
    )
    scene.add(bond_mesh)

    # Add camera
    camera = pyrender.PerspectiveCamera(yfov=np.pi / 3.0)

    # Create camera transformation matrix to look at the box center
    camera_translation = np.array(camera_position)
    look_at_target = box_center

    # Calculate camera orientation to look at the box center
    forward = look_at_target - camera_translation
    forward = forward / np.linalg.norm(forward)

    # Create up vector (world up is typically [0, 0, 1] or [0, 1, 0])
    world_up = np.array([0, 0, 1])

    # If forward is parallel to world_up, use a different up vector
    if abs(np.dot(forward, world_up)) > 0.99:
        world_up = np.array([0, 1, 0])

    # Calculate right and up vectors for camera orientation
    right = np.cross(forward, world_up)
    right = right / np.linalg.norm(right)
    up = np.cross(right, forward)
    up = up / np.linalg.norm(up)

    # Create camera transformation matrix
    camera_matrix = np.eye(4)
    camera_matrix[:3, 0] = right
    camera_matrix[:3, 1] = up
    camera_matrix[:3, 2] = -forward  # Camera looks down negative z-axis
    camera_matrix[:3, 3] = camera_translation

    camera_node = pyrender.Node(camera=camera, matrix=camera_matrix)
    scene.add_node(camera_node)

    # Add better lighting for non-interactive mode
    if not interactive:
        directional_intensity = (
            10.0 if style == "metallic" else 1.5 if style == "matte" else 3.0
        )
        ambient_intensity = (
            5.0 if style == "metallic" else 2.5 if style == "matte" else 1.0
        )

        # Add directional light
        directional_light = pyrender.DirectionalLight(
            color=np.ones(3), intensity=directional_intensity
        )
        # light_pose = np.eye(4)
        # light_pose[:3, 3] = camera_translation
        scene.add(directional_light, pose=camera_matrix)  # light_pose)
        # Add ambient light
        ambient_light = pyrender.DirectionalLight(
            color=np.ones(3), intensity=ambient_intensity
        )
        ambient_pose = np.eye(4)
        ambient_pose[:3, 3] = box_center + np.array([0, 0, max_dimension])
        scene.add(ambient_light, pose=ambient_pose)

    if interactive:
        # Show interactive preview with improved performance materials
        print("Interactive mode: Use mouse to navigate the camera.")
        print("Controls:")
        print("  - Mouse drag: Rotate camera")
        print("  - Mouse wheel: Zoom in/out")
        print("  - Mouse middle button: Pan camera")
        print("  - Q/ESC: Quit and save the image with current camera position")
        print("  - F: Toggle fullscreen")
        print("  - A: Toggle axes display")
        print("  - S: Save screenshot")
        print("  - Z: Reset camera to default view")
        print("Note: Using simplified materials for better interactive performance.")

        # Create viewer for interactive preview with improved lighting
        pyrender.Viewer(scene, use_raymond_lighting=True, viewport_size=(1200, 1200))

        # After closing the viewer, get the final camera position
        # Note: The viewer modifies the camera node in place
        final_camera_position = camera_node.translation
        print(f"Final camera position: {final_camera_position}")

    # Render scene with higher resolution
    resolution = 1920 * 4
    r = pyrender.OffscreenRenderer(resolution, resolution)
    try:
        render_result = r.render(scene)
        if isinstance(render_result, tuple):
            color, depth = render_result
        else:
            color = render_result
    except Exception as e:
        print(f"Rendering error: {e}")
        return

    # Save image
    import imageio

    if color is not None:
        imageio.imwrite(filename, color)
        print(f"Image saved to: {filename}")
    else:
        print("Error: Failed to render image")

    # Clean up renderer
    r.delete()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize LAMMPS structure file in 3D")
    parser.add_argument("input_file", help="Path to LAMMPS structure file")
    parser.add_argument("output_file", help="Path to save the rendered image")
    parser.add_argument(
        "--interactive",
        "-i",
        action="store_true",
        help="Show interactive preview to set camera position",
    )
    parser.add_argument(
        "--style",
        "-s",
        choices=["glassy", "metallic", "matte"],
        default="glassy",
        help="Material style for atoms and bonds (default: glassy)",
    )
    parser.add_argument(
        "--camera-x", type=float, default=0.0, help="Initial camera X position"
    )
    parser.add_argument(
        "--camera-y", type=float, default=0.0, help="Initial camera Y position"
    )
    parser.add_argument(
        "--camera-z", type=float, default=10.0, help="Initial camera Z position"
    )

    args = parser.parse_args()

    camera_position = (args.camera_x, args.camera_y, args.camera_z)
    print("Reading input file:", args.input_file)
    universe = read_data_file(args.input_file)
    print(
        f"Visualizing universe with {universe.get_nr_of_atoms()} atoms in 3D with style {args.style}..."
    )
    visualize_universe_3d(
        universe, args.output_file, camera_position, args.interactive, args.style
    )
