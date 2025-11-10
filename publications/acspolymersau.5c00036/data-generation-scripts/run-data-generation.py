#!/usr/bin/env python

import os
import time

import numpy as np
import pandas as pd
import pylimer_tools_cpp as ptc
from pylimer_tools.calc.miller_macosko_theory import (predict_maximum_p)
from pylimer_tools.generate_network import generate_structure
from pylimer_tools.io.bead_spring_parameter_provider import (
    ParameterType, get_parameters_for_polymer)
from tqdm import tqdm

base_path = os.path.join(os.path.dirname(__file__), "..")

n_chains_tot = int(1e5)
n_rep = 5

params = get_parameters_for_polymer("PDMS", ParameterType.GAUSSIAN)
r02_slope = params.get("R02")
r02_slope_magnitude = r02_slope.to(params.get("distance_units") ** 2).magnitude
kbt = params.get("T") * params.get("kb")
gamma_conversion_factor = (
    (kbt / ((params.get("distance_units")) ** 3)).to("MPa").magnitude
)

input_data = pd.read_csv(
    os.path.join(
        base_path, "input-data", "combined-all-experimental-systems.csv"
    )
)

output_path = os.path.join(
    base_path, "output-data", "force-balance-results.csv"
)

try:
    output_data = pd.read_csv(output_path)

    # join input and output data
    combined_data = pd.merge(input_data, output_data, how="outer", on=["system_id"])

    # find groups where we have less than 10 replications in the force balance data
    to_redo = combined_data.groupby("system_id").filter(lambda x: len(x) < 10)
except FileNotFoundError:
    output_data = pd.DataFrame()

    to_redo = input_data.copy()

# for these systems, we want to re-do the force balance calculations
for system_id in tqdm(to_redo["system_id"].unique(), desc="Processing systems"):
    total_start_time = time.time()

    system_data = input_data[input_data["system_id"] == system_id].iloc[
        0
    ]  # Get first row for this system

    # use `r` and `b2` to compute the number of chains we need
    n_chains_bi = n_chains_tot
    n_chains_mono = int(
        2 * (n_chains_bi / system_data["b2"] - n_chains_bi)
        if not pd.isna(system_data["b2"])
        else 0
    )
    n_chains_xlink = int(
        (((n_chains_bi) * 2 + n_chains_mono) * system_data["r"]) / (system_data["f"])
    )

    data_for_this_universe = {
        "system_id": system_id,
        "fb_n_rep": n_rep,
        "n_chains_bi": n_chains_bi,
        "n_chains_mono": n_chains_mono,
        "n_chains_xlink": n_chains_xlink,
        "using_p_or_wsol": (
            "p" if pd.isna(system_data["w_sol_experimental"]) else "wsol"
        ),
        "target_w_sol": system_data["w_sol_experimental"],
        "target_p": system_data["p"],
    }

    # Generate the universe
    universe_start_time = time.time()
    try:
        universe = generate_structure(
            params=params,
            target_p=(
                None
                if not pd.isna(system_data["w_sol_experimental"])
                else min(
                    system_data["p"],
                    predict_maximum_p(
                        r=system_data["r"],
                        f=int(system_data["f"]),
                        b2=system_data["b2"] if not pd.isna(system_data["b2"]) else 1.0,
                    ),
                )
            ),
            target_wsol=(
                None
                if pd.isna(system_data["w_sol_experimental"])
                else system_data["w_sol_experimental"]
            ),
            n_chains_crosslinkers=n_chains_xlink,
            target_f=int(system_data["f"]),
            n_solvent_chains=0,  # Default, not specified in CSV
            n_beads_per_solvent_chain=0,  # Default, not specified in CSV
            n_chains_1=n_chains_bi,
            n_beads_per_chain_1=int(system_data["n_beads_bifunctional"]),
            n_chains_2=0,  # Default, not specified in CSV
            n_beads_per_chain_2=0,  # Default, not specified in CSV
            n_mono_chains=n_chains_mono,
            n_mono_beads_per_chain=(
                0
                if pd.isna(system_data.get("n_beads_monofunctional"))
                else int(system_data["n_beads_monofunctional"])
            ),
            n_beads_per_xlink=(
                1
                if pd.isna(system_data.get("n_beads_crosslink"))
                else int(system_data["n_beads_crosslink"])
            ),
            remove_wsol=system_data["wsol_status"] == "removed",
            functionalize_discrete=False,  # Default
            disable_primary_loops=False,  # Default
            disable_secondary_loops=False,  # Default
            z_score_std_mult=3.0,  # Default
        )
    except Exception as e:
        print(f"Error generating universe for system {system_id}: {e}")
        continue
    universe_generation_time = time.time() - universe_start_time

    data_for_universe = {
        "fb_n_entanglements_initial": [],
        "fb_n_entanglements_final": [],
        "fb_n_steps": [],
        "fb_modulus [MPa]": [],
        "fb_wsol": [],
        "fb_wdang": [],
        "fb_time_per_rep": [],
    }

    # Use the structure with force balance
    fb_total_start_time = time.time()
    for rep in range(n_rep):
        fb_rep_start_time = time.time()
        fb = ptc.MEHPForceBalance2(
            universe=universe,
            crosslinker_type=2,
            nr_of_entanglements_to_sample=int(
                params.get_entanglement_density() * universe.get_volume()
            ),
            upper_sampling_cutoff=params.get_sampling_cutoff(),
            same_strand_cutoff=0,
            entanglements_as_springs=False,
        )
        data_for_universe["fb_n_entanglements_initial"].append(
            fb.get_nr_of_extra_atoms()
        )
        fb.run_force_relaxation(
            simplification_mode=ptc.StructureSimplificationMode.INACTIVE_THEN_X2F
        )
        data_for_universe["fb_n_entanglements_final"].append(fb.get_nr_of_extra_atoms())
        data_for_universe["fb_n_steps"].append(fb.get_nr_of_iterations())
        data_for_universe["fb_modulus [MPa]"].append(
            gamma_conversion_factor
            * np.sum(fb.get_gamma_factors(r02_slope_magnitude))
            / universe.get_volume()
        )
        data_for_universe["fb_wsol"].append(fb.get_soluble_weight_fraction())
        data_for_universe["fb_wdang"].append(fb.get_dangling_weight_fraction())
        data_for_universe["fb_time_per_rep"].append(time.time() - fb_rep_start_time)

    fb_total_time = time.time() - fb_total_start_time

    for key, val in data_for_universe.items():
        data_for_this_universe[key] = np.mean(val)
        data_for_this_universe[key + "_std"] = np.std(val)

    # Also measure some phantom quantities
    phantom_start_time = time.time()
    phantom_fb = ptc.MEHPForceBalance2(universe=universe)
    phantom_fb.run_force_relaxation()
    phantom_time = time.time() - phantom_start_time

    data_for_this_universe["phantom_fb_n_steps"] = phantom_fb.get_nr_of_iterations()
    data_for_this_universe["phantom_fb_modulus [MPa]"] = (
        gamma_conversion_factor
        * np.sum(phantom_fb.get_gamma_factors(r02_slope_magnitude))
        / universe.get_volume()
    )
    data_for_this_universe["phantom_fb_wsol"] = phantom_fb.get_soluble_weight_fraction()
    data_for_this_universe["phantom_fb_wdang"] = (
        phantom_fb.get_dangling_weight_fraction()
    )

    # Add timing information
    total_time = time.time() - total_start_time
    data_for_this_universe["timing_universe_generation"] = universe_generation_time
    data_for_this_universe["timing_fb_total"] = fb_total_time
    data_for_this_universe["timing_phantom"] = phantom_time
    data_for_this_universe["timing_total"] = total_time

    # Also add some metadata
    data_for_this_universe["pylimertools_v"] = ptc.version_information()
    data_for_this_universe["processed_at"] = pd.Timestamp.now()

    # add the new data to the output data, save it
    new_row = pd.DataFrame([data_for_this_universe])

    # Append to output data
    output_data = pd.concat([output_data, new_row], ignore_index=True)

    # Save the updated data
    output_data.to_csv(output_path, index=False)

    print(f"Completed system {system_id}, saved results to {output_path}")
