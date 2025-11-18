#!/usr/bin/env bash

cd "$(dirname "$0")" || exit

# Variables to easily switch between different targets
OUTFOLDER="output"
SCRATCH="/cluster/scratch/betim"
SCRATCH_OUTFOLDER="$SCRATCH/doctorate/uniaxial-deformation/hexagonal-lattice/$OUTFOLDER"
LAMMP_IN_FOLDER="input"
INPUT_FILE_TO_DIRECTION="input_raw/uniaxial_fast_tension_2d_phantom_cont.in"
mkdir -p $OUTFOLDER
mkdir -p $LAMMP_IN_FOLDER
EULER_NUM_CORES=27
CONT_NR=0
EULER_RUN_FILE="startUniaxialCont$CONT_NR""Euler.sh"

mkdir -p $OUTFOLDER
mkdir -p $LAMMP_IN_FOLDER

echo "#!/usr/bin/env bash" >$EULER_RUN_FILE
{
  echo "module load gcc/8.2.0"
  echo "module load openmpi/4.0.2"
  echo "export OMP_DYNAMIC=true"
  echo "export OMP_NUM_THREADS=1"
  echo "mkdir -p $OUTFOLDER"
} >>$EULER_RUN_FILE

EULER_CMD_PREFIX="bsub -G es_gusev -R \"span[ptile=$EULER_NUM_CORES]\" -n $EULER_NUM_CORES"

for CHAINFILE in structure/*0000.structure.out; do
  # for CHAINFILE in structure/M9999*.out; do

  CHAINNAME=$(basename "$CHAINFILE")

  for TEMP in 150 300 450; do
    for ERATE in 2 6 10; do
      for DIR_COMBI in xy yx; do

        INPUT_DIR1="${DIR_COMBI:0:1}"
        INPUT_DIR2="${DIR_COMBI:1:1}"

        echo "$INPUT_DIR1 $INPUT_DIR2 $ERATE $TEMP"

        TARGET_INPUT="$LAMMP_IN_FOLDER/$INPUT_DIR1/uniaxial_tension_cont_$ERATE""_T_$TEMP""_$CHAINNAME.in"
        mkdir -p "$LAMMP_IN_FOLDER/$INPUT_DIR1/"

        RESTART_INPUT_FILE1="./$OUTFOLDER/$INPUT_DIR1/deformed_restart_1_$ERATE""e-2_T_$TEMP""_$CHAINNAME.out"
        RESTART_INPUT_FILE2="./$OUTFOLDER/$INPUT_DIR1/deformed_restart_2_$ERATE""e-2_T_$TEMP""_$CHAINNAME.out"
        # find the newer of the two restart files to use afterwards
        unset -v RESTART_INPUT_FILE
        for file in $RESTART_INPUT_FILE1 $RESTART_INPUT_FILE2; do
          [[ $file -nt $RESTART_INPUT_FILE ]] && RESTART_INPUT_FILE=$file
        done

        if [ -z ${RESTART_INPUT_FILE+x} ]; then
          RESTART_INPUT_FILE=$RESTART_INPUT_FILE1
          echo "Restart file set to $RESTART_INPUT_FILE, but does not exist."
          continue 
        fi

        cp "$INPUT_FILE_TO_DIRECTION" "$TARGET_INPUT"
        ESCAPED_RESTARTFILE=$(printf '%s\n' "$RESTART_INPUT_FILE" | sed -e 's/[\/&]/\\&/g')
        ESCAPED_CHAINFILE=$(printf '%s\n' "$CHAINFILE" | sed -e 's/[\/&]/\\&/g')
        ESCAPED_CHAINNAME=$(printf '%s\n' "$CHAINNAME" | sed -e 's/[\/&]/\\&/g')
        ESCAPED_SCRATCH_OUT=$(printf '%s\n' "$SCRATCH_OUTFOLDER" | sed -e 's/[\/&]/\\&/g')

        # Replace input, output & direction
        sed -i "s/RESTART_INPUT_FILE/$ESCAPED_RESTARTFILE/gi" "$TARGET_INPUT"
        sed -i "s/TEMPERATURE/$TEMP/gi" "$TARGET_INPUT"
        sed -i "s/INPUT_ERATE/$ERATE""e-2/gi" "$TARGET_INPUT"
        sed -i "s/INPUT_FILE/$ESCAPED_CHAINFILE/gi" "$TARGET_INPUT"
        sed -i "s/INPUT_NAME/$ESCAPED_CHAINNAME/gi" "$TARGET_INPUT"
        sed -i "s/SCRATCH_OUT/$ESCAPED_SCRATCH_OUT/gi" "$TARGET_INPUT"
        sed -i "s/INPUT_DIR1/$INPUT_DIR1/gi" "$TARGET_INPUT"
        sed -i "s/INPUT_DIR2/$INPUT_DIR2/gi" "$TARGET_INPUT"

        # Write to run files
        {
          echo "mkdir -p $OUTFOLDER/$INPUT_DIR1/"
          echo "mkdir -p $SCRATCH_OUTFOLDER/$INPUT_DIR1/"
          echo "$EULER_CMD_PREFIX -W 120:00 -N -J \"$INPUT_DIR1-$ERATE-$TEMP-uniaxial-$CHAINNAME\" -o \"$OUTFOLDER/$INPUT_DIR1/uniaxial_tension_$ERATE""_T_$TEMP""_$CHAINNAME-run.oout\" \"mpirun /cluster/home/betim/masters/lammps/build/lmp_mpi -suffix hybrid omp opt -in \\\"$TARGET_INPUT\\\" -l \\\"$SCRATCH_OUTFOLDER/$INPUT_DIR1/uniaxial_tension_$ERATE""_T_$TEMP""_$CHAINNAME-log.out\\\" >> \\\"$OUTFOLDER/$INPUT_DIR1/uniaxial_tension_$ERATE""_T_$TEMP""_$CHAINNAME-run.out\\\"\";"
        } >>"$EULER_RUN_FILE"
      done
    done
  done

  echo "Processed input-file: $CHAINFILE"
done

chmod u+x "$EULER_RUN_FILE"
