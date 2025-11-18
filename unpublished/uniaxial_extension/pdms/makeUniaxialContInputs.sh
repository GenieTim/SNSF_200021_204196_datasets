#!/usr/bin/env bash

cd "$(dirname "$0")" || exit

# Variables to easily switch between different targets
OUTFOLDER="output"
SCRATCH="/cluster/scratch/betim"
SCRATCH_OUTFOLDER="$SCRATCH/doctorate/strain-relaxation/uniaxial-deformation/$OUTFOLDER"
LAMMP_IN_FOLDER="input"
RUN_TYPE_DIR="tension"
# RUN_TYPE_DIR="compression"
RUN_TYPE="uniaxial_$RUN_TYPE_DIR"
INPUT_FILE_TO_DIRECTION="input_raw/$RUN_TYPE""_cont.in"
mkdir -p $OUTFOLDER
mkdir -p $LAMMP_IN_FOLDER
EULER_NUM_CORES=48
CONT_NR=0
EULER_RUN_FILE="start_$RUN_TYPE""_cont$CONT_NR""_euler.sh"

echo "#!/usr/bin/env bash" >$EULER_RUN_FILE
{
  echo "module load gcc/8.2.0"
  echo "module load openmpi/4.0.2"
  echo "export OMP_DYNAMIC=true"
  echo "export OMP_NUM_THREADS=1"
  echo "mkdir -p $OUTFOLDER"
} >>$EULER_RUN_FILE

EULER_CMD_PREFIX="bsub -G es_gusev -R \"span[ptile=$EULER_NUM_CORES]\" -n $EULER_NUM_CORES"

# for CHAINFILE in $(gls -v1 molecules/*/crosslinking_*off*.out molecules/*/*/crosslinking_*off*.out); do
for CHAINFILE in structure/*.out; do

  CHAINNAME=$(basename "$CHAINFILE")

  for DIR_COMBI in xyz yzx zxy; do

    INPUT_DIR1="${DIR_COMBI:0:1}"
    INPUT_DIR2="${DIR_COMBI:1:1}"
    INPUT_DIR3="${DIR_COMBI:2:3}"

    echo "$INPUT_DIR1 $INPUT_DIR2 $INPUT_DIR3"

    TARGET_INPUT="$LAMMP_IN_FOLDER/$INPUT_DIR1/$RUN_TYPE""_$CHAINNAME""_cont$CONT_NR.in"
    RESTART_INPUT_FILE1=""
    RESTART_INPUT_FILE2=""
    if ((CONT_NR == 0)); then
      # since I accidently only made one of the two restart files unique:
      # RESTART_INPUT_FILE1="./$OUTFOLDER/$INPUT_DIR1/deformed_restart_1_$CHAINNAME.out"
      RESTART_INPUT_FILE2="./$OUTFOLDER/$INPUT_DIR1/deformed_$RUN_TYPE_DIR""_restart_2_$CHAINNAME.out"
      RESTART_INPUT_FILE1=$RESTART_INPUT_FILE2
    else
      CONT_NR_M1=$((CONT_NR - 1))
      RESTART_INPUT_FILE1="./$OUTFOLDER/$INPUT_DIR1/deformed_restart_1_$CHAINNAME""_cont$CONT_NR_M1.out"
      RESTART_INPUT_FILE2="./$OUTFOLDER/$INPUT_DIR1/deformed_restart_2_$CHAINNAME""_cont$CONT_NR_M1.out"
    fi

    # find the newer of the two restart files to use afterwards
    unset -v RESTART_INPUT_FILE
    for file in $RESTART_INPUT_FILE1 $RESTART_INPUT_FILE2; do
      [[ $file -nt $RESTART_INPUT_FILE ]] && RESTART_INPUT_FILE=$file
    done

    if [ -z ${RESTART_INPUT_FILE+x} ]; then
      RESTART_INPUT_FILE=$RESTART_INPUT_FILE1
      echo "Restart file set to $RESTART_INPUT_FILE, but does not exist."
    fi

    mkdir -p "inputs/$INPUT_DIR1/"
    cp "$INPUT_FILE_TO_DIRECTION" "$TARGET_INPUT"

    ESCAPED_RESTARTFILE=$(printf '%s\n' "$RESTART_INPUT_FILE" | sed -e 's/[\/&]/\\&/g')
    ESCAPED_CHAINNAME=$(printf '%s\n' "$CHAINNAME" | sed -e 's/[\/&]/\\&/g')
    ESCAPED_SCRATCH_OUT=$(printf '%s\n' "$SCRATCH_OUTFOLDER" | sed -e 's/[\/&]/\\&/g')

    # Replace input, output & direction
    sed -i "s/RESTART_INPUT_FILE/$ESCAPED_RESTARTFILE/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_NAME/$ESCAPED_CHAINNAME""_cont$CONT_NR/gi" "$TARGET_INPUT"
    sed -i "s/SCRATCH_OUT/$ESCAPED_SCRATCH_OUT/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_DIR1/$INPUT_DIR1/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_DIR2/$INPUT_DIR2/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_DIR3/$INPUT_DIR3/gi" "$TARGET_INPUT"

    # Write to run files
    {
      echo "mkdir -p $OUTFOLDER/$INPUT_DIR1/"
      echo "mkdir -p $SCRATCH_OUTFOLDER/$INPUT_DIR1/"
      echo "$EULER_CMD_PREFIX -W 120:00 -N -B -J \"$INPUT_DIR1-$RUN_TYPE-$CHAINNAME\" -o \"$OUTFOLDER/$INPUT_DIR1/$RUN_TYPE""_$CHAINNAME""_cont$CONT_NR-run.out\" mpirun /cluster/home/betim/masters/lammps/build/lmp_mpi -suffix omp -in \"$TARGET_INPUT\" -l \"$SCRATCH_OUTFOLDER/$INPUT_DIR1/$RUN_TYPE""_$CHAINNAME""_cont$CONT_NR-log.out\";"
    } >>"$EULER_RUN_FILE"
  done
done

chmod u+x "$EULER_RUN_FILE"
