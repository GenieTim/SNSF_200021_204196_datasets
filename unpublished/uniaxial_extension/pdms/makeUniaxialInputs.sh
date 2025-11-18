#!/usr/bin/env bash

cd "$(dirname "$0")" || exit

# Variables to easily switch between different targets
OUTFOLDER="output"
SCRATCH="/cluster/scratch/betim"
SCRATCH_OUTFOLDER="$SCRATCH/doctorate/strain-relaxation/uniaxial-deformation/$OUTFOLDER"
LAMMP_IN_FOLDER="input"
RUN_TYPE="uniaxial_tension"
# RUN_TYPE="uniaxial_compression"
INPUT_FILE_TO_DIRECTION="input_raw/$RUN_TYPE.in"
mkdir -p $OUTFOLDER
mkdir -p $LAMMP_IN_FOLDER
EULER_NUM_CORES=48
EULER_RUN_FILE="start_$RUN_TYPE""_Euler.sh"

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

for CHAINFILE in structure/*.out; do

  CHAINNAME=$(basename "$CHAINFILE")

  for DIR_COMBI in xyz yzx zxy; do

    INPUT_DIR1="${DIR_COMBI:0:1}"
    INPUT_DIR2="${DIR_COMBI:1:1}"
    INPUT_DIR3="${DIR_COMBI:2:3}"

    echo "$INPUT_DIR1 $INPUT_DIR2 $INPUT_DIR3"

    TARGET_INPUT="$LAMMP_IN_FOLDER/$INPUT_DIR1/$RUN_TYPE""_$CHAINNAME.in"
    mkdir -p "$LAMMP_IN_FOLDER/$INPUT_DIR1/"
    cp "$INPUT_FILE_TO_DIRECTION" "$TARGET_INPUT"

    ESCAPED_CHAINFILE=$(printf '%s\n' "$CHAINFILE" | sed -e 's/[\/&]/\\&/g')
    ESCAPED_CHAINNAME=$(printf '%s\n' "$CHAINNAME" | sed -e 's/[\/&]/\\&/g')
    ESCAPED_SCRATCH_OUT=$(printf '%s\n' "$SCRATCH_OUTFOLDER" | sed -e 's/[\/&]/\\&/g')

    # Replace input, output & direction
    sed -i "s/INPUT_FILE/$ESCAPED_CHAINFILE/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_NAME/$ESCAPED_CHAINNAME/gi" "$TARGET_INPUT"
    sed -i "s/SCRATCH_OUT/$ESCAPED_SCRATCH_OUT/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_DIR1/$INPUT_DIR1/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_DIR2/$INPUT_DIR2/gi" "$TARGET_INPUT"
    sed -i "s/INPUT_DIR3/$INPUT_DIR3/gi" "$TARGET_INPUT"

    # Write to run files
    { 
      echo "mkdir -p $OUTFOLDER/$INPUT_DIR1/" 
      echo "mkdir -p $SCRATCH_OUTFOLDER/$INPUT_DIR1/"
      echo "$EULER_CMD_PREFIX -W 120:00 -N -B -J \"$INPUT_DIR1-$RUN_TYPE-$CHAINNAME\" -o \"$OUTFOLDER/$INPUT_DIR1/$RUN_TYPE""_$CHAINNAME-run.out\" mpirun /cluster/home/betim/masters/lammps/build/lmp_mpi -suffix omp -in \"$TARGET_INPUT\" -l \"$SCRATCH_OUTFOLDER/$INPUT_DIR1/$RUN_TYPE""_$CHAINNAME-log.out\";"
    } >> "$EULER_RUN_FILE"
  done

  echo "Processed input-file: $CHAINFILE"
done

chmod u+x "$EULER_RUN_FILE"
