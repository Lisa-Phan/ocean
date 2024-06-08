#!/bin/bash
#SBATCH -o %x.out      # Name of stdout output file (%j expands to job ID)
#SBATCH -e %x.err      # Name of stderr output file (%j expands to job ID)
#SBATCH -A CHE23010
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=dhp563@my.utexas.edu


eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate ash_reinstall

source /home1/09069/dhp563/orca_env_config.sh

#CDTools to work with tmp files at remote /tmp directory
export CDTools=/scratch/tacc/apps/CDTools
export PATH=${PATH}:${CDTools}/bin

#input file is first argument, output file is inputfile name with extension replace with .out
INPUT_FILE=$1
JOB_NAME=${INPUT_FILE%.*}
OUTPUT_FILE=${INPUT_FILE%.*}.out
COLLECTING_DIR=$(pwd)

#distribute files, with inputdir being the current dir
distribute.bash $(pwd)

wait

#run job
python $INPUT_FILE > $OUTPUT_FILE

wait

#collect file
collect.bash /tmp/outputdir $COLLECTING_DIR
