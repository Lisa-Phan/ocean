#!/bin/bash

export ORCA=/work/09069/dhp563/ls6/orcaa/all_three/orca

module use /scratch/projects/compilers/nvhpc_23.7/Linux_x86_64/23.7/comm_libs/12.2/hpcx/hpcx-2.15/modulefiles/
module use /scratch/projects/compilers/nvhpc_23.7/modulefiles/
module load nvhpc-hpcx-cuda12
module load hpcx
module unload intel impi
module list

#get custom functions to work with orca
source /scratch/09069/dhp563/workflow_code_mcpb_cp2k/anotherpb/script/minibash/func.sh

#$ORCA water > out_again
