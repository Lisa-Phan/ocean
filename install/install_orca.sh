#!/bin/bash
#===========================================
# Note of commands used in installing orca
#===========================================

#go on development node
idev -t 2:00:00
module load gcc
cdw

#install libpciaccess
wget https://www.x.org/pub/individual/lib/libpciaccess-0.17.tar.gz
tar -xzvf libpciaccess-0.17.tar.gz
cd libcpiacces-0.17
configure
make

#install openmpi
#https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html#binary-packages
#wget the link from this site, https://www-lb.open-mpi.org/software/ompi/v4.1/ in a directory at $HOME
cdw
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz
tar -xvzf openmpi-4.1.5.tar.gz 
cd openmpi-4.1.5
./configure --prefix=$HOME/openmpi/v3/openmpi-4.1.1 --enable-shared --enable-static --disable-builtin-atomics --disable-mpi-fortran
make -j4 all 
make -j4 install

#modify path variable by changing .bashrc
export PATH=/work/09069/dhp563/ls6/openmpi/v2/openmpi-4.1.5/bin/:$PATH
export LD_LIBRARY_PATH=/work/09069/dhp563/ls6/openmpi/v2/openmpi-4.1.5/lib/:$LD_LIBRARY_PATH
#### Test if the installation works with the example script included

#restart and refresh
#then source.bashrc
source .bashrc
#check $PATH

#install orca
#this requires a registration step and is therefore a bit cumbersome...
# orca's main website https://orcaforum.kofo.mpg.de/app.php/portal
#### Download and transfer .xz file onto TACC
#### These files are a bit big. the combined content of the three part download can be found at 
# /work/09069/dhp563/ls6/orcaa/all_three
# orca works fine out-of-the-box, but parallelization requires openmpi

#### TEST SCRIPT: named 'water' #### 
! B3LYP def2-SVP Opt

%PAL NPROCS 6 END

*xyz 0 1
O        0.000000000      0.000000000      0.000000000
H        0.000000000      0.759337000      0.596043000
H        0.000000000     -0.759337000      0.596043000

*
#### END OF TEST SCRIPT ####


#get full path to orca executable
#place this under .bashrc custom variable
ORCA=/work/09069/dhp563/ls6/orcaa/all_three/orca

#restart terminal 
#test run orca
idev -t 2:00:00
#unpublished module that tacc support recommend using
module use /scratch/projects/compilers/nvhpc_23.7/modulefiles
module load nvhpc
