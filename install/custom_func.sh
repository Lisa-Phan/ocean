#!/bin/bash

save_code () { cp "$1" /work/09069/dhp563/ls6/code_lib/ ; }

disk_usage () { du -hs -- * .[^.]* | sort -h ; }

refresh () { 
rm "$1"; vim "$1" 
}

spin_check() {
    grep -n -i 'mulliken population' -A $1 *out | grep -i 'fe'
}

clean_orca () { rm *.{tmp,txt,opt,gbw,X,xyz} ; }

cq () {
    echo " gpu-a100-dev: $(squeue -p gpu-a100-dev | wc -l)"
    echo " gpu-a100-small: $(squeue -p gpu-a100-small | wc -l)"
    echo " development: $(squeue -p development | wc -l)"
    echo " vm-small: $(squeue -p vm-small | wc -l)"
}


## Sample file format AF-O15552-F1-model_v4.pdb
## first var is input file id, second var is directory of downloaded pdb location
download_AF () {
	HOST='https://alphafold.ebi.ac.uk/files'
	for i in `cat $1` ; do wget "$HOST/AF-${i}-F1-model_v4.pdb" -P $2 ; sleep 0.1 ; done
}

createGBWDirectory () {
#6/6/24
#Mini helper function to setup gbw directory to restart a halted NEB job

    NUM_IMAGES=$1  # first argument is the number of images, counting from 0 to NUM_IMAGES

    if [[ -d "current_images_gbw" ]]; then
        echo "gbw directory already exists"
        return
    else
        mkdir current_images_gbw
    fi

    # Use a C-style for loop to properly iterate from 0 to NUM_IMAGES
    for ((image_index=0; image_index<=$((NUM_IMAGES)); image_index++)); do
        echo $image_index
        cp Pooljob_image_$image_index/*.gbw current_images_gbw/current_image$image_index.gbw
    done
}





#-------- Variables --------#
ROSETTA_BIN=/work/09069/dhp563/ls6/rosetta_bin/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin
ORCA=/work/09069/dhp563/ls6/orcaa/all_three/orca
MOBIL=/work2/projects/utprojections/safegraph_data
PRODIGAL=/home1/09069/dhp563/prodigal/Prodigal/prodigal


#-------- Custom cd function --------#
alias amber='cd /scratch/09069/dhp563/amber_test'
alias mc='cd /work/09069/dhp563/mambaforge/envs/ambertools/lib/python3.11/site-packages/pymsmt/mcpb/'
alias whale='cd /work/09069/dhp563/ls6/orcaa/'
alias inter='cd /scratch/09069/dhp563/guess_intermediate/'
alias methane='cd /scratch/09069/dhp563/guess_intermediate/methane_6YD0/gaussian_job/'
alias introbio='cd /scratch/09069/dhp563/Intro_to_BioML/IntroBioMLProject/IntroBioMLProject/'

#-------- Custom copy ---------#
alias orc='cp $HOME/cdtools_remote.sh ./'
alias cash='cp /home1/09069/dhp563/cdtools_remote_ash.sh ./'


# -------- Custom file view function ---------#
alias mcview='vim /work/09069/dhp563/mambaforge/envs/ambertools/bin/MCPB.py'


# -------- Custom idev --------------- #
alias norm='idev -t 2:00:00 -A CHE23010'
alias small='idev -t 2:00:00 -A CHE23010 -p vm-small'
alias gpu='idev -t 2:00:00 -A CHE23010 -p gpu-a100-dev'
