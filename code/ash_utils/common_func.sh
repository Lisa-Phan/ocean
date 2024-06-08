#!/bin/bash

#num images = number of images + 2, counting from 0
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
        echo copying .gbw file for image $image_index
        cp Pooljob_image_$image_index/*.gbw current_images_gbw/current_image$image_index.gbw
    done
}

