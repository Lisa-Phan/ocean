#!/bin/bash
#Lisa P
#2025-02-05
#helper script to plot orca orbital with orca_plot
#input: 
# Positional argument 1: orbital file, which are .gbw, .loc, etc
# Positional argument 2: string representing orbital indices, e.g., "1-10, 12, 15-20" for orbitals 1 to 10, 12, and 15 to 20
# Positional argument 3: string representing alpha (0) or beta (1) orbitals
# String input for printf would need to be modified for other use case
# example: ./plot_orca_orbital.sh h2o.gbw "1-5" 0 to plot alpha orbital 1 to 5 of h2o.gbw



# Function to expand a string of numbers and ranges into a list of numbers
# Take a string like "1-5,7,10-12" and return "1 2 3 4 5 7 10 11 12"
expand_numbers() {
    local input="$1"
    local result=()
    
    IFS=',' read -ra parts <<< "$input"  # Split by commas

    for part in "${parts[@]}"; do
        if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            # If it's a range (e.g., 0-10)
            for ((i=${BASH_REMATCH[1]}; i<=${BASH_REMATCH[2]}; i++)); do
                result+=("$i")
            done
        else
            # If it's a single number (e.g., 12)
            result+=("$part")
        fi
    done

    echo "${result[@]}"  # Print the expanded list
}

#template string used was written for local spin analysis
#Positional argument 1 is the gbw file
#Positional argument 2 is the string representing orbital indices
#Positional argument 3 is the string representing alpha (0) or beta (1)
#iterate over the string and call the orca_plot command

get_orbital_iterate() {

    local orbital_file="$1"
    local orbital_index_string="$2"
    local operator="$3"
    echo "Orbital file: $orbital_file"
    orbitals=($(expand_numbers "$orbital_index_string"))

    for orbital in "${orbitals[@]}"; do
        echo Plotting orbital $orbital for operator $operator
        printf "1\n1\n2\n$orbital\n3\n$operator\n5\n7\n11\n12\n" | $ORCA_DIR/orca_plot $orbital_file -i
        wait
    done
}

get_orbital_iterate "$1" "$2" "$3" > orbital_by_indices.out 2> orbital_by_indices.err
