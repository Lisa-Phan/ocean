"""
4/26/24

Script to set up broken symmetry wavefunction for NEB using 
initial interpolated images from reactant to product 

Input 
an *initial_path_trj.xyz file, one image for each step
index of atom to flip spin on
"""
import os
import textwrap
from itertools import zip_longest
import subprocess


XYZ_FILE = r"/scratch/09069/dhp563/ash/test/orca_neb/orca_read_gbw_long_0_11/knarr_current.xyz"
ORCA_SUB_FILE = r"/home1/09069/dhp563/cdtools_remote.sh"


FLIP_ATOM = 28
CHARGE = 0
MULTIPLICITY = 11
JOB_NAME = 'NEB_broken_sym'
NPROCS = 30
IMAGE_COUNT = 10

#ORCA job template for single point broken symmetry calculation


TEMPLATE = textwrap.dedent(f"""\
! BP86 DEF2-SVP RIJCOSX KDIIS TightSCF UKS

%pal nprocs {NPROCS} end
 
%basis
    newgto Fe "def2-TZVP" end
    newgto O "def2-TZVP" end
    newgto N "def2-TZVP" end
end

%scf 
    FlipSpin {FLIP_ATOM} #Flip spin one iron center
    FinalMs 0
    convergence tight
    CNVZerner True
    end

* xyz {CHARGE} {MULTIPLICITY}
                           """)



def output_writer(block: str,  name: str, template = TEMPLATE) -> None:
    """
    Combine template and xyz block to output file
    """
    with open(name, 'w') as file:
        file.write(template)
        file.write(''.join(block))
        file.write("*\n")
        file.close()

def read_write_xyz_file(filename: str) -> None:
    """
    Read input file and return a list of coordinate blocks
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
        atom_count = int(lines[0].strip())
        xyz_block_size = atom_count + 2
        
        #reset file pointer
        file.seek(0)
        counter = 0 #image counter
        for next_n_lines in zip_longest(*[file] * xyz_block_size):
            coordinate_block = ''
            for index, line in enumerate(next_n_lines):                
                if index in [0, 1]: #atom count and comment line, continue
                    continue
                else:
                    #clean up 
                    atom, x, y ,z = line.split()
                    #append to text block and return
                    coordinate_block += f"{atom} {x} {y} {z}\n"
            #write to file
            output_writer(coordinate_block, f'{JOB_NAME}_im{counter}.inp' )
            counter += 1  


#TODO: Currently job submission with setup_orca_jobs doesn't work because of environment variable issues
#Temp fix is to run command directly from terminal
COMMAND_BLOCK = textwrap.dedent("""
for i in {1..9}; 
    do cd $i; 
    sbatch -t 1:00:00 -p gpu-a100-small -N 1 -n 30 -J $i cdtools_remote.sh NEB_broken_sym_im$i.inp ;
    cd ../; 
    done
    """)                        
                                                            

def setup_orca_jobs(num_images: int) -> None:
    """
    Copy files and set up working directory with subprocess
    Each image has one directory, 
    current_dir|
               |_restart_image_sp| 
                                 |_0
                                 |_1
                                 |...
    """
    current_script_path = os.path.realpath(__file__)
    current_script_dir  = os.path.dirname(current_script_path)

    for i in range(num_images):
        mkdir = ['mkdir', '-p', f'restart_image_sp/{i}']
        move_file = ['mv', f'{JOB_NAME}_im{i}.inp', f'restart_image_sp/{i}']
        copy_orca_subfile = ['cp', ORCA_SUB_FILE, f'restart_image_sp/{i}']
        run_orca = [ 'sbatch', '-t', '1:30:00',
                     '-N', '1',
                     '-n', f'{NPROCS}',
                     '-p', 'normal',
                     '-J', f'image_{i}_{JOB_NAME}', 
                     f"{current_script_dir}/restart_image_sp/{i}/cdtools_remote.sh", 
                     f'{current_script_dir}/restart_image_sp/{i}/{JOB_NAME}_im{i}.inp']

        subprocess.run(mkdir)
        subprocess.run(move_file)
        subprocess.run(copy_orca_subfile)
        subprocess.run(run_orca)

if __name__ == '__main__':
    read_write_xyz_file(XYZ_FILE)
    setup_orca_jobs(IMAGE_COUNT)

    

            



            

