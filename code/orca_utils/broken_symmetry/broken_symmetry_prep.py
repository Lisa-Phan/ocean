"""
Orca broken symmetry helper script

Takes in a single point calculation input 
set up directory for broken symmetry job

take one input file,
create an opt directory, single point directory, and antiferromagnetic singlet directory

the rationale is that the single point wavefunction from single point calc
can be used either for the antiferromagnetic singlet or the high spin normal optimization state

ideally, when antiferromagnetic case work, then the regulat high spin opt would not be needed
but, sometimes there are stability issues, and thus getting an optimized geometry with 
high spin state first will help

single point ---- > high spin opt ---> singlet antiferro opt
and
single point ---- > singlet antiferro opt

input: 
    single point calculation input file
    job_name


    single point spin multiplicity strip
    high spin spin multiplicity strip
    singlet spin multiplicity strip

output: 
    three directory, copied orca submission template over and run job
    file names will follow the format: 
    sp/job_name_sp.inp
    opt_after_sp/job_name_opt.inp
    singlet_after_sp/job_name_singlet.inp
    singlet_after_opt/job_name_singlet_opt.inp

"""
import sys
import os


################ POSITIONAL VARS #####################

SP_INPUT = sys.argv[1]
JOB_NAME = sys.argv[2]

# full system charge
CHARGE = sys.argv[3]

# full system high spin multiplicity
HIGH_MULT = sys.argv[4]

# full system "antiferromagnetic singlet" multiplicity
SINGLET_MULT = sys.argv[5]

WORK_PLACE = sys.argv[6]

################# LOCAL VARS #######################
#POD Orca path
POD_ORCA = '/stor/work/YiLu/orca/orca_three_parts/orca'

#time for single point calc 
TIME_SP = '03:00:00'

#nodes for single point calc
NODES_SP = 10

#partition for single point calc
PARTITION_SP = 'gpu-a100-small'

def create_dirs():

    # create directories
    os.mkdir('sp')
    os.mkdir('highspin_after_sp')
    os.mkdir('singlet_after_sp')
    os.mkdir('singlet_after_highspin')

def create_input_files(single_point_input, output, gbw_name, charge, spin):
    """
    Take in the single point input file as a template
    create new opt input by reading in the gbw

    gbw_name: the name of the gbw file to read in, without .gbw extension
    charge: charge of the system
    spin: spin multiplicity of the system
    output: the name of the output file
    single_point_input: the name of the single point input file
    
    """
    # create input files
    with open(single_point_input, 'r') as f:
        sp_lines = f.readlines()
    
        # this is the high spin input optimization
        with open(output, 'w') as output:
            for index, line in enumerate(sp_lines):
                
                #add the OPT and MOREAD keywords to the header
                if index == 0:
                    output.write('! OPT MOREAD ' + line.strip('!') + '\n')
                    output.write('\n')
                    output.write('%moinp "{}.gbw"\n'.format(gbw_name))
                elif "* xyz" in line: 
                    #change the charge multiplicity strip
                    output.write("* xyz " + charge + ' ' + spin + '\n')
                else:
                    output.write(line)


def main():
    create_dirs()

    #create input for highspin opt, read sp gbw
    create_input_files(SP_INPUT, 
                       os.path.join('highspin_after_sp', f'{JOB_NAME}_highspin_after_sp.inp'), 
                       gbw_name = JOB_NAME + '_sp', 
                       charge = CHARGE, 
                       spin = HIGH_MULT)
    
    #create input for singlet opt, read sp gbw
    create_input_files(SP_INPUT, 
                       os.path.join('singlet_after_sp', f'{JOB_NAME}_singlet_after_sp.inp'),
                       gbw_name = JOB_NAME + '_sp', 
                       charge = CHARGE, 
                       spin = SINGLET_MULT)
    
    #create input for singlet opt after high spin opt
    create_input_files(SP_INPUT, 
                       os.path.join('singlet_after_highspin', f'{JOB_NAME}_singlet_after_highspin.inp'),
                       gbw_name = JOB_NAME + '_highspin_after_sp', 
                       charge = CHARGE, 
                       spin = SINGLET_MULT)
    
    #move sp input to sp directory
    os.rename(SP_INPUT, os.path.join('sp', f'{JOB_NAME}_sp.inp'))

    #copy orca submission template to each directory
    if WORK_PLACE == 'POD': 
        #submit sp calc with tmux, orca_sub is custom function
        os.system(f'tmux new-session -d -s {JOB_NAME} "{POD_ORCA} sp/{JOB_NAME}_sp.inp > sp/{JOB_NAME}_sp.out"')
    
    
    elif WORK_PLACE == 'TACC':
        #copy template file over and submit, redundant for now
        os.system('cp $TEMP/orca_submission_temp.sh sp/')
        os.system('cp $TEMP/orca_submission_temp.sh highspin_after_sp/')
        os.system('cp $TEMP/orca_submission_temp.sh singlet_after_sp/')
        os.system('cp $TEMP/orca_submission_temp.sh singlet_after_highspin/')

        #submit single point calc
        os.system('sbatch -t {} -p {} -n {} -J {} sp/orca_submission_temp.sh sp/{}_sp.inp'.format(TIME_SP, 
                                                                                                  PARTITION_SP, 
                                                                                                  NODES_SP, 
                                                                                                  JOB_NAME, 
                                                                                                  JOB_NAME))

    
        

    
if __name__ == "__main__":
    main()
