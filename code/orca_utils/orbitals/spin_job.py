"""
1/30/25
Script for generating ORCA file to use for spin population analysis
IN_PROGRESS UNTESTED

INPUT: 
a directory of ASH result
a new directory for result

Script will do: 

1. make new directory
2. copy .gbw file to new directory, renaming it to former_orca.gbw
3. using orca.inp geometry block, write a new inp file to run localized spin analysis
based on give header

"""
import textwrap
import sys
import os



def get_spin_calc_block(block_type: str, gbw_inp: str, job_name: str) -> str:
    """
    return string header for ORCA input file
    
    block_tye: str, either DO_SPIN_CALC or READ_LOC
    gbw_inp: str, path to gbw file
    job_name: str, name of job, used to infer loc file name
    """
    DO_SPIN_CALC_BLOCK = textwrap.dedent(
        f"""! MOREAD noiter printmos printbasis UKS BP86 DEF2-SVP UCO UNO
        %PAL nprocs 4 end
        
        ! Engrad
        %pointcharges "orca.pc"
        
        %moinp {gbw_inp}
        
        %loc 
            locmet pm
        end
        """)
    
    READ_LOCALIZED_MOL_ORB_BLOCK = textwrap.dedent(
        f"""! MOREAD noiter printmos printbasis UKS BP86 DEF2-SVP KDIIS 
        %PAL nprocs 1 end                                 
        %moinp "{job_name}.loc" """)
    
    if block_type == "DO_SPIN_CALC":
        return DO_SPIN_CALC_BLOCK
    elif block_type == "READ_LOC":
        return READ_LOCALIZED_MOL_ORB_BLOCK
    else:
        raise ValueError("block_type must be either DO_SPIN_CALC or READ_LOC")



def main():
    ASH_JOB_DIR = sys.argv[1]
    NEW_JOB_DIR = sys.argv[2]
    #check if directory exists
    if not os.path.exists(ASH_JOB_DIR):
        raise FileNotFoundError(f"{ASH_JOB_DIR} does not exist")
    if not os.path.exists(NEW_JOB_DIR):
        os.makedirs(NEW_JOB_DIR)
    
    #copy gbw file to new directory
    gbw_file = os.path.join(ASH_JOB_DIR, "orca.gbw")
    os.system(f"cp {gbw_file} {NEW_JOB_DIR}/former_orca.gbw")
    
    #copy point charge file to new directory
    pc_file = os.path.join(ASH_JOB_DIR, "orca.pc")
    os.system(f"cp {pc_file} {NEW_JOB_DIR}/orca.pc")

    
    #extract xyz block from orca.inp

def test():
    get_spin_calc_block("DO_SPIN_CALC", "orca.gbw", "test")

if __name__ == "__main__":
    test()
    
