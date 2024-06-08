"""
6/2/24

Template code from ASH tutorial for
NEB calculation with selective active region
Adapt code to re-run original NEB calculation with ferromagnetic coupling

"""

import textwrap 
from ash import *

##################################
# GLOBAL VARIABLES
################################

RELAXED_REACTANT = r"/scratch/09069/dhp563/ash/test/orca_neb/pmul_neb/pmul_file_collection_troubleshoot/broken_sym_NEB/reactant.xyz"
RELAXED_PRODUCT = r"/scratch/09069/dhp563/ash/test/orca_neb/pmul_neb/pmul_file_collection_troubleshoot/broken_sym_NEB/product.xyz"


reactant = Fragment(xyzfile=RELAXED_REACTANT, conncalc=False, charge = 0, mult = 1)
product = Fragment(xyzfile=RELAXED_PRODUCT, conncalc=False, charge = 0, mult = 1)
#NEB-CI job. Final saddlepoint structure stored in new object "Saddlepoint"


#create Orca theory from simple input and blocks
input = "!BP86 def2-SVP TightSCF"
block = textwrap.dedent("""
    %basis
        newgto Fe "def2-TZVP" end
        newgto O "def2-TZVP" end
        newgto N "def2-TZVP" end
    end
                        
    %scf
        CNVZerner True
    end
                                            
    """)

theory_obj = ORCATheory(orcasimpleinput = input, 
                        orcablocks = block, 
                        numcores = 1, 
                        orcadir = r"/work/09069/dhp563/ls6/orcaa/all_three")


Saddlepoint = interfaces.interface_knarr.NEB(reactant = reactant, 
                                             product = product, 
                                             theory = theory_obj, 
                                             images=8, 
                                             CI=True,  
                                             runmode = 'parallel', 
                                             numcores = 10, 
                                             restart_file = r"/scratch/09069/dhp563/ash/test/orca_neb/pmul_neb/pmul_file_collection_troubleshoot/initial_input_NEB/NEB_longer_Fe_MEP_trj.xyz",
                                             mofilesdir = r"/scratch/09069/dhp563/ash/test/orca_neb/pmul_neb/pmul_file_collection_troubleshoot/renamed_gbw_file_0_1")

Saddlepoint.print_system(filename='saddlepoint_PQ.ygg')
