"""
12/24/2024
Lisa P. 

Indexing atom for QMMM for structure of 6YDI and 3DHI
"""
################
# IMPORTS
################
from ash import *

################
# RUN SPECS
################
NUM_CORES = 8

#Change this to the new orca
ORCA_DIR = r"/stor/home/dhp563/orca_6_0_0/orca_6_0_0"

################
# PDB and XML files
################

#snapshot chosen from MD for 6YDI
#SOLVATED_PDB = r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_6YDI/step3_6YDI_NVT_lastframe_imaged_del_wat.pdb"

#snapshot chosen from MD for 3DHI
SOLVATED_PDB = r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_3DHI/last_snapshot_3DHI_NVT_imaged.pdb"

#xml containing information on the metal ions and water
SPECIAL_RES_XML = r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_6YDI/fe_o2_nowat.xml"

#indices built using prep_qmmm_active_site.py, 6 angstrom C-C around diiron center
ACTIVE_SITE_INDEX = r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_3DHI/3DHI_active_site_indices.txt"

EXTRA_BASIS_INDEX = r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_3DHI/3DHI_active_site_indices_extrabasis.txt"
################
# READ ATOM INDICES
################

#====== for 3DHI =======#
ATOM_INDICES = read_intlist_from_file(ACTIVE_SITE_INDEX)
print("Atom indices for QM region: ", ATOM_INDICES)

EXTRABASIS_ATOM_INDICES = read_intlist_from_file(EXTRA_BASIS_INDEX)
print("Atom indices for extra basis: ", EXTRABASIS_ATOM_INDICES)

###########################
# C. ORCA and QM/MM SPEC
###########################
# 1. Schulz paper setup, cheaper BP86 that is generally good
input_line = "! BP86 def2-SVP tightscf uks"
extrabasis = 'def2-TZVP'

# 2. Ragnar's recommended r2SCAN, proven reliable for spin coupled iron system like FeS
# input_line = "! r2SCAN def2-SVP tightscf uks"
# extrabasis = 'def2-TZVP'


#Charge and mult of QM-region  
#-4 charge from carboxylates, 
#+4 from 2Fe(II)
charge = 0
mult = 9 #per revision, Fe(II)-Fe(II) is ferromagnetically coupled, 4 unpaired electron per iron, 8 total

###########################
# MAIN
###########################

frag = Fragment(pdbfile = SOLVATED_PDB)

omm = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", SPECIAL_RES_XML], 
                            pdbfile = SOLVATED_PDB, 
                            numcores = NUM_CORES, 
                            autoconstraints = 'None',
                            rigidwater = False)

qmregion = ATOM_INDICES
extrabasisatoms = EXTRABASIS_ATOM_INDICES

#QM theory. Using xTB as a test, can change to ORCA
qm = ORCATheory(orcadir = ORCA_DIR,
                 numcores = NUM_CORES, 
                 orcasimpleinput = input_line,
                 brokensym = True,
                 extrabasisatoms = extrabasisatoms,
                 extrabasis = extrabasis,
                 HSmult = 11, 
                 atomstoflip = [15831])


#Combining all value list to a big list
# primary_atoms = [atom_index for atom_key in PRIMARY_SPHERE for atom_index in PRIMARY_SPHERE[atom_key]]
# secondary_atoms = [atom_index for atom_key in SECOND_SPHERE for atom_index in SECOND_SPHERE[atom_key]]


#Defining QM/MM Theory
qm_mm = QMMMTheory(qm_theory=qm, 
                   mm_theory=omm,
                   qmatoms=qmregion,
                   fragment=frag,
                   qm_charge=charge,
                   qm_mult=mult)

#Run MD or Opt
#MolecularDynamics(theory=qm_mm, fragment=frag, simulation_steps=100, temperature=100)
Optimizer(theory=qm_mm,
          maxiter=1000,
          fragment=frag,
          ActiveRegion=True,
          actatoms=qmregion)
