"""
5/24/2024

QM/MM space builder helper
Script to print out indices of atoms to use in constructing QM/MM space

Input: 
    - pdb file created from input prmtop and inpcrd file by ambpdb
    - index/indices of atoms 
    - search radius
    - rules trim option (current options: -c (c-calpha))
    
Construct PDB file from Amber prmtop and inpcrd file
Find residues within given radius from given atom indices
Build QM/MM space from the residues found
"""

import biotite
import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArray, Atom, distance, array

#==============================================================================
# Input and Global Variables
#==============================================================================

#INPCRD = r"/scratch/09069/dhp563/ash/test/amber/prmtop_inpcrd/prmtop_inpcrd_relaxed/6YD0_solv_tip3.inpcrd"
#PRMTOP = r"/scratch/09069/dhp563/ash/test/amber/prmtop_inpcrd/prmtop_inpcrd_relaxed/noflag.prmtop"
PDB= r"/scratch/09069/dhp563/ash/test/amber/prmtop_inpcrd/prmtop_inpcrd_relaxed/6YD0_solv_tip3.pdb"
ATOM_INDICES = [8155, 8154]  #iron indices from 6YD0
SPECIAL_RESIDUE = ['HH1', 'HH2', 'HH3']
DISTANCE = 4.0

#==============================================================================
# Functions
#==============================================================================

def get_pdb(pdb: str) -> AtomArray:
    """
    Read pdb file and return AtomArray
    """
    pdb_file = PDBFile.read(pdb).get_structure()[0]
    return pdb_file

def get_residue_index_within_distance(pdb: AtomArray, spec_atom_array: AtomArray, search_radius: float) -> list:
    """
    Iterate over pdb and get residues indices within given distance of speficied atom
    """
    nearby_atoms = []
    for spec_atom in spec_atom_array: 
        imm_nearby_atoms = [atom for atom in pdb if distance(atom, spec_atom) <= search_radius]
        nearby_atoms.extend(imm_nearby_atoms)
    
    #return non-duplicate residue indices
    return list(set([atom.res_id for atom in nearby_atoms]))

def get_atoms_from_residue_index(pdb: AtomArray, residue_index: list) -> AtomArray:
    """
    Read PDB file and return an AtomArray of atoms from the specified residue indices
    """
    return array([atom for atom in pdb if atom.res_id in residue_index])

def trim_atom_array(residues: AtomArray, rules: str = 'Ca-Cb', keep_glycine = False, keep_box_water = False) -> AtomArray:
    """
    Trim AtomArray based on rules
    
    C-Ca: remove all backbone atoms (C, O, N) to create a QM space selection with carbon-carbon boundary
    this doesn't work with glycine
    """
    if not keep_box_water:
        residues = array([atom for atom in residues if atom.res_name not in ['WAT', 'HOH']])

    if not keep_glycine:
        residues = array([atom for atom in residues if atom.res_name != 'GLY'])

    if rules == 'Ca-Cb':
        no_backbone_residues = array([atom for atom in residues if atom.atom_name not in ['C', 'O', 'N', 'CA', 'H', 'HA']])
        return no_backbone_residues
    
    else:
        raise ValueError('Invalid rule')
    
def print_atom_indices_by_element(atom_array: AtomArray) -> list:
    """
    Print atom indices of specified element
    """
    #get all elements in array
    elements = set([atom.elements for atom in atom_array])
    
    element_dict_indices = {}

    for element in elements:
        element_dict_indices[element] = [atom.atom_index for atom in atom_array if atom.element == element]
        print(f"Indices of {element}: \n")
        print(element_dict_indices[element])

    return element_dict_indices
        

#==============================================================================
# Main
#==============================================================================

if __name__ == '__main__':
    print('Reading pdb file... \n') 
    pdb = get_pdb(PDB)

    
    #create atom indices annotation
    pdb.set_annotation('atom_index', np.arange(0, len(pdb)))

    #check special residue
    special_residue = array([atom for atom in pdb if atom.res_name in SPECIAL_RESIDUE])
    print("Special residue: \n")
    print(special_residue)

    #get atoms from specified indices
    spec_atoms = pdb[ATOM_INDICES]
    print("Atom from specified indices: \n")
    print(spec_atoms)

    print(f"Getting residues within {DISTANCE} of specified atom indices... \n")
    #get residues within given distance of specified atom indices
    residue_index = get_residue_index_within_distance(pdb, spec_atoms, DISTANCE)

    #get active site atoms
    active_site_atoms = get_atoms_from_residue_index(pdb, residue_index)

    #trim atom to create qm space
    trimmed_array = trim_atom_array(active_site_atoms, 
                                    rules = 'Ca-Cb', 
                                    keep_glycine = False)

    print(f"Active site atoms: /n ")
    print(trimmed_array)

    #get complete active site by combining special residue and trimmed array
    complete_active_site = special_residue + trimmed_array

    print(f"Complete active site: \n")
    print(complete_active_site)  

    print("Atom indices to use for QM space: \n")
    print(list(complete_active_site.get_annotation('atom_index')))    

    print_atom_indices_by_element(complete_active_site)
