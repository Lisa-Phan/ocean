"""
1/3/2025

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

Rewrite to better handle multichain structures and updated biotite 1.0.1 syntax


"""
import sys
import biotite
import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArray, Atom, distance, array, infer_elements, get_residue_positions

#==============================================================================
# Input and Global Variables
#==============================================================================

#INPCRD = r"/scratch/09069/dhp563/ash/test/amber/prmtop_inpcrd/prmtop_inpcrd_relaxed/6YD0_solv_tip3.inpcrd"
#PRMTOP = r"/scratch/09069/dhp563/ash/test/amber/prmtop_inpcrd/prmtop_inpcrd_relaxed/noflag.prmtop"

#chain index file
#temp fix for dependency in get_chains only working properly when everything is sorted nicely chain-wise
#command to generate file:
# sed -i '/TER/d' $PDB
# cut -c21,22 $PDB > chain_index.txt
# sed -i '/^\s*$/d' chain_index.txt
CHAIN_INDEX_FILE = r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_6YDI/chain_index.txt"

#PDB file to use
PDB= r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_6YDI/step3_6YDI_NVT_lastframe_imaged.pdb"

# A. TRIM JOB
#residues that bypass the implemented selection rule
TRIM_RESIDUE = []
ATOM_INDICES = []


# B. LINKER JOB
LINKER_RESIDUE = ['GLU_99_A', 'ARG_230_A', 'ARG_131_A', 
                'ALA_195_A', 'ILE_130_A', 'GLN_125_A', 
                'HIS_132_A', 'LEU_229_A', 'THR_198_A', 
                'GLU_228_A', 'GLU_129_A', 'ASP_128_A', 
                'HIS_231_A', 'GLY_193_A', 'GLU_194_A']

#catalytically relevant water, bypass selection rule
SPECIAL_WATER = ['HOH_514_A']
SPECIAL_RESIDUE = ['FE_512_A', 'FE_513_A']

#search for however many angstroms around the specified atom indices
DISTANCE = 3.0

#output index file, to use in ASH QMMM indexing
INDEX_FILE_NAME = r"/stor/scratch/YiLu/dhp563/ash/sandbox/qmmm_fullcomplex_6YDI/6YDI_active_site_linker"

#selection rule for indexing special atoms with extra basis set, including 
# elements specified by letters and last character as number for distance
EXTRA_BASIS_SELECTION = 'N_O_3'

#==============================================================================
# Functions
#==============================================================================

def read_chain_index(file: str) -> list:
    """
    Read chain index file and return a list of chain indices
    """
    with open(file, 'r') as f:
        chain_index = f.readlines()
    #remove all empty lines
    chain_index = [line.strip() for line in chain_index if line.strip()]
    return chain_index

def get_pdb(pdb: str) -> AtomArray:
    """
    Read pdb file and return AtomArray
    Add a few annotations
        element: guess element type from atom name
        atom_index: counting from 0, assign each atom an index based on their ordering in the pdb file
    """
    pdb_structure = PDBFile.read(pdb).get_structure()[0]
    pdb_structure.set_annotation('element', infer_elements(pdb_structure))
    #create atom indices annotation
    pdb_structure.set_annotation('atom_index', np.arange(0, len(pdb_structure)))
    pdb_structure.set_annotation('chain', read_chain_index(CHAIN_INDEX_FILE))
    pdb_structure.set_annotation('res_index', get_residue_positions(pdb_structure, pdb_structure.atom_index))
    return pdb_structure

def get_residue_index_within_distance(pdb: AtomArray, 
                                           spec_atom_array: AtomArray, 
                                           search_radius: float) -> dict:
    """
    Iterate over pdb and get residues indices within given distance of speficied atom
    Return two things:
    A dictionary of chain: residue index
    A list of atoms to use for extrabasis, which are coordinating nitrogen or oxygen 
    """
    nearby_atoms = []
    for spec_atom in spec_atom_array: 
        imm_nearby_atoms = [atom for atom in pdb if distance(atom, spec_atom) <= search_radius]
        nearby_atoms.extend(imm_nearby_atoms)
    
    chain_res_dict = {}

    for atom in nearby_atoms:
        current_chain = atom.chain
        try:
            chain_res_dict[current_chain].append(atom.res_index)
        except KeyError:
            chain_res_dict[current_chain] = [atom.res_index] 

    for chain in chain_res_dict:
        chain_res_dict[chain] = list(set(chain_res_dict[chain]))

    return chain_res_dict

def get_extra_basis_atoms(qmatoms: AtomArray, spec_atoms:AtomArray ,selection: str) -> AtomArray:
    """
    Get atoms to use for extra basis set
    """
    #infer element and distance from selection string
    search_radius = selection.split('_')[-1]
    element = selection.split('_')[:-1]
    
    print("Searching for atoms with element: ", element, " and distance: ", search_radius)  

    extra_basis_atoms = []
    for spec_atom in spec_atoms:
        extra_basis_atoms.append(spec_atom)
        for atom in qmatoms:
            if atom.element in element:
                if distance(atom, spec_atom) <= float(search_radius):
                    if atom not in extra_basis_atoms:
                        extra_basis_atoms.append(atom)
    return array(extra_basis_atoms)

def get_special_residue(pdb: AtomArray, special_residue: list) -> AtomArray:
    """
    Get atoms from special residue
    """
    selected = []
    for residue in special_residue:
        #iterate over structure, each time adding atoms from the residue to the selected array
        res_name, res_num, chain = residue.split('_')
        special_residue = [atom for atom in pdb if 
                           atom.res_name == res_name and 
                           atom.res_id == int(res_num) and 
                           atom.chain == chain]
        if len(special_residue) == 0:
            raise ValueError(f"Special residue {residue} not found in pdb")    
        selected += special_residue

    return array(selected)

def get_atoms_from_chain_residue_index(pdb: AtomArray, residue_index: dict) -> AtomArray:
    """
    Read PDB file and return an AtomArray of atoms from the chain and specified residue indices
    
    pdb: the atom array to search over
    residue_index: a dictionary of chain: indices
    """
    atoms_in_index = []
    for chain in residue_index:
        for res_index in residue_index[chain]:
            atoms_in_index.extend([atom for atom in pdb if atom.res_index == res_index and atom.chain == chain])
    return array(atoms_in_index)

def manual_sort_continuous(iterable):
    """
    stack overflow code for sorting continuous intervals
    """
    iterable = iter(iterable)
    start = end = next(iterable)
    for index in iterable:
        if index != end + 1:
            yield (start, end)
            start = index
        end = index
    yield (start, end)


def trim_atom_linker_rule(residue_array : AtomArray, 
                          linker_residue = LINKER_RESIDUE, 
                          special_residue = SPECIAL_RESIDUE, 
                          special_water = SPECIAL_WATER) -> AtomArray:
    """
    For given residues, include indices with linker atoms
    """
    #make a chain : residue index dictionary from linker_residue
    chain_res_dict = {}
    for atom in linker_residue:
        res_name, res_num, chain = atom.split('_')
        try:
            chain_res_dict[chain].append(int(res_num))
        except KeyError:
            chain_res_dict[chain] = [int(res_num)]

    print('chain_res_dict', chain_res_dict)

    #create a list of all necessary residues
    combined_residue_list = linker_residue + special_residue + special_water

    for chain in chain_res_dict:
        res_ind_to_search = chain_res_dict[chain]
        #sort all residues of the same chain
        residues = iter(sorted(res_ind_to_search))
        start_res, end_res = zip(*manual_sort_continuous(residues))
        print('start res:', start_res)
        print('end res:', end_res)
        atom_indices = []

        #pick out C and O atoms for residues at n - 1 
        for start in start_res:
            print('start:', start)
            tail_backbone_index = int(start) - 1
            tail_atoms = [atom for atom in residue_array if atom.res_id == tail_backbone_index and atom.atom_name in ['C', 'O'] and atom.chain == chain]
            print('tail atoms:', tail_atoms)
            atom_indices.extend(tail_atoms)    

        #include all mentioned atoms in linker and special residue
    for residue in combined_residue_list:
        res_name, res_num, chain = residue.split('_')
        atom_indices.extend([atom for atom in residue_array if atom.res_name == res_name and atom.res_id == int(res_num) and atom.chain == chain])
        
    return array(atom_indices)


def trim_atom_array(residues: AtomArray, rules: str = 'Ca-Cb', keep_glycine = False, keep_box_water = False) -> AtomArray:
    """
    Trim AtomArray based on rules
    
    C-Ca: remove all backbone atoms (C, O, N) to create a QM space selection with carbon-carbon boundary
    this doesn't work with glycine
    keep_box_water argument assumes bulk water has different naming than active site water
    """
    if not keep_box_water:
        residues = array([atom for atom in residues if atom.res_name not in ['WAT', 'HOH']])

    if not keep_glycine:
        residues = array([atom for atom in residues if atom.res_name != 'GLY'])

    if rules == 'Ca-Cb':
        no_backbone_residues = array([atom for atom in residues if atom.atom_name not in ['C', 'O', 'N', 'CA', 'H', 'HA']])
        return no_backbone_residues
    
    elif rules == 'linker':
        #Special trim rule based on Schulz et al's model builing
        #For specified residue ranges, cut the larger atom at the C-N bond
        #and the smaller atom at the C-C bond
        #PDB residues naturally enda th te C-N bond, so we only need account for the extra C-C bond at the end 
        return trim_atom_linker_rule(residues)
    
    else:
        raise ValueError('Invalid rule')
    
def print_atom_indices_by_element(atom_array: AtomArray) -> list:
    """
    Print atom indices of specified element
    """
    #get all elements in array
    elements = set([atom.element for atom in atom_array])
    
    element_dict_indices = {}

    for element in elements:
        element_dict_indices[element] = [atom.atom_index for atom in atom_array if atom.element == element]
        print(f"Indices of {element}: \n")
        print(element_dict_indices[element])

    return element_dict_indices

def write_atom_indices_to_file(atom_indices_dict: dict, file: str):
    """
    Write atom indices to file
    """
    with open(file, 'w') as f:
        for element in atom_indices_dict:
            for atom_index in atom_indices_dict[element]:
                f.write(f"{atom_index}")
                f.write('\n')

#==============================================================================
# Main
#==============================================================================

def run_qm_space_builder_trim():
    print('Reading pdb file... \n') 
    pdb = get_pdb(PDB)


    #check special residue
    if len (SPECIAL_RESIDUE) > 0:
        special_residue = get_special_residue(pdb, SPECIAL_RESIDUE)
        print("Special residue: \n")
        print(special_residue)
    else:
        special_residue = None

    #get atoms from specified indices
    spec_atoms = pdb[ATOM_INDICES]
    print("Atom from specified indices: ")
    print(spec_atoms)
    print(spec_atoms.res_id)

    print(f"Getting residues within {DISTANCE} of specified atom indices... ")
    #get residues within given distance of specified atom indices
    residue_index = get_residue_index_within_distance(pdb, spec_atoms, DISTANCE)

    #get active site atoms
    active_site_atoms = get_atoms_from_chain_residue_index(pdb, residue_index)

    #get complete active site by combining special residue and active site atoms
    if len(special_residue) > 0:
        active_site_atoms = active_site_atoms + special_residue

    #trim atom to create qm space
    trimmed_array = trim_atom_array(active_site_atoms, 
                                    rules = 'Ca-Cb', 
                                    keep_glycine = False)

    #add water to the trimmed array and active site atoms
    if len(SPECIAL_WATER) > 0:
        special_water = get_special_residue(pdb, SPECIAL_WATER)
        trimmed_array = trimmed_array + special_water
        active_site_atoms = active_site_atoms + special_water

    print(f"Active site atoms, applied trimming rule:")
    print(trimmed_array)

    #get extrabasis atoms
    extrabasis_atoms = get_extra_basis_atoms(active_site_atoms, spec_atoms, EXTRA_BASIS_SELECTION)

    print(f"Complete active site, include all parts of residue:")
    print(active_site_atoms)  

    #write atom indices to text file
    # print("Printing QM atom indices")
    element_dict = print_atom_indices_by_element(trimmed_array)
    
    # print("Printing extrabasis atom indices")
    extrabasis_element_dict = print_atom_indices_by_element(extrabasis_atoms)
    
    write_atom_indices_to_file(element_dict, INDEX_FILE_NAME + '.txt')
    write_atom_indices_to_file(extrabasis_element_dict, INDEX_FILE_NAME + '_extrabasis.txt')

def run_qm_space_builder_linker():
    pdb = get_pdb(PDB)
    long_trim = trim_atom_linker_rule(pdb, LINKER_RESIDUE)
    print(long_trim)


if __name__ == '__main__':
    job_type = sys.argv[1]
    if job_type == 'trim':
        run_qm_space_builder_trim()
    elif job_type == 'linker':
        run_qm_space_builder_linker()
    else: 
        raise ValueError(f'Invalid job type {job_type}')
