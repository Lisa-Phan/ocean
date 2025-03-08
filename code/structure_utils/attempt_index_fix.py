"""
3/7/2025
Fix atom indexing with hex string

In the situation where atom indexing 
in PDB exceeds 99999, indexing number
does not increase and cause issues

Code attempt to patch this by making 
atom indexing turns in to hex string
after 99999

attempt to fix pymol's
auto reindexing issue where water atoms are 
grouped incorrectly. This is too hard to fix, the 
better solution is to not do read-write with water on pymol in the 
first place. 

"""
import numpy as np
import biotite
from biotite.structure.io.pdb import PDBFile
from biotite.structure import array

TEST_PDB = r"/stor/scratch/YiLu/dhp563/ash/qmmm/reduced_o2/6YDI_waterbridge/lastsnapshot_6YDI_reduced_waterbridge_short_dioxygen_nearsite.pdb"

PDB_TEMPLATE = (
    "{ATOM:<6}{ATOM_INDEX:>5} {ATOM_NAME:^4} {RESIDUE_NAME:<3} {CHAIN_ID:1}{RESIDUE_NUMBER:>4}    "
    "{X:>8.3f}{Y:>8.3f}{Z:>8.3f}{OCCUPANCY:>6.2f}{TEMP_FACTOR:>6.2f}          {ELEMENT:>2}"
)

OUTFILE = "test.pdb"

def hex_filter_rule(array) -> np.ndarray:
    """
    array: an array of integers
    return: an array of integers where 
    integers greater than 99999 are turned
    into hex string, VMD style indexing 

    100000 -> A0000
    100001 -> A0001
    etc 

    """
    new_array = np.zeros(len(array), dtype=object)
    for index, value in enumerate(array):
        if value > 99999:
            new_array[index] = hex(value - 100000 + 0xA0000)[2:].upper()
        else:
            new_array[index] = value
    return new_array

def get_atom_annotation(atom) -> dict:
    """
    Read atom information and return a dictionary
    atom: an Atom object from Biotite
    """
    field_dictionary = {
        "ATOM": "ATOM",
        "ATOM_INDEX": atom.atom_id,
        "ATOM_NAME": atom.atom_name,
        "RESIDUE_NAME": atom.res_name,
        "CHAIN_ID": atom.chain_id,
        "RESIDUE_NUMBER": atom.res_id,
        "X": atom.coord[0],
        "Y": atom.coord[1],
        "Z": atom.coord[2],
        "OCCUPANCY": 1.0,
        "TEMP_FACTOR": 0.0,
        "ELEMENT": atom.element
    }
    return field_dictionary

def read_pdb(pdb_file_path):
    """
    Read PDB file and return PDBFile object
    """
    pdb_file = PDBFile.read(pdb_file_path)
    pdb_structure = pdb_file.get_structure()[0]
    pdb_structure.set_annotation("atom_id", hex_filter_rule(np.arange(0, len(pdb_structure))))
    return pdb_structure

# def fix_water_indexing_issue(pdb_structure):
#     """
#     Chain H of structure has water atoms numbering being weird, where 5 water
#     molecules are assigned the same residue index. Fix by reindexing
    
#     reassign residue index 15 atoms at a time
#     For instance
#     O 0
#     O 0
#     O 0
#     O 0
#     O 0
#     H 0 
#     H 0


#     """
#     #first infer the number of times to repeat values in next_series_of_indices
#     #by counting the number of times O appears in the chain consecutively
#     chain_h = pdb_structure[pdb_structure.chain_id == "H"]
#     element_string = chain_h.get_annotation('element')
#     print(element_string) 
    
#     water_molecules = []

#     for index, atom in enumerate(chain_h):
#         if atom.atom_name == "O":
#             inferred_H1 = chain_h[index + 6]
#             inferred_H2 = chain_h[index + 11]
#             if inferred_H1.atom_name == "H1" and inferred_H2.atom_name == "H2":
#                 print("Found a water molecule")
#                 #add water molecule to list
#                 water_molecules.append([atom, inferred_H1, inferred_H2])
#             else:
#                 inferred_H1 = chain_h[index + 7]
#                 inferred_H2 = chain_h[index + 1]
#                 if inferred_H1.atom_name == "H1" and inferred_H2.atom_name == "H2":
#                     print("Found a water molecule")
#                     #add water molecule to list
#                     water_molecules.append([atom, inferred_H1, inferred_H2])
#                 else: 
#                     print(f'inferred atom1 {inferred_H1.atom_name} and atom2 {inferred_H2.atom_name} are not H1 and H2')
#     #unlist water molecules
#     water_molecules = [atom for sublist in water_molecules for atom in sublist]
#     print(water_molecules)

    # current_residue_index = 0 
    # next_series_of_indices = list(np.arange(current_residue_index, 
    #                                  current_residue_index + 6))*3  #repeat 3 times

    # #fix residue index, 15 atoms at a time
    # #create a copy of chain_h, modify this copy
    # fix_residue_id = []
    
    # for atom in chain_h:
    #     id = next_15_indices.pop(0)
    #     fix_residue_id.append(id)
    #     if len(next_15_indices) == 0:
    #         current_residue_index += 6
    #         next_15_indices = list(np.arange(current_residue_index, 
    #                                  current_residue_index + 6))*3
    
    # chain_h.set_annotation("res_id", np.array(fix_residue_id))
    # pdb_structure[pdb_structure.chain_id == "H"] = chain_h
    # print(pdb_structure)
    # return pdb_structure
    
def custom_write_pdb(field_dictionary):
    out_line = PDB_TEMPLATE.format(**field_dictionary)
    return out_line

def main():
    pdb_structure = read_pdb(TEST_PDB)
    with open (OUTFILE, "w") as f:
        for atom in pdb_structure:
            f.write(custom_write_pdb(get_atom_annotation(atom)) + "\n")

if __name__ == "__main__":
    main()
