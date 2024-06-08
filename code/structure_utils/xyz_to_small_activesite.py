"""
exhaustive search distance
"""
import numpy as np
from typing import Callable
from biotite.structure import distance, stack, check_duplicate_atoms
from biotite.structure import AtomArray, AtomArrayStack
from biotite.structure.io.xyz import *

def get_neighbor(structure: AtomArray, target: AtomArray, radius: float) -> AtomArray:
    """
    Find all atoms in atomArray that are within radius of specified target atom(s)
    """
    neighbors = AtomArray(0)
    for atom in target:
        distance_set = distance(structure, atom)
        neighbors += structure[distance_set <= radius]
    #check duplicate return the indices of duplicated atoms
    neighbors[~check_duplicate_atoms(neighbors)]    
    return neighbors

def get_active_site(structures: AtomArrayStack, radius: float, get_target: Callable) -> AtomArrayStack:
    """
    Find all atoms in atomArray that are within radius of specified target atom(s)
    target is a function that takes a structure and return the target atoms as an AtomArray
    """
    structure_length, num_atoms = len(structures), len(structures[0])

    neighbors_stack = []
    for index, structure in enumerate(structures):
        site_atoms = get_target(structure)
        neighbors = get_neighbor(structure, site_atoms, radius)
        neighbors_stack.append(neighbors)
    return neighbors_stack

def get_target(structure: AtomArray) -> AtomArray:
    """
    Get the target atoms for a given structure
    """
    return structure[np.isin(structure.element, ['Fe', 'FE'])]

def write_atom_list(atom_list: list, outfile: str, header: set):
    """
    Write a list of atoms to a file

    header is a set of two lists, coming from the original xyz
    first list is for the number of atoms in the xyz
    second list is the energy line

    """
    with open(outfile, 'w') as xyz_outfile:
        for index, atom_array in enumerate(atom_list):
            xyz_outfile.write(f"{header[0][index]} \n")
            xyz_outfile.write(f"{header[1][index]} \n") 
            for atom in atom_array:
                line = " " + str(atom.element)
                line += " {:>11.{}f}".format(atom.coord[0], 6)
                line += "    {:>11.{}f}".format(atom.coord[1], 6)
                line += "    {:>11.{}f}".format(atom.coord[2], 6)
                line += " " + "\n"
                xyz_outfile.write(line)
    xyz_outfile.close()

def test():
    xyz = XYZFile()
    xyz.read('../redo_NPT.06/NPT-pos-1.xyz')
    structure_stack = xyz.get_structure()
    active_site = get_active_site(structure_stack, 5, get_target)
    header = xyz.get_header()
    #print(header)
    #write to outfile
    write_atom_list(active_site, 'big_active_site.xyz', header)
    
if __name__ == "__main__":
    test()
    


