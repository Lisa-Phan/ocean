"""
Mini script that sanitize input coordinate orca file

Only reads the coordinate block sandwiched by 
*xyz and *.

Remove duplicated atoms, with first occurence kept.
"""
import sys

def get_xyz_lines(file):
    """
    return *xyz charge multiplicity line
    """
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('*xyx') or line.startswith('* xyx'):
                return line
                
    
def read_xyz(file):
    """
    Read only the lines that are between *xyz and *.
    """
    with open(file, 'r') as f:
        #print everything until the next *
        for line in f:
            if line.startswith('*xyx') or line.startswith('* xyx'):
                yield line
                line = next(f)
                while not line.startswith('*'):
                    yield line
                    line = next(f)
                
def get_atom_dictionary(text):
    """
    Get the atom and coordinate from the text
    """
    atom_dictionary = {}
    for idx, line in enumerate(text):
        atom, x, y, z = line.strip("\n").split()
        atom_dictionary[idx] = {atom, (x, y, z)}
    return atom_dictionary

def check_coord_dupes(atom_coord):
    """
    remove dupes 
    """
    tracker = []
    for key, value in atom_coord.items():
        if value not in tracker:
            tracker.append(value)
        else:
            print('Duplicate found at index {}'.format(key))
            continue
    return tracker

def print_non_dupe(non_dupe: list):
    """
    print to stdout the non-duplicate atoms
    """
    for value in non_dupe:
        atom, coord = value
        print(atom)
        print(coord)
        # try:
        #     x, y, z = coord
        # except ValueError:
        #     print('ValueError: {}'.format(coord))
        #     continue
        # print('  {}    {}    {}    {}\n'.format(atom, x, y, z))

file = sys.argv[1]
#read_xyz returns a generator
text = list(read_xyz(file))

#read text into atom and coordinate
atom_coord = get_atom_dictionary(text)

non_dup_ind = check_coord_dupes(atom_coord)
print_non_dupe(non_dup_ind)

