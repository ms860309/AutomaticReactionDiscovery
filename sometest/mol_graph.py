import networkx as nx
import numpy as np
from scipy.spatial import distance_matrix
import itertools
import os
from networkx.algorithms import isomorphism
import signal

def xyz_file_to_atoms(filename):
    """
    From an .xyz file get a list of atoms

    Arguments:
        filename (str): .xyz filename

    Returns:
        (list(autode.atoms.Atom)): Atoms
    """

    atoms = []

    if not os.path.exists(filename):
        raise XYZfileDidNotExist

    if not filename.endswith('.xyz'):
        raise XYZfileWrongFormat

    # Open the file that exists and should(!) be in the correct format
    with open(filename, 'r') as xyz_file:

        try:
            # First item in an xyz file is the number of atoms
            n_atoms = int(xyz_file.readline().split()[0])

        except IndexError:
            raise XYZfileWrongFormat

        # XYZ lines should be the following 2 + n_atoms lines
        xyz_lines = xyz_file.readlines()[1:n_atoms+1]

        for line in xyz_lines:

            try:
                atom_label, x, y, z = line.split()[:4]
                atoms.append(Atom(atomic_symbol=atom_label, x=x, y=y, z=z))

            except (IndexError, TypeError, ValueError):
                raise XYZfileWrongFormat

    return atoms


def atoms_to_xyz_file(atoms, filename, title_line='', append=False):
    """
    Print a standard .xyz file from a list of atoms

    Arguments:
        atoms (list(autode.atoms.Atom)): 
        filename (str):

    Keyword Arguments:
        title_line (str):
        append (bool):
    """

    with open(filename, 'a' if append else 'w') as xyz_file:
        print(len(atoms), title_line, sep='\n', file=xyz_file)

        for atom in atoms:
            x, y, z = atom.coord
            print(f'{atom.label:<3}{x:^10.5f}{y:^10.5f}{z:^10.5f}',
                  file=xyz_file)
    return None

class Atom:

    def __repr__(self):
        x, y, z = self.coord
        return f'[{self.label}, {x:.4f}, {y:.4f}, {z:.4f}]'

    def translate(self, vec):
        """
        Translate this atom by a vector

         Arguments:
             vec (np.ndarray): Shape = (3,)
          """
        self.coord += vec
        return None

    def rotate(self, axis, theta, origin=None):
        """Rotate this atom theta radians around an axis given an origin

        Arguments:
            axis (np.ndarray): Axis to rotate in. shape = (3,)
            theta (float): Angle in radians (float)

        Keyword Arguments:
            origin (np.ndarray): Rotate about this origin. shape = (3,)
                                 if no origin is specified then the atom
                                 is rotated without translation.
        """
        # If specified shift so that the origin is at (0, 0, 0), apply the
        # rotation, and shift back
        if origin is not None:
            self.translate(vec=-origin)

        # Normalise the axis
        axis = np.asarray(axis)
        axis = axis / np.linalg.norm(axis)

        # Compute the 3D rotation matrix using
        # https://en.wikipedia.org/wiki/Euler–Rodrigues_formula
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        rot_matrix = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                               [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                               [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

        # Apply the rotation
        self.coord = np.matmul(rot_matrix, self.coord)

        if origin is not None:
            self.translate(vec=origin)

        return None

    def __init__(self, atomic_symbol, x=0.0, y=0.0, z=0.0):
        """
        Atom class. Centered at the origin by default

        Arguments:
            atomic_symbol (str): Symbol of an element e.g. 'C'

        Keyword Arguments:
            x (float or str): x coordinate in 3D space (Å)
            y (float or str): y coordinate in 3D space (Å)
            z (float or str): z coordinate in 3D space (Å)
        """
        assert atomic_symbol in elements

        self.label = atomic_symbol
        self.coord = np.array([float(x), float(y), float(z)])


elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
            'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
            'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
            'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
            'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La',
            'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
            'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
            'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
            'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md',
            'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
            'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

# A set of reasonable valances for anionic/neutral/cationic atoms
valid_valances = {'H': [0, 1],
                  'B': [3, 4],
                  'C': [2, 3, 4],
                  'N': [2, 3, 4],
                  'O': [1, 2, 3],
                  'F': [0, 1],
                  'Si': [2, 3, 4],
                  'P': [2, 3, 4, 5, 6],
                  'S': [2, 3, 4, 5, 6],
                  'Cl': [0, 1, 2, 3, 4],
                  'Br': [0, 1, 2, 3, 4],
                  'I': [0, 1, 2, 3, 4, 5, 6],
                  'Rh': [0, 1, 2, 3, 4, 5, 6]
                  }

# masses from https://ciaaw.org/atomic-masses.htm
atomic_weights = {'H': 1.01, 'He': 4.0, 'Li': 6.94, 'Be': 9.01, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.0, 'F': 19.0, 'Ne': 20.18, 'Na': 22.99, 'Mg': 24.3,
                  'Al': 26.98, 'Si': 28.08, 'P': 30.97, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.95, 'K': 39.1, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.0,
                  'Mn': 54.94, 'Fe': 55.84, 'Co': 58.93, 'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38, 'Ga': 69.72, 'Ge': 72.63, 'As': 74.92, 'Se': 78.97, 'Br': 79.9, 'Kr': 83.8,
                  'Rb': 85.47, 'Sr': 87.62, 'Y': 88.91, 'Zr': 91.22, 'Nb': 92.91, 'Mo': 95.95, 'Tc': 97.91, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87,
                  'Cd': 112.41, 'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.9, 'Xe': 131.29, 'Cs': 132.91, 'Ba': 137.33, 'La': 138.91, 'Ce': 140.12,
                  'Pr': 140.91, 'Nd': 144.24, 'Pm': 144.91, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.5, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93,
                  'Yb': 173.04, 'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59,
                  'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.04, 'Pa': 231.04,
                  'U': 238.03, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0, 'Cm': 247.0, 'Bk': 247.0, 'Cf': 251.0, 'Es': 252.0, 'Fm': 257.0, 'Md': 258.0, 'No': 259.0, 'Lr': 262.0,
                  'Rf': 267.0, 'Db': 268.0, 'Sg': 271.0, 'Bh': 274.0, 'Hs': 269.0, 'Mt': 276.0, 'Ds': 281.0, 'Rg': 281.0, 'Cn': 285.0, 'Nh': 286.0, 'Fl': 289.0, 'Mc': 288.0,
                  'Lv': 293.0, 'Ts': 294.0, 'Og': 294.0}

# vdw radii from https://books.google.no/books?id=bNDMBQAAQBAJ
vdw_radii = {'H': 1.1, 'He': 1.4, 'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73, 'Al': 1.84,
             'Si': 2.1, 'P': 1.8, 'S': 1.8, 'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ca': 2.31, 'Sc': 2.15, 'Ti': 2.11, 'V': 2.07, 'Cr': 2.06, 'Mn': 2.05, 'Fe': 2.04,
             'Co': 2.0, 'Ni': 1.97, 'Cu': 1.96, 'Zn': 2.01, 'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.9, 'Br': 1.85, 'Kr': 2.02, 'Rb': 3.03, 'Sr': 2.49, 'Y': 2.32,
             'Zr': 2.23, 'Nb': 2.18, 'Mo': 2.17, 'Tc': 2.16, 'Ru': 2.13, 'Rh': 2.1, 'Pd': 2.1, 'Ag': 2.11, 'Cd': 2.18, 'In': 1.93, 'Sn': 2.17, 'Sb': 2.06, 'Te': 2.06,
             'I': 1.98, 'Xe': 2.16, 'Cs': 3.43, 'Ba': 2.68, 'La': 2.43, 'Ce': 2.42, 'Pr': 2.4, 'Nd': 2.39, 'Pm': 2.38, 'Sm': 2.36, 'Eu': 2.35, 'Gd': 2.34, 'Tb': 2.33,
             'Dy': 2.31, 'Ho': 2.3, 'Er': 2.29, 'Tm': 2.27, 'Yb': 2.26, 'Lu': 2.24, 'Hf': 2.23, 'Ta': 2.22, 'W': 2.18, 'Re': 2.16, 'Os': 2.16, 'Ir': 2.13, 'Pt': 2.13,
             'Au': 2.14, 'Hg': 2.23, 'Tl': 1.96, 'Pb': 2.02, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02, 'Rn': 2.2, 'Fr': 3.48, 'Ra': 2.83, 'Ac': 2.47, 'Th': 2.45, 'Pa': 2.43,
             'U': 2.41, 'Np': 2.39, 'Pu': 2.43, 'Am': 2.44, 'Cm': 2.45, 'Bk': 2.44, 'Cf': 2.45, 'Es': 2.45, 'Fm': 2.45, 'Md': 2.46, 'No': 2.46, 'Lr': 2.46}

pi_valencies = {'B': [1, 2], 'N': [1, 2], 'O': [1], 'C': [1, 2, 3], 'P': [1, 2, 3, 4], 'S': [1, 3, 4, 5],
                'Si': [1, 2, 3]}

metals = ['Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga',
          'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce',
          'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
          'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm',
          'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl',
          'Mc', 'Lv']


def get_maximal_valance(atom_label):
    """Get the maximum valance of an atom

    Arguments:
        atom_label (str): atom label e.g. C or Pd

    Returns:
        (int): maximal valence of the atom
    """

    if atom_label in valid_valances.keys():
        return valid_valances[atom_label][-1]
    else:
        return 6


def get_atomic_weight(atom_label):
    """Get the atomic weight of an atom

    Arguments:
        atom_label (str): atom label e.g. C or Pd

    Returns:
        (float): atomic weight of the atom
    """

    if atom_label in atomic_weights.keys():
        return atomic_weights[atom_label]
    else:
        return 70


def get_vdw_radius(atom_label):
    """Get the van der waal's radius of an atom
    Arguments:
        atom_label (str): atom label e.g. C or Pd
    Returns:
        (float): van der waal's radius of the atom
    """
    if atom_label in vdw_radii.keys():
        return vdw_radii[atom_label]
    else:
        return 2.3


def is_pi_atom(atom_label, valency):
    """
    Determine if an atom is a 'π-atom' i.e. is unsaturated and is a first or
    second row element

    Arguments;
        atom_label (str):
        valency (int):

    Returns:
        (bool)
    """

    if atom_label not in pi_valencies.keys():
        return False

    if valency in pi_valencies[atom_label]:
        return True

    return False

class XYZfileWrongFormat(Exception):
    pass

class XYZfileDidNotExist(Exception):
    pass

class Species:
    def __init__(self, atoms):
        self.atoms = atoms
        self.n_atoms = 0 if atoms is None else len(atoms)
        self.graph = None

    def get_coordinates(self):
        """Return a np.ndarray of size n_atoms x 3 containing the xyz
        coordinates of the molecule"""
        return np.array([atom.coord for atom in self.atoms])

def make_graph(species, rel_tolerance=0.25, bond_list=None,
               allow_invalid_valancies=False):
    """
    Make the molecular graph from the 'bonds' determined on a distance criteria
    or a smiles parser object. All attributes default to false

    Nodes attributes:
        (0) atom_label: Atomic symbol of this atom
        (1) stereo: Is this atom part of some stereochemistry e.g. R/S or E/Z

    Edge attributes:
        (1) pi: Is this bond a pi bond. If it is then there should be no
                rotation the bond axis in conformer generation
        (2) active: Is this bond being made/broken
                   (applies only to TransitionState objects)

    Arguments:
        species (autode.species.Species):

    Keyword Arguments:
        rel_tolerance (float):
        bond_list (list(tuple)):
        allow_invalid_valancies (bool):
    """

    graph = nx.Graph()

    # Add the atoms to the graph all are initially assumed not to be
    # stereocenters
    for i in range(species.n_atoms):
        graph.add_node(i, atom_label=species.atoms[i].label, stereo=False)

    # If bonds are specified then add edges to the graph and return
    if bond_list is not None:
        [graph.add_edge(bond[0], bond[1], pi=False, active=False) for bond in bond_list]
        species.graph = graph
        return None

    else:
        # Loop over the unique pairs of atoms and add 'bonds'
        coordinates = species.get_coordinates()
        dist_mat = distance_matrix(coordinates, coordinates)

        for i in get_atom_ids_sorted_type(species):

            # Iterate through the closest atoms to atom i
            for j in np.argsort(dist_mat[i]):

                if i == j:
                    # Don't bond atoms to themselves
                    continue

                avg_bond_length = get_avg_bond_length(atom_i_label=species.atoms[i].label,
                                                      atom_j_label=species.atoms[j].label)

                # If the distance between atoms i and j are less or equal to 1.2x average length add a 'bond'
                if dist_mat[i, j] <= avg_bond_length * (1.0 + rel_tolerance) and (i, j) not in graph.edges:
                    graph.add_edge(i, j, pi=False, active=False)

    species.graph = graph
    set_graph_attributes(species)

    if not allow_invalid_valancies:
        remove_bonds_invalid_valancies(species)

    return None

def get_bond_type_list(graph):
    """Finds the types (i.e CH) of all the bonds in a molecular graph
    Arguments:
        graph (nx.Graph): Molecular graph
    Returns:
        bond_list_dict (dict): key = bond type, value = list of bonds of this type
    """
    bond_list_dict = {}
    atom_types = set()

    for _, atom_label in graph.nodes.data('atom_label'):
        atom_types.add(atom_label)

    ordered_atom_labels = sorted(atom_types)

    for index, atom_label in enumerate(ordered_atom_labels):
        for i in range(index, len(ordered_atom_labels)):
            key = atom_label + ordered_atom_labels[i]
            bond_list_dict[key] = []

    for bond in graph.edges:
        atom_i_label = graph.nodes[bond[0]]['atom_label']
        atom_j_label = graph.nodes[bond[1]]['atom_label']
        key1, key2 = atom_i_label + atom_j_label, atom_j_label + atom_i_label

        if key1 in bond_list_dict.keys():
            bond_list_dict[key1].append(bond)
        elif key2 in bond_list_dict.keys():
            bond_list_dict[key2].append(bond)

    return bond_list_dict

def get_fbonds(graph, key):
    """Get all the possible forming bonds of a certain type
    Arguments:
        graph (nx.Graph): graph object of a molecule
        key (str): string representing the bond type to be examined
    Returns:
        list: list of bonds that can be made of this type
    """
    possible_fbonds = []
    bonds = list(graph.edges)
    for atom_i in graph.nodes:
        for atom_j in graph.nodes:
            if atom_i < atom_j:
                if not (atom_i, atom_j) in bonds and not (atom_j, atom_i) in bonds:
                    bond = (atom_i, atom_j)
                    atom_i_label = graph.nodes[bond[0]]['atom_label']
                    atom_j_label = graph.nodes[bond[1]]['atom_label']
                    key1, key2 = atom_i_label + atom_j_label, atom_j_label + atom_i_label
                    if key1 == key or key2 == key:
                        possible_fbonds.append(bond)

    return possible_fbonds

def get_bond_rearrangs(reactant, product, name):
    """For a reactant and product (complex) find the set of breaking and
    forming bonds that will turn reactants into products. This works by
    determining the types of bonds that have been made/broken (i.e CH) and
    then only considering rearrangements involving those bonds.

    Arguments:
        reactant (autode.complex.ReactantComplex):
        product (autode.complex.ProductComplex):
        name (str): Name of the reaction

    Returns:
        list: list of bond rearrang objects linking reacs and prods
    """

    if os.path.exists(f'{name}_bond_rearrangs.txt'):
        return get_bond_rearrangs_from_file(f'{name}_bond_rearrangs.txt')

    if is_isomorphic(reactant.graph, product.graph) and product.n_atoms > 3:
        return None

    possible_bond_rearrs = []

    reac_bond_dict = get_bond_type_list(reactant.graph)
    prod_bond_dict = get_bond_type_list(product.graph)

    # list of bonds where this type of bond (e.g C-H) has less bonds in
    # products than reactants
    all_possible_bbonds = []

    # list of bonds that can be formed of this bond type. This is only used
    # if there is only one type of bbond, so can be overwritten for each new
    # type of bbond
    bbond_atom_type_fbonds = None

    # list of bonds where this type of bond (e.g C-H) has more bonds in
    #  products than reactants
    all_possible_fbonds = []

    # list of bonds that can be broken of this bond type. This is only used
    # if there is only one type of fbond, so can be overwritten for each new
    # type of fbond
    fbond_atom_type_bbonds = None

    # list of bonds where this type of bond (e.g C-H) has the same number of
    # bonds in products and reactants
    possible_bbond_and_fbonds = []

    for reac_key, reac_bonds in reac_bond_dict.items():
        prod_bonds = prod_bond_dict[reac_key]
        possible_fbonds = get_fbonds(reactant.graph, reac_key)
        if len(prod_bonds) < len(reac_bonds):
            all_possible_bbonds.append(reac_bonds)
            bbond_atom_type_fbonds = possible_fbonds
        elif len(prod_bonds) > len(reac_bonds):
            all_possible_fbonds.append(possible_fbonds)
            fbond_atom_type_bbonds = reac_bonds
        else:
            if len(reac_bonds) != 0:
                possible_bbond_and_fbonds.append([reac_bonds, possible_fbonds])
    print(possible_bbond_and_fbonds)
    print('------')
    """
    # The change in the number of bonds is > 0 as in the reaction
    # initialisation reacs/prods are swapped if this is < 0
    delta_n_bonds = reactant.graph.number_of_edges() - product.graph.number_of_edges()
    if delta_n_bonds == 0:
        funcs = [get_fbonds_bbonds_1b1f, get_fbonds_bbonds_2b2f]
    elif delta_n_bonds == 1:
        funcs = [get_fbonds_bbonds_1b, get_fbonds_bbonds_2b1f]
    elif delta_n_bonds == 2:
        funcs = [get_fbonds_bbonds_2b]
    else:
        return None
    
    for func in funcs:
        possible_bond_rearrs = func(reactant, product, possible_bond_rearrs,
                                    all_possible_bbonds,
                                    all_possible_fbonds,
                                    possible_bbond_and_fbonds,
                                    bbond_atom_type_fbonds,
                                    fbond_atom_type_bbonds)

        if len(possible_bond_rearrs) > 0:
            # This function will return with the first bond rearrangement
            # that leads to products

            n_bond_rearrangs = len(possible_bond_rearrs)
            if n_bond_rearrangs > 1:
                logger.info(f'Multiple *{n_bond_rearrangs}* possible bond '
                            f'breaking/makings are possible')
                possible_bond_rearrs = strip_equivalent_bond_rearrangs(reactant, possible_bond_rearrs)

            save_bond_rearrangs_to_file(possible_bond_rearrs,
                                        filename=f'{name}_bond_rearrangs.txt')

            return possible_bond_rearrs

    return None
    """
def get_atom_ids_sorted_type(species):
    """
    Get a list of atom ids sorted by increasing atomic weight, useful for when
     a molecular graph depends on the order
    of atoms in what will be considered bonded

    Arguments:
        species (autode.species.Species):

    Returns:
        (list(int)):
    """
    return sorted(list(range(species.n_atoms)), key=lambda i: get_atomic_weight(atom_label=species.atoms[i].label))

def get_avg_bond_length(atom_i_label, atom_j_label):
    """Get the average bond length from either a molecule and a bond or two
    atom labels (e.g. atom_i_label = 'C'
    atom_j_label = 'H')
    Keyword Arguments:
        atom_i_label (str): atom label e.g 'C'
        atom_j_label (str): atom label e.g 'C'
    Returns:
        float: avg bond length of the bond in Å
    """
    key1, key2 = atom_i_label + atom_j_label, atom_j_label + atom_i_label

    if key1 in avg_bond_lengths.keys():
        return avg_bond_lengths[key1]
    elif key2 in avg_bond_lengths.keys():
        return avg_bond_lengths[key2]
    else:
        return 1.5

def set_graph_attributes(species):
    """
    For a molecular species set the π bonds and stereocentres in the molecular graph.

    Arguments:
        species (autode.species.Species):
    """

    for bond in species.graph.edges:
        atom_i, atom_j = bond

        if all([is_pi_atom(atom_label=species.atoms[atom].label, valency=species.graph.degree[atom]) for atom in bond]):
            # TODO fix this for alternating single and double bonds, currently all shown as pi
            species.graph.edges[atom_i, atom_j]['pi'] = True

    # List of atom indexes that are rings in the species
    rings = find_cycles(species.graph)

    for (i, j) in species.graph.edges:

        if species.graph.edges[(i, j)]['pi'] is False:
            continue

        if any(i in ring for ring in rings):
            # The ring should define the stereochemistry of this pi bond
            continue

        if is_chiral_pi_bond(species, bond=(i, j)):
            species.graph.nodes[i]['stereo'] = True
            species.graph.nodes[j]['stereo'] = True

    for i in range(species.n_atoms):
        if is_chiral_atom(species, atom_index=i):
            species.graph.nodes[i]['stereo'] = True

    return None

def find_cycles(graph):
    """Finds all the cycles in a graph

    Arguments:
        graph (nx.Graph): the molecular graph

    Returns:
        list(list): each list has the atoms in a cycle
    """
    return nx.cycle_basis(graph)

def is_chiral_atom(species, atom_index):
    """Determine if an atom is chiral, by seeing if any of the bonded groups
    are the same"""
    neighbours = list(species.graph.neighbors(atom_index))

    if len(neighbours) != 4:
        return False

    graphs = []
    for neighbour in neighbours:
        graph = species.graph.copy()
        graph.remove_edge(atom_index, neighbour)
        split_subgraphs = get_separate_subgraphs(graph)
        graphs.append([subgraph for subgraph in split_subgraphs if neighbour in list(subgraph.nodes())][0])

    for graph1, graph2 in itertools.combinations(graphs, 2):
        if is_isomorphic(graph1, graph2, ignore_active_bonds=True):
            return False

    return True

def get_separate_subgraphs(graph):
    """Find all the unconnected graphs in a graph

    Arguments:
        graph (nx.Graph): graph

    Returns:
        list: list of graphs separate graphs
    """
    return [graph.subgraph(c).copy() for c in nx.connected_components(graph)]

def is_isomorphic(graph1, graph2, ignore_active_bonds=False, timeout=5):
    """Check whether two NX graphs are isomorphic. Contains a timeout because
    the gm.is_isomorphic() method occasionally gets stuck

    Arguments:
        graph1 (nx.Graph): graph 1
        graph2 (nx.Graph): graph 2

    Keyword Arguments:
        ignore_active_bonds (bool):
        timeout (float): Timeout in seconds

    Returns:
        (bool): if the graphs are isomorphic
    """

    if ignore_active_bonds:
        graph1, graph2 = get_graphs_ignoring_active_edges(graph1, graph2)

    if not isomorphism.faster_could_be_isomorphic(graph1, graph2):
        return False

    # Always match on atom types
    node_match = isomorphism.categorical_node_match('atom_label', 'C')

    if ignore_active_bonds:
        gm = isomorphism.GraphMatcher(graph1, graph2,
                                      node_match=node_match)

    else:
        # Also match on edges
        edge_match = isomorphism.categorical_edge_match('active', False)
        gm = isomorphism.GraphMatcher(graph1, graph2,
                                      node_match=node_match,
                                      edge_match=edge_match)

    # NX can hang here for not very large graphs, so kill after a timeout

    def handler(signum, frame):
        raise TimeoutError

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(int(timeout))
    try:
        result = gm.is_isomorphic()
        # Cancel the timer
        signal.alarm(0)
        return result

    except TimeoutError:
        return False

def get_graphs_ignoring_active_edges(graph1, graph2):
    """
    Remove any active edges that are in either graph1 or graph2 from both graphs
    Arguments:
        graph1 (nx.Graph):
        graph2 (nx.Graph):

    Returns:
        (tuple(nx.Graph))
    """
    graph1_no_ae, graph2_no_ae = graph1.copy(), graph2.copy()

    # Iterate through the pairs removing any active edges from both ga and gb
    for (ga, gb) in [(graph1_no_ae, graph2_no_ae), (graph2_no_ae, graph1_no_ae)]:

        for (i, j) in [edge for edge in ga.edges if ga.edges[edge]['active'] is True]:
            ga.remove_edge(i, j)

            if (i, j) in gb.edges:
                gb.remove_edge(i, j)

    return graph1_no_ae, graph2_no_ae

def remove_bonds_invalid_valancies(species):
    """
    Remove invalid valencies for atoms that exceed their maximum valencies e.g.
    H should have no more than 1 'bond'

    Arguments:
        species (autode.species.Species):
    """

    for i in species.graph.nodes:

        max_valance = get_maximal_valance(atom_label=species.atoms[i].label)
        neighbours = list(species.graph.neighbors(i))

        if len(neighbours) <= max_valance:
            # All is well
            continue

        # Get the atom indexes sorted by the closest to atom i
        closest_atoms = sorted(neighbours, key=lambda j: species.get_distance(i, j))

        # Delete all the bonds to atom(s) j that are above the maximal valance
        for j in closest_atoms[max_valance:]:
            species.graph.remove_edge(i, j)

    return None

avg_bond_lengths = {'LiNi': 2.7070695291559197, 'BiO': 2.4298898886030154, 'OCe': 2.488400221718809, 'UAs': 3.0832212197773314, 'ClHf': 2.4410942125986943, 'AlAl': 2.6497629851622038, 'RuZr': 2.9322616006579074, 'NMn': 2.1825246259337505, 'CoP': 2.212743178999564, 'ReSn': 2.76496140321245, 'HDy': 2.3369442313247566, 'ReSb': 2.7600883493906703, 'SPb': 2.8404475127884847, 'HgTe': 2.7829486101514727, 'ReSe': 2.5200811821293843, 'SIr': 2.365427951760859, 'SnSe': 2.587833445421715, 'RuSn': 2.665534715621561, 'CuCuB': 2.625194962612407, 'ZnNi': 2.385386968013931, 'SnSn': 2.9557148767199877, 'CW': 2.1537652648432757, 'PbBr': 3.0364502625041947, 'GePt': 2.4302410821648883, 'RuSe': 2.5010691628423873, 'HAs': 1.3821823502708392, 'HfHf': 3.238548960971258, 'ZrZr': 3.5033533268020207, 'NaGa': 3.235677503004142, 'PdRu': 2.7964154581439007, 'CPt': 2.0458786557983455, 'OYb': 2.3210598544872436, 'HfCo': 2.598229427398477, 'CdDO': 2.5357603660043044, 'CuNi': 2.5876013957174835, 'ClLi': 2.434495721224297, 'SNa': 2.900766971378299, 'MoW': 3.0219761431733434, 'ErSb': 3.2408414088973103, 'ReAg': 2.86078494248277, 'ClAl': 2.1581370198560235, 'NPb': 2.598903956152097, 'SiC': 1.8706641953661955, 'OSm': 2.4532509021736635, 'HfP': 2.651584617810598, 'IrGe': 2.4953169411715947, 'NbCl': 2.4502796936460514, 'ClLu': 2.6147551904085646, 'NPt': 2.0571937272126757, 'FFe': 1.9725602652283045, 'TaO': 2.015913792977548, 'TiCl': 2.3409141913191305, 'TaI': 2.956659601603536, 'ReAs': 2.53723737234476, 'MoK': 3.6136518878313493, 'WRe': 3.0143077198252866, 'AgAg': 2.972633675366598, 'CoCr': 2.8296368299927708, 'HU': 2.397434733507184, 'OCu': 2.0751221088156697, 'PtRu': 2.86385151676487, 'ONi': 2.060416531647872, 'SeLa': 3.2335324850815, 'CdZn': 3.4635061112499175, 'AgAs': 2.4723414885141795, 'YbP': 2.875395388837503, 'MnMo': 2.856953732689341, 'LuN': 2.4151042406446295, 'PtRh': 2.737304349631165, 'FZr': 2.0352206872306686, 'CoCo': 2.53733109826256, 'LiLi': 2.6142783550548208, 'CoCl': 2.343959337048827, 'BYb': 2.7714463338931106, "O'C'": 1.367004011868391, 'VAu': 2.7232344217987174, 'IrCl': 2.419047449309889, 'SbS': 2.5499231287434205, 'PuN': 2.5724468731493872, 'PuO': 2.2618262650268433, 'AsAs': 2.4058277867142013, 'TeGe': 2.5945191925146154, 'BaCl': 3.1862781201587045, 'ZrCl': 2.4797172982935995, 'ZrCo': 2.670307934004338, 'PuC': 2.6443646997361556, 'SbC': 2.132113482191541, 'BaCu': 3.739825960000002, 'SbF': 1.8511473126286953, 'LiPt': 2.6299896548221047, 'KBr': 3.441093450468664, 'SbK': 3.870437612525152, 'TaS': 2.4816509559298807, 'SbI': 3.0472108766656607, 'SbN': 2.2434570766140953, 'SbO': 2.0691517648663407, 'CeCO': 2.4001812539842544, 'CMo': 2.2145079857893064, 'CMn': 1.9837872493003639, 'GeFe': 2.3984152697585848, 'CMg': 2.26881101386439, 'ReS': 2.3858463285215215, 'AsHf': 2.8814942152910685, 'NaAO': 2.688639106175858, 'SmS': 2.8778599065454715, 'MoF': 2.1633530212354493, 'AuOs': 2.835653698986364, 'SmN': 2.5557739832992024, 'ReI': 2.796770328862999, 'ReN': 2.121275303287943, 'PCl': 2.003020159979613, 'BEr': 2.744166093706119, 'SmF': 2.2996860345765193, 'RhSe': 2.462592722831039, 'SmC': 2.7882308599200307, 'RhSb': 2.683056046050556, 'CrCl': 2.364697606106688, 'SbMo': 2.8395955061906952, 'TlPd': 2.905245558526089, 'HI': 2.218978184776487, 'PdCr': 2.7627734125361263, 'HoSb': 3.2449838096481174, 'NiCl': 2.3523327661768025, 'CBe': 1.8820985438653695, 'HLu': 2.22180857593825, 'CBi': 2.238168922584931, 'CdRu': 2.7750087384283115, 'HgCo': 2.6226667288972534, 'ClSb': 2.48183638149146, 'PbPb': 3.212310373964439, 'INi': 2.644211380214665, 'AgSb': 2.7360373456422886, 'AgSe': 2.67830139116289, 'CoGa': 2.4057684904854333, 'ClSn': 2.453330636168858, 'NOs': 2.0634811390709555, 'AgO': 2.4470096563328387, 'AgN': 2.264469252451163, 'NFe': 2.0576788465239626, 'HgMo': 2.849141743931727, 'AgI': 2.880027859736813, 'HCu': 1.80505421099909, 'BLu': 2.8731649340045298, 'AgF': 2.6110684522307284, 'AgC': 2.320949348027622, 'RhMo': 2.85967019685819, 'CrO': 1.9522848316088155, 'FTa': 1.8867118115801058, 'OsSn': 2.694857547698872, 'SnMn': 2.635944388984265, 'HCd': 2.3583388643139798, 'HCe': 2.396134565068835, 'NbBr': 2.688552474282795, 'SSi': 2.1295079288620147, 'PCr': 2.4078451699059693, 'BLa': 3.130318432197145, 'FTl': 2.973751369738073, 'HCl': 1.2686612421405876, 'AgP': 2.4364427875613246, 'FP': 1.5725881208061836, 'FS': 1.5847562619640267, 'FCu': 2.2431351549399645, 'FW': 1.8730211905778162, 'FeRh': 2.739129307945338, 'CuCl': 2.378832803343361, 'SiSn': 2.627672464352516, 'FCa': 2.2770791929975207, 'FB': 1.364609269178264, 'FC': 1.3316802921574742, 'CuCu': 2.652880598891229, 'BSi': 2.0242076598656356, 'FF': 1.0347156484970932, 'BTl': 2.396510365596341, 'FH': 1.0903010775564081, 'FCl': 1.8532837692019963, 'IrBr': 2.5224187182729723, 'BSb': 2.3014573604163604, 'FO': 0.8979349326194215, 'TlBr': 2.5530629671173117, 'AsS': 2.2918446692969283, 'AuB': 2.2924803915857854, 'AsSi': 2.350065081892182, 'HfI': 2.898584433680038, 'OsAs': 2.4820208378788355, 'AuF': 1.3281604360040509, 'AuI': 2.649211887597211, 'PtNi': 2.609069510498008, 'RbH': 2.9227548513033317, 'CrGa': 2.4291796102413055, 'AuO': 2.0463584824661876, 'SMo': 2.4035314394056835, 'VBr': 2.4347393369443098, 'AuS': 2.3289994941915406, 'BrTh': 2.8396056464290056, 'CTa': 2.3884669227791617, 'CTb': 2.683026004627555, 'ZnRe': 2.7109328347449093, 'CTm': 2.61005709818862, 'PdNi': 2.944097151387155, 'BrTa': 2.5920261406399754, 'BrSi': 2.262645045023553, 'RuRu': 2.8278265792870823, 'HgCu': 2.8637118016934187, 'BrSn': 2.5930231206051335, 'SeAs': 2.450691322603182, 'SiW': 2.527771135949235, 'BrSb': 2.8289569507108543, 'BrSe': 2.633502127126532, 'InS': 2.4998384271276177, 'NTi': 2.0616330263531255, 'NTh': 2.4953421703412246, 'InO': 2.220924487962098, 'NTl': 2.6390296930679957, 'ORh': 2.0860036442612793, 'SeAl': 2.443103861150297, 'NTa': 2.0604958458936653, 'NTc': 2.1331457848150808, 'ORe': 1.8767156811358139, 'NTe': 2.2027217446732887, 'AgFe': 2.671497480493031, 'MoCu': 2.707050820515104, 'YY': 3.5806514664404228, 'SY': 2.754699784944037, 'RhCu': 2.633560812549602, 'BeN': 1.7367525583582504, 'BeO': 1.6197606533183304, 'CoACl': 2.462164820110209, 'SS': 2.0024507383729566, 'ClAg': 2.666974415984924, 'YbS': 2.773166300259441, 'SU': 2.846280344069924, 'OAlA': 1.883416542924304, 'SH': 1.1807795612025116, 'SO': 1.4521232933366812, 'NiCo': 2.55475208745558, 'YbN': 2.41908778014619, 'AlSi': 2.456978304786579, 'SB': 1.8974426823584578, 'MoCo': 2.732725039315575, 'YbC': 2.634894457938566, 'AlSb': 2.7536885566843163, 'HgMn': 2.6022174064756785, 'YFe': 3.443652001424321, 'NHg': 2.4079411984296253, 'NHf': 2.2108752002269103, 'YN': 2.4237832829417374, 'RhP': 2.2941801296650364, 'RhS': 2.3642651147949656, 'YB': 2.6179614080661007, 'LiBi': 2.92362838364096, 'RhW': 2.839354358753573, 'YF': 2.309328978968292, 'HgFe': 2.5811726710161405, 'GeCl': 2.223475300321501, 'GeCo': 2.3772474773607737, 'HgAs': 2.6037019294254784, 'CSc': 2.4037886425209085, 'HTi': 1.9389801096462909, 'ONp': 2.3070950862472355, 'RhB': 2.217250213640156, 'ODy': 2.3867500584189862, 'NSr': 2.693928634037443, 'HIn': 1.6800186373968955, 'RhF': 2.174973993362162, 'ZnS': 2.365484999519886, 'LiF': 1.964131494980243, 'ZrS': 2.6374205499668273, 'BAs': 2.0810635000367412, 'LiB': 2.3919452118899542, 'ZrP': 2.724243416627099, 'IBa': 3.4404989522808136, 'LiN': 2.0776507286792265, 'LiO': 1.9878128104848398, 'LiH': 2.0164875577680186, 'HgS': 2.555412991705643, 'CAu': 2.040513968886866, 'SiCl': 2.077049767436305, 'AuSe': 2.449947154672638, 'CPb': 2.3735934124758886, 'LiP': 2.610876232456476, 'ZrB': 2.595887002045729, 'HgI': 2.774745322705845, 'HgF': 2.5630806418705085, 'AuSn': 2.7354079201752692, 'IBr': 2.787010496324369, 'HgC': 2.1012303408331285, 'ZrI': 2.903537156856002, 'ZrH': 1.9547161110956406, 'AlN': 1.9574655915756556, 'AlO': 1.8298592897984274, 'PtAg': 2.8606880130135317, 'BHo': 2.666511997349767, 'AlI': 2.533411672794378, 'AlF': 1.8124488000702828, 'NCr': 2.049695839422264, 'CFe': 2.0062412966181116, 'LaCl': 2.8466572126432634, 'TlTe': 2.9549369369658938, 'ScCl': 2.502595673103318, 'PtAs': 2.3595753635138252, 'RhNi': 2.628674910751867, 'SnIr': 2.6151250454659407, 'NCl': 1.737938963856322, 'ThCl': 2.791227723983537, 'CGa': 2.058738937772583, 'NCo': 2.044057968181483, 'PNa': 2.899943295473412, 'OGd': 2.4196802769229513, 'TiP': 2.5810028432777967, 'NCd': 2.3356338878872838, 'NCe': 2.59914597572961, 'AlP': 2.4432491451428655, 'AgACl': 3.024914959543951, 'IrIr': 2.7484351960124775, 'NZr': 2.2208753195555735, 'UF': 2.3014472280740277, 'WCu': 2.684999424093159, 'RhAs': 2.4686645351792134, 'STc': 2.377630278250582, 'STm': 2.8701089163718843, 'PtCl': 2.331436212064129, 'LiC': 2.340462307932625, 'STh': 2.9760108558104297, 'II': 2.8594237302823537, 'IrP': 2.316959015976385, 'BiMo': 2.9304071841576538, 'SeSe': 2.399349541751457, 'OsBr': 2.5180146620043886, 'SeSi': 2.2524180144152046, 'RhAl': 2.4269662662743268, 'HGd': 2.188752781429149, 'NaSi': 2.991378752378181, 'SeSm': 2.8947338725798475, 'NZn': 2.0901726850158466, 'WCo': 2.7893783834639616, 'RbCl': 3.430980518407841, 'PRe': 2.4415480348714906, 'TiBr': 2.5459459037733407, 'MgF': 1.9868403262248449, 'CuOs': 2.6810148978025192, 'SAl': 2.275105678696092, 'BrC': 1.9005783459901855, 'TiTi': 2.951455327905903, 'TbS': 2.8564230255143186, 'MgH': 2.0141991267066204, 'MgO': 2.072057970605808, 'MgN': 2.115454422884608, 'DyS': 2.857605784340484, 'MgS': 2.6014177702382475, 'TbH': 2.224429928709167, 'DyN': 2.525022789021828, 'MgP': 2.631409117233076, 'CsAg': 3.856563111086396, 'TbO': 2.4038480843091867, 'TbN': 2.5416554472695987, 'WTi': 3.04752739345006, 'LiAl': 2.8208560131723606, 'WSn': 2.76519323299985, 'CNa': 2.830512159373349, 'DyB': 2.72541428182172, 'DyC': 2.6658013318072933, 'SeMo': 2.509359556597272, 'CdF': 2.4696858966160873, 'CdCl': 2.59907561535338, 'BAg': 2.6879310292462484, 'TmP': 2.8487081641672227, 'CdCd': 1.3067759383988198, 'SiCo': 2.205473834234736, 'CoB': 2.101426610727582, 'WPt': 2.7595589121225026, 'FNd': 2.3219093264334294, 'TiS': 2.4397971423064284, 'HgPt': 2.8780664339623074, "C'C'": 1.4471672875078008, 'NCoA': 2.0109540921567306, 'PtSe': 2.4575007725449924, 'CCa': 2.6992626912076596, 'FNi': 2.0035477036535747, 'RuBr': 2.5499438303619484, 'PdAs': 2.380916865059004, 'KAg': 3.9100151393218, 'IPd': 2.634522066982884, 'SbCu': 2.5917538588636644, 'AuAg': 2.911627690350485, 'LiS': 2.5114259300401827, 'SnC': 2.1690738376969567, 'NiGa': 2.1700262347821075, 'DW': 1.8306514979173771, 'HoO': 2.3718022435900643, 'LuRu': 3.0233723227976625, 'PdAg': 2.8772554295273864, 'GeMo': 3.0126856882224615, 'IPt': 2.671523027643532, 'GaTe': 2.5262597386426333, 'AuAu': 2.906411448306194, 'HAl': 1.6178942581635314, "N'C'": 1.3578461583327843, 'CsAs': 3.8225261266391692, 'AuAs': 2.4676892887086734, 'NiOs': 2.5828114539636404, 'CHf': 2.495127831184204, 'PrCl': 2.844305925332444, 'CaS': 2.946036329488045, 'BrP': 2.216751809185485, 'CaP': 2.954280371026261, 'FeGa': 2.4360820341618554, 'TeBr': 2.7529448680238335, 'CHo': 2.573397856668976, 'BaF': 2.7618407600411166, 'LiSi': 2.6951808993158934, 'TcCl': 2.3993200091327345, 'CaH': 2.1993966090666635, 'CuSi': 2.408825892630392, 'CaN': 2.4899762347653125, 'LiSn': 2.7848469930735207, 'ErC': 2.6005051171913824, 'HRe': 1.7798468259987001, 'BaO': 2.8174660779161216, 'ErN': 2.470254070657434, 'TiC': 2.3721360262355855, 'PdAN': 1.9940048915116648, 'IY': 3.13881546101777, 'CuRu': 2.665542233203348, 'SiGe': 2.3708078291153893, 'BiPt': 2.747136998483202, 'CaSi': 3.0641540443217012, 'SiRu': 2.386489312246205, 'BiS': 2.8310205702231097, 'BMo': 2.3897452364545275, 'PdBr': 2.4763117176785276, 'MoMo': 2.6646892011395273, 'NLa': 2.6942173649395276, 'HMn': 1.8080228806185128, 'NdCl': 2.8043807643073473, 'CuRb': 3.4714636373158907, 'EuSe': 3.0682388838514045, 'TlTl': 3.1684778314103976, 'EuSi': 3.1936532233274497, 'ClRh': 2.398267216126521, 'NBa': 2.9275464298084715, 'GaS': 2.324136200479328, 'ReRe': 2.6882177251610986, 'CsO': 3.2449239444942695, 'CsN': 3.2874161056524116, 'BiNi': 2.7816421652225407, 'NGa': 1.9975437978384538, 'LuW': 2.9297607436364252, 'OCa': 2.4274403193341434, 'YSb': 3.26010526420154, 'CsF': 3.313462072840273, 'IrTe': 2.638104920249451, 'YSe': 2.792118072445514, 'WFe': 2.787789267328108, 'ReRu': 2.921122407016875, 'IFe': 2.632912663042025, 'GeNi': 2.282450784691833, 'BrY': 2.768940561279656, 'LuI': 3.2301510290383373, 'NbCo': 2.992327656096426, 'KSi': 3.413761600262787, 'CsP': 3.768471830308635, 'RhBr': 2.553446208644088, 'GaCl': 2.1901062644348257, 'AuCl': 2.29980815775298, 'RuCl': 2.4129185505712765, 'CZr': 2.5130305803644535, 'AlB': 2.2081624575331484, 'TeC': 2.1278388549891907, 'CdI': 2.823089518850739, 'SSn': 2.5004527449934435, 'BTa': 2.381557480297976, 'PdTi': 2.8155345721513423, 'GaNd': 3.2198533022938136, 'SiNi': 2.2263166777144954, 'PdTe': 2.6071841209832414, 'BHf': 2.528059445144819, 'BTi': 2.4079862103548413, 'BNa': 2.802397839666889, 'ZnC': 2.081255410634223, 'BrMg': 2.5410832081277936, 'PtRe': 2.952482333710576, 'HfGe': 2.670431282452255, 'ZnF': 2.050169101508828, 'HfF': 2.1353818890690883, 'PdW': 2.848471023959414, 'PtTe': 2.6068058134682213, 'ZnI': 2.558819406635391, 'ZnH': 1.6103290804861352, 'PdSe': 2.440641851159273, 'IK': 3.5757095416214435, 'WOs': 2.9011041318094795, 'MoRu': 2.8509602295119802, 'TlC': 2.7435123416812965, 'PdN': 2.0460254932871367, 'PdK': 3.237563346276828, 'TaAs': 2.7653335062090507, 'NiIr': 2.7795657821581736, 'PbCl': 2.9050078511924182, 'CdP': 2.5473892447037656, 'TlI': 2.7589043808197506, 'BiGa': 2.7062416806107263, 'TlO': 2.7173246794809973, "H'C": 1.0570248273641463, 'MnCr': 2.847297116775744, 'SHo': 2.8868743303718842, 'SmAs': 3.068654859071786, 'NC': 1.3819259523536414, 'CeC': 2.713582402523283, 'SEr': 2.8489252548645267, "C'Cl'": 1.734234744007109, 'CeEO': 2.400181253984254, 'BOs': 2.1897849338011586, 'SiRh': 2.319829479360674, 'CdS': 2.618152657850771, 'SbP': 2.601624904516032, 'CrTe': 2.661371628013206, 'PdFe': 2.667709496900076, 'CoBCl': 2.3669885935382675, 'YbGa': 3.1295315844110814, 'TeO': 1.999760532296311, 'CeS': 2.963281586983438, 'ClB': 1.8145465263264116, 'AuBr': 2.423997989322173, 'NiBr': 2.426788653134236, 'GdCl': 2.7498031750747094, 'GaGa': 2.6106059680253964, 'ZnO': 2.0383064710160714, 'IC': 2.104470795537292, 'NbTe': 2.959032890513408, 'SeIr': 2.5122213595386507, 'EuO': 2.4332539421167643, 'ClO': 1.4082742022499197, 'SbGa': 2.719111160487302, 'PtZn': 2.7055888948083835, 'EuC': 2.852224536054547, 'ClC': 1.737175476844588, 'AgMo': 2.919017003649842, 'ISi': 2.4703467016996177, 'NaO': 2.4342540084597286, 'ISn': 3.0789096416947115, 'ISm': 3.2448235218271373, 'AgMn': 2.882929339014945, 'ISc': 2.8529567261074047, 'TeMo': 2.735911542912268, 'NSi': 1.7507875024918487, 'ISe': 2.9693749129827287, 'FeS': 2.2797436787819754, 'CoOs': 2.6221512630952883, 'SeS': 2.290440559765073, 'TeFe': 2.5698476813193, 'SeP': 2.170645761441978, 'NbS': 2.48775026484812, 'ITe': 3.008228271198096, 'OMn': 2.098639533111047, 'SrCl': 2.930053179018782, 'SeF': 1.812164500734087, 'SeC': 1.8986939056539431, 'SeB': 2.0851571977471326, 'RbRb': 2.586215690236642, 'AsBe': 2.1761615289833, 'SeO': 1.6771480180695022, 'SeN': 1.8242048027686264, 'DyDO': 2.1836621696954794, 'SeK': 3.3158805773501667, 'YSi': 3.1857297760988477, 'SeH': 1.2658911486738011, 'GeGa': 2.520195264517661, 'AsIn': 2.6555805550174196, 'SrSi': 3.208931782697763, 'GeGe': 2.5827290331147283, 'KF': 2.88512708679889, 'AuLi': 2.8700498834742656, 'KK': 4.1573977560358975, 'HFe': 1.588749489503702, 'OCoA': 2.140079185875381, 'KH': 2.7575056818810975, 'RbS': 3.430472892045786, 'KN': 2.9234663410588895, 'IMo': 2.8297101564111555, 'RbP': 3.6014081869250303, 'RbO': 3.006683449285733, 'RbN': 3.0894588740035345, 'KP': 3.4481581477550147, 'EuS': 2.945767643884223, 'FRu': 2.163634760030494, 'HPr': 2.3080143016652204, 'MgSi': 2.647818145945258, 'BrBe': 2.1114972752617103, 'InCl': 2.504250227450158, 'PtOs': 2.8221437761822044, 'BIn': 2.34952564243821, 'TeH': 1.5504904583633161, 'NiN': 2.042167598137738, 'AsV': 2.5358639025485554, 'FeOs': 2.695603094299921, 'DD': 0.49821199999999966, 'DB': 1.1982963251471692, 'POs': 2.363230068974513, 'NpCl': 2.685569246207976, 'AsC': 1.936537579414117, 'LaBr': 2.944427617508985, 'NiY': 2.6487728149383574, 'AsF': 1.6912386745633576, 'NiS': 2.2354522759098985, 'NiP': 2.2070973698690133, 'UCl': 2.682394260309698, 'BIr': 2.2349683862315257, 'TlS': 2.8927251220982835, 'SmSi': 3.197916767342875, 'GaSi': 2.4276068350164524, 'SbPd': 2.6269130658366913, 'WBr': 2.5802416663391194, 'BPd': 2.12064221651576, 'ITi': 2.667927049445489, 'TeTh': 3.226352165270648, 'SbPt': 2.5720244277242124, "Br'C'": 1.9386273716669462, 'BPt': 2.235278192286815, 'PaF': 2.1945202429228385, 'ClMo': 2.4389006037115384, "H'N'": 0.921667680049344, 'NbSe': 2.5761470974671177, 'BRu': 2.245226515278222, 'AlD': 1.5969332856475835, 'CrCr': 2.5407153018092594, 'BiCl': 2.7295624506819527, 'NbP': 2.615262944180665, 'FeSe': 2.4060389117677095, 'OTm': 2.349269366988559, 'ClNa': 2.9245008400562296, 'OsOs': 2.8610563315186353, 'GaI': 2.5734319225661975, 'CoF': 2.0417769875323395, 'IrAg': 2.7162870403422916, 'WK': 3.121295714340231, 'AgHg': 3.0445289853106567, "OC'": 1.2942134085660113, 'ReOs': 2.640639485589761, 'IBe': 2.2439051616465227, 'IrAu': 2.7993504231587965, 'CdFe': 2.6091086006446615, 'HGe': 1.5416055090155396, 'PtPt': 2.7990195721454416, 'NiMn': 2.668525291818046, 'NiMo': 2.723261395159539, 'GaH': 1.3969224989437439, 'CRu': 2.1123675977415015, 'VCl': 2.347459367093525, 'WCl': 2.4204780658826395, 'TlCl': 2.607149175308074, 'PtPb': 2.9586130131228794, 'OsSb': 2.6788542583812816, 'ZnSi': 2.6063782480017794, 'MnGe': 2.3213053067536875, 'ThI': 3.0789700812400285, 'TlP': 3.005140032773963, 'ZnSe': 2.4377732346766274, 'ThF': 2.3351592953900284, 'PNd': 2.8799139584418607, 'SeTh': 2.888423211572957, 'ThC': 2.795215946922777, 'SLa': 3.0300625651485333, 'BN': 1.5334570038254807, 'NbN': 2.1173909116949163, 'HSc': 2.0214095762626836, 'PY': 2.9179286961099202, 'BiI': 3.1035879118616485, 'BiF': 2.345757210378398, 'WN': 2.0924336809571473, 'MoP': 2.4813925022507046, 'NiFe': 2.510998948578899, 'AuFe': 2.7110899459731184, 'WI': 2.847627323149393, 'HSm': 2.243249009147614, 'HSr': 2.580044966637364, 'PZn': 2.378386672567416, 'MoO': 1.9803544844571817, 'MoN': 2.1318890000236403, 'PdBi': 2.7769322811544734, 'HoCl': 2.655878115625434, 'DCu': 1.9061705765389758, 'InF': 2.0740344160873945, 'InMo': 2.750332596668137, 'BrHg': 2.6442248255920306, "O'H'": 0.8929053336369367, 'NS': 1.620409424588893, 'PAs': 2.317872012166414, 'BGd': 2.758589129541411, 'BGe': 2.1361973260204583, 'ZrAs': 2.8804683815791363, 'TeAs': 3.033404807694192, 'InFe': 2.623687716362218, 'TcBr': 2.681827561691311, 'MoAMoA': 2.7171365655800983, 'ClK': 3.209937233370751, 'ErCl': 2.610302747592545, 'CIn': 2.1773407277339505, 'TeAl': 2.5777912022728695, 'CrAl': 2.3760885734272836, 'HgHg': 3.189909956809123, 'TeAg': 2.775359339541472, 'TbBr': 2.7724647871151604, 'ZrAO': 2.2156217561680336, 'BrMn': 2.602347973518951, 'EuTi': 3.3023514517744452, 'NIn': 2.270872921035983, 'HgAu': 2.9194651112929266, 'WW': 2.7314951191293035, 'RhTi': 2.3378948311871093, 'FeFe': 2.6130506106325053, 'DyP': 2.980810650611608, 'BGa': 2.178336234918593, 'ReBi': 2.8316051531600817, 'SrGe': 3.1724670342310546, 'NNa': 2.5427956416139703, 'MgMg': 2.8465888639647106, 'NIr': 2.089648537884844, 'CaCl': 2.744593466018559, 'GaP': 2.4003931040962443, 'MnCl': 2.5028004410909093, 'RbSi': 3.6403468627715703, 'TeSn': 2.7375105702495848, 'SnCr': 2.6130168561697578, 'TiSe': 2.4941375980137894, 'AuPt': 2.8054244602224316, 'CdAO': 2.458171212075381, 'TeSe': 2.778961226509462, 'SnGe': 2.8090945830638665, 'MoIr': 2.8654958338967735, 'BrPr': 2.9210754240877566, 'TiSn': 2.84620656539249, 'CeCl': 2.736850483070955, 'AuPb': 3.5759139304534413, 'BNi': 2.119779988430656, 'AsMg': 2.7572427102455865, 'RuCo': 2.642915526763196, 'AuPd': 2.850523575920352, 'LaIn': 3.404890713914852, 'SnCo': 2.5448449658655297, 'NiTi': 2.6035377113307723, 'PCe': 3.0301027339823587, 'AsMo': 2.588155754395591, 'AsMn': 2.3881836918863546, 'IAs': 2.778472827770457, 'WTe': 2.7602377343203184, 'YbYb': 1.5744591720000018, 'RuH': 1.7164385412642682, 'RuI': 2.74035407192368, 'RuN': 2.073528839431745, 'AgRh': 2.98094229046807, 'CuAS': 2.3557632321120496, 'KO': 2.827877328159513, 'SiOs': 2.379325696902409, 'MoSn': 2.7510303915503167, 'MoSi': 2.5361748637476658, 'DSb': 1.5165224424816284, 'RbSe': 3.501440844857324, 'HN': 0.8921467459940361, 'SiFe': 2.32098397757047, 'CuAO': 2.0429110940227395, 'CuAN': 2.0242934319187778, 'BNb': 2.5356319474389153, 'RuP': 2.330875390331582, 'RuW': 2.8748765221693167, 'AsN': 1.9281510570996176, 'PS': 1.9940386329872093, 'PP': 2.1560540558895984, 'OsCl': 2.3920823689089676, 'PV': 2.491768193357048, 'PW': 2.5023297422467143, 'NFeA': 2.296223600445284, 'PU': 3.035422429324163, 'CuAs': 2.568678103360258, 'ClANb': 2.467968118670316, 'CuAu': 2.912666802837814, 'FePt': 2.719528631381874, 'ZnPd': 2.5018953315574097, 'PC': 1.8262919562735758, 'AsAl': 2.452763582423827, 'CRe': 2.0183928670427447, 'BrZn': 2.384847867594012, 'TmI': 3.079573108894576, 'CRh': 2.111344200104029, 'CuACl': 2.5502244000000007, 'WS': 2.338820957690799, 'OPt': 2.0473685896609504, 'PSi': 2.2545511059775016, 'TlGe': 3.0837356162053426, 'ReCo': 2.648760502443849, 'OPr': 2.5054502614356067, 'CrSe': 2.487799682066019, 'YbFe': 3.0729412069529425, 'BrBr': 2.2275316441453024, 'CdEO': 2.6699790316182757, 'CrSi': 2.3634311151541, 'PuCl': 2.627863559831191, 'NiAS': 2.1625014472403636, 'OPd': 2.0547925075521776, 'TcO': 1.8331664730846837, 'NiAO': 2.0446899341852087, 'ZrSe': 2.6244516597525447, 'OPa': 2.3559660325647394, 'CLu': 2.5912611509039243, 'FeBi': 2.6376236892078335, 'SeGd': 2.973658338452849, 'ReCu': 2.7120079714735237, 'NdI': 3.1727432630223036, 'SrC': 2.8878734665229393, 'NdO': 2.4824773649433207, 'HgSi': 2.4503866648810693, 'NiAu': 3.287853255360337, 'RhZn': 2.6222190944071055, 'HgSe': 2.6824404126901933, 'ClCl': 1.4284798744447706, 'SrO': 2.6109437195371035, 'CoAs': 2.3668713513907864, 'NiAs': 2.3617173068794566, 'SrS': 3.07429131577701, 'SrP': 3.1379275590619797, "N'C": 1.473930709315362, 'RhPd': 2.757960974510443, 'RhZr': 2.9020701705998473, 'NdS': 2.9483572812083674, 'AgCu': 2.804525301668806, 'GeS': 2.240594576917445, 'DN': 1.0002402057332835, 'GeP': 2.368451947059199, 'GeW': 2.4077642570509945, 'CAl': 2.0023846502380325, 'DO': 0.9428463313734669, 'SiZr': 2.9240165891092538, 'NbZn': 2.5407375000000005, 'NNiA': 1.9546918936736817, 'RhRu': 2.8900603863282166, 'OLa': 2.561284694720194, 'CuTi': 2.951505169945368, 'AsCr': 2.465175691480388, 'AsCl': 2.3658070361991217, 'OLu': 2.308598019855343, 'SbH': 1.633286664961096, 'CoLi': 2.778393705922493, 'GeK': 3.607906864683807, 'AsCd': 2.6351978853034006, 'SmSm': 3.4728911906176445, 'GeO': 1.7813985755094819, 'GeN': 1.9521721328566612, 'LaI': 3.2019764366969374, 'InBr': 2.6018185986710716, 'COs': 1.9589066303707947, 'HNb': 1.7538817112351432, "H'C'": 1.0285237491926498, 'HNd': 2.3435860595375586, 'AuMo': 2.8225346199505004, 'AuMn': 2.7393035511546038, 'HNi': 1.6315833486379205, 'FSi': 1.6624071924331534, 'BrCr': 2.598616744009055, 'BrCu': 2.502668810801956, 'PEr': 2.761842077463879, 'CsCl': 3.485355984586813, 'CsCs': 3.3573255344256006, 'FSr': 2.9270995947050924, 'HfH': 2.057125671064434, 'HLa': 2.366113857201757, 'BrCo': 2.4283595756360206, 'HPt': 1.738276535437267, 'BCu': 2.1739127360593815, 'OsH': 1.7725802903524328, "O'B": 1.4355308689027027, 'DyCl': 2.6923552201779786, 'LaBi': 3.4337072186204862, 'BCr': 2.228695212048293, 'MnSi': 2.325682369947051, 'CK': 3.1893036316567107, 'PLu': 2.7848464941449294, 'FeBr': 2.43415518163969, 'PtGa': 2.369164447064784, 'MnPb': 3.079965781764246, 'CC': 1.4278606086920171, 'GaSe': 2.4195104260480007, 'IrRu': 2.7887385117289436, 'IGe': 2.5336086260928306, 'CD': 1.0061070992836214, 'NiNi': 2.660995741603897, 'CY': 2.6418389296825118, 'ReF': 1.9937956644526862, 'CuSe': 2.456321400288201, 'NbMn': 3.1756004039628927, 'PLa': 2.9504649170840445, 'CS': 1.7545906237540243, 'IrRe': 2.9156905024662936, 'LaC': 2.815000526005823, 'CaGe': 3.044632839047009, 'PtTl': 3.098262863559045, 'CuSn': 2.799118113604529, 'RhCo': 2.552133672053358, 'RbAs': 3.6493409399203856, 'HEr': 2.2278645596182747, 'WAg': 2.9708904509009937, 'DRe': 2.0129770191027263, 'AlBr': 2.318483707632079, 'NdN': 2.5930873250541016, 'IEu': 3.2443953208893244, 'ZnBO': 1.9260739739061694, 'BeCl': 2.0088845447988386, 'CrC': 2.070622591407702, 'CrF': 1.9189603989162227, 'SRu': 2.347673896232662, 'MnI': 2.790141081016998, 'PrC': 2.786590055601066, 'WAs': 2.5876886820544636, 'PAu': 2.2896676288755278, 'FOs': 1.9783024604183, 'BiBi': 3.046725638746145, 'CuMn': 2.7321457049412428, 'SrI': 3.294067529055998, 'ClMg': 2.4566027099533208, 'PPd': 2.2917954731601653, 'ErTe': 3.015731497051462, 'MoZn': 2.972685772578208, 'CRb': 3.3533850016717435, 'BaC': 3.1916324231753643, 'SeDy': 2.979113068962466, 'SnPt': 2.646031216456008, 'SOs': 2.408917162499986, 'TaCl': 2.418382774515689, 'SeNd': 3.1123216308872315, 'SmTe': 3.285868710839062, 'CoBN': 1.9529043358574516, 'BiBr': 2.893977780657314, 'LiCu': 2.490242362745092, 'HfBr': 2.504250763256304, 'SnPd': 2.8140763876320887, 'HYb': 2.128456548261408, 'WSb': 2.8135318559511977, 'SeCe': 3.2103689251833987, 'BrGd': 2.8109429103021575, 'IYb': 3.0691898249106093, 'GeSe': 2.3603898357671875, 'AsO': 1.7332007410698766, 'NaH': 2.4863615476478373, 'NaI': 3.355087963350627, 'GaBr': 2.450684901609612, 'SeCl': 2.456812224395904, 'SeCo': 2.3585747821468575, 'GeSb': 2.654577876206654, 'ScSc': 1.693722394303888, 'IrI': 2.7227700939519277, 'LaSb': 3.347386479008939, 'TiO': 1.9540161448034254, 'OsGe': 2.494591195538242, 'IrO': 2.1241678570757676, 'FeMo': 2.765380684871657, 'TiF': 1.935869581032637, 'LiAs': 2.6509039685670372, 'IrF': 2.00064069782114, 'PdHg': 3.013470929362246, 'RuO': 2.0808557314850673, 'SbAu': 2.6139642879702363, 'PdS': 2.3164993475310482, 'PrN': 2.6378928410773304, 'SiSi': 2.333120192077337, 'VF': 2.0025491317270943, 'CdC': 2.2862026950838907, 'AuCo': 2.7269736171525882, 'IP': 2.4371528512098823, 'IS': 2.655698099201238, 'CdO': 2.3640565585786337, 'SmBr': 2.95766406762764, 'InPd': 2.6144501680005563, 'IV': 2.711715275989715, 'PdGe': 2.382429470550691, 'SrCu': 3.6178715440000007, 'VS': 2.362532068550082, 'TmTe': 3.0041837506819173, 'AuCr': 2.7364235952976736, 'VV': 2.251505188520253, 'IN': 2.245471524825764, 'PdPd': 2.819957604538408, 'ClFe': 2.272219683120052, 'IB': 2.1727332432332442, 'TbTb': 3.5313312014386304, 'TeSi': 2.448078893959033, 'EuGa': 3.3239707272141943, 'SiV': 2.5351500703345464, 'HgOs': 2.6857396311781496, 'CoS': 2.2892767595478007, 'SnO': 2.1868370314452723, 'CsS': 3.6035052082476287, 'UTe': 3.168840296970683, 'InIn': 2.860215501762649, 'AgS': 2.5262264024638936, 'WSe': 2.4573477009058915, 'HCo': 1.6190908456460569, 'IRh': 2.70584161608011, 'HoN': 2.4867997024268416, 'CoC': 1.9491530430992414, 'PtTi': 2.76522749261721, 'GeZn': 2.3817424013307638, 'BrDy': 2.713760052126478, 'RuSb': 2.5773039727717753, 'HPb': 1.6019427005920979, 'CuS': 2.3196400605269685, 'CoO': 2.0677373237784495, 'HPd': 1.660638726757391, 'LiZr': 2.918254536793423, 'PrS': 2.9740852581366237, 'BB': 1.7774699182943026, 'CuTe': 2.6357145298244866, 'InSi': 2.5905654714068183, 'OHf': 2.1419968555491993, 'OHg': 2.6126935884555444, 'ThTh': 3.9460972251098703, 'PbSe': 2.7697042914016428, 'BK': 3.2463244500753254, 'GaO': 1.9407729967481644, 'BV': 2.3666924604134296, 'BW': 2.3458572534248763, 'BP': 1.945299661925178, 'DYb': 2.0678145291617493, 'ClTb': 2.7053296682040355, 'BrAs': 2.649334933041739, 'GaF': 1.9790109725159994, 'FeAs': 2.3904306133944653, 'CuPt': 2.768228875878326, 'CdSe': 2.7047965520834736, 'ThH': 2.4036330516146815, 'TeCl': 2.5733312190033573, 'RePd': 2.8256381813683964, 'FeAO': 2.041850202280819, 'HBe': 1.4460738345243147, 'OsRu': 2.8992151196141105, 'TeRe': 2.6653356060991147, 'SnRh': 2.6115067906423843, 'HBa': 2.761770328414377, 'CCs': 3.5133284185716915, 'BFe': 2.130917525521103, 'OsRh': 2.7325751988047062, 'InNi': 2.3042346232336426, 'NAu': 2.062847212480521, 'MnSe': 2.453687709909738, 'NbAs': 2.6270542290899774, 'HfS': 2.638270891477117, 'MnSb': 2.4920323534526383, 'HBr': 1.6399934171815664, 'ReMo': 2.9696533784158365, 'OO': 1.366154117806917, 'FeTa': 2.9639828435248283, 'PtCo': 2.6009149802305007, 'CrRu': 2.7616608182333824, 'BaGe': 3.4142224604860796, 'MnTe': 2.6601277171364224, 'TaZn': 2.702070202430332, 'PbEO': 2.5434937929290355, 'ZnFe': 2.4529837270911403, 'YbSi': 3.101443374327582, 'SnP': 2.5976962697528543, 'YbSe': 2.7604614271640573, 'NBi': 2.505030520934481, 'PtCr': 2.5251836334557423, 'IrHg': 2.5738320252620506, 'SbSb': 2.802680911333311, 'ON': 1.258434482640632, 'PbP': 2.715035502790811, 'SmP': 3.0272079026226417, 'CCu': 2.0343671036932784, 'TeTe': 2.8055992871049424, 'IrOs': 2.857131315770538, 'OC': 1.3189485764079245, 'CuCs': 3.9400890067265504, 'MnMn': 2.754770161042098, 'RhGa': 2.4880204728457915, 'CuCr': 2.7691574542932353, 'WMn': 3.1349806301390886, 'OY': 2.352992208428384, 'PbB': 2.4092985756007863, 'OW': 1.9611397596509963, 'FeRu': 2.7300356785462463, 'PbO': 2.623598335224935, 'HgNb': 2.790107922609078, 'HSi': 1.430642116605243, 'MoPb': 2.5495315889143177, 'GeF': 1.7741034411618908, 'SiPd': 2.3406570033860965, 'OTh': 2.4673425489449947, 'SbSi': 2.5686626529398753, 'BrLi': 2.5637982323354827, 'InSe': 2.6230330533262465, 'YbAl': 3.0136499911929375, 'PTh': 3.0283849999875194, 'IO': 2.0595984596278045, 'SbSn': 2.913667440859153, 'HP': 1.2894155398572211, 'BrRe': 2.5563675398710677, 'BrRb': 3.434812600128606, 'HW': 1.7549184557854658, 'PTa': 2.5418840961757767, 'InGa': 2.6571565000040525, 'TmSe': 2.975830377906687, 'HH': 0.7414, 'NRh': 2.0848430691878663, 'SmGa': 3.2912064011752973, 'HO': 0.8674578130425501, 'ZnZn': 2.2778086052719817, 'HB': 1.1075833040744982, 'SbSe': 2.7310487828817545, 'PtCd': 2.848252054144855, 'TbMn': 3.410803246132555, 'ZnCl': 2.254890064907885, 'SCr': 2.2976858376894174, 'TmN': 2.4650454616816235, 'ZrCu': 2.657093583395445, 'MoCr': 3.1069983964812864, 'HV': 1.7647753693446877, 'PdC': 2.0971946363908947, 'SeOs': 2.4989919260534004, 'CB': 1.6603229435721027, 'PtBr': 2.4965162833705947, 'HgRe': 2.748063547849207, 'CaI': 3.145951057474339, 'PbF': 2.276231596583891, 'FN': 1.4164453513944077, 'BaP': 3.272577245186612, 'MoBr': 2.5993800783258756, 'BSc': 2.416839639399264, 'HgRu': 2.761791364982894, 'CoFe': 2.5689338397676984, 'HRh': 1.717171961586181, 'KCe': 3.5966969299104297, 'NiCr': 2.714584819451809, 'ScO': 2.1093421970878667, 'ScN': 2.160112812438677, 'HTc': 1.7195701527692873, 'RhGe': 2.466565999916851, 'HTa': 1.8263096487373998, 'NbSb': 2.8648024012222635, 'TeRh': 2.6272256576663793, 'MnS': 2.468302310483927, 'NbSi': 2.6142208030735414, 'SiAu': 2.337338356343304, 'AuGe': 2.4065524901023805, 'OsI': 2.7603345938181256, 'IrPt': 2.6810729703926226, 'ZrBr': 2.6380218769752823, 'PdCl': 2.341846490023714, 'PdCo': 2.617267673479155, 'NpH': 2.1790892002553224, 'NiSn': 2.6669024060829956, 'MnB': 2.134621615003593, 'VCu': 2.6694459003326396, 'PtS': 2.309725772178123, 'ScP': 2.8138919275213743, 'MnF': 1.8366755459781807, 'PtP': 2.288529614718863, 'AsGa': 2.473482805693036, 'FPd': 2.0242914235951828, 'CeBr': 2.905610546965658, 'TcTc': 2.209833573094656, 'PtMo': 2.7595658805427985, 'PtMn': 2.701297129859302, 'ICo': 2.60210287876188, 'BO': 1.4383841448699741, 'LiSe': 2.618845228165374, 'PBe': 2.2165807157746453, 'SnFe': 2.6269513523609405, 'UO': 2.241559461856216, 'UN': 2.4562895632572403, 'UI': 3.0836117963190075, 'PdCN': 1.9940048915116648, 'OP': 1.5439212600607006, 'PBi': 2.634553641232749, 'AuK': 3.6805, 'UC': 2.7384208591925163, 'UB': 2.5669506678750764, 'ScS': 2.678926725654115, 'KS': 3.402341845456659, 'ZnCo': 3.038519021789625, 'TlSn': 2.5572693922217065, 'SnNa': 3.145821302893083, 'AmO': 2.4115088275062897, 'SCl': 2.0830058653533965, 'TlSe': 3.1739068475869754, 'LaLa': 2.829192387773285, 'PMn': 2.320452667753737, 'CGe': 1.9734039758337243, 'CGd': 2.654349995035471, 'VC': 2.246765414398232, 'RhIr': 2.7934392444080487, 'FeMn': 2.7184543549981437, 'HY': 2.219192245678811, 'TaIr': 2.5554529530482566, 'FeCr': 2.789046571763547, 'AmS': 2.9215373685372836, 'ILi': 2.7859594602278195, 'OZr': 2.1606514111110204, 'SiIr': 2.3317922401531193, 'CrRh': 2.6939610289946576, 'SrSr': 0.7335446187048932, 'CsSe': 3.6741559895451408, 'YCl': 2.6707856359285924, 'NN': 1.330033935948112, 'CuP': 2.2645219444700997, 'SbNa': 3.264259481527289, 'FeO': 2.0274719191168713, 'CuFe': 2.7460524750111515, 'HAu': 1.6059511150607075, 'FeV': 2.800544813190718, 'TlAs': 2.762235584787372, 'HIr': 1.5870143970375838, 'FeP': 2.267795233825141, 'SrSe': 3.4923507376705603, 'NP': 1.64217971831557, 'CuN': 2.0278504125653027, 'NEu': 2.6020031263149725, 'CuK': 3.640395692368059, "S'S": 2.018293605521786, 'NV': 2.1103381712349356, 'DSi': 1.6034684787562132, 'AlBi': 3.0047009028390748, 'NiRu': 2.5765066144273767, 'LiFe': 2.465107370524253, 'LuOs': 3.039382169422943, 'ZnB': 2.426485294994012, 'AuW': 2.7558885169332727, "C'N": 1.477161550193356, 'PbDO': 2.8671791861165157, 'PTc': 2.417624388341136, "C'C": 1.496568392644728, 'CTc': 1.9793217953001914, 'SmB': 2.845753595951782, 'TeF': 1.8986546637431971, 'GeNa': 3.0972966113278764, 'PbAs': 2.8142808504111225, 'CoTe': 2.534643134036995, 'NbO': 2.012097014683972, 'WIr': 2.8292739378632423, 'YbCl': 2.6483417090945385, 'TeP': 2.4641931841422857, 'IrCu': 2.8320672566591076, 'BRe': 2.232990357070685, 'TeS': 2.61218338202344, 'NbF': 1.9896177225113034, 'VZn': 2.6447190553809197, 'SiPt': 2.3441149024036325, 'TlAu': 3.0406064363772076, 'NiSe': 2.3243977826515905, 'RuY': 3.0801167184510576, 'MnAl': 2.393184738260873, 'CNi': 2.000529859640795, 'CuBN': 1.9162655200280156, 'CuBO': 2.0102368460050157, 'NaF': 2.354232880002539, 'PbAO': 2.5296901442768607, 'BrCd': 2.7167330155619593, 'NdC': 2.778936569405223, 'OMoA': 1.9840180201893027, 'FEu': 2.2901996187792437, 'CuCuA': 2.625194962612406, 'RuTe': 2.70196271164699, 'HgP': 2.4982307527369354, 'PdMo': 2.8532858954245146, 'PdMn': 2.690400108143116, 'RuTi': 2.970010556713201, 'DyI': 3.198165189854756, 'CrH': 1.7481685377172924, 'ONiD': 2.044689934185208, 'USe': 2.953318875573285, 'SnF': 2.1754286103055325, 'SBr': 2.077748072379574, 'NpN': 2.61085148278171, 'SnB': 2.3532215724828234, 'ClRe': 2.395630522912228, 'LuSb': 3.232531029868632, 'InOs': 2.714185512294792, 'SnN': 2.283129893118019, 'SiTa': 2.705266727547982, 'LaF': 2.7913874730974384, 'HAg': 2.3946701840897844, 'NpF': 2.3237400527370276, 'SnK': 3.833953476759634, 'SnH': 1.7855694679968206, 'ClEu': 2.736817626469801, 'MoFO': 1.944543010804509, 'AuTe': 2.762099508775193, 'MoV': 3.195637103181653, 'SmCl': 2.765633456308135, 'BrNd': 2.86656399594838, 'HC': 0.9661353550249995, 'ReMn': 3.0202719044856976, 'RuAs': 2.4527478657734143, 'RhRh': 2.7608586400186748, 'BaBa': 4.129486096503967, "HO'": 0.8418044800576917, 'SbIn': 2.8528990740006415, 'ClHg': 2.5653789379130254, 'BiIn': 2.963422884591158, 'FPt': 2.0951017084468093, 'SAgA': 2.5061889074168247, 'AsZn': 2.4491663879014567, 'NiDN': 2.0144755446730023, 'BiSn': 3.056186129488765, 'CIr': 2.120433121775423, 'MgFe': 2.617695179260555, 'CdMn': 2.682625327547832, 'OEr': 2.3504500247873876, 'SbNi': 2.574306799711436, 'PtPd': 2.8163665028307534, 'BiSe': 2.9157963791147963, 'SiO': 1.6308808035010696, 'SbFe': 2.4807199112048903, 'TeIn': 2.775521272985103, 'GeBr': 2.42767767637931, 'NbCu': 2.791592502109571, 'CeI': 3.1640940392225985, 'GdS': 2.8860570102954908, 'PtD': 1.6804424680000007, 'HgRh': 2.683694182921498, 'AgRu': 2.8331366331627943, 'InP': 2.706573728814413, 'BrU': 2.815333937615591, 'GeRu': 2.4267829769754092, 'PbI': 3.2302743599466854, 'HgGe': 2.7554270773232354, 'BrO': 1.628063827390518, 'BrN': 2.055150724988527, 'LuS': 2.7348472241410446, 'LuTe': 2.9934878798936313, 'HMo': 1.7944723983635615, 'OOs': 2.0156460305860606, 'CoMn': 2.590284898303259, 'TaSe': 2.8017969138802776, 'RuGa': 2.376290189068315, 'CaAs': 3.0738493126243043, 'BrB': 1.9637960589107368, 'GdN': 2.5438536795049265, 'NaBr': 3.0132736732014007, 'VO': 1.8810400502279732, 'AuRu': 2.8180689682210467, 'NaBO': 2.337386839179595, 'SiBi': 2.651392984933389, 'NbNb': 2.932778593577897, 'ErEr': 3.540101472857107, 'ScF': 2.0770747424587666, 'ICl': 2.5535690923632424, 'LaTi': 3.502310262831787, 'IrFe': 2.8414910593510387, 'TaTa': 2.876880910851504, 'NaNa': 1.8750472892443093, 'ThB': 2.5564250810121734, 'AgBr': 2.7290231196494696, 'AuRe': 2.88913951592911, 'SBa': 3.35995859854814, 'CNb': 2.339840065560297, 'NdB': 2.888237619938327, 'ICs': 3.924340242444734, 'ICr': 2.7090468522482323, 'HCs': 3.0788147650609368, 'SnAs': 2.730918123955366, 'AuRh': 2.7957678776342045, 'ICu': 2.6691189889543376, 'InI': 2.714599492673125}

reactant = Species(xyz_file_to_atoms('../script/reactant.xyz'))
product = Species(xyz_file_to_atoms('../script/reactant.xyz'))

bonds = '../script/bonds.txt'
with open(bonds, 'r') as f:
    lines = f.read()
reactant_bonds = [(i[0]-1, i[1]-1) for i in eval(lines)]
product_bonds = [(1,2,1),(1,3,1),(1,4,1),(1,5,1),(1,13,1),(2,3,1),(2,4,1),(2,6,1),(3,5,1),(3,6,1),(3,10,1),(4,5,1),(4,6,1),(4,21,1),(5,6,1),(7,14,1),(8,15,1),(9,10,1),(9,13,1),(10,11,1),(10,16,1),(11,12,1),(11,17,1),(11,20,1),(12,13,1),(12,18,1),(13,14,1),(14,15,1),(14,19,1)]
product_bonds = [(i[0]-1, i[1]-1) for i in product_bonds]

make_graph(reactant, bond_list= reactant_bonds)
make_graph(product, bond_list= product_bonds)

a = get_bond_rearrangs(reactant, product, 'a')