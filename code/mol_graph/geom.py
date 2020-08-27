#third party
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import distance_matrix

def get_neighbour_list(species, atom_i):
    """Calculate a neighbour list from atom i as a list of atom labels
    Arguments:
        atom_i (int): index of the atom
        species (autode.species.Species):
    Returns:
        (list(int)): list of atom ids in ascending distance away from atom_i
    """
    coords = species.get_coordinates()
    distance_vector = cdist(np.array([coords[atom_i]]), coords)[0]

    dists_and_atom_labels = {}
    for atom_j, dist in enumerate(distance_vector):
        dists_and_atom_labels[dist] = species.atoms[atom_j].label

    atom_label_neighbour_list = []
    for dist, atom_label in sorted(dists_and_atom_labels.items()):
        atom_label_neighbour_list.append(atom_label)

    return atom_label_neighbour_list