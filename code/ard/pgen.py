#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains the :class:`Generate` for generating product structures.
"""

#third party
import pybel
import networkx as nx

# local application imports
import props
import constants
import gen3D
import numpy as np
from graph import make_graph, is_isomorphic
from bond_rearrangement import get_bond_rearrangs, BondRearrangement
from atoms import Atom
from species import Species

ELEMENT_TABLE = props.ElementData()
###############################################################################

class StructureError(Exception):
    """
    An exception class for errors that occur while generating structures.
    """
    pass

###############################################################################

class Generate(object):
    """
    Generation of product structures.
    The attributes are:

    =============== ========================== ================================
    Attribute       Type                       Description
    =============== ========================== ================================
    `reac_mol`      :class:`gen3D.Molecule`    A molecule object for the reactant structure
    `reac_smi`      ``str``                    Canonical SMILES of the reactant structure
    `atoms`         ``tuple``                  A tuple containing the atomic numbers of reactant/product structures
    `prod_mols`     ``list``                   A list of :class:`gen3D.Molecule` product structures
    =============== ========================== ================================

    Note: Bonds are represented as
          (beginAtomIdx, endAtomIdx, bondOrder)
    """

    def __init__(self, reac_mol, reactant_inchikey, reactant_graph, bond_dissociation_cutoff):
        self.reac_mol = reac_mol
        self.reactant_inchikey = reactant_inchikey
        self.reac_mol_graph = reactant_graph
        self.bond_dissociation_cutoff = float(bond_dissociation_cutoff)
        self.atoms = None
        self.prod_mols = []
        self.add_bonds = []
        self.break_bonds = []
        self.initialize()
        coords = [atom.coords for atom in self.reac_mol]
        self.reactant_coords = [np.array(coords).reshape(len(self.atoms), 3)]
        self.constraint = []
        for idx, i in enumerate(self.atoms):
            if i == 28:
                self.constraint.append(idx)
                
    def initialize(self):
        """
        Set the canonical SMILES for the reactant and extract the atomic
        numbers.
        """
        #self.reac_smi = self.reac_mol.write('can').split()[0]
        self.atoms = tuple(atom.atomicnum for atom in self.reac_mol)

    def DoU(self, products_bonds):
        """
        Degrees of Unsaturation, DoU is also known as Double Bond Equivalent
        DoU = (2*C+2+N-X-H)/2
        C is the number of carbons
        N is the number of nitrogens
        X is the number of halogens (F,Cl,Br,I)
        H is the number of hydrogens
        Oxygen and Sulfer are not included in the formula because saturatuin is unaffected by these elements.
        So if Oxygen or Sulfer in atoms, pass DoU filtration.
        Now only consider double bonds and triple bonds, ignoring rings.
        """
        if 8 not in self.atoms and 16 not in self.atoms:
            count_atoms = dict((a,self.atoms.count(a)) for a in self.atoms)
            try:
                nH_atoms = int(count_atoms[1])
            except:
                nH_atoms = 0
            try:
                nC_atoms = int(count_atoms[6]) 
            except:
                nC_atoms = 0
            try:
                nN_atoms = int(count_atoms[7])
            except:
                nN_atoms = 0
            try:
                nF_atoms = int(count_atoms[9])
            except:
                nF_atoms = 0
            try:
                nCl_atoms = int(count_atoms[17])
            except:
                nCl_atoms = 0
            try:
                nBr_atoms = int(count_atoms[35])
            except:
                nBr_atoms = 0

            DoU = (nC_atoms*2 + 2 + nN_atoms -nF_atoms - nCl_atoms - nBr_atoms - nH_atoms)/2
            bond_order =[]
            for i in products_bonds:
                double = [j[2] for j in i]
                bond_order.append(double)
            nDouble_bonds = []
            for order in bond_order:
                nDouble = dict((o,order.count(o)) for o in order)
                nDouble_bonds.append(nDouble)

            products_bonds = list(products_bonds)
            del_idx = []
            for idx,order in enumerate(nDouble_bonds):
                try:
                    if order[3] != 0:
                        order[2] += 2*order[3]
                        if order[2] > DoU:
                            del_idx.append(idx)
                except:
                    pass
                try:
                    if order[2] > DoU:
                        del_idx.append(idx)
                except:
                    pass
            del_idx = set(del_idx)
            products_bonds_idx = [i for i in range(len(products_bonds))]
            products_bonds_idx = set(products_bonds_idx)
            index = list(products_bonds_idx - del_idx)
            products_bonds = set([products_bonds[i] for i in index])
        else:
            products_bonds = products_bonds
        return products_bonds

    def generateProducts(self, nbreak=2, nform=2):
        """
        Generate all possible products from the reactant under the constraints
        of breaking a maximum of `nbreak` and forming a maximum of `nform`
        bonds.
        """
        if nbreak > 3 or nform > 3:
            raise Exception('Breaking/forming bonds is limited to a maximum of 3')

        # Extract bonds as an unmutable sequence (indices are made compatible with atom list)
        reactant_bonds = tuple(sorted(
            [(bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1, bond.GetBondOrder())
             for bond in pybel.ob.OBMolBondIter(self.reac_mol.OBMol)]
        ))
        # Extract valences as a mutable sequence
        reactant_valences = [atom.OBAtom.BOSum() for atom in self.reac_mol]
        # Initialize set for storing bonds of products
        # A set is used to ensure that no duplicate products are added
        products_bonds = set()

        # Generate all possibilities for forming bonds
        natoms = len(self.atoms)
        bonds_form_all = [(atom1_idx, atom2_idx, 1)
                          for atom1_idx in range(natoms - 1)
                          for atom2_idx in range(atom1_idx + 1, natoms)]
        # Generate products
        bf_combinations = ((0, 1), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3))
        for bf in bf_combinations:
            if bf[0] <= nbreak and bf[1] <= nform:
                self._generateProductsHelper(
                    bf[0],
                    bf[1],
                    products_bonds,
                    reactant_bonds,
                    reactant_valences,
                    bonds_form_all
                )

        if products_bonds:
            #Filter the products_bonds which doesn't follow DoU rule
            #products_bonds = self.DoU(products_bonds)
            for bonds in products_bonds:
                #for SSM calculation
                break_bonds = []
                form_bonds = []
                for i in list(set(bonds) ^ set(reactant_bonds)):
                    if i not in bonds:
                        # which is breaked
                        break_bonds.append(i)
                    else:
                        # which is formed
                        form_bonds.append(i)
                # deal with double bonds
                for i in form_bonds:
                    if i[2] > 1 and i in bonds:
                        for j in break_bonds:
                            if i[0] == j[0] and i[1] == j[1]:
                                break_bonds.remove(j)
                # deal with double bonds
                for i in break_bonds:
                    if i[2] > 1 and i in reactant_bonds:
                        for j in form_bonds:
                            if i[0] == j[0] and i[1] == j[1]:
                                form_bonds.remove(j)

                if self.check_bond_type(bonds):
                    for i in form_bonds:
                        if i[2] == 2:
                            form_bonds.remove(i)
                    for i in break_bonds:
                        if i[2] == 2:
                            break_bonds.remove(i)

                    #dist = self.check_bond_length(self.reactant_coords, form_bonds)

                    self.add_bonds.append(form_bonds)
                    self.break_bonds.append(break_bonds)
                    mol = gen3D.makeMolFromAtomsAndBonds(self.atoms, bonds, spin=self.reac_mol.spin)
                    mol.setCoordsFromMol(self.reac_mol)
                    if mol.write('inchiKey').strip() not in self.reactant_inchikey:
                        if self.check_bond_dissociation_energy(bonds, break_bonds):
                            self.prod_mols.append(mol)

    def check_bond_type(self, bonds):
        
        bond_type = {}
        
        for i in range(len(self.atoms)):
            num = 0
            for j in bonds:
                if j[0] == i or j[1] == i:
                    num += j[2]
            bond_type[i] = num
        tmp = []
        for idx, i in enumerate(self.atoms):
            if i == 6:
                if bond_type[idx] > 4:   # use !=  or  >   need test
                    tmp.append(idx)
        if len(tmp) == 0:
            return 1
        else:
            return 0
    
    def check_bond_dissociation_energy(self, bond_list, bbond_list):
        energy = 0.0

        reactant_graph = self.reac_mol_graph
        atoms = []
        for i in range(reactant_graph.n_atoms):
            atoms.append(Atom(atomic_symbol=reactant_graph.atoms[i].label))
        product = Species(atoms)

        graph = nx.Graph()

        for i in range(reactant_graph.n_atoms):
            graph.add_node(i, atom_label=reactant_graph.atoms[i].label, stereo=False)

        if bond_list is not None:
            [graph.add_edge(bond[0], bond[1], pi=False, active=False) for bond in bond_list]
            product.graph = graph

        for break_bond in bbond_list:
            first_atom = reactant_graph.atoms[break_bond[0]].label
            second_atom = reactant_graph.atoms[break_bond[1]].label
            
            try:
                energy += props.bond_dissociation_energy[first_atom + second_atom]
            except:
                energy += props.bond_dissociation_energy[second_atom + first_atom]
        
        if energy/constants.cal_to_J >= self.bond_dissociation_cutoff:
            return False
        else:
            return True

    def check_bond_length(self, coords, add_bonds):
        """
        Use reactant coordinate to check if the add bonds's bond length is too long.
        Return a 'list of distance'.
        """
        dist = []
        for bond in add_bonds:
            coord_vect_1 = coords[0][bond[0]]
            coord_vect_2 = coords[0][bond[1]]
            diff = coord_vect_1 - coord_vect_2
            dist.append(np.linalg.norm(diff))
        return dist
        
    def _generateProductsHelper(self, nbreak, nform, products, bonds, valences, bonds_form_all, bonds_broken=None):
        """
        Generate products recursively given the number of bonds that should be
        broken and formed, a set for storing the products, a sequence of atoms,
        of bonds, and of valences. `bonds_form_all` should contain a tuple of
        tuples of bonds that contains all possibilities for forming bonds.
        Nothing is returned, but formed products are added to `products`.
        """
        if bonds_broken is None:
            bonds_broken = []
        if nbreak == 0 and nform == 0:
            # If no more bonds are to be changed, then add product (base case)
            if len(bonds_broken) == 2:
                if not (bonds_broken[0][0] in self.constraint and bonds_broken[0][1] in self.constraint):
                    if not (bonds_broken[1][0] in self.constraint and bonds_broken[1][1] in self.constraint):
                        products.add((tuple(sorted(bonds))))
            elif len(bonds_broken) == 3:
                if not (bonds_broken[0][0] in self.constraint and bonds_broken[0][1] in self.constraint):
                    if not (bonds_broken[1][0] in self.constraint and bonds_broken[1][1] in self.constraint):
                        if not (bonds_broken[2][0] in self.constraint and bonds_broken[2][1] in self.constraint):
                            products.add((tuple(sorted(bonds))))
            else:
                products.add((tuple(sorted(bonds))))
        if nbreak > 0:
            # Break bond
            for bond_break_idx, bond_break in enumerate(bonds):
                valences_break = self.changeValences(valences, bond_break, -1)
                bonds_break = self.breakBond(bonds, bond_break_idx)

                # Keep track of bonds that have been broken
                if bond_break_idx == 0:
                    bonds_broken.append(bond_break)
                else:
                    bonds_broken[-1] = bond_break

                # Call function recursively to break next bond
                self._generateProductsHelper(
                    nbreak - 1,
                    nform,
                    products,
                    bonds_break,
                    valences_break,
                    bonds_form_all,
                    bonds_broken
                )

            # Remove last bond that has been broken after loop terminates
            del bonds_broken[-1]
        elif nform > 0:
            # Form bond
            for bond_form in bonds_form_all:
                # Do not add bond if it has previously been broken
                if bond_form[:2] in [bond[:2] for bond in bonds_broken]:
                    continue
                if bond_form[0] in self.constraint and bond_form[1] in self.constraint:
                    continue
                # Form new bond and catch exception if it violates constraints
                try:
                    valences_form = self.changeValences(valences, bond_form, 1)
                    bonds_form = self.formBond(bonds, bond_form)
                except StructureError:
                    continue

                # Call function recursively to form next bond
                self._generateProductsHelper(
                    nbreak,
                    nform - 1,
                    products,
                    bonds_form,
                    valences_form,
                    bonds_form_all,
                    bonds_broken
                )

    @staticmethod
    def breakBond(bonds, break_idx):
        """
        Break a bond given a tuple of tuples of bonds in the form:
            (atom index 1, atom index 2, bond type)
        The bond at index `break_idx` is broken. A tuple of tuples of updated
        bonds is returned.
        """
        # Break double or triple bond
        if bonds[break_idx][2] > 1:
            return bonds[:break_idx] + (bonds[break_idx][:2] + (bonds[break_idx][2] - 1,),) + bonds[(break_idx + 1):]

        # Break single bond
        return bonds[:break_idx] + bonds[(break_idx + 1):]

    @staticmethod
    def formBond(bonds, new_bond):
        """
        Form a bond given a tuple of tuples of bonds in the form:
            (atom index 1, atom index 2, bond type)
        If a bond already exists in `bonds`, the bond order is incremented. If
        a bond between atoms that are not yet connected is to be formed, then a
        new bond is added to the tuple of bonds. The bond to be added is
        specified in `new_bond`. A tuple of tuples of updated bonds is
        returned.
        """
        # Ensure that only one bond is added at a time
        assert new_bond[2] == 1

        try:  # Check if bond exists as single bond
            idx = bonds.index(new_bond)
        except ValueError:
            try:  # Check if bond exists as double bond
                idx = bonds.index(new_bond[:2] + (2,))
            except ValueError:
                try:  # Check if bond exists as triple bond
                    idx = bonds.index(new_bond[:2] + (3,))
                except ValueError:  # Add new bond if it does not exist yet
                    return bonds + (new_bond,)
                else:  # Raise exception if trying to exceed triple bond
                    raise StructureError('Bond type cannot be higher than triple bond for bond {}'.format(bonds[idx]))
            else:  # Return bonds with double bond increased to triple bond
                return bonds[:idx] + (bonds[idx][:2] + (3,),) + bonds[(idx + 1):]
        else:  # Return bonds with single bond increased to double bond
            return bonds[:idx] + (bonds[idx][:2] + (2,),) + bonds[(idx + 1):]

    def changeValences(self, valences, bond, inc):
        """
        Update the valences corresponding to each atom in `self.atoms` given
        the current valences in `valences`, the bond that is being affected,
        and the increment (typically 1 or -1). The maximum valences for each
        atom type are respected and an error is raised if they would be
        violated. A sequence of updated valences is returned.
        """
        # Create copy of current valences
        valences_temp = valences[:]

        # Check for invalid operation
        if inc < 0 and (valences_temp[bond[0]] < inc or valences_temp[bond[1]] < inc):
            raise Exception('Cannot decrease valence below zero-valence')

        # Change valences of both atoms participating in bond
        valences_temp[bond[0]] += inc
        valences_temp[bond[1]] += inc

        # Check if maximum valences are exceeded
        element0 = ELEMENT_TABLE.from_atomic_number(self.atoms[bond[0]])
        element1 = ELEMENT_TABLE.from_atomic_number(self.atoms[bond[1]])
        if valences_temp[bond[0]] > element0.max_bonds:
            raise StructureError('Maximum valence on atom {} exceeded'.format(bond[0]))
        if valences_temp[bond[1]] > element1.max_bonds:
            raise StructureError('Maximum valence on atom {} exceeded'.format(bond[1]))

        # Return valid valences
        return valences_temp

