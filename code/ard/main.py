#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains the :class:`ARD` for running an automatic reaction discovery. This
includes filtering reactions, generating 3D geometries, and running transition
state searches.
"""

# standard library imports
import os
import time

#third party
import pybel
import openbabel as ob
import multiprocessing as mp

# local application imports
import constants
import gen3D
from gen3D import Molecule
import util
from quantum import QuantumError
from node import Node
from pgen import Generate
from network import Network
from input_output import xyz_file_to_atoms
from species import Species
from graph import make_graph

###############################################################################

class ARD(object):
    """
    Automatic reaction discovery class. Filters reactions based on estimated
    thermo of reactant and products, generates force field 3D geometries, and
    runs transition state searches.
    The attributes are:

    =============== ======================== ==================================
    Attribute       Type                     Description
    =============== ======================== ==================================
    `reac_smi`      ``str``                  A valid SMILES string describing the reactant structure
    `nbreak`        ``int``                  The maximum number of bonds that may be broken
    `nform`         ``int``                  The maximum number of bonds that may be formed
    `dh_cutoff`     ``float``                Heat of reaction cutoff (kcal/mol) for reactions that are too endothermic
    `forcefield`    ``str``                  The force field for 3D geometry generation
    `distance`      ``float``                The initial distance between molecules
    `output_dir`    ``str``                  The path to the output directory
    =============== ======================== ==================================

    """

    def __init__(self, reac_smi, imaginarybond=0, nbreak=3, nform=3, dh_cutoff=20.0, theory_low=None,
                 forcefield='uff', distance=3.5, output_dir='', **kwargs):
        self.reac_smi = reac_smi
        self.reactant_graph = kwargs['graph']
        self.forcefield = forcefield

    def executeXYZ(self, **kwargs):
        
        reac_mol = self.reac_smi
        reactant_graph = self.reactant_graph
        #reac_mol.gen3D(forcefield=self.forcefield)
        network = Network(reac_mol, reactant_graph, forcefield = self.forcefield, **kwargs)
        network.genNetwork(reac_mol)
    
###############################################################################

def readInput(input_file):
    # Allowed keywords
    keys = ('reac_smi', 'imaginarybond', 'nbreak', 'nform', 'dh_cutoff', 'dh_cutoff_method', 'manual_bonds', 'graph', 'bond_dissociation_cutoff')
    # Read all data from file
    with open(input_file, 'r') as f:
        input_data = f.read().splitlines()
    # Create dictionary
    input_dict = {}
    # Extract remaining keywords and values
    for line in input_data:
        if line != '' and not line.strip().startswith('#'):
            key = line.split()[0].lower()
            if key not in keys:
                continue
            if line.split()[1] == '=':
                input_dict[key] = line.split()[2]
            else:
                input_dict[key] = line.split()[1]
                
    return input_dict

def extract_bonds(bonds):
    with open(bonds, 'r') as f:
        lines = f.read()
    lines = eval(lines)
    return lines
    
def readXYZ(xyz, bonds = None):
    mol = next(pybel.readfile('xyz', xyz))
    if bonds:
        m = Molecule(pybel.ob.OBMol())
        OBMol = m.OBMol
        for i in mol:
            a = pybel.ob.OBAtom()
            a.SetAtomicNum(i.atomicnum)
            a.SetVector(i.coords[0], i.coords[1], i.coords[2])
            OBMol.AddAtom(a)
        for bond in bonds:
            OBMol.AddBond(bond[0], bond[1], bond[2])
        mol_obj = gen3D.Molecule(OBMol)
    else:
        mol_obj = gen3D.Molecule(mol.OBMol)

    reactant_graph = Species(xyz_file_to_atoms(xyz))
    reactant_bonds = [(i[0]-1, i[1]-1) for i in bonds]
    make_graph(reactant_graph, bond_list= reactant_bonds)

    """
    #debuging
    bs = [(1,2,1),(1,3,1),(1,4,1),(1,5,1),(1,13,1),(2,3,1),(2,4,1),(2,6,1),(3,5,1),(3,6,1),(3,10,1),(4,5,1),(4,6,1),(4,21,1),(5,6,1),(7,14,1),(8,15,1),(9,10,1),(9,13,1),(10,11,1),(10,16,1),(11,12,1),(11,17,1),(11,20,1),(12,13,1),(12,18,1),(13,14,1),(14,15,1),(14,19,1)]
    atoms = tuple(atom.atomicnum for atom in mol_obj)
    mol = gen3D.makeMolFromAtomsAndBonds(atoms, bs)
    mol.setCoordsFromMol(mol_obj)

    Hatom = gen3D.readstring('smi', '[H]')
    ff = pybel.ob.OBForceField.FindForceField('UFF')

    mol_obj.separateMol()
    if len(mol_obj.mols) > 1:
        mol_obj.mergeMols()
    mol.separateMol()
    if len(mol.mols) > 1:
        mol.mergeMols()

    gen3D.make3DandOpt(mol_obj, 'UFF', make3D = False)
    gen3D.make3DandOpt(mol, 'UFF', make3D = False)

    arrange3D = gen3D.Arrange3D(mol_obj, mol)
    msg = arrange3D.arrangeIn3D()

    ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
    mol_obj.gen3D(make3D=False)
    ff.Setup(Hatom.OBMol)
    mol.gen3D(make3D=False)
    ff.Setup(Hatom.OBMol)

    print(mol.toNode())
    raise
    """
    return mol_obj, reactant_graph