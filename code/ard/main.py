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
        self.imaginarybond = int(imaginarybond)
        self.nbreak = int(nbreak)
        self.nform = int(nform)
        self.dh_cutoff = float(dh_cutoff)
        self.theory_low = theory_low
        self.forcefield = forcefield
        self.distance = float(distance)
        self.output_dir = output_dir

    def executeXYZ(self, **kwargs):
        
        reac_mol = self.reac_smi
        #reac_mol.gen3D(forcefield=self.forcefield)
        network = Network(reac_mol, forcefield = self.forcefield, **kwargs)
        network.genNetwork(reac_mol)
    
###############################################################################

def readInput(input_file):
    # Allowed keywords
    keys = ('reac_smi', 'imaginarybond', 'nbreak', 'nform', 'dh_cutoff', 'dh_cutoff_method', 'manual_bonds')
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
    return mol_obj