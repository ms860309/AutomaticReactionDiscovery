#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Contains the :class:`ARD` for running an automatic reaction discovery. This
includes filtering reactions, generating 3D geometries, and running transition
state searches.
"""

from __future__ import print_function

import os
import time

import pybel
from rmgpy import settings
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.molecule import Molecule

import constants
import gen3D
import util
from quantum import QuantumError
from node import Node
from pgen import Generate
from network import Network

from imaginary import Imaginary
from filter_rule import _filter
import openbabel as ob
import multiprocessing as mp
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
        reac_mol.gen3D(forcefield=self.forcefield)
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
    
def readXYZ(xyz, reactant_bonds = None, product_bonds = None):
    mol = next(pybel.readfile('xyz', xyz))
    mol = gen3D.Molecule(mol.OBMol, reactant_bonds, product_bonds)
    return mol
