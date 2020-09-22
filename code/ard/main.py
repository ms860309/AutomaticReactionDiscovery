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
from openbabel import openbabel as ob
from openbabel import pybel as pb
import multiprocessing as mp
import logging
import networkx as nx

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
    `reactant`      ``str``                  Reactant smiles string or reactant openbabel object
    `nbreak`        ``int``                  The maximum number of bonds that may be broken
    `nform`         ``int``                  The maximum number of bonds that may be formed
    `dh_cutoff`     ``float``                Heat of reaction cutoff (kcal/mol) for reactions that are too endothermic
    `forcefield`    ``str``                  The force field for 3D geometry generation
    `distance`      ``float``                The initial distance between molecules
    `output_dir`    ``str``                  The path to the output directory
    =============== ======================== ==================================

    """

    def __init__(self, reactant, forcefield='uff', distance=3.5, output_dir='', **kwargs):
        self.reactant = reactant
        self.reactant_graph = kwargs['graph']
        self.forcefield = forcefield
        self.nbreak = kwargs['nbreak']
        self.nform = kwargs['nform']
        self.dh_cutoff = float(kwargs['dh_cutoff'])
        self.bond_dissociation_cutoff = float(kwargs['bond_dissociation_cutoff'])
        if kwargs['use_inchi_key'] == '1':
            self.use_inchi_key = True
        else:
            self.use_inchi_key = False
        log_level = logging.INFO
        self.logger = util.initializeLog(log_level, os.path.join(os.getcwd(), 'ARD.log'), logname='main')
        self.initialize()

    def initialize(self):
        self.logger.info('######################################################################')
        self.logger.info('#################### AUTOMATIC REACTION DISCOVERY ####################')
        self.logger.info('######################################################################')
        self.logger.info('Reactant Geometry: \n{}'.format(str(self.reactant.toNode())))
        self.logger.info('Maximum number of bonds to be broken: ' + str(self.nbreak))
        self.logger.info('Maximum number of bonds to be formed: ' + str(self.nform))
        self.logger.info('Heat of reaction cutoff: {} kcal/mol'.format(self.dh_cutoff))
        self.logger.info('Bond dissociation energy cutoff: {} kcal/mol'.format(self.bond_dissociation_cutoff))
        self.logger.info('Force field for 3D structure generation: {}'.format(self.forcefield))
        self.logger.info('######################################################################\n')

    def executeXYZ(self, **kwargs):
        self.logger.info('ARD initiated on {} \n'.format(time.asctime()))
        start_time = time.time()
        
        reac_mol = self.reactant
        network = Network(reac_mol, self.reactant_graph, forcefield = self.forcefield, logger = self.logger, **kwargs)
        network.genNetwork(reac_mol, self.use_inchi_key, nbreak = kwargs['nbreak'], nform = kwargs['nform'])
        self.finalize(start_time)

    def finalize(self, start_time):
        """
        Finalize the job.
        """
        self.logger.info('\nARD terminated on ' + time.asctime())
        self.logger.info('Total ARD run time: {:.1f} s'.format(time.time() - start_time))

###############################################################################

def readInput(input_file):
    # Allowed keywords
    keys = ('reactant', 'nbreak', 'nform', 'dh_cutoff', 'dh_cutoff_method', 
            'manual_bonds', 'bond_dissociation_cutoff', 'constraint', 'use_inchi_key')
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

def extract_constraint_index(constraint):
    with open(constraint, 'r') as f:
        lines = f.read()
    lines = eval(lines)
    return lines

def readXYZ(xyz, bonds = None):
    # extract molecule information from xyz
    mol = next(pb.readfile('xyz', xyz))
    # Manually give bond information 
    # (Because in metal system the bond information detect by openbabel usually have some problem)
    if bonds:
        m = Molecule(pb.ob.OBMol())
        obmol = m.OBMol

        obmol.BeginModify()
        for atom in mol:
            coords = [coord for coord in atom.coords]
            atomno = atom.atomicnum
            obatom = ob.OBAtom()
            obatom.thisown = 0
            obatom.SetAtomicNum(atomno)
            obatom.SetVector(*coords)
            obmol.AddAtom(obatom)
            del obatom

        for bond in bonds:
            obmol.AddBond(bond[0], bond[1], bond[2])
        #obmol.ConnectTheDots()
        #obmol.PerceiveBondOrders()
        #obmol.SetTotalSpinMultiplicity(1)
        obmol.SetTotalCharge(int(mol.charge))
        obmol.Center()
        obmol.EndModify()
        mol_obj = gen3D.Molecule(obmol)

        reactant_graph = Species(xyz_file_to_atoms(xyz))
        reactant_bonds = [(i[0]-1, i[1]-1) for i in bonds]
        make_graph(reactant_graph, bond_list= reactant_bonds)
    else:
        mol_obj = gen3D.Molecule(mol.OBMol)
        reactant_graph = Species(xyz_file_to_atoms(xyz))
        reactant_bonds = tuple(sorted(
            [(bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1)
             for bond in pb.ob.OBMolBondIter(mol.OBMol)]
        ))
        make_graph(reactant_graph, bond_list= reactant_bonds)
    return mol_obj, reactant_graph