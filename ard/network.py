import os
import time
import psutil

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
from mopac import mopac
import copy
import shutil
import time
#database
from connect import db

info = psutil.virtual_memory()

class Network(object):

    def __init__(self, reac_mol, forcefield, **kwargs):
        self.reac_mol = reac_mol
        self.nbreak = int(kwargs['nbreak'])
        self.nform = int(kwargs['nform'])
        self.forcefield = forcefield
        self.dh_cutoff = float(kwargs['dh_cutoff'])
        self.output_dir = kwargs['output_dir']
        self.reactions = {}
        self.network_prod_mols = []
        self.add_bonds = []
        self.ard_path = kwargs['ard_path']
        self.generations = kwargs['generations']
        self.method = kwargs["dh_cutoff_method"]
        

    def genNetwork(self, mol_object, **kwargs):
        """
        Execute the automatic reaction discovery procedure.
        """
        collection = db['qm_calculate_center']
        initial_pool = db['pool']
        statistics = db['statistics']
        targets = list(initial_pool.find({}, {'reactant_inchi_key':1}))
        #Add all reactant to a list for pgen filter isomorphic
        inchi_key_list = [i['reactant_inchi_key'] for i in targets]
        gen = Generate(mol_object, inchi_key_list)
        gen.generateProducts(nbreak=self.nbreak, nform=self.nform)
        prod_mols = gen.prod_mols
        add_bonds = gen.add_bonds
        break_bonds = gen.break_bonds

        # Load thermo database and choose which libraries to search
        thermo_db = ThermoDatabase()
        thermo_db.load(os.path.join(settings['database.directory'], 'thermo'))
        thermo_db.libraryOrder = ['primaryThermoLibrary', 'NISTThermoLibrary', 'thermo_DFT_CCSDTF12_BAC',
                                    'CBS_QB3_1dHR', 'DFT_QCI_thermo', 'BurkeH2O2', 'GRI-Mech3.0-N', ]
        # Filter reactions based on standard heat of reaction
        if self.method == "mopac":
            H298_reactant = mopac(mol_object, self.forcefield)
            H298_reac = H298_reactant.mopac_get_H298(mol_object)
            prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_mopac(H298_reac, mol)]
        else:
            H298_reac = self.reac_mol.getH298(thermo_db)
            prod_mols_filtered = [mol for mol in prod_mols if self.filterThreshold(H298_reac, mol, thermo_db)]

        #check product isomorphic and filter them
        if self.nbreak == 3 and self.nform == 3:
            gen_2.generateProducts(nbreak=2, nform=2)
            prod_mols_2 = gen_2.prod_mols
            #prod_mols_filtered_2 after filter by delta H
            if self.method == "mopac":
                prod_mols_filtered_2 = [mol for mol in prod_mols if self.filter_dh_mopac(H298_reac, mol)]
            else:
                prod_mols_filtered_2 = [mol for mol in prod_mols_2 if self.filterThreshold(H298_reac, mol, thermo_db)]
            #prod_mols_filtered_2 after filter by isomorphic
            prod_mols_filtered_2 = self.unique_key_filterIsomorphic_itself(prod_mols_filtered_2)
            prod_mols_filtered = self.unique_key_filterIsomorphic(prod_mols_filtered, prod_mols_filtered_2)
            prod_mols_filtered = self.unique_key_filterIsomorphic_itself(prod_mols_filtered)
            prod_mols_filtered += prod_mols_filtered_2 
        else:
            prod_mols_filtered = self.unique_key_filterIsomorphic_itself(prod_mols_filtered)

        # initial round add all prod to self.network
        reactant_name = mol_object.toRMGMolecule().to_inchi_key()
        prod_mols_filtered = self.unique_key_filterIsomorphic(reactant_name, prod_mols_filtered)

        for idx, mol in enumerate(prod_mols_filtered):
            index = prod_mols.index(mol)
            self.network_prod_mols.append(mol)
            # gen geo return path
            dir_path = self.gen_geometry(mol_object, mol, add_bonds[index], break_bonds[index])
            product_name = mol.toRMGMolecule().to_inchi_key()
            rxn_idx = '{}_{}'.format(reactant_name, idx+1)
            self.reactions[rxn_idx] = [reactant_name, product_name]
            collection.insert_one({
                                    rxn_idx: [reactant_name, product_name], 
                                   'Reactant SMILES':mol_object.write('can').strip(), 
                                   'reactant_inchi_key':reactant_name, 
                                   'product_inchi_key':product_name, 
                                   'Product SMILES':mol.write('can').strip(), 
                                   'path':dir_path, 
                                   'ssm_status':'job_unrun', 
                                   'generations':self.generations
                                   }
                                  )
        statistics.insert_one({'Reactant SMILES':mol_object.write('can').strip(), 'reactant_inchi_key':reactant_name, 'add how many products':len(prod_mols_filtered)})


        
    def filterThreshold(self, H298_reac, prod_mol, thermo_db):
        """
        Filter threshold based on standard enthalpies of formation of reactants
        and products. Returns `True` if the heat of reaction is less than
        `self.dh_cutoff`, `False` otherwise.
        """
        H298_prod = prod_mol.getH298(thermo_db)
        dH = H298_prod - H298_reac
        if dH < self.dh_cutoff:
            return 1
        return 0
    
    def filter_dh_mopac(self, H298_reac, prod_mol):
        H298_product = mopac(prod_mol, self.forcefield)
        H298_prod = H298_product.mopac_get_H298(prod_mol)
        dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            return 1
        return 0


    def unique_key_filterIsomorphic(self, reactant_name, compare):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        collect = db['qm_calculate_center']
        collection = db['reactions']
        targets = list(collect.find({'ts_status':'job_success'}))
        base_unique = [i['product_inchi_key'] for i in targets]
        compare_unique = [mol.toRMGMolecule().to_inchi_key() for mol in compare]
        isomorphic_idx = [compare_unique.index(i) for i in set(compare_unique) - set(base_unique)]
        result = [compare[i] for i in isomorphic_idx]
        
        product_pool = db['same_product']
        same_unique_key = list(set(compare_unique) & set(base_unique))

        for i in same_unique_key:
            reactant_target = list(collection.find({'reactant_inchi_key':reactant_name}))
            number = len(reactant_target)
            path_target = list(collection.find({'product_inchi_key':i}))
            reactions_name = '{}_{}'.format(reactant_name, number+1)
            collection.insert_one({
                                   reactions_name:[reactions_name, i],
                                   'reactant_smi':reactant_target[0]['Reactant SMILES'],
                                   'product_smi':path_target[0]['Product SMILES'],
                                   'path':path_target[0]['path'],
                                   'generations':self.generations})

        return result

    def unique_key_filterIsomorphic_itself(self, base):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        
        base_unique = [mol.toRMGMolecule().to_inchi_key() for mol in base]
        result = [base[base_unique.index(i)] for i in set(base_unique)]
        
        return result


    def gen_geometry(self, reactant_mol, network_prod_mol, add_bonds, break_bonds, **kwargs):
        #database
        rxn = db['qm_calculate_center']
        # These two lines are required so that new coordinates are
        # generated for each new product. Otherwise, Open Babel tries to
        # use the coordinates of the previous molecule if it is isomorphic
        # to the current one, even if it has different atom indices
        # participating in the bonds. a hydrogen atom is chosen
        # arbitrarily, since it will never be the same as any of the
        # product structures.
        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
        # Generate 3D geometries
        reactant_mol.gen3D(forcefield=self.forcefield, make3D=False)
        network_prod_mol.gen3D(forcefield=self.forcefield, make3D=False)

        reactant_mol_copy, network_prod_mol_copy= reactant_mol.copy(), network_prod_mol.copy()
        try:
            arrange3D = gen3D.Arrange3D(reactant_mol, network_prod_mol)
            msg = arrange3D.arrangeIn3D()
            if msg != '':
                print(msg)
        except:
            reactant_mol, network_prod_mol = reactant_mol_copy, network_prod_mol_copy

        ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
        reactant_mol.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)
        network_prod_mol.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)

        reactant = reactant_mol.toNode()
        product = network_prod_mol.toNode()
        subdir = os.path.join(os.path.dirname(self.ard_path), 'reactions')
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        dirname = network_prod_mol.toRMGMolecule().to_inchi_key()
        targets = list(rxn.find({'product_inchi_key':dirname}))
        dirname = '{}_{}'.format(dirname, len(targets)+1)
        output_dir = util.makeOutputSubdirectory(subdir, dirname)
        kwargs['output_dir'] = output_dir
        self.makeInputFile(reactant, product, **kwargs)
        self.makeCalEnergyFile(product, **kwargs)
        self.makeDrawFile(reactant, 'reactant.xyz', **kwargs)
        self.makeDrawFile(product, 'product.xyz', **kwargs)
        self.makeisomerFile(add_bonds, break_bonds, **kwargs)
        return output_dir

        
    @staticmethod
    def makeInputFile(reactant, product, **kwargs):
        """
        Create input file for TS search and return path to file.
        """
        path = os.path.join(kwargs['output_dir'], 'input.xyz')
        nreac_atoms = len(reactant.getListOfAtoms())
        nproduct_atoms = len(product.getListOfAtoms())

        with open(path, 'w') as f:
            f.write('{}\n\n{}\n{}\n\n{}\n'.format(nreac_atoms, reactant, nproduct_atoms, product))

        return path

    @staticmethod
    def makeCalEnergyFile(_input, exchange = 'b3lyp', basis = '6-31g*', spin = 0, multiplicity = 1, **kwargs):
        """
        Create input file for energy calculation.
        """
        path = os.path.join(kwargs['output_dir'], 'energy.in')
        with open(path, 'w') as f:
            f.write('$rem\njobtype = freq\nexchange = {}\nbasis = {}\n$end\n'.format(exchange, basis))
            f.write('\n$molecule\n{} {}\n{}\n$end'.format(spin, multiplicity, _input))

    @staticmethod
    def makeDrawFile(_input, filename = 'draw.xyz', **kwargs):
        """
        Create input file for network drawing.
        """
        path = os.path.join(kwargs['output_dir'], filename)
        ninput_atoms = len(_input.getListOfAtoms())

        with open(path, 'w') as f:
            f.write('{}\n\n{}'.format(ninput_atoms, _input))

    @staticmethod
    def makeisomerFile(add_bonds, break_bonds, **kwargs):
        """
        Create input file(add which bonds) for Single ended String Method (SSM) calculation.
        only for break 2 form 2 if more then need modify
        """
        try:
            path = os.path.join(kwargs['output_dir'], 'add_bonds.txt')

            with open(path, 'w') as f:
                for i in add_bonds:
                    f.write('ADD {} {}\n'.format(i[0]+1,i[1]+1))
                for i in break_bonds:
                    f.write('BREAK {} {}\n'.format(i[0]+1,i[1]+1))
        except:
            print("maybe list out of range, check add bond")

