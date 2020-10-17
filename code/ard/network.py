# standard library imports
import os
import time
import psutil
import shutil

#third party
from openbabel import pybel
from rmgpy import settings
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.molecule import Molecule
from subprocess import Popen, PIPE
import numpy as np

# local application imports
import constants
import gen3D
import util
from quantum import QuantumError
from node import Node
from pgen import Generate
import mopac
from mopac import Mopac


#database
from connect import db

info = psutil.virtual_memory()

class Network(object):

    def __init__(self, reac_mol, reactant_graph, logger, **kwargs):
        self.reac_mol = reac_mol
        self.reactant_graph = reactant_graph
        self.logger = logger
        self.forcefield = kwargs['forcefield']
        self.constraintff_alg = kwargs['constraintff_alg']
        self.mopac_method = kwargs['mopac_method']
        self.dh_cutoff = float(kwargs['dh_cutoff'])
        self.bond_dissociation_cutoff = kwargs['bond_dissociation_cutoff']
        self.ard_path = kwargs['ard_path']
        self.generations = kwargs['generations']
        self.method = kwargs["dh_cutoff_method"]
        self.binding_mode_energy_cutoff = float(kwargs['binding_mode_energy_cutoff'])
        self.constraint = kwargs['constraint_index']
        self.binding_cutoff_select = kwargs['binding_cutoff_select']
        self.count = 0
        self.moac_reac = None

    def genNetwork(self, mol_object, use_inchi_key, nbreak, nform):
        """
        Execute the automatic reaction discovery procedure.
        """
        # Database
        qm_collection = db['qm_calculate_center']
        pool_collection = db['pool']
        statistics_collection = db['statistics']
        targets = list(pool_collection.find({'generations': 1}))

        # Add all reactant to a list for pgen filter isomorphic
        initial_reactant_inchi_key = targets[0]['reactant_inchi_key']

        # Reactant information
        reactant_key = mol_object.write('inchiKey').strip()  # inchikey

        # Generate all possible products
        gen = Generate(mol_object, initial_reactant_inchi_key, self.reactant_graph, self.bond_dissociation_cutoff, use_inchi_key, self.constraint)
        self.logger.info('Generating all possible products...')
        gen.generateProducts(nbreak = int(nbreak), nform =  int(nform))
        prod_mols = gen.prod_mols
        self.logger.info('{} possible products generated\n'.format(len(prod_mols)))
        add_bonds = gen.add_bonds
        break_bonds = gen.break_bonds
        # Filter reactions based on standard heat of reaction  delta H
        if self.method == "mopac":
            self.logger.info('Now use {} to filter the delta H of reactions....\n'.format(self.method))
            reac_mol_copy = mol_object.copy()
            if self.generations == 1:
                H298_reac = self.get_mopac_H298(mol_object)
                update_field = {'reactant_energy':H298_reac}
                pool_collection.update_one(targets[0], {"$set": update_field}, True)
                prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_mopac(mol_object, reac_mol_copy, mol, add_bonds[prod_mols.index(mol)], break_bonds[prod_mols.index(mol)], self.logger, len(prod_mols))]
            else:
                query = [{'$match':{'reactant_inchi_key':reactant_key}},
                        {'$group':{'_id':'$reactant_inchi_key', 'reactant_mopac_hf':{'$min':'$reactant_mopac_hf'}}}]
                H298_reac = list(qm_collection.aggregate(query))[0]['reactant_mopac_hf']
                prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_mopac(mol_object, reac_mol_copy, mol, add_bonds[prod_mols.index(mol)], self.logger, len(prod_mols), H298_reac)]
        else:
            self.logger.info('Now use {} to filter the delta H of reactions....\n'.format(self.method))
            # Load thermo database and choose which libraries to search
            thermo_db = ThermoDatabase()
            thermo_db.load(os.path.join(settings['database.directory'], 'thermo'))
            thermo_db.libraryOrder = ['primaryThermoLibrary', 'NISTThermoLibrary', 'thermo_DFT_CCSDTF12_BAC',
                                        'CBS_QB3_1dHR', 'DFT_QCI_thermo', 'BurkeH2O2', 'GRI-Mech3.0-N', ]
            if self.generations == 1:
                H298_reac = self.reac_mol.getH298(thermo_db)
                update_field = {'reactant_energy':H298_reac}
                pool_collection.update_one(targets[0], {"$set": update_field}, True)
            elif self.generations > 1:
                H298_reac = targets[0]['reactant_energy']
            prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_rmg(H298_reac, mol, thermo_db)]

        self.logger.info('After delta H filter {} product remain.\n'.format(len(prod_mols_filtered)))
        
        if self.method != "mopac":
            self.logger.info('Generate geometry........\n')
            for mol in prod_mols_filtered:
                index = prod_mols.index(mol)
                # Generate geometry and return path
                dir_path = self.gen_geometry(mol_object, mol, add_bonds[index], break_bonds[index])
                product_name = mol.write('inchiKey').strip()
                self.logger.info('\nReactant inchi key: {}\nProduct inchi key: {}\nDirectory path: {}\n'.format(reactant_key, product_name, dir_path))
                qm_collection.insert_one({
                                    'reaction': [reactant_key, product_name], 
                                    'Reactant SMILES':mol_object.write('can').split()[0], 
                                    'reactant_inchi_key':reactant_key, 
                                    'product_inchi_key':product_name, 
                                    'Product SMILES':mol.write('can').split()[0], 
                                    'path':dir_path, 
                                    'ssm_status':'job_unrun',
                                    'generations':self.generations
                                    })
            # Generate geometry and insert to database
            statistics_collection.insert_one({
                'Reactant SMILES':mol_object.write('can').split()[0], 
                'reactant_inchi_key':reactant_key, 
                'add how many products':len(prod_mols_filtered),
                'generations': self.generations})
        else:
            if self.binding_cutoff_select == 'lowest':
                self.logger.info('Now find the lowest binding mode energy(mopac)....')
                query = [{'$match':{'reactant_inchi_key':reactant_key}},
                        {'$group':{'_id':'$reactant_inchi_key', 'reactant_mopac_hf':{'$min':'$reactant_mopac_hf'}}}]
                minimum_energy = list(qm_collection.aggregate(query))[0]['reactant_mopac_hf']
                self.logger.info('The lowest energy of binding mode is {}'.format(minimum_energy))
            else:
                # Assume user is using starting reactant as reference
                self.logger.info('Now the binding mode energy reference will use the starting reactant mopac HF...')
                minimum_energy = H298_reac
                self.logger.info('The lowest energy of binding mode is {}'.format(minimum_energy))

            ssm_target_query = {'$and': 
                    [{ 'reactant_inchi_key':reactant_key},
                    {'reactant_mopac_hf':
                        {'$lte':minimum_energy + self.binding_mode_energy_cutoff}}
                    ]}
            ssm_unrun_targets = list(qm_collection.find(ssm_target_query))
            self.logger.info('There are {} products remain after binding energy filter'.format(len(ssm_unrun_targets)))
            for unrun_job in ssm_unrun_targets:
                update_field = {"ssm_status":"job_unrun"}
                qm_collection.update_one(unrun_job, {"$set": update_field}, True)
            delete_query = {'$and':
                    [{'reactant_inchi_key':reactant_key},
                    {'reactant_mopac_hf':
                        {'$gt':minimum_energy + self.binding_mode_energy_cutoff}}
                    ]}
            delete_targets = list(qm_collection.find(delete_query))
            for target in delete_targets:
                qm_collection.delete_one(target)

            # Generate geometry and insert to database
            statistics_collection.insert_one({
                'Reactant SMILES':mol_object.write('can').split()[0], 
                'reactant_inchi_key':reactant_key, 
                'add how many products':len(ssm_unrun_targets),
                'generations': self.generations})

    def filter_dh_rmg(self, H298_reac, prod_mol, thermo_db):
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

    def filter_dh_mopac(self, reac_obj, reac_mol_copy, prod_mol, form_bonds, break_bonds, logger, total_prod_num, refH = None):
        self.count += 1
        mopac_object = Mopac(reac_obj, prod_mol, self.mopac_method, self.forcefield, self.constraintff_alg, form_bonds, logger, total_prod_num, self.count, self.constraint)
        H298_reac, H298_prod, reactant, product = mopac_object.mopac_get_H298(reac_mol_copy)

        if H298_prod == False or H298_reac == False:
            return 0
        self.logger.info('Product energy calculate by mopac is {} kcal/mol and reactant is {} kcal/mol'.format(H298_prod, H298_reac))
        if refH:
            self.logger.info('In the {} generations, reactant hf use {} instead.'.format(self.generations, refH))
            dH = H298_prod - refH
        else:
            dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            self.logger.info('Delta H is {}, smaller than threshold'.format(dH))
            self.logger.info('Finished {}/{}'.format(self.count, total_prod_num))

            qm_collection = db['qm_calculate_center']
            dir_path = self.output(reactant, product, form_bonds, break_bonds, prod_mol)
            reactant_key = reac_obj.write('inchiKey').strip()
            product_name = prod_mol.write('inchiKey').strip()
            self.logger.info('\nReactant inchi key: {}\nProduct inchi key: {}\nDirectory path: {}\n'.format(reactant_key, product_name, dir_path))
            qm_collection.insert_one({
                                'reaction': [reactant_key, product_name], 
                                'Reactant SMILES':reac_obj.write('can').split()[0], 
                                'reactant_inchi_key':reactant_key, 
                                'product_inchi_key':product_name, 
                                'Product SMILES':prod_mol.write('can').split()[0],
                                'reactant_mopac_hf':H298_reac,
                                'product_mopac_hf':H298_prod,
                                'path':dir_path, 
                                'generations':self.generations
                                })
            return 1
        self.logger.info('Delta H is {}, greater than threshold'.format(dH))
        self.logger.info('Finished {}/{}\n'.format(self.count, total_prod_num))
        return 0

    def get_mopac_H298(self, mol_object, charge = 0, multiplicity = 'SINGLET'):
        tmpdir = os.path.join(os.path.dirname(os.getcwd()), 'tmp')
        reactant_path = os.path.join(tmpdir, 'reactant.mop')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)

        reac_geo = str(mol_object.toNode()).splitlines()
        reactant_geometry = []
        for i in reac_geo:
            i_list = i.split()
            atom = i_list[0] + " "
            k = i_list[1:] + [""]
            l = " 0 ".join(k)
            out = atom + l
            reactant_geometry.append(out)
        reactant_geometry = "\n".join(reactant_geometry)

        with open(reactant_path, 'w') as f:
            f.write("NOSYM 1SCF CHARGE={} {} {}\n\n".format(charge, multiplicity, self.mopac_method))
            f.write("\n{}".format(reactant_geometry))
        mopac.runMopac(tmpdir, 'reactant.mop')
        mol_hf = mopac.getHeatofFormation(tmpdir, 'reactant.out')
        return float(mol_hf)

    def unique_key_filterIsomorphic_itself(self, base):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        base_unique = [mol.write('inchiKey').strip() for mol in base]
        #base_unique = [mol.toRMGMolecule().to_inchi_key() for mol in base]
        result = [base[base_unique.index(i)] for i in set(base_unique)]
        return result

    def gen_geometry(self, reactant_mol, product_mol, add_bonds, break_bonds, **kwargs):
        # Database
        qm_collection = db['qm_calculate_center']

        reactant_mol_copy = reactant_mol.copy()
        # Initial optimization
        if self.constraint == None:
            Hatom = gen3D.readstring('smi', '[H]')
            ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
        else:
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)

        # Arrange
        # If arrange error can use try
        arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, self.constraint)
        msg = arrange3D.arrangeIn3D()
        if msg != '':
            self.logger.info(msg)

        # After arrange to prevent openbabel use the previous product coordinates if it is isomorphic
        # to the current one, even if it has different atom indices participating in the bonds.
        if self.constraint == None:
            ff.Setup(Hatom.OBMol)
            reactant_mol.gen3D(make3D=False)
            ff.Setup(Hatom.OBMol)
            product_mol.gen3D(make3D=False)
            ff.Setup(Hatom.OBMol)
        else:
            reactant_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)
            product_mol.gen3D(self.constraint, forcefield=self.forcefield, method = self.constraintff_alg, make3D=False)

        reactant = reactant_mol.toNode()
        product = product_mol.toNode()
        self.logger.info('Reactant and product geometry is :\n{}\n****\n{}'.format(str(reactant), str(product)))
        subdir = os.path.join(os.path.dirname(self.ard_path), 'reactions')
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        b_dirname = product_mol.write('inchiKey').strip()
        targets = list(qm_collection.find({'product_inchi_key':b_dirname}))
        dirname = self.dir_check(subdir, b_dirname, len(targets) + 1)

        output_dir = util.makeOutputSubdirectory(subdir, dirname)
        kwargs['output_dir'] = output_dir
        #self.makeInputFile(reactant, product, **kwargs)
        self.makeDrawFile(reactant, 'reactant.xyz', **kwargs)
        self.makeDrawFile(product, 'product.xyz', **kwargs)
        self.makeisomerFile(add_bonds, break_bonds, **kwargs)

        reactant_mol.setCoordsFromMol(reactant_mol_copy)
        return output_dir

    def output(self, reactant_mol, product_mol, add_bonds, break_bonds, prod_mol, **kwargs):
        # Database
        qm_collection = db['qm_calculate_center']

        subdir = os.path.join(os.path.dirname(self.ard_path), 'reactions')
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        b_dirname = prod_mol.write('inchiKey').strip()
        targets = list(qm_collection.find({'product_inchi_key':b_dirname}))
        dirname = self.dir_check(subdir, b_dirname, len(targets) + 1)

        output_dir = util.makeOutputSubdirectory(subdir, dirname)
        kwargs['output_dir'] = output_dir
        #self.makeInputFile(reactant, product, **kwargs)
        self.makeDrawFile(reactant_mol, 'reactant.xyz', **kwargs)
        self.makeDrawFile(product_mol, 'product.xyz', **kwargs)
        self.makeisomerFile(add_bonds, break_bonds, **kwargs)
        return output_dir

    @staticmethod
    def dir_check(subdir, b_dirname, num):
        """
        When parallely run job, the dir is constructed but data is not on database yet
        """
        check = False
        number = num
        while check == False:
            new_name = '{}_{}'.format(b_dirname, number)
            if os.path.exists(os.path.join(subdir, new_name)):
                number += 1
            else:
                check = True
        return new_name

    @staticmethod
    def makeInputFile(reactant, product, **kwargs):
        """
        Create input file for TS search and return path to file.
        """
        path = os.path.join(kwargs['output_dir'], 'de_ssm_input.xyz')
        nreac_atoms = len(reactant.getListOfAtoms())
        nproduct_atoms = len(product.getListOfAtoms())

        with open(path, 'w') as f:
            f.write('{}\n\n{}\n{}\n\n{}\n'.format(nreac_atoms, reactant, nproduct_atoms, product))
        return path

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
        path = os.path.join(kwargs['output_dir'], 'add_bonds.txt')

        with open(path, 'w') as f:
            if len(add_bonds) != 0:
                for i in add_bonds:
                    f.write('ADD {} {}\n'.format(i[0]+1,i[1]+1))
            if len(break_bonds) != 0:
                for i in break_bonds:
                    f.write('BREAK {} {}\n'.format(i[0]+1,i[1]+1))

