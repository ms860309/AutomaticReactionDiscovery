# standard library imports
import os
import time
import psutil
import copy
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
from mopac import Mopac


#database
from connect import db

info = psutil.virtual_memory()

class Network(object):

    def __init__(self, reac_mol, reactant_graph, forcefield, logger, **kwargs):
        self.reac_mol = reac_mol
        self.reactant_graph = reactant_graph
        self.forcefield = forcefield
        self.logger = logger
        self.dh_cutoff = float(kwargs['dh_cutoff'])
        self.bond_dissociation_cutoff = kwargs['bond_dissociation_cutoff']
        self.ard_path = kwargs['ard_path']
        self.generations = kwargs['generations']
        self.method = kwargs["dh_cutoff_method"]
        self.constraint = kwargs['constraint_index']
        self.count = 0

    def genNetwork(self, mol_object, nbreak, nform):
        """
        Execute the automatic reaction discovery procedure.
        """
        # Database
        qm_collection = db['qm_calculate_center']
        pool_collection = db['pool']
        statistics_collection = db['statistics']
        targets = list(pool_collection.find({}))

        # Add all reactant to a list for pgen filter isomorphic
        inchi_key_list = [i['reactant_inchi_key'] for i in targets]

        # Generate all possible products
        gen = Generate(mol_object, inchi_key_list, self.reactant_graph, self.bond_dissociation_cutoff, self.constraint)
        self.logger.info('Generating all possible products...')
        gen.generateProducts(nbreak = int(nbreak), nform =  int(nform))
        prod_mols = gen.prod_mols
        self.logger.info('{} possible products generated\n'.format(len(prod_mols)))
        add_bonds = gen.add_bonds
        break_bonds = gen.break_bonds

        # Filter reactions based on standard heat of reaction  delta H
        if self.method == "mopac":
            self.logger.info('Now use {} to filter the delta H of reactions....\n'.format(self.method))
            prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_mopac(mol_object, mol, add_bonds[prod_mols.index(mol)], self.logger, len(prod_mols))]
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

        self.logger.info('After delta H filter {} product remain.\n'.format(len(prod_mols)))
        # Reactant information
        reactant_key = mol_object.write('inchiKey').strip()  # inchikey
        reactant_smi = mol_object.write('can').split()[0]    # smiles

        # Check isomorphic with products in database
        prod_mols_filtered = self.unique_key_filterIsomorphic(reactant_key, reactant_smi, prod_mols_filtered, add_bonds, break_bonds)

        # Generate geometry and insert to database
        statistics_collection.insert_one({
            'Reactant SMILES':mol_object.write('can').split()[0], 
            'reactant_inchi_key':reactant_key, 
            'add how many products':len(prod_mols_filtered)})

        self.logger.info('Generate geometry........\n')
        for mol in prod_mols_filtered:
            index = prod_mols.index(mol)
            # Generate geometry and return path
            dir_path = self.gen_geometry(mol_object, mol, add_bonds[index], break_bonds[index])
            product_name = mol.write('inchiKey').strip()
            self.logger.info('Reactant inchi key: {}\nProduct inchi key: {}\n Directory path: {}'.format(reactant_key, product_name, dir_path))
            qm_collection.insert_one({
                                   'reaction': [reactant_key, product_name], 
                                   'Reactant SMILES':mol_object.write('can').split()[0], 
                                   'reactant_inchi_key':reactant_key, 
                                   'product_inchi_key':product_name, 
                                   'Product SMILES':mol.write('can').split()[0], 
                                   'path':dir_path, 
                                   'ssm_status':'job_unrun', 
                                   'generations':self.generations
                                   }
                                  )

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
    
    def filter_dh_mopac(self, reac_obj, prod_mol, form_bonds, logger, total_prod_num):
        self.count += 1
        mopac_object = Mopac(self.forcefield, form_bonds, logger, total_prod_num, self.count, self.constraint)
        H298_reac, H298_prod = mopac_object.mopac_get_H298(reac_obj, prod_mol)

        if H298_prod == False and H298_reac == False:
            return 0

        self.logger.info('Product energy calculate by mopac is {} kcal/mol and reactant is {} kcal/mol'.format(H298_prod, H298_reac))
        dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            self.logger.info('Delta H is {}, smaller than threshold'.format(dH))
            self.logger.info('Finished {}/{}'.format(self.count, total_prod_num))
            return 1
        self.logger.info('Delta H is {}, greater than threshold'.format(dH))
        self.logger.info('Finished {}/{}\n'.format(self.count, total_prod_num))
        return 0

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

    def unique_key_filterIsomorphic(self, reactant_key, reactant_smi, compare, add_bonds, break_bonds):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        qm_collection = db['qm_calculate_center']
        reactions_collection = db['reactions']
        targets = list(qm_collection.find({'ts_status':'job_success'}))
        base_unique = [i['product_inchi_key'] for i in targets]
        compare_unique = [mol.write('inchiKey').strip() for mol in compare]
        isomorphic_idx =[]
        same = {}
        for idx, i in enumerate(compare_unique):
            if i not in base_unique:
                isomorphic_idx.append(idx)
            else:
                if i not in same.keys():
                    same[i] = [idx]
                else:
                    same[i] += [idx]

        result = [compare[i] for i in isomorphic_idx]
        
        same_unique_key = list(set(compare_unique) & set(base_unique))
        for idx, i in enumerate(same_unique_key):
            for j in same[i]:
                dc = ''
                for k in add_bonds[j]:
                    dc += 'ADD {} {}\n'.format(k[0]+1,k[1]+1)
                for l in break_bonds[j]:
                    dc += 'BREAK {} {}\n'.format(l[0]+1,l[1]+1)
        
                path_target = list(qm_collection.find({'product_inchi_key':i}))
                check_1 = [reactant_key, i]
                check_duplicate_1 = list(reactions_collection.find({'reaction':check_1}))
                
                if len(check_duplicate_1) == 0 and reactant_key != i:
                    reactions_collection.insert_one({
                                        'reaction':[reactant_key, i],
                                        'reactant_smi':reactant_smi,
                                        'product_smi':path_target[0]['Product SMILES'],
                                        'path':path_target[0]['path'],
                                        'generations':self.generations,
                                        'for_debug':'from same',
                                        'unique': 'waiting for check',
                                        'driving_coordinate':dc,
                                        'ssm':'job_unrun'})
                elif len(check_duplicate_1) > 0 and reactant_key != i:
                    # aleady have
                    reactions_collection.insert_one({
                                        'reaction':[reactant_key, i],
                                        'reactant_smi':reactant_smi,
                                        'product_smi':path_target[0]['Product SMILES'],
                                        'path':path_target[0]['path'],
                                        'generations':self.generations,
                                        'for_debug':'from same',
                                        'unique': 'waiting for check',
                                        'driving_coordinate':dc,
                                        'ssm':'job_unrun'})
                else:
                    # aleady have
                    reactions_collection.insert_one({
                                        'reaction':[reactant_key, i],
                                        'reactant_smi':reactant_smi,
                                        'product_smi':path_target[0]['Product SMILES'],
                                        'path':path_target[0]['path'],
                                        'generations':self.generations,
                                        'for_debug':'from same',
                                        'unique': 'reactant equal to product'})

        return result

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

        reactant_mol.separateMol()
        if len(reactant_mol.mols) > 1:
            reactant_mol.mergeMols()
        product_mol.separateMol()
        if len(product_mol.mols) > 1:
            product_mol.mergeMols()

        # Initial optimisation
        if self.constraint == None:
            Hatom = gen3D.readstring('smi', '[H]')
            ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
            ff.Setup(Hatom.OBMol)
            reactant_mol.gen3D(make3D=False)
            ff.Setup(Hatom.OBMol)
            product_mol.gen3D(make3D=False)
            ff.Setup(Hatom.OBMol)
        else:
            gen3D.constraint_force_field(reactant_mol.OBMol, self.constraint)
            gen3D.constraint_force_field(product_mol.OBMol, self.constraint)

        # Arrange
        reactant_mol_copy = reactant_mol.copy()
        arrange3D = gen3D.Arrange3D(reactant_mol, product_mol, self.constraint)
        msg = arrange3D.arrangeIn3D()
        if msg != '':
            self.logger.info(msg)

        # After arrange to prevent openbabel use the previous product coordinates if it is isomorphic
        # to the current one, even if it has different atom indices participating in the bonds.
        if self.constraint == None:
            Hatom = gen3D.readstring('smi', '[H]')
            ff = pybel.ob.OBForceField.FindForceField(self.forcefield)
            ff.Setup(Hatom.OBMol)
            reactant_mol.gen3D(make3D=False)
            ff.Setup(Hatom.OBMol)
            product_mol.gen3D(make3D=False)
            ff.Setup(Hatom.OBMol)
        else:
            gen3D.constraint_force_field(reactant_mol.OBMol, self.constraint)
            gen3D.constraint_force_field(product_mol.OBMol, self.constraint)

        reactant = reactant_mol.toNode()
        product = product_mol.toNode()
        self.logger.info('Reactant and product geometry is :\n{}****\n{}'.format(str(reactant), str(product)))
        subdir = os.path.join(os.path.dirname(self.ard_path), 'reactions')
        if not os.path.exists(subdir):
            os.mkdir(subdir)
        b_dirname = product_mol.write('inchiKey').strip()
        targets = list(qm_collection.find({'product_inchi_key':b_dirname}))
        dirname = self.dir_check(subdir, b_dirname, len(targets) + 1)

        output_dir = util.makeOutputSubdirectory(subdir, dirname)
        kwargs['output_dir'] = output_dir
        self.makeInputFile(reactant, product, **kwargs)
        self.makeCalEnergyFile(product, **kwargs)
        self.makeDrawFile(reactant, 'reactant.xyz', **kwargs)
        self.makeDrawFile(product, 'product.xyz', **kwargs)
        self.makeisomerFile(add_bonds, break_bonds, **kwargs)

        reactant_mol.setCoordsFromMol(reactant_mol_copy)
        return output_dir

    def dir_check(self, subdir, b_dirname, num):
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
    def makeCalEnergyFile(_input, spin = 0, multiplicity = 1, **kwargs):
        """
        Create input file for energy calculation.
        """
        path = os.path.join(kwargs['output_dir'], 'reactant_energy.in')
        with open(path, 'w') as f:
            f.write('$molecule\n{} {}\n{}\n$end\n'.format(spin, multiplicity, _input))
            f.write('$rem\n')
            f.write('JOBTYPE OPT\n')
            f.write('METHOD B97-D3\n')
            f.write('DFT_D D3_BJ\n')
            f.write('BASIS def2-mSVP\n')
            f.write('SCF_ALGORITHM DIIS\n')
            f.write('MAX_SCF_CYCLES 150\n')
            f.write('SCF_CONVERGENCE 8\n')
            f.write('SYM_IGNORE TRUE\n')
            f.write('SYMMETRY FALSE\n')
            f.write('GEOM_OPT_MAX_CYCLES 150\n')
            f.write('GEOM_OPT_TOL_GRADIENT 100\n')
            f.write('GEOM_OPT_TOL_DISPLACEMENT 400\n')
            f.write('GEOM_OPT_TOL_ENERGY 33\n')
            f.write('WAVEFUNCTION_ANALYSIS FALSE\n')
            f.write('$end\n@@@\n')
            f.write('$molecule\nread\n$end\n')
            f.write('JOBTYPE FREQ\n')
            f.write('METHOD B97-D3\n')
            f.write('DFT_D D3_BJ\n')
            f.write('BASIS def2-mSVP\n')
            f.write('SCF_GUESS READ\n')
            f.write('SCF_ALGORITHM DIIS\n')
            f.write('MAX_SCF_CYCLES 150\n')
            f.write('SCF_CONVERGENCE 8\n')
            f.write('SYM_IGNORE TRUE\n')
            f.write('SYMMETRY FALSE\n')
            f.write('WAVEFUNCTION_ANALYSIS FALSE\n')
            f.write('$end')

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

