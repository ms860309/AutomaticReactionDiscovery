import logging
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
from mopac import mopac
import copy
import shutil

class Network(object):

    def __init__(self, reac_mol, forcefield, **kwargs):
        self.reac_mol = reac_mol
        self.nbreak = int(kwargs['nbreak'])
        self.nform = int(kwargs['nform'])
        self.forcefield = forcefield
        self.dh_cutoff = float(kwargs['dh_cutoff'])
        self.output_dir = kwargs['output_dir']
        log_level = logging.INFO
        self.logger = util.initializeLog(log_level, os.path.join(self.output_dir, 'ARD.log'), logname='main')
        self.network_log = util.initializeLog(log_level, os.path.join(self.output_dir, 'NETWORK.log'), logname='sec')
        self.reactions = {}
        self.network_prod_mols = []
        self.add_bonds = []
        self.pre_product = []
        self.reactant_list = []
        self.nround = -1
        self.first_num = 0
        self.tmp = 0
        self._count = 0
        self.generation = 1
        self.rxn_num = -1
        self.method = kwargs["dh_cutoff_method"]

    def genNetwork(self, mol_object, _round = 1, **kwargs):
        """
        Execute the automatic reaction discovery procedure.
        """
        self.logger.info('Generating all possible products...')
        #Add all reactant to a list for pgen filter isomorphic
        if mol_object not in self.reactant_list:
            self.reactant_list.append(mol_object)
        gen = Generate(mol_object, self.reactant_list)
        gen_2 = Generate(mol_object, self.reactant_list)
        gen.generateProducts(nbreak=self.nbreak, nform=self.nform)
        prod_mols = gen.prod_mols
        add_bonds = gen.add_bonds
        self.logger.info('{} possible products generated\n'.format(len(prod_mols)))

        # Load thermo database and choose which libraries to search
        thermo_db = ThermoDatabase()
        thermo_db.load(os.path.join(settings['database.directory'], 'thermo'))
        thermo_db.libraryOrder = ['primaryThermoLibrary', 'NISTThermoLibrary', 'thermo_DFT_CCSDTF12_BAC',
                                    'CBS_QB3_1dHR', 'DFT_QCI_thermo', 'BurkeH2O2', 'GRI-Mech3.0-N', ]
        # Filter reactions based on standard heat of reaction
        if self.method == "mopac":
            H298_reactant = mopac(mol_object, self.forcefield)
            H298_reac = H298_reactant.mopac_get_H298(mol_object)
            self.logger.info('Filtering reactions...')
            prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_mopac(H298_reac, mol, **kwargs)]
        else:
            H298_reac = mol_object.getH298(thermo_db)
            self.logger.info('Filtering reactions...')
            prod_mols_filtered = [mol for mol in prod_mols if self.filterThreshold(H298_reac, mol, thermo_db, **kwargs)]

        self.logger.info('{} products remaining\n'.format(len(prod_mols_filtered)))
        #check product isomorphic and filter them
        if self.nbreak == 3 and self.nform == 3:
            gen_2.generateProducts(nbreak=2, nform=2)
            prod_mols_2 = gen_2.prod_mols
            self.logger.info('{} possible products generated with nbreak = 2 and nbreak = 2\n'.format(len(prod_mols_2)))
            self.logger.info('Filtering reactions...')
            #prod_mols_filtered_2 after filter by delta H
            if self.method == "mopac":
                prod_mols_filtered_2 = [mol for mol in prod_mols if self.filter_dh_mopac(H298_reac, mol, **kwargs)]
            else:
                prod_mols_filtered_2 = [mol for mol in prod_mols_2 if self.filterThreshold(H298_reac, mol, thermo_db, **kwargs)]
            self.logger.info('{} products remaining with nbreak = 2 and nbreak = 2 after filter by delta H\n'.format(len(prod_mols_filtered_2)))
            #prod_mols_filtered_2 after filter by isomorphic

            #prod_mols_filtered_2 = self.filterIsomorphic_itself(prod_mols_filtered_2, prod_mols_filtered_2)
            prod_mols_filtered_2 = self.unique_key_filterIsomorphic_itself(prod_mols_filtered_2, prod_mols_filtered_2)
            self.logger.info('{} products remaining with nbreak = 2 and nbreak = 2 after isomorphic screening\n'.format(len(prod_mols_filtered_2)))
            #prod_mols_filtered = self.filterIsomorphic(prod_mols_filtered, prod_mols_filtered_2)#
            prod_mols_filtered = self.unique_key_filterIsomorphic(prod_mols_filtered, prod_mols_filtered_2)
            self.logger.info('{} products remaining after isomorphic screening by nbreak=2 and nform=2\n'.format(len(prod_mols_filtered)))
            #prod_mols_filtered = self.filterIsomorphic_itself(prod_mols_filtered, prod_mols_filtered)
            prod_mols_filtered = self.unique_key_filterIsomorphic_itself(prod_mols_filtered, prod_mols_filtered)
            self.logger.info('{} products remaining with nbreak = 3 and nbreak = 3 after isomorphic screening\n'.format(len(prod_mols_filtered)))
            prod_mols_filtered += prod_mols_filtered_2 
            self.logger.info('{} total products remaining(break=2, form=2 + break=3, form=3)\n'.format(len(prod_mols_filtered)))
        else:
            #prod_mols_filtered = self.filterIsomorphic_itself(prod_mols_filtered, prod_mols_filtered)
            prod_mols_filtered = self.unique_key_filterIsomorphic_itself(prod_mols_filtered, prod_mols_filtered)
            self.logger.info('{} products remaining after isomorphic screening\n'.format(len(prod_mols_filtered)))

        # append product_mol to network
        pre_products = []
        # initial round add all prod to self.network
        self.nround += _round
        if self.nround == 0:
            self.network_log.info("Here is the {} generation\n".format(self.generation))
            self.network_log.info("There are {} rounds need to generate possible products at next generation\n".format(len(prod_mols_filtered)))
            self.first_num = len(prod_mols_filtered)
            self.network_log.info("starting generate geometry\n")
            num = 0
            for idx, mol in enumerate(prod_mols_filtered):
                index = prod_mols.index(mol)
                self.add_bonds.append(add_bonds[index])
                num += 1
                self.network_prod_mols.append(mol)
                self.gen_geometry(mol_object, mol, add_bonds[index])
                rxn_idx = 'reaction{}'.format(idx)
                rxn_num = '{:05d}'.format(self.rxn_num)
                self.reactions[rxn_idx] = ['00000', rxn_num]
            self.recurrently_gen(self.network_prod_mols, 0)
        else:
            for mol in prod_mols_filtered:
                pre_products.append(mol)
            filtered = []
            filtered = self.unique_key_filterIsomorphic(self.network_prod_mols, pre_products)
            det = len(filtered)
            for mol in filtered:
                self.pre_product.append(mol)
                self.network_prod_mols.append(mol)
            self.network_log.info("Add {} new products into network\n".format(det))
            self.network_log.info("Now network have {} products\n".format(len(self.network_prod_mols)))
            if det != 0:
                self.network_log.info("starting generate geometry\n")
                index = self.network_prod_mols.index(mol_object)
                rxn_idx = 'reaction{}'.format(index)
                tmp_list = self.reactions[rxn_idx]
                num = 0                
                for idx, prod_mol in enumerate(filtered):
                    index = prod_mols.index(mol)
                    self.add_bonds.append(add_bonds[index])
                    num += 1
                    a = copy.deepcopy(tmp_list)
                    rxn_idx = 'reaction{}'.format(self.rxn_num)
                    rxn_num = '{:05d}'.format(self.rxn_num + 1)
                    a.append(rxn_num)
                    self.reactions[rxn_idx] = a
                    self.gen_geometry(mol_object, mol, add_bonds[index])
    
    def recurrently_gen (self, prod_mols_filtered, num):
        """
        For generate recyclely
        First gen geometry and then recurrently gen
        """
        self.pre_product = []
        for mol in prod_mols_filtered[num:]:
            self.genNetwork(mol)    
            self._count += 1
            self.network_log.info("Finished {}/{} rounds at {} generations\n".format(self._count, self.first_num+self.tmp, self.generation))
            self.checker(self.nround, self.first_num)
    
    def checker(self, count, check):
        if count == check:
            if len(self.pre_product) !=0:
                self.generation += 1
                self.network_log.info("Here is the {} generation\n".format(self.generation))
                self.network_log.info("There are {} rounds need to generate possible products at this generation\n".format(len(self.pre_product)))
                self.tmp = len(self.pre_product)
                self.recurrently_gen(self.network_prod_mols, len(self.network_prod_mols) - self.tmp) 
            elif len(self.pre_product) == 0:
                self.network_log.info("Network is finished.\nself.reactions = {}".format(self.reactions))
        elif count == check + self.tmp:
            if len(self.pre_product) !=0:
                self.generation += 1
                self.network_log.info("Here is the {} generation\n".format(self.generation))
                self.network_log.info("There are {} rounds need to generate possible products at this generation\n".format(len(self.pre_product)))
                self.tmp += len(self.pre_product)
                self.recurrently_gen(self.network_prod_mols, len(self.network_prod_mols) - len(self.pre_product)) 
            elif len(self.pre_product) == 0:
                self.network_log.info("Network is finished.\nself.reactions = {}".format(self.reactions))


    def finalize(self, start_time):
        """
        Finalize the job.
        """
        self.logger.info('\nARD terminated on ' + time.asctime())
        self.logger.info('Total ARD run time: {:.1f} s'.format(time.time() - start_time))

    def filterThreshold(self, H298_reac, prod_mol, thermo_db, **kwargs):
        """
        Filter threshold based on standard enthalpies of formation of reactants
        and products. Returns `True` if the heat of reaction is less than
        `self.dh_cutoff`, `False` otherwise.
        """
        H298_prod = prod_mol.getH298(thermo_db)
        dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            return True
        return False
    
    def filter_dh_mopac(self, H298_reac, prod_mol, **kwargs):
        H298_product = mopac(prod_mol, self.forcefield)
        H298_prod = H298_product.mopac_get_H298(prod_mol)
        dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            return True
        return False

    def filterIsomorphic(self, base, compare):

        isomorphic_idx = []
        for idx_1, i in enumerate(compare):
            prod_rmg_mol = i.toRMGMolecule()
            for j in base:
                tmp = j.toRMGMolecule()
                if prod_rmg_mol.isIsomorphic(tmp):
                    isomorphic_idx.append(idx_1)
        isomorphic = set(isomorphic_idx)
        prod_mols_filtered_idx = set([i for i in range(len(compare))])
        index = list(prod_mols_filtered_idx - isomorphic)
        result = [compare[i] for i in index]
        return result

    def filterIsomorphic_itself(self, base, compare):
        isomorphic_idx = []
        for idx_1, i in enumerate(compare):
            prod_rmg_mol = i.toRMGMolecule()
            for j in base[idx_1+1:]:
                tmp = j.toRMGMolecule()
                if prod_rmg_mol.isIsomorphic(tmp):
                    isomorphic_idx.append(idx_1)
        isomorphic = set(isomorphic_idx)
        result = [compare[i] for i in isomorphic]
        return result

    def unique_key_filterIsomorphic(self, base, compare):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        isomorphic_idx = []
        base_unique = []
        compare_unique = []
        result = []
        for i in base:
            mol = i.toRMGMolecule()
            #unique = mol.toAugmentedInChIKey()
            unique = mol.toInChIKey()
            base_unique.append(unique)
        for i in compare:
            mol = i.toRMGMolecule()
            #unique = mol.toAugmentedInChIKey()
            unique = mol.toInChIKey()
            compare_unique.append(unique)
        for idx_1, i in enumerate(compare_unique):
            if i not in base_unique:
                isomorphic_idx.append(idx_1)
        for i in isomorphic_idx:
            result.append(compare[i])
        return result

    def unique_key_filterIsomorphic_itself(self, base, compare):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        isomorphic_idx = []
        base_unique = []
        compare_unique = []
        result = []
        for i in base:
            mol = i.toRMGMolecule()
            #unique = mol.toAugmentedInChIKey()
            unique = mol.toInChIKey()
            base_unique.append(unique)
        for i in compare:
            mol = i.toRMGMolecule()
            #unique = mol.toAugmentedInChIKey()
            unique = mol.toInChIKey()
            compare_unique.append(unique)
        for idx_1, i in enumerate(compare_unique):
            if i not in base_unique[idx_1+1:]:
                isomorphic_idx.append(idx_1)
        for i in isomorphic_idx:
            result.append(compare[i])
        return result


    def gen_geometry(self, reactant_mol, network_prod_mol, add_bonds, **kwargs):
        start_time = time.time()
        self.rxn_num += 1
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

        arrange3D = gen3D.Arrange3D(reactant_mol, network_prod_mol)
        msg = arrange3D.arrangeIn3D()
        if msg != '':
            self.logger.info(msg)

        ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
        reactant_mol.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)
        network_prod_mol.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)

        reactant = reactant_mol.toNode()
        product = network_prod_mol.toNode()
        subdir = os.path.join(os.path.abspath(os.pardir), 'reactions')
        if self.rxn_num == 0:
            if os.path.exists(subdir):
                shutil.rmtree(subdir)
            rxn_dir = os.mkdir(subdir)
            rxn_num = '{:05d}'.format(self.rxn_num)
            output_dir = util.makeOutputSubdirectory(subdir, rxn_num)
            kwargs['output_dir'] = output_dir
            kwargs['name'] = rxn_num
            self.makeCalEnergyFile(reactant, **kwargs)
            self.makeDrawFile(reactant, filename = 'reactant.xyz', **kwargs)
            self.rxn_num += 1
            rxn_num = '{:05d}'.format(self.rxn_num)
            output_dir = util.makeOutputSubdirectory(subdir, rxn_num)
            kwargs['output_dir'] = output_dir
            kwargs['name'] = rxn_num
            self.logger.info('Product {}: {}\n{}\n****\n{}\n'.format(self.rxn_num, product.toSMILES(), reactant, product))
            self.makeInputFile(reactant, product, **kwargs)
            self.makeCalEnergyFile(product, **kwargs)
            self.makeDrawFile(reactant, filename = 'reactant.xyz', **kwargs)
            self.makeDrawFile(product, filename = 'product.xyz', **kwargs)
            self.makeAdd_bondsFile(add_bonds, **kwargs)
            self.finalize(start_time)
        else:
            rxn_num = '{:05d}'.format(self.rxn_num)
            output_dir = util.makeOutputSubdirectory(subdir, rxn_num)
            kwargs['output_dir'] = output_dir
            kwargs['name'] = rxn_num
            self.logger.info('Product {}: {}\n{}\n****\n{}\n'.format(self.rxn_num, product.toSMILES(), reactant, product))
            self.makeInputFile(reactant, product, **kwargs)
            self.makeCalEnergyFile(product, **kwargs)
            self.makeDrawFile(reactant, 'reactant.xyz', **kwargs)
            self.makeDrawFile(product, 'product.xyz', **kwargs)
            self.makeAdd_bondsFile(add_bonds, **kwargs)
            self.finalize(start_time)


    def gen_network_input(network):
        """
        Generate network input file
        """
        path = os.path.join(kwargs['output_dir'], 'network_input.in')
        with open(path, 'w') as f:
            f.write('{}'.format(reactions))

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
    def makeAdd_bondsFile(_input, **kwargs):
        """
        Create input file(add which bonds) for Single ended String Method (SSM) calculation.
        """
        try:
            path = os.path.join(kwargs['output_dir'], 'add_bonds.txt')
            first_bond = _input[0]
            second_bond = _input[1]

            with open(path, 'w') as f:
                f.write('ADD {} {}\nADD {} {}'.format(first_bond[0]+1, first_bond[1]+1, second_bond[0]+1, second_bond[1]+1))
        except:
            print("maybe list out of range, check add bond")

    def logHeader(self):
        """
        Output a log file header.
        """
        self.logger.info('######################################################################')
        self.logger.info('#################### AUTOMATIC REACTION DISCOVERY ####################')
        self.logger.info('######################################################################')
        self.logger.info('Reactant SMILES: ' + self.reac_mol)
        self.logger.info('Maximum number of bonds to be broken: ' + str(self.nbreak))
        self.logger.info('Maximum number of bonds to be formed: ' + str(self.nform))
        self.logger.info('Heat of reaction cutoff: {:.1f} kcal/mol'.format(self.dh_cutoff))
        self.logger.info('Force field for 3D structure generation: ' + self.forcefield)
        self.logger.info('######################################################################\n')
