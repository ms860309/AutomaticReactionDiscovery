import logging
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
        log_level = logging.INFO
        self.logger = util.initializeLog(log_level, os.path.join(self.output_dir, 'ARD.log'), logname='main')
        self.network_log = util.initializeLog(log_level, os.path.join(self.output_dir, 'NETWORK.log'), logname='sec')
        self.reactions = {}
        self.network_prod_mols = []
        self.add_bonds = []
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
        start_time = time.time()
        self.logHeader(mol_object)
        #Add all reactant to a list for pgen filter isomorphic
        if mol_object not in self.reactant_list:
            self.reactant_list.append(mol_object)
        gen = Generate(mol_object, self.reactant_list)
        self.logger.info('Generating all possible products...')
        gen.generateProducts(nbreak=self.nbreak, nform=self.nform)
        prod_mols = gen.prod_mols
        self.logger.info('{} possible products generated\n'.format(len(prod_mols)))
        self.reactant_list += prod_mols
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
            prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_mopac(H298_reac, mol, **kwargs)]
        else:
            H298_reac = self.reac_mol.getH298(thermo_db)
            self.logger.info('Filtering reactions..., dH cutoff')
            prod_mols_filtered = [mol for mol in prod_mols if self.filterThreshold(H298_reac, mol, thermo_db, **kwargs)]
            self.logger.info('{} products remaining\n'.format(len(prod_mols_filtered)))

        #check product isomorphic and filter them
        if self.nbreak == 3 and self.nform == 3:
            gen_2.generateProducts(nbreak=2, nform=2)
            prod_mols_2 = gen_2.prod_mols
            #prod_mols_filtered_2 after filter by delta H
            if self.method == "mopac":
                prod_mols_filtered_2 = [mol for mol in prod_mols if self.filter_dh_mopac(H298_reac, mol, **kwargs)]
            else:
                prod_mols_filtered_2 = [mol for mol in prod_mols_2 if self.filterThreshold(H298_reac, mol, thermo_db, **kwargs)]
            #prod_mols_filtered_2 after filter by isomorphic
            prod_mols_filtered_2 = self.unique_key_filterIsomorphic_itself(prod_mols_filtered_2)
            prod_mols_filtered = self.unique_key_filterIsomorphic(prod_mols_filtered, prod_mols_filtered_2)
            prod_mols_filtered = self.unique_key_filterIsomorphic_itself(prod_mols_filtered)
            prod_mols_filtered += prod_mols_filtered_2 
        else:
            self.logger.info('Filtering reactions..., isomorphic')
            prod_mols_filtered = self.unique_key_filterIsomorphic_itself(prod_mols_filtered)
            self.logger.info('{} products remaining\n'.format(len(prod_mols_filtered)))

        # append product_mol to network
        pre_products = []
        # initial round add all prod to self.network
        self.nround += _round
        if self.nround == 0:
            self.logger.info('Feasible products:\n')
            self.network_logHeader(mol_object)
            self.network_log.info("-----------------------------------------\n")
            self.network_log.info("\nGeneration : {}\n".format(self.generation))
            self.network_log.info("-----------------------------------------\n")
            self.network_log.info("Need {} Rounds\n".format(len(prod_mols_filtered)))
            self.network_log.info("-----------------------------------------\n")
            self.first_num = len(prod_mols_filtered)
            collection = db['reactions']
            for idx, mol in enumerate(prod_mols_filtered):
                data = {}
                index = prod_mols.index(mol)
                self.network_prod_mols.append(mol)
                self.gen_geometry(mol_object, mol, add_bonds[index], break_bonds[index])
                self.logger.info('Product SMILES {} : {}'.format(idx+1, mol.write('can').strip()))
                rxn_idx = 'reaction{}'.format(idx)
                rxn_num = '{:05d}'.format(self.rxn_num)
                self.reactions[rxn_idx] = ['00000', rxn_num]
                data[rxn_idx] = ['00000', rxn_num]
                data['for_ssm_check'] = rxn_num
                data['Reactant SMILES'] = mol_object.write('can').strip()
                data['Product SMILES'] = mol.write('can').strip()
                collection.insert_one(data)
            self.finalize(start_time)
            self.logger.info('\n\nStart nerwork exploring......')
            self.first_generation(self.network_prod_mols)
        else:
            self.logger.info('Feasible products:\n')
            for mol in prod_mols_filtered:
                pre_products.append(mol)
            filtered = []
            same = []
            filtered, same, same_key = self.unique_key_filterIsomorphic(self.network_prod_mols, pre_products)
            det = len(filtered)
            for mol in filtered:
                self.pre_product.append(mol)
                self.network_prod_mols.append(mol)
            self.network_log.info('Product SMILES : {} ---> Generate {} possible products.'.format(mol_object.write('can').strip(), det))
            self.network_log.info("Total products : {}".format(len(self.network_prod_mols)))
            if det != 0:
                collection = db['reactions']
                index = self.network_prod_mols.index(mol_object)
                rxn_idx = 'reaction{}'.format(index)
                tmp_list = self.reactions[rxn_idx]          
                for j, mol in enumerate(filtered):
                    data = {}
                    idx = prod_mols.index(mol)
                    a = copy.deepcopy(tmp_list)
                    rxn_idx = 'reaction{}'.format(self.rxn_num)
                    rxn_num = '{:05d}'.format(self.rxn_num + 1)
                    a.append(rxn_num)
                    self.reactions[rxn_idx] = a
                    self.gen_geometry(mol_object, mol, add_bonds[idx], break_bonds[idx])
                    self.logger.info('Product SMILES {} : {}'.format(j+1, mol.write('can').strip()))
                    data[rxn_idx] = a
                    data['for_ssm_check'] = rxn_num
                    data['Reactant SMILES'] = mol_object.write('can').strip()
                    data['Product SMILES'] = mol.write('can').strip()
                    collection.insert_one(data)
                self.finalize(start_time)
                self.logger.info('\n\nStart nerwork exploring......')
                for key, mol in enumerate(same):
                    data = {}
                    idx = prod_mols.index(mol)
                    self.add_bonds.append(add_bonds[idx])
                    self.break_bonds.append(break_bonds[idx])     
                    a = copy.deepcopy(tmp_list)
                    rxn_idx = 'reaction{}'.format(self.rxn_num+key)
                    rxn_num = '{:05d}'.format(same_key[key]+1)
                    a.append(rxn_num)
                    self.reactions[rxn_idx] = a 
                    data[rxn_idx] = a
                    data['for_ssm_check'] = rxn_num
                    data['Reactant SMILES'] = mol_object.write('can').strip()
                    data['Product SMILES'] = mol.write('can').strip()
                    collection.insert_one(data)
            """
            else:
                # filter = 0
                index = self.network_prod_mols.index(mol_object)
                rxn_idx = 'reaction{}'.format(index)
                tmp_list = self.reactions[rxn_idx]
                for key, same_prod in enumerate(same):
                    idx = prod_mols.index(same_prod)
                    self.add_bonds.append(add_bonds[idx])
                    self.break_bonds.append(break_bonds[idx])     
                    a = copy.deepcopy(tmp_list)
                    rxn_idx = 'reaction{}'.format(self.rxn_num)
                    rxn_num = '{:05d}'.format(same_key[key]+1)
                    a.append(rxn_num)
                    self.reactions[rxn_idx] = a
                    self.rxn_num += 1 
            """

    def recurrently_gen (self, prod_mols_filtered, num):
        """
        For generate recyclely
        First gen geometry and then recurrently gen
        """
        self.genNetwork(prod_mols_filtered[num])
        self._count += 1
        self.network_log.info("Finished {}/{} rounds at {} generations\n".format(self._count, self.first_num+self.tmp, self.generation))
        self.network_log.info("-----------------------------------------\n")
        self.checker(self.nround, self.first_num)
    
    def checker(self, count, check):
        if count == check:
            if len(self.pre_product) !=0:
                self.recurrently_checker(count, check)
            elif len(self.pre_product) == 0:
                self.network_log.info("Network is finished.\nself.reactions = {}".format(self.reactions))
        elif count == check + self.tmp:
            if len(self.pre_product) !=0:
                self.recurrently_checker(count, check)
            elif len(self.pre_product) == 0:
                self.network_log.info("Network is finished.\nself.reactions = {}".format(self.reactions))

    def recurrently_checker(self, count, check):
        """
        Use database to check number,
        if success + fail == target number:
            recurrently_gen
        else:
            wait 5min (300s)
            back to checker
        """
        collect = db['moleculues']
        query = {"ssm_status":
            {"$in":
                ["job_success", "job_fail"]
            }
        }
        targets = list(collect.find(query))
        if len(targets) == count:
            self.generation += 1
            self.network_log.info("-----------------------------------------\n")
            self.network_log.info("Generation : {}\n".format(self.generation))
            self.network_log.info("-----------------------------------------\n")
            self.network_log.info("Need {} Rounds\n".format(len(self.pre_product)))
            self.network_log.info("-----------------------------------------\n")
            self.tmp = len(self.pre_product)
            self.pre_product = []
            req = {'ssm_status':'job_success'}
            targets = list(collect.find(query))
            for target in targets:
                number = int(target['dir'])
                self.recurrently_gen(self.network_prod_mols, number)
        else:
            time.sleep(300)
            self.checker(self.nround, self.first_num)
    
    def first_generation(self, products):
        collect = db['moleculues']
        query = {"ssm_status":
            {"$in":
                ["job_success", "job_fail"]
            }
        }
        targets = list(collect.find(query))
        if len(targets) == len(products):
            req = {'ssm_status':'job_success'}
            targets = list(collect.find(query))
            for target in targets:
                number = int(target['dir'])
                self.recurrently_gen(self.network_prod_mols, number)
        else:
            time.sleep(300)
            self.first_generation(products)
        
    
    def filter_dh_mopac(self, H298_reac, prod_mol, **kwargs):
        H298_product = mopac(prod_mol, self.forcefield)
        H298_prod = H298_product.mopac_get_H298(prod_mol)
        dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            return True
        return False


    def unique_key_filterIsomorphic(self, base, compare):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        
        base_unique = [mol.toRMGMolecule().to_inchi_key() for mol in base]
        compare_unique = [mol.toRMGMolecule().to_inchi_key() for mol in compare]
        isomorphic_idx = [compare_unique.index(i) for i in set(compare_unique) - set(base_unique)]
        result = [compare[i] for i in isomorphic_idx]
        same_idx = [compare_unique.index(i) for i in set(compare_unique) & set(base_unique)]
        same_key = [base_unique.index(i) for i in set(compare_unique) & set(base_unique)]
        same = [compare[i] for i in same_idx]
        
        return result, same, same_key

    def unique_key_filterIsomorphic_itself(self, base):
        """
        Convert rmg molecule into inchi key(unique key) and check isomorphic
        """
        
        base_unique = [mol.toRMGMolecule().to_inchi_key() for mol in base]
        result = [base[base_unique.index(i)] for i in set(base_unique)]
        
        return result


    def gen_geometry(self, reactant_mol, network_prod_mol, add_bonds, break_bonds, **kwargs):
        start_time = time.time()
        self.rxn_num += 1
        #database
        collect = db['molecules']
        
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
                self.logger.info(msg)
        except:
            reactant_mol, network_prod_mol = reactant_mol_copy, network_prod_mol_copy
            
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
            os.mkdir(subdir)
            rxn_num = '{:05d}'.format(self.rxn_num)
            output_dir = util.makeOutputSubdirectory(subdir, rxn_num)
            kwargs['output_dir'] = output_dir
            kwargs['name'] = rxn_num
            self.makeCalEnergyFile(reactant, **kwargs)
            self.makeDrawFile(reactant, filename = 'reactant.xyz', **kwargs)

            collect.insert_one({'dir':rxn_num, 'path' : output_dir, "energy_status":"job_unrun", "ssm_status":"job_unrun"})

            self.rxn_num += 1
            rxn_num = '{:05d}'.format(self.rxn_num)
            output_dir = util.makeOutputSubdirectory(subdir, rxn_num)
            kwargs['output_dir'] = output_dir
            kwargs['name'] = rxn_num
            self.makeInputFile(reactant, product, **kwargs)
            self.makeCalEnergyFile(product, **kwargs)
            self.makeDrawFile(reactant, filename = 'reactant.xyz', **kwargs)
            self.makeDrawFile(product, filename = 'product.xyz', **kwargs)
            self.makeisomerFile(add_bonds, break_bonds, **kwargs)

            collect.insert_one({'dir':rxn_num, 'path' : output_dir, "energy_status":"job_unrun", "ssm_status":"job_unrun"})

        else:
            rxn_num = '{:05d}'.format(self.rxn_num)
            output_dir = util.makeOutputSubdirectory(subdir, rxn_num)
            kwargs['output_dir'] = output_dir
            kwargs['name'] = rxn_num
            self.makeInputFile(reactant, product, **kwargs)
            self.makeCalEnergyFile(product, **kwargs)
            self.makeDrawFile(reactant, 'reactant.xyz', **kwargs)
            self.makeDrawFile(product, 'product.xyz', **kwargs)
            self.makeisomerFile(add_bonds, break_bonds, **kwargs)

            collect.insert_one({'dir':rxn_num, 'path' : output_dir, "energy_status":"job_unrun", "ssm_status":"job_unrun", "ts_status":"job_unrun", "irc_status":"job_unrun"})

    def logHeader(self, mol_object):
        """
        Output a log file header.
        """
        self.logger.info('######################################################################')
        self.logger.info('#################### AUTOMATIC REACTION DISCOVERY ####################')
        self.logger.info('######################################################################')
        self.logger.info('Reactant SMILES: ' + mol_object.write('can').strip())
        self.logger.info('Maximum number of bonds to be broken: ' + str(self.nbreak))
        self.logger.info('Maximum number of bonds to be formed: ' + str(self.nform))
        self.logger.info('Heat of reaction cutoff: {:.1f} kcal/mol'.format(self.dh_cutoff))
        self.logger.info('######################################################################\n')

    def network_logHeader(self, mol_object):
        """
        Output a log file header.
        """
        self.network_log.info('######################################################################')
        self.network_log.info('#################### NETWORK EXPLORATION ####################')
        self.network_log.info('######################################################################')
        self.network_log.info('Reactant SMILES: ' + mol_object.write('can').strip())
        self.network_log.info('Maximum number of bonds to be broken: ' + str(self.nbreak))
        self.network_log.info('Maximum number of bonds to be formed: ' + str(self.nform))
        self.network_log.info('Heat of reaction cutoff: {:.1f} kcal/mol'.format(self.dh_cutoff))
        self.network_log.info('######################################################################\n')
        
    def finalize(self, start_time):
        """
        Finalize the job.
        """
        self.logger.info('\nTotal Memory : {}'.format(info.total))
        self.logger.info('Memory used : {}'.format(psutil.Process(os.getpid()).memory_info().rss))
        self.logger.info('Memory used percent : {}'.format(info.percent))
        self.logger.info('\nARD terminated on ' + time.asctime())
        self.logger.info('Total ARD run time: {:.1f} s'.format(time.time() - start_time))
        
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

