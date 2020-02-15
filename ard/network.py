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

class Network(object):

    def __init__(self, reac_mol, gen, gen_2, forcefield, **kwargs):
        self.reac_mol = reac_mol
        self.gen = gen
        self.gen_2 = gen_2
        self.nbreak = int(kwargs['nbreak'])
        self.nform = int(kwargs['nform'])
        self.forcefield = forcefield
        self.dh_cutoff = float(kwargs['dh_cutoff'])
        self.output_dir = kwargs['output_dir']
        log_level = logging.INFO
        self.logger = util.initializeLog(log_level, os.path.join(self.output_dir, 'ARD.log'), logname='main')
        self.network_prod_mols = []
        self.pre_product = []
        self.nround = -1
        self.next_num = 0
        self.first_num = 0
        self.method = kwargs["dh_cutoff_method"]

    def genNetwork(self, reac_mol, _round = 1, **kwargs):
        """
        Execute the automatic reaction discovery procedure.
        """

        self.logger.info('Generating all possible products...')
        self.gen.generateProducts(nbreak=self.nbreak, nform=self.nform)
        prod_mols = self.gen.prod_mols
        self.logger.info('{} possible products generated\n'.format(len(prod_mols)))

        # Load thermo database and choose which libraries to search
        thermo_db = ThermoDatabase()
        thermo_db.load(os.path.join(settings['database.directory'], 'thermo'))
        thermo_db.libraryOrder = ['primaryThermoLibrary', 'NISTThermoLibrary', 'thermo_DFT_CCSDTF12_BAC',
                                    'CBS_QB3_1dHR', 'DFT_QCI_thermo', 'BurkeH2O2', 'GRI-Mech3.0-N', ]

        # Filter reactions based on standard heat of reaction
        if self.method == "mopac":
            H298_reactant = mopac(reac_mol, self.forcefield)
            H298_reac = H298_reactant.makeInputFile(reac_mol)
            self.logger.info('Filtering reactions...')
            prod_mols_filtered = [mol for mol in prod_mols if self.filter_dh_mopac(H298_reac, mol, **kwargs)]
        else:
            H298_reac = reac_mol.getH298(thermo_db)
            self.logger.info('Filtering reactions...')
            prod_mols_filtered = [mol for mol in prod_mols if self.filterThreshold(H298_reac, mol, thermo_db, **kwargs)]

        self.logger.info('{} products remaining\n'.format(len(prod_mols_filtered)))
        #check product isomorphic and filter them
        if self.nbreak == 3 and self.nform == 3 :
            self.gen_2.generateProducts(nbreak=2, nform=2)
            prod_mols_2 = self.gen_2.prod_mols
            self.logger.info('{} possible products generated with nbreak = 2 and nbreak = 2\n'.format(len(prod_mols_2)))
            self.logger.info('Filtering reactions...')
            #prod_mols_filtered_2 after filter by delta H
            prod_mols_filtered_2 = [mol for mol in prod_mols_2 if self.filterThreshold(H298_reac, mol, thermo_db, **kwargs)]
            self.logger.info('{} products remaining with nbreak = 2 and nbreak = 2 after filter by delta H\n'.format(len(prod_mols_filtered_2)))
            #prod_mols_filtered_2 after filter by isomorphic
            prod_mols_filtered_2 = self.filterIsomorphic_itself(prod_mols_filtered_2, prod_mols_filtered_2)
            self.logger.info('{} products remaining with nbreak = 2 and nbreak = 2 after isomorphic screening\n'.format(len(prod_mols_filtered_2)))
            prod_mols_filtered = self.filterIsomorphic(prod_mols_filtered, prod_mols_filtered_2)
            self.logger.info('{} products remaining after isomorphic screening by nbreak=2 and nform=2\n'.format(len(prod_mols_filtered)))
            prod_mols_filtered = self.filterIsomorphic_itself(prod_mols_filtered, prod_mols_filtered)
            self.logger.info('{} products remaining with nbreak = 3 and nbreak = 3 after isomorphic screening\n'.format(len(prod_mols_filtered)))
            prod_mols_filtered += prod_mols_filtered_2 
            self.logger.info('{} total products remaining(break=2, form=2 + break=3, form=3)\n'.format(len(prod_mols_filtered)))
        else:
            prod_mols_filtered = self.filterIsomorphic_itself(prod_mols_filtered, prod_mols_filtered)
            self.logger.info('{} products remaining after isomorphic screening\n'.format(len(prod_mols_filtered)))

        # append product_mol to network
        pre_products = []
        # initial round add all prod to self.network
        self.nround += _round
        if self.nround == 0:
            print("Here is the next generation")
            print("There are {} rounds need to generate possible products at next generation".format(len(prod_mols_filtered)))
            self.first_num = len(prod_mols_filtered)
            for mol in prod_mols_filtered:
                mol.gen3D(forcefield=self.forcefield, make3D=False)
                self.network_prod_mols.append(mol)
            self.recurrently_gen(self.network_prod_mols, 0)
        else:
            for mol in prod_mols_filtered:
                mol.gen3D(forcefield=self.forcefield, make3D=False)
                pre_products.append(mol)
            filtered = []
            filtered = self.filterIsomorphic(self.network_prod_mols, pre_products)
            det = len(filtered)
            self.next_num += det
            for mol in filtered:
                self.pre_product.append(mol)
                self.network_prod_mols.append(mol)
            print("Add {} new products into network".format(det))
            print("Now network have {} products".format(len(self.network_prod_mols)))

    def recurrently_gen (self, prod_mols_filtered, num):
        """
        For generate recyclely
        """
        _count = 0
        self.pre_product = []
        for mol in prod_mols_filtered[num:]:
            self.genNetwork(mol)    
            _count += 1
            print("*Finished {} rounds".format(_count))
            self.checker(self.nround, self.first_num)
    
    def checker(self, count, check):
        if count == check:
            if len(self.pre_product) !=0:
                print("Here is the next generation")
                print("There are {} rounds need to generate possible products at first generation".format(self.next_num))
                parse = len(self.pre_product)
                self.recurrently_gen(self.network_prod_mols, parse)
            elif len(self.pre_product) == 0:
                print("starting generate geometry")
                self.gen_geometry(self.network_prod_mols)             
        elif count == check + self.next_num:
            if len(self.pre_product) !=0:
                print("Here is the next generation")
                print("There are {} rounds need to generate possible products at this generation".format(self.next_num))
                parse = len(self.pre_product)
                self.recurrently_gen(self.network_prod_mols, parse)
            elif len(self.pre_product) == 0:
                print("starting generate geometry")
                self.gen_geometry(self.network_prod_mols) 
        else:
            pass
            

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
        H298_prod = H298_product.makeInputFile(prod_mol)
        dH = H298_prod - H298_reac

        if dH < self.dh_cutoff:
            return True
        return False

    def filterIsomorphic(self, base, compare):

        isomorphic_idx = []
        for i in compare:
            prod_rmg_mol = i.toRMGMolecule()
            for idx_2, j in enumerate(base):
                compare = j.toRMGMolecule()
                if compare.isIsomorphic(prod_rmg_mol):
                    isomorphic_idx.append(idx_2)
        isomorphic = set(isomorphic_idx)
        prod_mols_filtered_idx = set([i for i in range(len(base))])
        index = list(prod_mols_filtered_idx - isomorphic)
        prod_mols_filtered = [base[i] for i in index]
        return prod_mols_filtered

    def filterIsomorphic_itself(self, base, compare):

        isomorphic_idx = []
        for idx_1, i in enumerate(compare):
            prod_rmg_mol = i.toRMGMolecule()
            for idx_2, j in enumerate(base[idx_1+1:]):
                compare = j.toRMGMolecule()
                if compare.isIsomorphic(prod_rmg_mol):
                    isomorphic_idx.append(idx_1 + 1 + idx_2)
        isomorphic = set(isomorphic_idx)
        prod_mols_filtered_idx = set([i for i in range(len(base))])
        index = list(prod_mols_filtered_idx - isomorphic)
        prod_mols_filtered = [base[i] for i in index]
        return prod_mols_filtered

    def gen_geometry(self, network_prod_mols, **kwargs):

        start_time = time.time()
        # These two lines are required so that new coordinates are
        # generated for each new product. Otherwise, Open Babel tries to
        # use the coordinates of the previous molecule if it is isomorphic
        # to the current one, even if it has different atom indices
        # participating in the bonds. a hydrogen atom is chosen
        # arbitrarily, since it will never be the same as any of the
        # product structures.
        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField(self.forcefield)

        reac_mol_copy = self.reac_mol.copy()
        # Generate 3D geometries
        if network_prod_mols:
            self.logger.info('Feasible products:\n')
            rxn_dir = util.makeOutputSubdirectory(self.output_dir, 'reactions')
            for rxn, mol in enumerate(network_prod_mols):

                arrange3D = gen3D.Arrange3D(self.reac_mol, mol)
                msg = arrange3D.arrangeIn3D()
                if msg != '':
                    self.logger.info(msg)

                ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
                self.reac_mol.gen3D(make3D=False)
                ff.Setup(Hatom.OBMol)
                mol.gen3D(make3D=False)
                ff.Setup(Hatom.OBMol)

                reactant = self.reac_mol.toNode()
                product = mol.toNode()

                rxn_num = '{:04d}'.format(rxn)
                output_dir = util.makeOutputSubdirectory(rxn_dir, rxn_num)
                kwargs['output_dir'] = output_dir
                kwargs['name'] = rxn_num

                self.logger.info('Product {}: {}\n{}\n****\n{}\n'.format(rxn, product.toSMILES(), reactant, product))
                self.makeInputFile(reactant, product, **kwargs)

                self.reac_mol.setCoordsFromMol(reac_mol_copy)
        else:
            self.logger.info('No feasible products found')
        # Finalize
        self.finalize(start_time)


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