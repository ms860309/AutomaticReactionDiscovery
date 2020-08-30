# standard library imports
import os
import shutil
import time
import psutil

#third party
import logging
from subprocess import Popen, PIPE
import difflib
import pybel
import numpy as np

# local application imports
from node import Node
import gen3D
import util


class MopacError(Exception):
    """
    An exception class for errors that occur during mopac calculations.
    """
    pass

class Mopac(object):

    def __init__(self, forcefield, form_bonds):
        self.forcefield = forcefield
        self.form_bonds = form_bonds
        self.constraint = [0,1,2,3,4,5]  # for debug


        log_level = logging.INFO
        process = psutil.Process(os.getpid())
        #self.logger = util.initializeLog(log_level, os.path.join(os.getcwd(), 'ARD.log'), logname='main')
        #self.logger.info('\nARD initiated on ' + time.asctime() + '\n')
        #self.logger.info('memory usage: {}'.format(process.memory_percent()))

    def mopac_get_H298(self, reac_obj, prod_obj, charge = 0, multiplicity = 'SINGLET', method = 'PM7'):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """

        tmpdir = os.path.join(os.path.dirname(os.getcwd()), 'tmp')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        reactant_path = os.path.join(tmpdir, 'reactant.mop')
        product_path = os.path.join(tmpdir, 'product.mop')
        reac_geo, prod_geo = self.genInput(reac_obj, prod_obj)
        
        if reac_geo == False and prod_geo == False:
            return float(0.0), float(999999.0)
        else:
            with open(reactant_path, 'w') as f:
                f.write("CHARGE={} {} {}\n\n".format(charge, multiplicity, method))
                f.write("\n{}".format(reac_geo))
            start_time = time.time()
            self.runMopac(tmpdir, 'reactant.mop')
            reactant = self.getHeatofFormation(tmpdir, 'reactant.out')

            with open(product_path, 'w') as f:
                f.write("CHARGE={} {} {}\n\n".format(charge, multiplicity, method))
                f.write("\n{}".format(prod_geo))
            self.runMopac(tmpdir, 'product.mop')
            product = self.getHeatofFormation(tmpdir, 'product.out')
            self.finalize(start_time, 'mopac')
            
            """
            info = psutil.virtual_memory()
            print("cpu numbers : {}".format(psutil.cpu_count()))
            print("total memory : {}".format(info.total))
            print("memory used : {}".format(psutil.Process(os.getpid()).memory_info().rss))
            print("memory used percent : {}".format(info.percent))
            """
            return float(reactant), float(product)
    
    def genInput(self, reac_mol, prod_obj, threshold = 4.0):
        #start_time = time.time()

        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField(self.forcefield)

        reac_mol.separateMol()
        if len(reac_mol.mols) > 1:
            reac_mol.mergeMols()
        prod_obj.separateMol()
        if len(prod_obj.mols) > 1:
            prod_obj.mergeMols()

        reac_mol_copy = reac_mol.copy()
        arrange3D = gen3D.Arrange3D(reac_mol, prod_obj)
        msg = arrange3D.arrangeIn3D()
        if msg != '':
            print(msg)

        ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
        gen3D.make3DandOpt(reac_mol, self.forcefield, make3D = False)
        ff.Setup(Hatom.OBMol)
        gen3D.make3DandOpt(prod_obj, self.forcefield, make3D = False)
        ff.Setup(Hatom.OBMol)
        
        # Check reactant expected forming bond length must smaller than 4 angstrom after arrange. Default = 4
        dist = self.check_bond_length(reac_mol, self.form_bonds) # return the maximum value in array

        if dist >= threshold:
            return False, False
        else:
            #self.logger.info('\nStructure:\n{}\n'.format(str(prod_obj.toNode())))
            geometry = str(prod_obj.toNode()).splitlines()
            output = []
            for idx, i in enumerate(geometry):
                i_list = i.split()
                atom = i_list[0] + " "
                k = i_list[1:] + [""]

                l = " 0 ".join(k)
                """
                if idx in self.constraint:
                    l = " 0 ".join(k)
                else:
                    l = " 1 ".join(k)
                """

                out = atom + l
                output.append(out)
            output = "\n".join(output)

            reactant_geo = str(reac_mol.toNode()).splitlines()
            reac_geo = []
            for idx, i in enumerate(reactant_geo):
                i_list = i.split()
                atom = i_list[0] + " "
                k = i_list[1:] + [""]

                l = " 0 ".join(k)
                """
                if idx in self.constraint:
                    l = " 0 ".join(k)
                else:
                    l = " 1 ".join(k)
                """

                out = atom + l
                reac_geo.append(out)
            reac_geo = "\n".join(reac_geo)

            reac_mol.setCoordsFromMol(reac_mol_copy)
            #self.finalize(start_time, 'arrange')
            
            return reac_geo, output
    
        
    def finalize(self, start_time, jobname):
        """
        Finalize the job.
        """
        #self.logger.info('Total {} run time: {:.1f} s'.format(jobname, time.time() - start_time))

    def getHeatofFormation(self, tmpdir, target = 'reactant.out'):
        """
        if Error return False, which HF may be 0.0
        """
        input_path = os.path.join(tmpdir, target)
        with open(input_path, 'r') as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.strip().startswith('FINAL HEAT OF FORMATION'):
                break
        string = lines[idx].split()
        if string[0] == 'FINAL':
            HeatofFormation = string[5]
        else:
            HeatofFormation = False

        return HeatofFormation

    def runMopac(self, tmpdir, target = 'reactant.mop'):
        input_path = os.path.join(tmpdir, target)
        p = Popen(['mopac', input_path])
        p.wait()

    def check_bond_length(self, product, add_bonds):
        """
        Use reactant coordinate to check if the add bonds's bond length is too long.
        Return a 'list of distance'.
        """
        coords = [atom.coords for atom in product]
        atoms = tuple(atom.atomicnum for atom in product)
        coords = [np.array(coords).reshape(len(atoms), 3)]

        dist = []
        for bond in add_bonds:
            coord_vect_1 = coords[0][bond[0]]
            coord_vect_2 = coords[0][bond[1]]
            diff = coord_vect_1 - coord_vect_2
            dist.append(np.linalg.norm(diff))
            
        if dist == []:
            dist = [0]
        return float(max(dist))