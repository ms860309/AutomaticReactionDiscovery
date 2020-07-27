from node import Node
import shutil
import os
from subprocess import Popen, PIPE
import gen3D
import pybel
import psutil
import time
import logging
import util
import pybel
import difflib

class MopacError(Exception):
    """
    An exception class for errors that occur during mopac calculations.
    """
    pass

class mopac(object):

    def __init__(self, forcefield, reactant_bonds, product_bonds):
        self.forcefield = forcefield
        self.reactant_bonds = reactant_bonds
        self.product_bonds = product_bonds
        log_level = logging.INFO
        process = psutil.Process(os.getpid())
        self.logger = util.initializeLog(log_level, os.path.join(os.getcwd(), 'ARD.log'), logname='main')
        self.logger.info('\nARD initiated on ' + time.asctime() + '\n')
        self.logger.info('memory usage: {}'.format(process.memory_percent()))

    def mopac_get_H298(self, reac_obj, InputFile, charge = 0, multiplicity = 'SINGLET', method = 'PM7'):
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
        reac_geo, geometry = self.genInput(reac_obj, InputFile)
        
        with open(reactant_path, 'w') as f:
            f.write("LARGE CHARGE={} {} {}\n\n".format(charge, multiplicity, method))
            f.write("\n{}".format(reac_geo))
        start_time = time.time()
        self.runMopac(tmpdir, 'reactant.mop')
        reactant = self.getHeatofFormation(tmpdir, 'reactant.out')

        with open(product_path, 'w') as f:
            f.write("LARGE CHARGE={} {} {}\n\n".format(charge, multiplicity, method))
            f.write("\n{}".format(geometry))
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
    
    def genInput(self, reac_mol, InputFile):
        start_time = time.time()

        Hatom = gen3D.readstring('smi', '[H]')
        ff = pybel.ob.OBForceField.FindForceField(self.forcefield)

        reac_mol.gen3D(forcefield=self.forcefield, make3D=False)
        InputFile.gen3D(forcefield=self.forcefield, make3D=False)
        reac_mol_copy, InputFile_copy= reac_mol.copy(), InputFile.copy()

        arrange3D = gen3D.Arrange3D(reac_mol, InputFile, self.reactant_bonds, self.product_bonds)
        msg = arrange3D.arrangeIn3D()
        if msg != '':
            print(msg)

        ff.Setup(Hatom.OBMol)  # Ensures that new coordinates are generated for next molecule (see above)
        reac_mol.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)
        InputFile.gen3D(make3D=False)
        ff.Setup(Hatom.OBMol)

        geometry = str(InputFile.toNode())

        #self.fast_bonds_filter(geometry)
        self.logger.info('\nStructure:\n{}\n'.format(str(geometry)))
        geometry = geometry.splitlines()

        output = []
        for i in geometry:
            i_list = i.split()
            atom = i_list[0] + " "
            k = i_list[1:] + [""]
            l = " 1 ".join(k)
            out = atom + l
            output.append(out)
        output = "\n".join(output)

        reactant_geo = str(reac_mol.toNode()).splitlines()
        reac_geo = []
        for i in reactant_geo:
            i_list = i.split()
            atom = i_list[0] + " "
            k = i_list[1:] + [""]
            l = " 1 ".join(k)
            out = atom + l
            reac_geo.append(out)
        reac_geo = "\n".join(reac_geo)

        reac_mol.setCoordsFromMol(reac_mol_copy)

        #self.finalize(start_time, 'arrange')
        
        return reac_geo, output
    
    def fast_bonds_filter(self, geometry):
        xyz_path = os.path.join(os.path.dirname(os.getcwd()), 'filter')
        if os.path.exists(xyz_path):
            shutil.rmtree(xyz_path)
        os.mkdir(xyz_path)
        prod_path = os.path.join(xyz_path, 'prod.xyz')
        with open(prod_path, 'w') as f:
            f.write('17\n\n')
            f.write(str(geometry))
        mol = next(pybel.readfile('xyz', prod_path))
        bonds_mol_1 = []
        for bond in pybel.ob.OBMolBondIter(mol.OBMol):
            a = bond.GetBeginAtomIdx() - 1
            b = bond.GetEndAtomIdx() - 1
            if a < b:
                bonds_mol_1.append((a, b))
            else:
                bonds_mol_1.append((b, a))
        bonds_mol_1 = sorted(bonds_mol_1)
        a = difflib.SequenceMatcher(None, self.reactant_bonds, bonds_mol_1).ratio()
        b = difflib.SequenceMatcher(None, self.product_bonds, bonds_mol_1).ratio()
        print(a*100)
        print(b*100)
        
    def finalize(self, start_time, jobname):
        """
        Finalize the job.
        """
        self.logger.info('Total {} run time: {:.1f} s'.format(jobname, time.time() - start_time))

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
