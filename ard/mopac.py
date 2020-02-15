from node import Node
import shutil
import os
from subprocess import Popen, PIPE

class mopac(object):

    def __init__(self, reac_mol):
        self.reac_mol = reac_mol

    def makeInputFile(self, InputFile, charge = 0, multiplicity = 'SINGLET', method = 'PM6'):
        """
        Create a directory folder called "tmp" for mopac calculation
        Create a input file called "input.mop" for mopac calculation
        """
        path = os.path.abspath("")
        tmpdir = os.path.join(path, 'tmp')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        input_path = os.path.join(tmpdir, 'input.mop')

        geometry = self.genInput(InputFile)
        with open(input_path, 'w') as f:
            f.write(" AUX LARGE CHARGE={} {} {}".format(charge, multiplicity, method))
            f.write("\n{}".format(geometry))

        self.runMopac(tmpdir)
        result = self.getHeatofFormation(tmpdir)
        return float(result)
    
    def genInput(self, InputFile):

        geometry = InputFile.toNode()
        geometry = str(geometry)
        geometry = geometry.splitlines()
        output = []
        for i in geometry:
            i_list = i.split()
            atom = i_list[0]
            atom = atom + " "
            k = i_list[1:]
            l = " 1 ".join(k)
            out = atom + l
            output.append(out)
        output = "\n".join(output)

        return output

    def getHeatofFormation(self, tmpdir):
        input_path = os.path.join(tmpdir, "input.out")
        with open(input_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.strip().startswith('FINAL HEAT OF FORMATION'):
                    HeatofFormation = line.split()[5]
        return HeatofFormation

    def runMopac(self, tmpdir):
        input_path = os.path.join(tmpdir, "input.mop")
        p = Popen(['mopac', input_path])
        p.wait()