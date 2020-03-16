import os
from subprocess import Popen
import shutil
def getXYZfilepath():
    """
    return xyz path and dir path
    """
    xyz = []
    rxn_dirs = []
    abspath_pardir = os.path.abspath(os.pardir)
    rxn_path = os.path.join(abspath_pardir, 'reactions')
    dirs = os.listdir(rxn_path)
    dirs.remove('00000')
    for i in dirs:
        xyz_path = os.path.join(os.path.join(rxn_path, i), "reactant.xyz")
        xyz.append(xyz_path)
        rxn_dirs.append(os.path.join(rxn_path, i))
    return xyz, rxn_dirs

def getISOMERSfilepath():
    """
    return add_bonds path
    """
    add_bonds = []
    abspath_pardir = os.path.abspath(os.pardir)
    rxn_path = os.path.join(abspath_pardir, 'reactions')
    dirs = os.listdir(rxn_path)
    dirs.remove('00000')
    for i in dirs:
        add_bonds_path = os.path.join(os.path.join(rxn_path, i), "add_bonds.txt")
        add_bonds.append(add_bonds_path)
    return add_bonds

def getQstartpath():
    abspath_pardir = os.path.abspath(os.pardir)
    qstart_path = os.path.join(os.path.join(abspath_pardir, 'submmit_required'), 'qstart')
    return qstart_path

def main():
    xyz, rxn_dirs= getXYZfilepath()
    addbonds = getISOMERSfilepath()
    qstart_path = getQstartpath()
    gsm_py_path = os.path.join(os.path.abspath(""), "gsm.py")
    for i in range(len(xyz)):
        # change workspace dir reaction/num
        if os.path.exists(rxn_dirs[i]):
            os.chdir(rxn_dirs[i])
        # mkdir SSM
        SSM_path = os.path.join(rxn_dirs[i], 'SSM')
        if os.path.exists(SSM_path):
            shutil.rmtree(SSM_path)
        os.mkdir(SSM_path)
        # change workspace dir reaction/num/SSM
        if os.path.exists(SSM_path):
            os.chdir(SSM_path)
        try:
            cmd = 'python {} -xyzfile {} -isomers {}  -lot_inp_file {}'.format(gsm_py_path, xyz[i], addbonds[i], qstart_path)
            p = Popen([cmd], shell = True)
            p.wait()
        except:
            continue
if __name__=='__main__':
    main()