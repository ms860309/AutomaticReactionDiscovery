import os
from subprocess import Popen

def getXYZfilepath():
    """
    return xyz path and dir path
    """
    xyz = []
    dirs = []
    abspath_pardir = os.path.abspath(os.pardir)
    rxn_path = os.path.join(abspath_pardir, 'reactions')
    dirs = os.listdir(rxn_path)
    for i in dirs[1:]:
        xyz_path = os.path.join(os.path.join(rxn_path, i), "reactant.xyz")
        xyz.append(xyz_path)
        dirs.append(os.path.join(rxn_path, i))
    return xyz, dirs

def getISOMERSfilepath():
    """
    return add_bonds path
    """
    add_bonds = []
    abspath_pardir = os.path.abspath(os.pardir)
    rxn_path = os.path.join(abspath_pardir, 'reactions')
    dirs = os.listdir(rxn_path)
    for i in dirs[1:]:
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
    for i in range(len(xyz)):
        cd_rxn_dir = 'cd {}'.format(rxn_dirs[i])
        p = Popen([cd_rxn_dir], shell = True)
        p.wait()
        mkdir_SSM = 'mkdir SSM'
        p = Popen([mkdir_SSM], shell = True)
        p.wait()
        cd_SSM = 'cd SSM'
        p = Popen([cd_SSM], shell = True)
        p.wait()
        cmd = 'python gsm.py -xyzfile {} -isomers {}  -lot_inp_file {}'.format(xyz[i], addbonds[i], qstart_path)
        p = Popen([cmd], shell = True)
        p.wait()

if __name__=='__main__':
    main()