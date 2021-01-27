from connect import db
import subprocess
import os
from os import path
import shutil
from shutil import copyfile
import sys
sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'script'))
# gsm is in script dir for launch gsm

"""
Submmit energy calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
def select_calE_target():
    qm_collection = db['qm_calculate_center']
    query = {"energy_status":"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_energy_jobs(num = 100, level_of_theory='ORCA', ncpus = 4, mpiprocs = 1, ompthreads = 4):
    targets = select_calE_target()
    
    for target in targets[:num]:
        Energy_dir_path = path.join(target['path'], 'ENERGY')
        if os.path.exists(Energy_dir_path):
            shutil.rmtree(Energy_dir_path)
            os.mkdir(Energy_dir_path)
            os.chdir(Energy_dir_path)
        else:
            os.mkdir(Energy_dir_path)
            os.chdir(Energy_dir_path)
        if level_of_theory.upper() == 'QCHEM':
            subfile = create_qchem_energy_sub_file(target['path'], Energy_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        elif level_of_theory == 'ORCA':
            subfile = create_orca_energy_sub_file(target['path'], Energy_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_energy_status(target, job_id)

def create_qchem_energy_sub_file(dir_path, Energy_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    subfile = path.join(Energy_dir_path, 'cal_E.job')
    energy_input_file = path.join(Energy_dir_path, 'energy.in')
    energy_output_file = path.join(Energy_dir_path, 'energy.out')
    reactant_xyz_path = path.join(dir_path, 'reactant.xyz')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(Energy_dir_path))), 'config')
    energy_lot = path.join(base_dir_path, 'qchem_opt_freq.lot')

    with open(energy_lot) as f:
        config = [line.strip() for line in f]

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(Energy_dir_path)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, energy_input_file, energy_output_file)
    clean_scratch = 'rm -r $QCSCRATCH'

    with open(reactant_xyz_path, 'r') as f:
        lines = f.read().splitlines()
    with open(energy_input_file, 'w') as f2:
        for i, text in enumerate(config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                config[(i+1):(i+1)] = cblock
                break
        for line in config:
            f2.write(line + '\n')

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))
    return subfile

def create_orca_energy_sub_file(dir_path, Energy_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    subfile = path.join(Energy_dir_path, 'cal_E.job')
    energy_input_file = path.join(Energy_dir_path, 'energy.in')
    reactant_xyz_path = path.join(dir_path, 'reactant.xyz')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(Energy_dir_path))), 'config')
    energy_lot = path.join(base_dir_path, 'orca_opt_freq.lot')

    with open(energy_lot) as f:
        config = [line.strip() for line in f]

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/energy.in $QCSCRATCH/\ncd $QCSCRATCH\n'
    command = '$orcadir/orca $QCSCRATCH/energy.in >> $PBS_O_WORKDIR/energy.out'
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(reactant_xyz_path, 'r') as f1:
        lines = f1.read().splitlines()
    with open(energy_input_file, 'w') as f2:
        for line in config:
            f2.write(line + '\n')
        f2.write('%pal\nnprocs {}\nend\n'.format(ncpus))
        f2.write('\n%maxcore 1000\n\n')
        f2.write('*xyz 0 1\n')
        for line in lines[2:]:
            f2.write(line + '\n')
        f2.write('*')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch, command, clean_scratch))
    return subfile

def update_energy_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    update_field = {"energy_status":"job_launched", "energy_jobid":job_id}
    qm_collection.update_one(target, {"$set": update_field}, True)
    
"""
Submmit SSM calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_ssm_target():
    qm_collection = db['qm_calculate_center']
    query = {"ssm_status":"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_ssm_jobs(num = 100, level_of_theory='QCHEM', ncpus = 4, mpiprocs = 1, ompthreads = 4):
    targets = select_ssm_target()
    
    for target in targets[:num]:
        SSM_dir_path = path.join(target['path'], 'SSM/')
        os.mkdir(SSM_dir_path)
        os.chdir(SSM_dir_path)
        if level_of_theory.upper() == 'QCHEM':
            subfile = create_qchem_ssm_sub_file(target['path'], SSM_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        elif level_of_theory == 'ORCA':
            subfile = create_orca_ssm_sub_file(target['path'], SSM_dir_path)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_ssm_status(target, job_id)
        
def create_qchem_ssm_sub_file(dir_path, SSM_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    subfile = path.join(SSM_dir_path, 'qchem_ssm.job')
    xyz_file = path.join(dir_path, 'reactant.xyz')
    isomers = path.join(dir_path, 'add_bonds.txt')
    lot_inp_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'qchem_qstart')
    ssm_args = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'ssm_argument')
    frozen_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'frozen.txt')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(SSM_dir_path)
    nes1 = 'module load qchem'
    # activate conda env is necessary because gsm install on the environment
    nes2 = 'source ~/.bashrc\nconda activate rmg3'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    if os.path.exists(frozen_file):
        command = 'gsm -xyzfile {} -mode SE_GSM -package Orca -isomers {} -lot_inp_file {} -frozen_coord_idx_file {} '.format(xyz_file, isomers, lot_inp_file, frozen_file)
    else:
        command = 'gsm -xyzfile {} -mode SE_GSM -package Orca -isomers {} -lot_inp_file {} '.format(xyz_file, isomers, lot_inp_file)
    with open(ssm_args, 'r') as f:
        lines = f.read().splitlines()
    command = command + ' '.join(lines) + ' > status.log 2>&1 '
    clean_scratch = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, scratch, command, clean_scratch))
    return subfile

def create_orca_ssm_sub_file(dir_path, SSM_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    subfile = path.join(SSM_dir_path, 'orca_ssm.job')
    xyz_file = path.join(dir_path, 'reactant.xyz')
    isomers = path.join(dir_path, 'add_bonds.txt')
    lot_inp_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'orca_qstart')
    ssm_args = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'ssm_argument')
    frozen_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'frozen.txt')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(SSM_dir_path)
    # activate conda env is necessary because gsm install on the environment
    nes2 = 'source ~/.bashrc\nconda activate rmg3'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    if os.path.exists(frozen_file):
        command = 'gsm -xyzfile {} -mode SE_GSM -package Orca -isomers {} -lot_inp_file {} -frozen_coord_idx_file {} '.format(xyz_file, isomers, lot_inp_file, frozen_file)
    else:
        command = 'gsm -xyzfile {} -mode SE_GSM -package Orca -isomers {} -lot_inp_file {} '.format(xyz_file, isomers, lot_inp_file)
      
    with open(ssm_args, 'r') as f:
        lines = f.read().splitlines()
    command = command + ' '.join(lines) + ' > status.log 2>&1 '
    clean_scratch = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes2, scratch, command, clean_scratch))
    return subfile

def update_ssm_status(target, job_id):
    collect = db['qm_calculate_center']
    update_field = {"ssm_status":"job_launched", "ssm_jobid":job_id}
    collect.update_one(target, {"$set": update_field}, True)

"""
Submmit OPT calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_opt_target():
    qm_collection = db['qm_calculate_center']
    query = {"opt_status":"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_opt_jobs(num=100, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    targets = select_opt_target()
    for target in targets[:num]:
        dir_path = target['path']
        OPT_dir_path = path.join(dir_path, 'OPT/')
        if not os.path.exists(OPT_dir_path):
            os.mkdir(OPT_dir_path)
        os.chdir(OPT_dir_path)
        
        create_opt_input(dir_path, OPT_dir_path)
        subfile = create_opt_sub_file(dir_path, OPT_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_opt_status(target, job_id)

def create_opt_input(dir_path, OPT_dir_path):
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(OPT_dir_path)))), 'config')
    opt_lot = path.join(base_dir_path, 'opt_freq.lot')
    opt_input = path.join(OPT_dir_path, 'opt.in')
    reactant_path = path.join(dir_path, 'reactant.xyz')
    with open(opt_lot) as f:
        config = [line.strip() for line in f]
    with open(reactant_path) as f:
        lines = f.read().split('\n')
        geo = [line for line in lines[2:]]
    with open(opt_input, 'w') as f:
        f.write('$molecule\n{} {}\n'.format(0, 1))
        f.write('\n'.join(geo))
        f.write('\n')
        f.write('$end\n\n')
        for line in config:
            f.write(line + '\n')

def create_opt_sub_file(dir_path, OPT_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    subfile = path.join(OPT_dir_path, 'cal_opt.job')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(OPT_dir_path)
    nes1 = 'module load qchem'
    nes2 = 'export QCSCRATCH=/tmp/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4 = 'qchem -nt {} {} {}'.format(ncpus, 'opt.in', 'opt.out')
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5))
    return subfile
    
def update_opt_status(target, job_id):
    collect = db['qm_calculate_center']
    update_field = {"opt_status":"job_launched", "opt_jobid":job_id}
    collect.update_one(target, {"$set": update_field}, True)

"""
Submmit Low level OPT calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_low_opt_target():
    qm_collection = db['qm_calculate_center']
    query = {"low_opt_status":"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_low_opt_jobs(num=100, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    targets = select_low_opt_target()
    
    for target in targets[:num]:
        dir_path = target['path']
        OPT_dir_path = path.join(dir_path, 'OPT/')
        if not os.path.exists(OPT_dir_path):
            os.mkdir(OPT_dir_path)
        os.chdir(OPT_dir_path)
        
        create_low_opt_input(dir_path, OPT_dir_path)
        subfile = create_low_opt_sub_file(dir_path, OPT_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_low_opt_status(target, job_id)

def create_low_opt_input(dir_path, OPT_dir_path):
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(OPT_dir_path)))), 'config')
    low_opt_lot = path.join(base_dir_path, 'low_opt.lot')
    low_opt_input = path.join(OPT_dir_path, 'low_opt.in')
    reactant_path = path.join(dir_path, 'reactant.xyz')
    with open(low_opt_lot) as f:
        config = [line.strip() for line in f]
    with open(reactant_path) as f:
        lines = f.read().split('\n')
        geo = [line for line in lines[2:]]
    with open(low_opt_input, 'w') as f:
        f.write('$molecule\n{} {}\n'.format(0, 1))
        f.write('\n'.join(geo))
        f.write('\n')
        f.write('$end\n\n')
        for line in config:
            f.write(line + '\n')

def create_low_opt_sub_file(dir_path, OPT_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    subfile = path.join(OPT_dir_path, 'cal_low_opt.job')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(OPT_dir_path)
    nes1 = 'module load qchem'
    nes2 = 'export QCSCRATCH=/tmp/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4 = 'qchem -nt {} {} {}'.format(ncpus, 'low_opt.in', 'low_opt.out')
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5))
    return subfile
    
def update_low_opt_status(target, job_id):
    collect = db['qm_calculate_center']
    update_field = {"low_opt_status":"job_launched", "low_opt_jobid":job_id}
    collect.update_one(target, {"$set": update_field}, True)

"""
Submmit TS refine calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
Here refine is not to get a high level TS.
It's just for ORCA with xtb GFN2-xtb level of theory SSM. To get a better TS initial guess
"""

def select_ts_refine_target():
    qm_collection = db['qm_calculate_center']
    query = {"ts_refine_status":"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_ts_refine_jobs(num=100, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    targets = select_ts_refine_target()
    
    for target in targets[:num]:
        TS_dir_path = path.join(target['path'], 'TS/')
        os.mkdir(TS_dir_path)
        os.chdir(TS_dir_path)
            
        SSM_dir_path = path.join(target['path'], 'SSM/')
        subfile = create_ts_refine_sub_file(SSM_dir_path, TS_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_ts_refine_status(target, job_id)
        
def create_ts_refine_sub_file(SSM_dir_path, TS_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts_refine.in')
    #ts_output_file = path.join(TS_dir_path, 'ts_refine.out')
    subfile = path.join(TS_dir_path, 'orca_ts_refine.job')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(SSM_dir_path)))), 'config')
    ts_lot = path.join(base_dir_path, 'orca_refine_ts')

    with open(ts_lot) as f:
        config = [line.strip() for line in f]

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/ts_refine.in $QCSCRATCH/\ncd $QCSCRATCH\n'
    command = '$orcadir/orca $QCSCRATCH/ts_refine.in >> $PBS_O_WORKDIR/ts_refine.out'
    copy_the_refine_xyz = 'cp $QCSCRATCH/ts_refine.xyz $PBS_O_WORKDIR'
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(tsnode_path, 'r') as f1:
        lines = f1.read().splitlines()
    with open(ts_input_file, 'w') as f2:
        for line in config:
            f2.write(line + '\n')
        f2.write('%pal\nnprocs {}\nend\n'.format(ncpus))
        f2.write('\n%maxcore 1000\n\n')
        f2.write('*xyz 0 1\n')
        for line in lines[2:]:
            f2.write(line + '\n')
        f2.write('*')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch, command, copy_the_refine_xyz, clean_scratch))
        
    return subfile
    

def update_ts_refine_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    update_field = {"ts_refine_status":"job_launched", "ts_refine_jobid":job_id}
    qm_collection.update_one(target, {"$set": update_field}, True)

"""
Submmit TS calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_ts_target():
    qm_collection = db['qm_calculate_center']
    query = {"ts_status":"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_ts_jobs(num=100, level_of_theory='ORCA', ncpus = 4, mpiprocs = 1, ompthreads = 4):
    targets = select_ts_target()
    
    for target in targets[:num]:
        SSM_dir_path = path.join(target['path'], 'SSM/')
        TS_dir_path = path.join(target['path'], 'TS/')
        if os.path.exists(TS_dir_path):
            os.chdir(TS_dir_path)
        else:
            os.mkdir(TS_dir_path)
            os.chdir(TS_dir_path)

        if level_of_theory.upper() == 'QCHEM':
            subfile = create_qchem_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        elif level_of_theory == 'ORCA':
            subfile = create_orca_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_ts_status(target, job_id)
        
def create_qchem_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    refine_path = path.join(TS_dir_path, 'ts_refine.xyz')
    if os.path.exists(refine_path):
        tsnode_path = refine_path
    else:
        tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts.in')
    ts_output_file = path.join(TS_dir_path, 'ts.out')
    subfile = path.join(TS_dir_path, 'qchem_ts.job')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(SSM_dir_path)))), 'config')
    ts_lot = path.join(base_dir_path, 'qchem_freq_ts_freq.lot')

    with open(ts_lot) as f:
        config = [line.strip() for line in f]

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(TS_dir_path)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, ts_input_file, ts_output_file)
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(tsnode_path, 'r') as f1:
        lines = f1.read().splitlines()
    with open(ts_input_file, 'w') as f2:
        for i, text in enumerate(config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                config[(i+1):(i+1)] = cblock
                break
        for line in config:
            f2.write(line + '\n')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command, clean_scratch))
    return subfile

def create_orca_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    refine_path = path.join(TS_dir_path, 'ts_refine.xyz')
    if os.path.exists(refine_path):
        tsnode_path = refine_path
    else:
        tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts_geo.in')
    subfile = path.join(TS_dir_path, 'orca_ts.job')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(SSM_dir_path)))), 'config')
    ts_lot = path.join(base_dir_path, 'orca_freq_ts_freq.lot')

    with open(ts_lot) as f:
        config = [line.strip() for line in f]

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/ts_geo.in $QCSCRATCH/\ncd $QCSCRATCH\n'
    command = '$orcadir/orca $QCSCRATCH/ts_geo.in >> $PBS_O_WORKDIR/ts_geo.out'
    copy_the_refine_xyz = 'cp $QCSCRATCH/ts_geo.xyz $PBS_O_WORKDIR'
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(tsnode_path, 'r') as f1:
        lines = f1.read().splitlines()
    with open(ts_input_file, 'w') as f2:
        for line in config:
            f2.write(line + '\n')
        f2.write('%pal\nnprocs {}\nend\n'.format(ncpus))
        f2.write('\n%maxcore 1000\n\n')
        f2.write('*xyz 0 1\n')
        for line in lines[2:]:
            f2.write(line + '\n')
        f2.write('*')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch, command, copy_the_refine_xyz, clean_scratch))
    return subfile

def update_ts_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    update_field = {"ts_status":"job_launched", "ts_jobid":job_id}
    qm_collection.update_one(target, {"$set": update_field}, True)
    
"""
Submmit IRC(Intrinsic Reaction Coordinate, in Qchem called 'rpath') calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_irc_target():
    qm_collection = db['qm_calculate_center']
    query = {"irc_status":"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_irc_jobs(num = 100, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    targets = select_irc_target()
    for target in targets[:num]:
        IRC_dir_path = path.join(target['path'], 'IRC/')
        os.mkdir(IRC_dir_path)
        os.chdir(IRC_dir_path)

        TS_dir_path = path.join(target['path'], 'TS/')
        subfile = create_irc_sub_file(TS_dir_path, IRC_dir_path, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_status(target, job_id)
        
def create_irc_sub_file(TS_dir_path, IRC_dir_path, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    ts_geo_path = path.join(TS_dir_path, 'ts_geo.xyz')
    irc_input_file = path.join(IRC_dir_path, 'pysisyphus_irc.yaml')
    subfile = path.join(IRC_dir_path, 'irc.job')
    new_ts_geo_path = path.join(IRC_dir_path, 'ts_geo.xyz')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(IRC_dir_path)))), 'config')
    irc_lot = path.join(base_dir_path, 'pysisyphus_irc.yaml')
    copyfile(irc_lot, irc_input_file)
    copyfile(ts_geo_path, new_ts_geo_path)
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(IRC_dir_path)
    nes = 'source ~/.bashrc\nconda activate rmg3\nexport TMPDIR=/tmp/$PBS_JOBID\nmkdir -p $TMPDIR\n'
    command = 'pysis pysisyphus_irc.yaml'
    deactivate = 'conda deactivate\nrm *.gbw\nrm -r $TMPDIR'

    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes, command, deactivate))
            
    return subfile
    

def update_irc_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    update_field = {'irc_status':'job_launched', 'irc_jobid':job_id}
    qm_collection.update_one(target, {"$set": update_field}, True)

"""
Submmit opt job which is from irc
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""

def select_irc_opt_target():
    
    qm_collection = db['qm_calculate_center']
    query = {'irc_opt_status':"job_unrun"}
    targets = list(qm_collection.find(query))
    return targets

def launch_irc_opt_jobs(num = 100, level_of_theory='ORCA', ncpus = 4, mpiprocs = 1, ompthreads = 4):
    
    targets = select_irc_opt_target()
    for target in targets[:num]:
        IRC_dir_path = path.join(target['path'], 'IRC/')
        os.chdir(IRC_dir_path)

        irc_equal = target['irc_equal']
        forward_end_output = os.path.join(IRC_dir_path, 'forward_end_opt.xyz')
        backward_end_output = os.path.join(IRC_dir_path, 'backward_end_opt.xyz')
        if os.path.exists(forward_end_output) and os.path.exists(backward_end_output):
            if irc_equal == 'forward equal to reactant and reverse equal to product':
                target_mol = backward_end_output
            elif irc_equal == 'reverse equal to reactant and forward equal to product':
                target_mol = forward_end_output
            elif irc_equal == 'forward equal to reactant but reverse does not equal to product':
                target_mol = backward_end_output
            elif irc_equal == 'reverse equal to reactant but forward does not equal to product':
                target_mol = forward_end_output
        else:
            if irc_equal == 'forward equal to reactant and reverse equal to product':
                target_mol = path.join(IRC_dir_path, 'finished_last.xyz')
            elif irc_equal == 'reverse equal to reactant and forward equal to product':
                target_mol = path.join(IRC_dir_path, 'finished_first.xyz')
            elif irc_equal == 'forward equal to reactant but reverse does not equal to product':
                target_mol = path.join(IRC_dir_path, 'finished_last.xyz')
            elif irc_equal == 'reverse equal to reactant but forward does not equal to product':
                target_mol = path.join(IRC_dir_path, 'finished_first.xyz')

        if level_of_theory.upper() == 'QCHEM':
            subfile = create_qchem_irc_opt_sub_file(IRC_dir_path, target_mol, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        elif level_of_theory == 'ORCA':
            subfile = create_orca_irc_opt_sub_file(IRC_dir_path, target_mol, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_opt_status(target, job_id)

def create_qchem_irc_opt_sub_file(irc_path, target_mol, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(irc_path)))), 'config')
    irc_opt_lot = path.join(base_dir_path, 'qchem_opt_freq.lot')
    irc_opt_input = path.join(irc_path, 'irc_opt.in')
    subfile = path.join(irc_path, 'irc_opt.job')

    with open(irc_opt_lot) as f:
        config = [line.strip() for line in f]
    with open(target_mol, 'r') as f1:
        lines = f1.read().splitlines()
    with open(irc_opt_input, 'w') as f2:
        for i, text in enumerate(config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                config[(i+1):(i+1)] = cblock
                break
        for line in config:
            f2.write(line + '\n')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(irc_path)
    nes1 = 'module load qchem'
    nes2 = 'export QCSCRATCH=/tmp/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4 = 'qchem -nt {} irc_opt.in irc_opt.out'.format(ncpus)
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5))
    return subfile

def create_orca_irc_opt_sub_file(irc_path, target_mol, ncpus = 4, mpiprocs = 1, ompthreads = 4):
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(irc_path)))), 'config')
    irc_opt_lot = path.join(base_dir_path, 'orca_opt_freq.lot')
    irc_opt_input = path.join(irc_path, 'irc_reactant.in')
    subfile = path.join(irc_path, 'irc_opt.job')

    with open(irc_opt_lot) as f:
        config = [line.strip() for line in f]

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe\nsource ~/.bashrc\n'.format(ncpus, mpiprocs, ompthreads)
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\ncp $PBS_O_WORKDIR/irc_reactant.in $QCSCRATCH/\ncd $QCSCRATCH\n'
    command = '$orcadir/orca $QCSCRATCH/irc_reactant.in >> $PBS_O_WORKDIR/irc_reactant.out'
    copy_the_refine_xyz = 'cp $QCSCRATCH/irc_reactant.xyz $PBS_O_WORKDIR'
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(target_mol, 'r') as f1:
        lines = f1.read().splitlines()
    with open(irc_opt_input, 'w') as f2:
        for line in config:
            f2.write(line + '\n')
        f2.write('%pal\nnprocs {}\nend\n'.format(ncpus))
        f2.write('\n%maxcore 1000\n\n')
        f2.write('*xyz 0 1\n')
        for line in lines[2:]:
            f2.write(line + '\n')
        f2.write('*')
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n'.format(shell, pbs_setting, scratch, command, copy_the_refine_xyz, clean_scratch))
    return subfile

def update_irc_opt_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    update_field = {'irc_opt_status':"job_launched", 'irc_opt_jobid':job_id}
    qm_collection.update_one(target, {"$set": update_field}, True)


def launch_jobs(num = 30, level_of_theory='ORCA', ncpus = 4, mpiprocs = 1, ompthreads = 4):
    launch_energy_jobs(num = num, level_of_theory='ORCA', ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
    launch_ssm_jobs(num = 50, level_of_theory='ORCA', ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
    launch_ts_refine_jobs(num = num, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
    launch_ts_jobs(num = num, level_of_theory='ORCA', ncpus = 8, mpiprocs = mpiprocs, ompthreads = 8)
    launch_irc_jobs(num = num, ncpus = 8, mpiprocs = mpiprocs, ompthreads = 8)
    launch_irc_opt_jobs(num = num, level_of_theory='ORCA', ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
    
    #launch_low_opt_jobs(num = num, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)
    #launch_opt_jobs(num = num, ncpus = ncpus, mpiprocs = mpiprocs, ompthreads = ompthreads)


launch_jobs(num = 30, level_of_theory='ORCA', ncpus = 4, mpiprocs = 1, ompthreads = 4)
