from connect import db
import subprocess
import os
from os import path
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
    collect = db['molecules']
    reg_query = {"energy_status":"job_unrun"}
    targets = list(collect.find(reg_query))

    selected_targets = []
    for target in targets:
        dir_path = target['path']
        selected_targets.append(dir_path)

    return selected_targets

def launch_energy_jobs():
    targets = select_calE_target()
    
    for target in targets:
        if os.path.exists(target):
            os.chdir(target)
        subfile = create_energy_sub_file(target)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode()
        # update status job_launched
        update_energy_status(target, job_id)
        
def create_energy_sub_file(path, ncpus = 1, mpiprocs = 1):
    subfile = os.path.join(path, 'calculate_energy.job')
    
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}'.format(ncpus, mpiprocs)
    target_path = 'cd {}'.format(path)
    nes1 = 'source ~/.bashrc\n conda activate rmg3'
    nes2 = 'export QCSCRATCH=/tmp/ypli/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4 = 'qchem -nt 2 energy.in energy.out'
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}').format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5)
    return subfile
    
def update_energy_status(target, job_id):
    collect = db['molecules']
    reg_query = {"path":target}
    update_field = {"energy_status":"job_launched", "energy_jobid":job_id}
    collect.update_one(reg_query, {"$set": update_field}, True)
    
"""
Submmit SSM calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_ssm_target():
    collect = db['molecules']
    reg_query = {"ssm_status":"job_unrun"}
    targets = list(collect.find(reg_query))

    selected_targets = []
    for target in targets:
        dir_path = target['path']
        selected_targets.append(dir_path)

    return selected_targets

def launch_ssm_jobs():
    targets = select_ssm_target()
    
    for target in targets:
        if os.path.exists(target):
            os.chdir(target)
            
        subfile = create_ssm_sub_file(target)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode()
        # update status job_launched
        update_ssm_status(target, job_id)
        
def create_ssm_sub_file(path, ncpus = 1, mpiprocs = 1):
    subfile = path.join(path, 'calculate_ssm.job')
    
    xyz_file = path.join(path, 'reactant.xyz')
    isomers = path.join(path, 'add_bonds.txt')
    lot_inp_file = path.join(path.join(path.dirname( path.dirname( path.abspath(__file__))),'submmit_required'), qstart)
    
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}'.format(ncpus, mpiprocs)
    target_path = 'cd {}'.format(path)
    nes1 = nes1 = 'source ~/.bashrc\n conda activate rmg3'
    nes2 = 'export QCSCRATCH=/tmp/ypli/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4 = 'gsm -xyzfile {} -mode SE_GSM -package QChem -isomers {} -lot_inp_file {}'.format(xyz_file, isomers, lot_inp_file)
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}').format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5)
    return subfile
    
def update_ssm_status(target, job_id):
    collect = db['molecules']
    reg_query = {"path":target}
    update_field = {"ssm_status":"job_launched", "ssm_jobid":job_id}
    collect.update_one(reg_query, {"$set": update_field}, True)
    
launch_energy_jobs()
launch_ssm_jobs()