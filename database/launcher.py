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
    reg_query = {"energy_status":"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    selected_targets = [target['path'] for target in targets]
    return selected_targets

def launch_energy_jobs():
    targets = select_calE_target()
    
    for target in targets:
        Energy_dir_path = path.join(target, 'ENERGY')
        if os.path.exists(Energy_dir_path):
            shutil.rmtree(Energy_dir_path)
            os.mkdir(Energy_dir_path)
            os.chdir(Energy_dir_path)
        else:
            os.mkdir(Energy_dir_path)
            os.chdir(Energy_dir_path)
            
        subfile = create_energy_sub_file(target, Energy_dir_path)
        
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_energy_status(target, job_id)

def create_energy_sub_file(dir_path, Energy_dir_path, ncpus = 8, mpiprocs = 1, ompthreads = 8):
    subfile = path.join(Energy_dir_path, 'cal_E.job')
    energy_input_file = path.join(Energy_dir_path, 'energy.in')
    energy_output_file = path.join(Energy_dir_path, 'energy.out')
    reactant_xyz_path = path.join(dir_path, 'reactant.xyz')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(Energy_dir_path))), 'config')
    energy_lot = path.join(base_dir_path, 'opt_sp.lot')

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
    
def update_energy_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    reg_query = {"path":target}
    update_field = {"energy_status":"job_launched", "energy_jobid":job_id}
    qm_collection.update_one(reg_query, {"$set": update_field}, True)
    
"""
Submmit SSM calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_ssm_target():
    qm_collection = db['qm_calculate_center']
    reg_query = {"ssm_status":"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    selected_targets = [target['path'] for target in targets]
    return selected_targets

def launch_ssm_jobs():
    targets = select_ssm_target()
    
    for target in targets:
        SSM_dir_path = path.join(target, 'SSM/')
        os.mkdir(SSM_dir_path)
        os.chdir(SSM_dir_path)
            
        subfile = create_ssm_sub_file(target, SSM_dir_path)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_ssm_status(target, job_id)
        
def create_ssm_sub_file(dir_path, SSM_dir_path, ncpus = 8, mpiprocs = 1, ompthreads = 8):
    subfile = path.join(SSM_dir_path, 'cal_ssm.job')
    xyz_file = path.join(dir_path, 'reactant.xyz')
    isomers = path.join(dir_path, 'add_bonds.txt')
    lot_inp_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'qstart')
    ssm_args = path.join(path.join(path.dirname(path.dirname(dir_path)), 'config'), 'ssm_argument')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(SSM_dir_path)
    nes1 = 'module load qchem'
    # activate conda env is necessary because gsm install on the environment
    nes2 = 'source ~/.bashrc\nconda activate rmg3'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'gsm -xyzfile {} -mode SE_GSM -package QChem -isomers {} -lot_inp_file {} '.format(xyz_file, isomers, lot_inp_file)    
    with open(ssm_args, 'r') as f:
        lines = f.read().splitlines()
    command = command + ' '.join(lines) + ' > status.log 2>&1 '
    clean_scratch = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, scratch, command, clean_scratch))
    return subfile
    
def update_ssm_status(target, job_id):
    collect = db['qm_calculate_center']
    reg_query = {"path":target}
    update_field = {"ssm_status":"job_launched", "ssm_jobid":job_id}
    collect.update_one(reg_query, {"$set": update_field}, True)

"""
Submmit OPT calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_opt_target():
    qm_collection = db['qm_calculate_center']
    reg_query = {"opt_status":"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    return targets

def launch_opt_jobs():
    targets = select_opt_target()
    
    for target in targets:
        dir_path = target['path']
        OPT_dir_path = path.join(dir_path, 'OPT/')
        if not os.path.exists(OPT_dir_path):
            os.mkdir(OPT_dir_path)
        os.chdir(OPT_dir_path)
        
        create_opt_input(dir_path, OPT_dir_path)
        subfile = create_opt_sub_file(dir_path, OPT_dir_path)
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
    opt_lot = path.join(base_dir_path, 'opt.lot')
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

def create_opt_sub_file(dir_path, OPT_dir_path, ncpus = 8, mpiprocs = 1, ompthreads = 8):
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
    reg_query = {"low_opt_status":"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    return targets

def launch_low_opt_jobs():
    targets = select_low_opt_target()
    
    for target in targets:
        dir_path = target['path']
        OPT_dir_path = path.join(dir_path, 'OPT/')
        if not os.path.exists(OPT_dir_path):
            os.mkdir(OPT_dir_path)
        os.chdir(OPT_dir_path)
        
        create_low_opt_input(dir_path, OPT_dir_path)
        subfile = create_low_opt_sub_file(dir_path, OPT_dir_path)
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

def create_low_opt_sub_file(dir_path, OPT_dir_path, ncpus = 8, mpiprocs = 1, ompthreads = 8):
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
Submmit TS calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_ts_target():
    qm_collection = db['qm_calculate_center']
    reg_query = {"ts_status":"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    selected_targets = [target['path'] for target in targets]
    return selected_targets

def launch_ts_jobs():
    targets = select_ts_target()
    
    for target in targets:
        TS_dir_path = path.join(target, 'TS/')
        os.mkdir(TS_dir_path)
        os.chdir(TS_dir_path)
            
        SSM_dir_path = path.join(target, 'SSM/')
        subfile = create_ts_sub_file(SSM_dir_path, TS_dir_path)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_ts_status(target, job_id)
        
def create_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus = 8, mpiprocs = 1, ompthreads = 8):
    tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts.in')
    ts_output_file = path.join(TS_dir_path, 'ts.out')
    subfile = path.join(TS_dir_path, 'cal_ts.job')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(SSM_dir_path)))), 'config')
    ts_lot = path.join(base_dir_path, 'freq_ts.lot')

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
    

def update_ts_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    reg_query = {"path":target}
    update_field = {"ts_status":"job_launched", "ts_jobid":job_id}
    qm_collection.update_one(reg_query, {"$set": update_field}, True)
    
"""
Submmit IRC(Intrinsic Reaction Coordinate, in Qchem called 'rpath') calculation job
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""
    
def select_irc_target():
    qm_collection = db['qm_calculate_center']
    reg_query = {"irc_status":"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    selected_targets = [target['path'] for target in targets]
    return selected_targets

def launch_irc_jobs():
    targets = select_irc_target()
    
    for target in targets:
        IRC_dir_path = path.join(target, 'IRC/')
        os.mkdir(IRC_dir_path)
        os.chdir(IRC_dir_path)
            
        TS_dir_path = path.join(target, 'TS/')
        subfile_1,  subfile_2= create_irc_sub_file(TS_dir_path, IRC_dir_path)
        cmd_1 = 'qsub {}'.format(subfile_1)
        process_1 = subprocess.Popen([cmd_1],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process_1.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_status(target, job_id, direction = 'forward')
        
        cmd_2 = 'qsub {}'.format(subfile_2)
        process_2 = subprocess.Popen([cmd_2],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process_2.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_status(target, job_id, direction = 'reverse')
        
def create_irc_sub_file(TS_dir_path, IRC_dir_path, ncpus = 8, mpiprocs = 1, ompthreads = 8):
    ts_geo_path = path.join(TS_dir_path, 'ts_geo.xyz')
    
    irc_forward_input_file = path.join(IRC_dir_path, 'irc_forward.in')
    irc_reverse_input_file = path.join(IRC_dir_path, 'irc_reverse.in')
    irc_forward_output_file = path.join(IRC_dir_path, 'irc_forward.out')
    irc_reverse_output_file = path.join(IRC_dir_path, 'irc_reverse.out')

    subfile_1 = path.join(IRC_dir_path, 'cal_irc_forward.job')
    subfile_2 = path.join(IRC_dir_path, 'cal_irc_reverse.job')

    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(IRC_dir_path)))), 'config')
    irc_forward_lot = path.join(base_dir_path, 'freq_irc_forward.lot')
    irc_reverse_lot = path.join(base_dir_path, 'freq_irc_reverse.lot')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(IRC_dir_path)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command_1 = 'qchem -nt {} {} {}'.format(ncpus, irc_forward_input_file, irc_forward_output_file)
    command_2 = 'qchem -nt {} {} {}'.format(ncpus, irc_reverse_input_file, irc_reverse_output_file)
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(irc_forward_lot) as f:
        forward_config = [line.strip() for line in f]
    with open(irc_reverse_lot) as f:
        reverse_config = [line.strip() for line in f]

    with open(ts_geo_path, 'r') as f1:
        lines = f1.read().splitlines()

    with open(irc_forward_input_file, 'w') as f2:
        for i, text in enumerate(forward_config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                forward_config[(i+1):(i+1)] = cblock
                break
        for line in forward_config:
            f2.write(line + '\n')
    with open(irc_reverse_input_file, 'w') as f2:
        for i, text in enumerate(reverse_config):
            if text.startswith('$molecule'):
                cblock = lines[2:]
                cblock.insert(0, '0  1')
                reverse_config[(i+1):(i+1)] = cblock
                break
        for line in reverse_config:
            f2.write(line + '\n')

    with open(subfile_1, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command_1, clean_scratch))
    with open(subfile_2, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, scratch, command_2, clean_scratch))
            
    return subfile_1, subfile_2
    

def update_irc_status(target, job_id, direction):
    qm_collection = db['qm_calculate_center']
    reg_query = {"path":target}
    irc_status = 'irc_{}_status'.format(str(direction))
    irc_jobid = 'irc_{}_jobid'.format(str(direction))
    update_field = {irc_status:"job_launched", irc_jobid:job_id}
    qm_collection.update_one(reg_query, {"$unset": {'irc_status':""}, "$set": update_field}, True)

"""
Submmit opt job which is from irc
1. select unrun job
2. push unrun job to qchem
3. update status "job_launched"
"""

def select_irc_opt_target(direction = 'forward'):
    
    qm_collection = db['qm_calculate_center']
    irc_opt_status = 'opt_{}_status'.format(direction)
    reg_query = {irc_opt_status:"job_unrun"}
    targets = list(qm_collection.find(reg_query))
    selected_targets = [target['path'] for target in targets]
    return selected_targets

def launch_irc_opt_jobs():
    
    targets = select_irc_opt_target(direction = 'forward')
    for target in targets:
        IRC_dir_path = path.join(target, 'IRC/')
        os.chdir(IRC_dir_path)
        subfile = create_irc_opt_sub_file(IRC_dir_path, direction = 'forward')
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_opt_status(target, job_id, direction = 'forward')
        
    targets = select_irc_opt_target(direction = 'reverse')
    for target in targets:
        IRC_dir_path = path.join(target, 'IRC/')
        os.chdir(IRC_dir_path)
        subfile = create_irc_opt_sub_file(IRC_dir_path, direction = 'reverse')
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_irc_opt_status(target, job_id, direction = 'reverse')


def update_irc_opt_status(target, job_id, direction):
    qm_collection = db['qm_calculate_center']
    reg_query = {"path":target}
    irc_opt_status = 'opt_{}_status'.format(direction)
    irc_opt_jobid = 'irc_{}_opt_jobid'.format(str(direction))
    update_field = {irc_opt_status:"opt_job_launched", irc_opt_jobid:job_id}
    qm_collection.update_one(reg_query, {"$set": update_field}, True)
    

def create_irc_opt_sub_file(irc_path, direction = 'forward', ncpus = 8, mpiprocs = 1, ompthreads = 8):
    job_name = 'irc_{}_opt.job'.format(direction)
    subfile = path.join(irc_path, job_name)
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(irc_path)
    nes1 = 'module load qchem'
    nes2 = 'export QCSCRATCH=/tmp/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    inputname = '{}_opt.in'.format(direction)
    outputname = '{}_opt.out'.format(direction)
    nes4 = 'qchem -nt {} {} {}'.format(ncpus, inputname, outputname)
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5))
    return subfile

#launch_energy_jobs()
#launch_ssm_jobs()
#launch_low_opt_jobs()
#launch_opt_jobs()
#launch_ts_jobs()
#launch_irc_jobs()
#launch_irc_opt_jobs()