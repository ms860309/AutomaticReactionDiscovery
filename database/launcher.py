from connect import db
import subprocess
import os
from os import path
import shutil
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
        Energy_dir_path = path.join(target, 'Energy')
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
        
def create_energy_sub_file(path, Energy_dir_path, ncpus = 1, mpiprocs = 1, ompthreads = 1):
    subfile = os.path.join(Energy_dir_path, 'cal_E.job')
    
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(Energy_dir_path)
    nes1 = 'module load qchem'
    nes2 = 'export QCSCRATCH=/tmp/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4 = 'qchem -nt 1 ../reactant_energy.in reactant_energy.out'
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5))
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
        
def create_ssm_sub_file(dir_path, SSM_dir_path, ncpus = 1, mpiprocs = 1, ompthreads = 1):
    subfile = path.join(SSM_dir_path, 'cal_ssm.job')
    xyz_file = path.join(dir_path, 'reactant.xyz')
    isomers = path.join(dir_path, 'add_bonds.txt')
    lot_inp_file = path.join(path.join(path.dirname(path.dirname(dir_path)), 'submmit_required'), 'qstart')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(SSM_dir_path)
    nes1 = 'module load qchem'
    # activate conda env is necessary because gsm install on the environment
    nes2 = 'source ~/.bashrc\nconda activate rmg3'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    coord_type = 'DLC'
    command = 'gsm -xyzfile {} -mode SE_GSM -package QChem -isomers {} -lot_inp_file {} -coordinate_type {} -max_gsm_iters 100 -max_opt_steps 30 -CONV_TOL 0.0005 -ADD_NODE_TOL 0.01 -DQMAG_MAX 0.8 -num_nodes 30 -conv_Ediff 300 >status.log'.format(xyz_file, isomers, lot_inp_file, coord_type)    
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
        
def create_ts_sub_file(SSM_dir_path, TS_dir_path, ncpus = 1, mpiprocs = 1, ompthreads = 1):
    tsnode_path = path.join(SSM_dir_path, 'TSnode.xyz')
    ts_input_file = path.join(TS_dir_path, 'ts.in')
    ts_output_file = path.join(TS_dir_path, 'ts.out')
    subfile = path.join(TS_dir_path, 'cal_ts.job')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(TS_dir_path)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command = 'qchem -nt {} {} {}'.format(ncpus, ts_input_file, ts_output_file)
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(tsnode_path, 'r') as f1:
        lines = f1.readlines()
        with open(ts_input_file, 'w') as f2:
            # geometry
            f2.write('$molecule\n')
            f2.write('0 1\n')
            for line in lines[2:]:
                f2.write(line)
            f2.write('$end\n')
            # jobtype
            f2.write('$rem\n')
            f2.write('JOBTYPE FREQ\n')
            f2.write('METHOD B97-D3\n')
            f2.write('DFT_D D3_BJ\n')
            f2.write('BASIS def2-mSVP\n')
            f2.write('SCF_ALGORITHM DIIS\n')
            f2.write('MAX_SCF_CYCLES 150\n')
            f2.write('SCF_CONVERGENCE 8\n')
            f2.write('SYM_IGNORE TRUE\n')
            f2.write('SYMMETRY FALSE\n')
            f2.write('WAVEFUNCTION_ANALYSIS FALSE\n')
            f2.write('$end\n')
            f2.write('\n@@@\n\n')
            f2.write('$molecule\n')
            f2.write('read\n')
            f2.write('$end\n')
            f2.write('$rem\n')
            f2.write('JOBTYPE TS\n')
            f2.write('METHOD B97-D3\n')
            f2.write('DFT_D D3_BJ\n')
            f2.write('BASIS def2-mSVP\n')
            f2.write('SCF_GUESS READ\n')
            f2.write('GEOM_OPT_HESSIAN READ\n')
            f2.write('SCF_ALGORITHM DIIS\n')
            f2.write('MAX_SCF_CYCLES 150\n')
            f2.write('SCF_CONVERGENCE 8\n')
            f2.write('SYM_IGNORE TRUE\n')
            f2.write('SYMMETRY FALSE\n')
            f2.write('GEOM_OPT_MAX_CYCLES 150\n')
            f2.write('GEOM_OPT_TOL_GRADIENT 100\n')
            f2.write('GEOM_OPT_TOL_DISPLACEMENT 400\n')
            f2.write('GEOM_OPT_TOL_ENERGY 33\n')
            f2.write('WAVEFUNCTION_ANALYSIS FALSE\n')
            f2.write('$end\n')
            
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
        
def create_irc_sub_file(TS_dir_path, IRC_dir_path, ncpus = 1, mpiprocs = 1, ompthreads = 1):
    ts_geo_path = path.join(TS_dir_path, 'ts_geo.xyz')
    
    irc_forward_input_file = path.join(IRC_dir_path, 'irc_forward.in')
    irc_reverse_input_file = path.join(IRC_dir_path, 'irc_reverse.in')
    irc_forward_output_file = path.join(IRC_dir_path, 'irc_forward.out')
    irc_reverse_output_file = path.join(IRC_dir_path, 'irc_reverse.out')
    
    subfile_1 = path.join(IRC_dir_path, 'cal_irc_forward.job')
    subfile_2 = path.join(IRC_dir_path, 'cal_irc_reverse.job')

    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(IRC_dir_path)
    nes1 = 'module load qchem'
    scratch = 'export QCSCRATCH=/tmp/$PBS_JOBID\nmkdir -p $QCSCRATCH\n'
    command_1 = 'qchem -nt {} {} {}'.format(ncpus, irc_forward_input_file, irc_forward_output_file)
    command_2 = 'qchem -nt {} {} {}'.format(ncpus, irc_reverse_input_file, irc_reverse_output_file)
    clean_scratch = 'rm -r $QCSCRATCH'
    
    with open(ts_geo_path, 'r') as f1:
        lines = f1.readlines()
        with open(irc_forward_input_file, 'w') as f2:
            # geometry
            f2.write('$molecule\n')
            f2.write('0 1\n')
            for line in lines[2:]:
                f2.write(line)
            f2.write('\n$end\n')
            # jobtype
            f2.write('\n$rem\n')
            f2.write('JOBTYPE FREQ\n')
            f2.write('METHOD B97-D3\n')
            f2.write('DFT_D D3_BJ\n')
            f2.write('BASIS def2-mSVP\n')
            f2.write('SCF_ALGORITHM DIIS\n')
            f2.write('MAX_SCF_CYCLES 150\n')
            f2.write('SCF_CONVERGENCE 8\n')
            f2.write('SYM_IGNORE TRUE\n')
            f2.write('SYMMETRY FALSE\n')
            f2.write('WAVEFUNCTION_ANALYSIS FALSE\n')
            f2.write('$end\n')
            f2.write('\n@@@\n\n')
            f2.write('$molecule\n\n')
            f2.write('read\n')
            f2.write('$end\n')
            f2.write('$rem\n')
            f2.write('JOBTYPE RPATH\n')
            f2.write('METHOD B97-D3\n')
            f2.write('DFT_D D3_BJ\n')
            f2.write('BASIS def2-mSVP\n')
            f2.write('SCF_GUESS READ\n')
            f2.write('rpath_direction 1\n')
            f2.write('rpath_max_cycles 60\n')
            f2.write('rpath_max_stepsize 80\n')
            f2.write('$end\n')
            
    with open(ts_geo_path, 'r') as f1:
        lines = f1.readlines()
        with open(irc_reverse_input_file, 'w') as f2:
            # geometry
            f2.write('$molecule\n')
            f2.write('0 1\n')
            for line in lines[2:]:
                f2.write(line)
            f2.write('\n$end\n')
            # jobtype
            f2.write('\n$rem\n')
            f2.write('JOBTYPE FREQ\n')
            f2.write('METHOD B97-D3\n')
            f2.write('DFT_D D3_BJ\n')
            f2.write('BASIS def2-mSVP\n')
            f2.write('SCF_ALGORITHM DIIS\n')
            f2.write('MAX_SCF_CYCLES 150\n')
            f2.write('SCF_CONVERGENCE 8\n')
            f2.write('SYM_IGNORE TRUE\n')
            f2.write('SYMMETRY FALSE\n')
            f2.write('WAVEFUNCTION_ANALYSIS FALSE\n')
            f2.write('$end\n')
            f2.write('\n@@@\n\n')
            f2.write('$molecule\n\n')
            f2.write('read\n')
            f2.write('$end\n')
            f2.write('$rem\n')
            f2.write('JOBTYPE RPATH\n')
            f2.write('METHOD B97-D3\n')
            f2.write('DFT_D D3_BJ\n')
            f2.write('BASIS def2-mSVP\n')
            f2.write('SCF_GUESS READ\n')
            f2.write('rpath_direction -1\n')
            f2.write('rpath_max_cycles 60\n')
            f2.write('rpath_max_stepsize 80\n')
            f2.write('$end\n')
            
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
    
launch_energy_jobs()
launch_ssm_jobs()
launch_ts_jobs()
launch_irc_jobs()