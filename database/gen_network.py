from connect import db
import subprocess
import os
from os import path
import shutil
import sys
sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'script'))

# Manual run 0th generations

def select_ard_target():
    collect = db['reactions']
    reg_query = {"ts_status":
                    {"$in": 
                        ["job_success"] 
                    }
                }
    targets = list(collect.find(reg_query))

    selected_targets = []
    gens = []
    for target in targets:
        dir_path = target['path']
        selected_targets.append(dir_path)
        gen = target['next_gen_num']
        gens.append(gen)
    zipped = zip(selected_targets, gens)
    return zipped

def launch_ard_jobs():
    targets = select_ard_target()
    script_path = os.path.join(os.path.abspath(os.pardir), 'script')
    for target in list(targets):
        dir_path, gen_num = target[0], target[1]
        
        if os.path.exists(script_path):
            os.chdir(script_path)
            
        subfile = create_ard_sub_file(dir_path, script_path, gen_num)
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        update_ard_status(target, job_id)

def create_ard_sub_file(dir_path, script_path, gen_num, ncpus = 1, mpiprocs = 1):
    subfile = path.join(dir_path, 'ard.job')
    product_xyz_path = path.join(dir_path, 'product.xyz')
    
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}'.format(ncpus, mpiprocs)
    target_path = 'cd {}'.format(script_path)
    nes1 = 'source ~/.bashrc'
    nes2 = 'conda activate rmg3'
    command = 'python ard.py {} {}'.format(product_xyz_path, gen_num)    
    deactivate = 'conda deactivate'
    
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, command, deactivate))
        
    return subfile

def update_ard_status(target, job_id):
    collect = db['reactions']
    reg_query = {"path":target}
    update_field = {"ard_status":"job_launched", "ard_jobid":job_id}
    collect.update_one(reg_query, {"$set": update_field}, True)