from connect import db
import subprocess
import os
from os import path
import shutil
import sys
sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'script'))
from openbabel import pybel

# Manual run 0th generations

def select_ard_target():
    # when lower generation finished return target else return []
    qm_collection = db['qm_calculate_center']
    reg_query = {'ard_status':
                    {"$in": 
                        ["job_unrun"]
                    }
                }
    targets = list(qm_collection.find(reg_query))
    return targets

def launch_ard_jobs():
    
    qm_collection = db['qm_calculate_center']
    pool_collection = db['pool']
    status_collection = db['status']

    if pool_collection.estimated_document_count() == 0:
        print('The ard not start')
        print('Starting ARD network exploring')
        script_path = path.join(path.dirname(path.dirname(path.abspath(__file__))), 'script')
        if os.path.exists(script_path):
            os.chdir(script_path)
        subfile = create_ard_sub_file(script_path, script_path, 1, 'reactant.xyz')
        # first reactant need to add to pool
        initial_reactant = next(pybel.readfile('xyz', path.join(script_path, 'reactant.xyz')))
        initial_reactant_inchi_key = initial_reactant.write('inchiKey').strip()
        pool_collection.insert_one({'reactant_inchi_key':initial_reactant_inchi_key})
        cmd = 'qsub {}'.format(subfile)
        process = subprocess.Popen([cmd],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, shell = True)
        stdout, stderr = process.communicate()
        # get job id from stdout, e.g., "106849.h81"
        job_id = stdout.decode().replace("\n", "")
        # update status job_launched
        print('ARD had launched')
        print('jobid is {}'.format(job_id))
        status_collection.insert_one({'status':'ARD had launched'})
    else:
        targets = select_ard_target()
        for target in targets:
            dir_path, gen_num, ard_ssm_equal, irc_equal = target['path'], target['generations'], target['ard_ssm_equal'], target['irc_equal']
            script_path = path.join(path.dirname(path.dirname(dir_path)), 'script')
            os.chdir(dir_path)

            if irc_equal == 'forward equal to reactant but reverse does not equal to product':
                irc_path = path.join(dir_path, 'IRC/')
                next_reactant = path.join(irc_path, 'reverse.xyz')
                subfile = create_ard_sub_file(dir_path, script_path, gen_num + 1, next_reactant)
            elif irc_equal == 'reverse equal to reactant but forward does not equal to product':
                irc_path = path.join(dir_path, 'IRC/')
                next_reactant = path.join(irc_path, 'forward.xyz')
                subfile = create_ard_sub_file(dir_path, script_path, gen_num + 1, next_reactant)
            else:
                if ard_ssm_equal == 'not_equal':
                    next_reactant = 'ssm_product.xyz'
                    subfile = create_ard_sub_file(dir_path, script_path, gen_num + 1, next_reactant)
                else:
                    next_reactant = 'product.xyz'
                    subfile = create_ard_sub_file(dir_path, script_path, gen_num + 1, next_reactant)

            cmd = 'qsub {}'.format(subfile)
            process = subprocess.Popen([cmd],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, shell = True)
            stdout, stderr = process.communicate()
            # get job id from stdout, e.g., "106849.h81"
            job_id = stdout.decode().replace("\n", "")
            # update status job_launched
            update_ard_status(target[0], job_id)

def create_ard_sub_file(dir_path, script_path, gen_num, next_reactant, ncpus = 8, mpiprocs = 1, ompthreads = 8):
    subfile = path.join(dir_path, 'ard.job')
    product_xyz_path = path.join(dir_path, next_reactant)
    ard_path = path.join(script_path, 'ard.py')
    input_path = path.join(script_path, 'input.txt')
    bonds_path = path.join(dir_path, 'bonds.txt')
    constraint = path.join(dir_path, 'constraint.txt')
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}:ompthreads={}\n#PBS -q workq\n#PBS -j oe'.format(ncpus, mpiprocs, ompthreads)
    target_path = 'cd {}'.format(script_path)
    nes1 = 'source ~/.bashrc'
    nes2 = 'conda activate rmg3'
    command = 'python {} {} {} -bonds {} -constraint {} -generations {}'.format(ard_path, input_path, product_xyz_path, bonds_path, constraint, gen_num)    
    deactivate = 'conda deactivate'
    
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}'.format(shell, pbs_setting, target_path, nes1, nes2, command, deactivate))
        
    return subfile

def update_ard_status(target, job_id):
    qm_collection = db['qm_calculate_center']
    reg_query = {"path":target}
    update_field = {"ard_status":"job_launched", "ard_jobid":job_id}
    qm_collection.update_one(reg_query, {"$set": update_field}, True)
    
launch_ard_jobs()