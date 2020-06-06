# connect to database

# search entries with status "job_launched" or "job_running"
# check the job pbs-status, e.g., qstat -f "job_id"
# 1. Invalid job id specified: means the job is off queue
#    1.1 check job folder if input.log file is there
#			Yes: job already finished
#			No:  job is aborted before even calculation start: job_aborted
#    1.2 if finished, check "Normal termination" key words should appear twice
#			Yes: job is converged
#			No:  job fails convergence (possibly need more wall-time): job_failed_convergence
#    1.3 if converged, check isomorphism between input molecule and output molecule
#			Yes: job is success: job_success
#			No:  job fails isomorphism (possibly need tweak on initial structure): job_failed_isomorphism
# 2. output string starts with e.g, "JobId=5037088"
#    2.1 check JobState
#		RUNNING: running: job_running
#		else:   still job_launched 

from connect import db
import subprocess
import os
from os import path
import sys
import rmgpy.molecule
from rmgpy.molecule.converter import from_ob_mol
import pybel
from collections import Counter

def select_energy_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    reg_query = {"energy_status":
                    {"$in":
                        ["job_launched", "job_running", "job_queueing"]
                    }
                }
    targets = list(qm_collection.find(reg_query))

    return targets

def check_energy_status(job_id):
    """
    This method checks pbs status of a job given job_id
    Returns off_queue or job_launched or job_running
    """

    commands = ['qstat', '-f', job_id]
    process = subprocess.Popen(commands,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if "Unknown Job Id" in stderr.decode():
        return "off_queue"

    # in pbs stdout is byte, so we need to decode it at first.
    stdout = stdout.decode().strip().split()
    idx = stdout.index('job_state')
    if stdout[idx+2] == 'R':
        return "job_running"
    elif stdout[idx+2] == 'Q':
        return 'job_queueing'
    else:
        return "job_launched"
    
def check_energy_content_status(rxn_path):
    
    energy_dir = path.join(rxn_path, "Energy")
    energy_path = path.join(energy_dir, "reactant_energy.out")
    if not path.exists(energy_path):
        return "job_aborted"
    else:
        try:
            with open(energy_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.strip().startswith('SCF   energy in the final basis set'):
                        SCF_Energy = line.split()[-1]
                        energy = float(SCF_Energy)
            return "job_success", energy
        except:
            return "job_fail", 0

def check_energy_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_energy_target()
    qm_collection = db['qm_calculate_center']
    
    # 2. check the job pbs_status
    for target in targets:
        job_id = target['energy_jobid']
        new_status = check_energy_status(job_id)
        if new_status == "off_queue":
                # 3. check job content
                new_status, energy = check_energy_content_status(target['path'])

                # 4. check with original status which
                # should be job_launched or job_running
                # if any difference update status
                orig_status = target['energy_status']
                if orig_status != new_status:
                    update_field = {
                                       'energy_status': new_status,
                                       'reactant_scf_energy':energy
                                   }
                    qm_collection.update_one(target, {"$set": update_field}, True)

def select_ssm_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    reg_query = {"ssm_status":
                    {"$in": 
                        ["job_launched", "job_running", "job_queueing"] 
                    }
                }
    targets = list(qm_collection.find(reg_query))

    return targets

def check_ssm_status(job_id):
    """
    This method checks pbs status of a job given job_id
    Returns off_queue or job_launched or job_running
    """
    commands = ['qstat', '-f', job_id]
    process = subprocess.Popen(commands,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if "Unknown Job Id" in stderr.decode():
        return "off_queue"

    # in pbs stdout is byte, so we need to decode it at first.
    stdout = stdout.decode().strip().split()
    idx = stdout.index('job_state')
    if stdout[idx+2] == 'R':
        return "job_running"
    elif stdout[idx+2] == 'Q':
        return "job_queueing"
    else:
        return "job_launched"
    
def check_ssm_content_status(target_path):

    ssm_path = path.join(target_path, 'SSM')
    ssm_pic_path = path.join(ssm_path, '0000_string.png')
    ssm_tsnode_path = path.join(ssm_path, 'TSnode.xyz')
    if not path.exists(ssm_pic_path) or not path.exists(ssm_tsnode_path):
        return "job_fail"
    else:
        if check_ssm(ssm_path) == 'Exiting early':
            return "Exiting early"
        elif check_ssm(ssm_path) == 'total dissociation':
            return 'total dissociation'
        else:
            generate_ssm_product_xyz(ssm_path)
            return "job_success"

def check_ssm(ssm_path):
    status_path = path.join(ssm_path, 'status.log')
    
    with open(status_path, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 1024, 0), 0)  # Read last 1 kB of file
        lines = f.readlines()
        
    for idx, i in enumerate(lines):
        if i.startswith('Finished GSM!'):
            break
    
    if lines[idx-1] == 'Exiting early\n':
        return 'Exiting early'
    elif lines[idx-2] == 'total dissociation\n':
        return 'total dissociation'
    else:
        return 0

def ard_prod_and_ssm_prod_checker(rxn_dir):
    # Use ssm product xyz to check whether ssm prod equal to ard product
    # If equal, insert the inchi key into products pool
    # If not equal, use ssm product as the product and insert inchi key into products pool
    # Next generation use the ssm product to generate
    qm_collection = db['qm_calculate_center']
    
    ard_prod_path = path.join(rxn_dir, 'product.xyz')
    ssm_prod_path = path.join(rxn_dir, 'ssm_product.xyz')
    OBMol_1 = readXYZ(ssm_prod_path)
    rmg_mol_1 = toRMGmol(OBMol_1)
    OBMol_2 = readXYZ(ard_prod_path)
    rmg_mol_2 = toRMGmol(OBMol_2)
    if rmg_mol_1.to_inchi_key() != rmg_mol_2.to_inchi_key():
        num = len(list(qm_collection.find({'product_inchi_key':rmg_mol_1.to_inchi_key()})))
        targets = list(qm_collection.find({'path':rxn_dir}))
        for i in targets:
            num += 1
            name = '{}_{}'.format(rmg_mol_1.to_inchi_key(), num)
            new_path = path.join(path.dirname(i['path']), name)
            os.rename(rxn_dir, new_path)
            prod_smi = OBMol_1.write('can').strip()
            update_field = {'product_inchi_key':rmg_mol_1.to_inchi_key(), 
                            'initial_dir_name':rxn_dir, 
                            'path':new_path, 
                            'ssm_status': 'job_success', 
                            "ts_status":"job_unrun", 
                            "energy_status":"job_unrun", 
                            'ard_ssm_equal':'not_equal',
                            'Product SMILES': prod_smi}
            qm_collection.update_one(i, {"$set": update_field}, True)
        return 'not_equal'
    else:
        return 'equal'

def toRMGmol(pyMol):
    rmg_mol = from_ob_mol(rmgpy.molecule.molecule.Molecule(), pyMol.OBMol)
    return rmg_mol

def readXYZ(xyz):
    mol = next(pybel.readfile('xyz', xyz))
    return mol

def check_ssm_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_ssm_target()

    qm_collection = db['qm_calculate_center']
    
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ssm_jobid']
        # 2. check the job pbs status
        new_status = check_ssm_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_ssm_content_status(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ssm_status']
        if orig_status != new_status:

            if new_status == 'job_success':
                equal = ard_prod_and_ssm_prod_checker(target['path'])
                if equal == 'equal':
                    update_field = {
                                    'ssm_status': new_status, "ts_status":"job_unrun", 'ard_ssm_equal':equal
                                }
                    qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                'ssm_status': new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)
	

def generate_ssm_product_xyz(target_path):
    opt_file = []
    path_list=os.listdir(target_path)
    path_list.sort()
    for filename in path_list:
        if filename.startswith('opt'):
            opt_file.append(filename)
    product_xyz = path.join(target_path, opt_file[-1])
    with open(product_xyz, 'r') as f:
        lines = f.readlines()
        for i in reversed(lines):
            a = i.split()
            if len(a) == 1:
                idx = lines.index(i)
                break
    parent_ssm_product_path  = path.join(path.dirname(target_path), 'ssm_product.xyz')
    with open(parent_ssm_product_path,'w') as q:
        q.write('{}\n{}'.format(lines[idx-1], ''.join(lines[idx+1:])))


def select_ts_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    reg_query = {"ts_status":
                    {"$in": 
                        ["job_launched", "job_running", "job_queueing"] 
                    }
                }
    targets = list(qm_collection.find(reg_query))

    return targets

def check_ts_status(job_id):
    """
    This method checks pbs status of a job given job_id
    Returns off_queue or job_launched or job_running
    """
    commands = ['qstat', '-f', job_id]
    process = subprocess.Popen(commands,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if "Unknown Job Id" in stderr.decode():
        return "off_queue"

    # in pbs stdout is byte, so we need to decode it at first.
    stdout = stdout.decode().strip().split()
    idx = stdout.index('job_state')
    if stdout[idx+2] == 'R':
        return "job_running"
    elif stdout[idx+2] == 'Q':
        return "job_queueing"
    else:
        return "job_launched"
    
def check_ts_content_status(target_path):

    ts_dir_path = path.join(target_path, 'TS')
    ts_out_path = path.join(ts_dir_path, 'ts.out')
    ts_geo_path = path.join(ts_dir_path, 'ts_geo.xyz')

    with open(ts_out_path, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 6144, 0), 0)  # Read last 6 kB of file
        lines = f.readlines()

    for idx, i in enumerate(lines):
        if i.startswith(' **  OPTIMIZATION CONVERGED  **\n'):
            break
    for idx2, i in enumerate(lines):
        if i.startswith('Z-matrix Print:\n'):
            break

    
    if lines[idx] != ' **  OPTIMIZATION CONVERGED  **\n':
        # here 0 is just let ts_energy do not error
        return "job_fail", 0
    else:
        # generate geometry xyz file
        geo = []
        for i in lines[idx+5:idx2-1]:
            atom = i.split()[1:]
            geo.append('  '.join(atom))
            
        with open(ts_geo_path, 'w') as f:
            f.write(str(len(geo)))
            f.write('\n\n')
            f.write('\n'.join(geo))
        # get ts energy
        ts_energy = lines[idx-4].split()[-1]

        return "job_success", float(ts_energy)
    
def insert_exact_rxn(reactant_inchi_key, product_inchi_key, reactant_smi, product_smi, path, generations, barrier):
    reactions_collection = db['reactions']

    reactions_collection.insert_one({
                        'reaction':[reactant_inchi_key, product_inchi_key],
                        'reactant_smi':reactant_smi,
                        'product_smi':product_smi,
                        'path':path,
                        'generations':generations,
                        'unique': 'new one',
                        'barrier_energy':barrier})

def check_ts_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_ts_target()

    qm_collection = db['qm_calculate_center']
    pool_collection = db['pool']
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ts_jobid']
        # 2. check the job pbs status
        new_status = check_ts_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, ts_energy= check_ts_content_status(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ts_status']
        next_gen_num = int(target['generations']) + 1
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                                'ts_status': new_status, 'ts_energy':ts_energy, 'irc_status':'job_unrun', 'next_gen_num':next_gen_num, "energy_status":"job_unrun"
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
                pool_collection.insert_one({'reactant_inchi_key':target['product_inchi_key']})
            else:
                update_field = {
                                'ts_status': new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)

def select_ts_barrier_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    query = {'$and': 
                    [
                    { "ts_status":
                        {"$in":
                        ['job_success']}},
                    {'energy_status':
                        {'$in':
                            ['job_success']}}
                    ]
                }
    targets = list(qm_collection.find(query))
    num =[]
    for idx, i in enumerate(targets):
        try:
            barrier_energy = i['barrier_energy']
        except:
            num.append(idx)
            
    targets = [targets[i] for i in num]
    
    return targets
    

def check_ts_barrier_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_ts_barrier_target()
    qm_collection = db['qm_calculate_center']
    
    for target in targets:
        reactant_energy = target['reactant_scf_energy']
        ts_energy = target['ts_energy']
        barrier_energy = (float(ts_energy) - float(reactant_energy)) * 627.5095
        update_field = {'barrier_energy':barrier_energy}
        qm_collection.update_one(target, {"$set": update_field}, True)

def check_ard_barrier_jobs(generations):
    # Select the lowest ts barrier reaction inserts to reaction collection and changes the ard status to job unrun.
    # check ard will start when ssm success = ts fail + ts success to make sure all reaction done.
    qm_collection = db['qm_calculate_center']
    ssm_query = {'$and':
                    [{"ssm_status":
                    {"$in":
                        ['job_success', 'job_running']
                    }
                },
                {
                    'generations':generations
                }]}
    ssm_success_number = len(list(qm_collection.find(ssm_query)))
    ts_query = {'$and':
                [{'$or': 
                    [
                    {'ts_status':
                        {"$in":
                            ['job_success', 'job_fail']}},
                    {'energy_status':
                        {'$in':
                            ['job_success', 'job_fail']}}
                    ]
                },
                {
                    'generations':generations
                }]}
    ts_fail_and_success_number = len(list(qm_collection.find(ts_query)))
    ard_query = {'ard_status':'job_unrun'}
    ard_number = len(list(qm_collection.find(ard_query)))
    if ssm_success_number == ts_fail_and_success_number and ard_number == 0:
        check_ts_barrier_jobs()
        
        targets = list(qm_collection.aggregate(
            [
                {"$match": {"$and":[
                    {
                        "ts_status": {"$in":['job_success']},
                        "energy_status": {"$in":['job_success']}
                    }
                ]}},
                {"$group":{"_id":"$reaction", 'barrier_energy':{'$min':"$barrier_energy"}}}
            ]
        ))
        for i in targets:
            reaction = i['_id']
            barrier = i['barrier_energy']
            query = {
                '$and':[
                    {
                        'reaction':reaction    
                    },
                    {
                        'barrier_energy':barrier
                    }
                ]
            }
            target = list(qm_collection.find(query)) # This is only have one result though it is a list (to visualize <pymongo.cursor.Cursor object>)
            qm_collection.update_one(target[0], {"$set": {'ard_status':'job_unrun'}}, True)
            insert_exact_rxn(target[0]['reactant_inchi_key'],
                            target[0]['product_inchi_key'],
                            target[0]['Reactant SMILES'],
                            target[0]['Product SMILES'],
                            target[0]['path'],
                            target[0]['generations'],
                            barrier)
         
check_energy_jobs()
check_ssm_jobs()
check_ts_jobs()

qm_collection = db['qm_calculate_center']
max_gen = qm_collection.find_one(sort=[("generations", -1)])
max_gen = max_gen['generations']

check_ard_barrier_jobs(max_gen)