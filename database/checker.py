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
from openbabel import pybel

"""
Energy check.
"""
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

def check_energy_job_status(job_id):
    """
    This method checks pbs status of a job given job_id
    Returns off_queue or job_launched or job_running
    """

    commands = ['qstat', '-f', job_id]
    process = subprocess.Popen(commands,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if 'Unknown Job Id' in stderr.decode():
        return 'off_queue'

    # in pbs stdout is byte, so we need to decode it at first.
    stdout = stdout.decode().strip().split()
    idx = stdout.index('job_state')
    if stdout[idx+2] == 'R':
        return 'job_running'
    elif stdout[idx+2] == 'Q':
        return 'job_queueing'
    else:
        return 'job_launched'

def check_energy_content(rxn_path):
    
    energy_dir = path.join(rxn_path, "ENERGY")
    energy_path = path.join(energy_dir, "energy.out")
    if not path.exists(energy_path):
        return "job_aborted"
    else:
        try:
            with open(energy_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.strip().startswith('SCF   energy in the final basis set'):
                        energy = float(line.split()[-1])
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
        new_status = check_energy_job_status(job_id)
        if new_status == "off_queue":
                # 3. check job content
                new_status, energy = check_energy_content(target['path'])

                # 4. check with original status which
                # should be job_launched or job_running
                # if any difference update status
                orig_status = target['energy_status']
                if orig_status != new_status:
                    update_field = {
                                       'energy_status': new_status,
                                       'reactant_scf_energy':energy, 
                                       'barrier': 'None'
                                   }
                    qm_collection.update_one(target, {"$set": update_field}, True)

"""
SSM check.
"""
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

def check_ssm_job_status(job_id):
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
    
def check_ssm_content(target_path):

    ssm_path = path.join(target_path, 'SSM')
    ssm_pic_path = path.join(ssm_path, '0000_string.png')
    ssm_tsnode_path = path.join(ssm_path, 'TSnode.xyz')

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

    if not path.exists(ssm_pic_path) or not path.exists(ssm_tsnode_path):
        return "job_fail"
    elif check_ssm(ssm_path) == 'Exiting early':
        return "Exiting early"
    elif check_ssm(ssm_path) == 'total dissociation':
        return 'total dissociation'
    else:
        generate_ssm_product_xyz(ssm_path)
        return "job_success"

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

def ard_prod_and_ssm_prod_checker(rxn_dir):
    # Use ssm product xyz to check whether ssm prod equal to ard product
    # If equal, insert the inchi key into products pool
    # If not equal, use ssm product as the product and insert inchi key into products pool
    # Next generation use the ssm product to generate
    qm_collection = db['qm_calculate_center']
    
    ard_prod_path = path.join(rxn_dir, 'product.xyz')
    ssm_prod_path = path.join(rxn_dir, 'ssm_product.xyz')
    pyMol_1 = xyz_to_pyMol(ssm_prod_path)
    pyMol_2 = xyz_to_pyMol(ard_prod_path)
    if pyMol_1.write('inchiKey').strip() != pyMol_2.write('inchiKey').strip():
        num = len(list(qm_collection.find({'product_inchi_key':pyMol_1.write('inchiKey').strip()})))
        targets = list(qm_collection.find({'path':rxn_dir}))
        for i in targets:
            if i['reactant_inchi_key'] == pyMol_1.write('inchiKey').strip():
                dirname = dir_check(path.dirname(i['path']), pyMol_1.write('inchiKey').strip(), num + 1)
                new_path = path.join(path.dirname(i['path']), dirname)
                os.rename(rxn_dir, new_path)
                prod_smi = pyMol_1.write('can').split()[0]
                update_field = {'product_inchi_key':pyMol_1.write('inchiKey').strip(), 
                                'initial_dir_name':rxn_dir, 
                                'path':new_path, 
                                'ssm_status': 'job_success',                                # Though the reactant equal to product but TS maybe fine. The SSM somehow do not use constrained optimize in product
                                'ard_ssm_equal':'ssm reactant equal to product',            # This is a special case but sometimes it will happen.
                                'Product SMILES': prod_smi}
                qm_collection.update_one(i, {"$set": update_field}, True)
            else:
                dirname = dir_check(path.dirname(i['path']), pyMol_1.write('inchiKey').strip(), num + 1)
                new_path = path.join(path.dirname(i['path']), dirname)
                os.rename(rxn_dir, new_path)
                prod_smi = pyMol_1.write('can').split()[0]
                update_field = {'product_inchi_key':pyMol_1.write('inchiKey').strip(), 
                                'initial_dir_name':rxn_dir, 
                                'path':new_path, 
                                'ssm_status': 'job_success', 
                                "ts_status":"job_unrun", 
                                'ard_ssm_equal':'not_equal',
                                'Product SMILES': prod_smi}
                qm_collection.update_one(i, {"$set": update_field}, True)
        return 'not_equal'
    else:
        return 'equal'
    
def dir_check(subdir, b_dirname, num):
    """
    When parallely run job, the dir is constructed but data is not on database yet
    """
    check = False
    number = num
    while check == False:
        new_name = '{}_{}'.format(b_dirname, number)
        if os.path.exists(os.path.join(subdir, new_name)):
            number += 1
        else:
            check = True
    
    return new_name

def xyz_to_pyMol(xyz):
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
        new_status = check_ssm_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_ssm_content(target['path'])

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
                # if not equal the 'ts_status': 'job_unrun' is update on ard_prod_and_ssm_prod_checker fuction
            else:
                update_field = {
                                'ssm_status': new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)

"""
TS check.
"""
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

def check_ts_job_status(job_id):
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
    
def check_ts_content(target_path):

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
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ts_jobid']
        # 2. check the job pbs status
        new_status = check_ts_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, ts_energy= check_ts_content(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ts_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                                'ts_status': new_status, 'ts_energy':ts_energy, 'irc_status':'job_unrun'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                'ts_status': new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)

"""
IRC check.
"""
def select_irc_target(direction = 'forward'):
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    irc_status = 'irc_{}_status'.format(direction)
    query = {irc_status:
                    {"$in": 
                        ["job_launched", "job_running", "job_queueing"] 
                    }
                }
    targets = list(qm_collection.find(query))

    return targets

def check_irc_status(job_id):
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
    
def check_irc_content(target_path, direction = 'forward'):
    reactant_path = path.join(target_path, 'reactant.xyz')
    irc_path = path.join(target_path, 'IRC/')
    base_dir_path = path.join(path.dirname(path.dirname(path.dirname(path.dirname(irc_path)))), 'config')
    opt_lot = path.join(base_dir_path, 'opt.lot')
    opt_name = '{}_opt.in'.format(direction)
    opt_in = path.join(irc_path, opt_name)
    irc_output = path.join(irc_path, 'irc_{}.out'.format(direction))

    with open(opt_lot) as f:
        config = [line.strip() for line in f]
    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])
    with open(irc_output, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 10240, 0), 0)  # Read last 10 kB of file
        lines = f.readlines()

    if lines[-2] == ' IRC backup failure\n' or lines[-2] == ' IRC failed final bisector step\n':

        with open(irc_output, 'r') as f:
            full_lines = f.readlines()

        # Sometimes irc success in forward(reverse) but fail in reverse(forward).
        # We wan't to catch the final structure to optimize. 
        # But if "convergence criterion reached" in first direction fail in second that will cause reactant equal to product.

        if '  IRC -- convergence criterion reached.\n' in full_lines:
            for idx, i in enumerate(full_lines):
                if i.startswith('  IRC -- convergence criterion reached.\n'):
                    break
            full_lines = full_lines[:idx]
            for idx2, j in enumerate(reversed(full_lines)):
                if j.startswith('             Standard Nuclear Orientation (Angstroms)\n'):
                    break
            geo = []
            for i in full_lines[-idx2 + 2 : -idx2 + 2 + atom_number]:
                atom = i.split()[1:]
                geo.append('  '.join(atom))
            with open(opt_in, 'w') as f:
                f.write('$molecule\n{} {}\n'.format(0, 1))
                f.write('\n'.join(geo))
                f.write('\n$end\n\n')
                for line in config:
                    f.write(line + '\n')
        else:
            for idx, i in enumerate(lines):
                if i.startswith('             Standard Nuclear Orientation (Angstroms)\n'):
                    break
            geo = []
            for i in lines[idx + 3 : idx + 3 + atom_number]:
                atom = i.split()[1:]
                geo.append('  '.join(atom))
                
            with open(opt_in, 'w') as f:
                f.write('$molecule\n{} {}\n'.format(0, 1))
                f.write('\n'.join(geo))
                f.write('\n$end\n\n')
                for line in config:
                    f.write(line + '\n')

        return 'need opt'
    elif lines[-5] == '        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n':
        return 'job_success'
    elif lines[-2] == ' Bad initial gradient\n':
        return 'Bad initial gradient'
    elif lines[-2] == ' IRC --- Failed line search\n':
        return 'Failed line search'
    elif lines[-6] == ' Error in gen_scfman\n':
        return 'Error in gen_scfman'
    else:
        return 'unknown fail information'

def generate_irc_product_xyz(target, direction='forward'):
    irc_path = path.join(target['path'], 'IRC')
    reactant_path = path.join(target['path'], 'reactant.xyz')
    output_name = 'irc_{}.out'.format(direction)
    output = path.join(irc_path, output_name)
    name = path.join(irc_path, '{}.xyz'.format(direction))

    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])

    with open(output, 'r') as f:
        full_lines = f.readlines()
    count = 1
    for idx, i in enumerate(full_lines):
        if i.startswith('  IRC -- convergence criterion reached.\n'):
            count += 1
            if count == 2:
                break
    full_lines = full_lines[:idx]
    for idx2, j in enumerate(reversed(full_lines)):
        if j.startswith('             Standard Nuclear Orientation (Angstroms)\n'):
            break
    geo = []
    for i in full_lines[-idx2 + 2 : -idx2 + 2 + atom_number]:
        atom = i.split()[1:]
        geo.append('  '.join(atom))
    with open(name, 'w') as f:
        f.write(str(atom_number))
        f.write('\n\n')
        f.write('\n'.join(geo))

def check_irc_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_irc_target(direction = 'forward')

    qm_collection = db['qm_calculate_center']
    # 2. check the job pbs status
    for target in targets:
        job_id = target['irc_forward_jobid']
        # 2. check the job pbs status
        new_status = check_irc_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_irc_content(target['path'], direction='forward')

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        irc_status = 'irc_{}_status'.format('forward')
        orig_status = target[irc_status]
        if orig_status != new_status:
            if new_status == 'job_success':
                generate_irc_product_xyz(target, direction='forward')
                update_field = {
                    irc_status:new_status, 'irc_equal':'waiting for check'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            elif new_status == 'need opt':
                opt_status = 'opt_{}_status'.format('forward')
                update_field = {
                    irc_status:new_status, opt_status:'job_unrun'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                irc_status: new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)
    
    # 1. select jobs to check
    targets = select_irc_target(direction = 'reverse')

    qm_collection = db['qm_calculate_center']
    # 2. check the job pbs status
    for target in targets:
        job_id = target['irc_reverse_jobid']
        # 2. check the job pbs status
        new_status = check_irc_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_irc_content(target['path'], direction='reverse')

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        irc_status = 'irc_{}_status'.format('reverse')
        orig_status = target[irc_status]
        if orig_status != new_status:
            if new_status == 'job_success':
                generate_irc_product_xyz(target, direction='reverse')
                update_field = {
                    irc_status:new_status, 'irc_equal':'waiting for check'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            elif new_status == 'need opt':
                opt_status = 'opt_{}_status'.format('reverse')
                update_field = {
                    irc_status:new_status, opt_status:'job_unrun'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                irc_status: new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)

"""
IRC  success check.
This is to check whether irc forward and reverse direction equal to expectation.
"""
def select_irc_equal_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    query = {'$and': 
                    [
                    { "irc_forward_status":
                        {"$in":
                        ['job_success', 'opt_success']}},
                    {'irc_reverse_status':
                        {'$in':
                            ['job_success', 'opt_success']}},
                    {'irc_equal':
                        {'$in':
                            ['waiting for check']}}
                    ]
                }
    targets = list(qm_collection.find(query))

    return targets

def check_irc_equal():

    targets = select_irc_equal_target()
    qm_collection = db['qm_calculate_center']
    
    for target in targets:
        if target['irc_forward_status'] in ['job_success', 'opt_success'] and target['irc_reverse_status'] in ['job_success', 'opt_success']:
            new_status = check_irc_equal_status(target)
            orig_status = target['irc_equal']
            if orig_status != new_status:
                if new_status == 'forward equal to reactant and reverse equal to product' or new_status == 'reverse equal to reactant and forward equal to product' or new_status == 'reverse equal to reactant but forward does not equal to product' or new_status == 'forward equal to reactant but reverse does not equal to product':
                    update_field = {
                                        'irc_equal': new_status, 'energy_status': 'job_unrun'
                                    }
                    qm_collection.update_one(target, {"$set": update_field}, True)
                else:
                    update_field = {
                                        'irc_equal': new_status
                                    }
                    qm_collection.update_one(target, {"$set": update_field}, True)
        elif target['opt_reverse_status'] == 'job_fail' or target['opt_forward_status'] == 'job_fail':
            update_field = {
                                'irc_equal': 'opt fail'
                            }
            qm_collection.update_one(target, {"$set": update_field}, True)

def check_irc_equal_status(target):

    irc_path = path.join(target['path'], 'IRC/')
    reactant_path = path.join(target['path'], 'reactant.xyz')
    product_path = path.join(target['path'], 'ssm_product.xyz')
    forward_output = path.join(irc_path, 'forward.xyz')
    reverse_output = path.join(irc_path, 'reverse.xyz')
    
    pyMol_1 = xyz_to_pyMol(reactant_path)
    pyMol_2 = xyz_to_pyMol(product_path)
    pyMol_3 = xyz_to_pyMol(forward_output)
    pyMol_4 = xyz_to_pyMol(reverse_output)
    
    if pyMol_3.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'forward equal to reverse'
    elif (pyMol_1.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip()) and (pyMol_1.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip()):
        return 'forward and reverse equal to reactant'
    elif (pyMol_2.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip()) and (pyMol_2.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip()):
        return 'forward and reverse equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'forward equal to reactant and reverse equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip():
        return 'reverse equal to reactant and forward equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() != pyMol_3.write('inchiKey').strip():
        return 'reverse equal to reactant but forward does not equal to product'
    elif pyMol_1.write('inchiKey').strip() != pyMol_4.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip():
        return 'reverse does not equal to reactant but forward equal to product'
    elif pyMol_1.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() != pyMol_4.write('inchiKey').strip():
        return 'forward equal to reactant but reverse does not equal to product'
    elif pyMol_1.write('inchiKey').strip() != pyMol_3.write('inchiKey').strip() and pyMol_2.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'forward does not equal to reactant but reverse equal to product'
    else:
        return 'unknown (Maybe both of them are not equal to reactant&product)'

"""
IRC opt check.
"""
def select_irc_opt_target(direction = 'forward'):
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    irc_opt_status = 'opt_{}_status'.format(direction)
    reg_query = {irc_opt_status:
                    {"$in":
                        ["opt_job_launched", "opt_job_running", "opt_job_queueing"]
                    }
                }
    targets = list(qm_collection.find(reg_query))

    return targets

def check_irc_opt_job_status(job_id):
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
        return "opt_job_running"
    elif stdout[idx+2] == 'Q':
        return 'opt_job_queueing'
    else:
        return "opt_job_launched"
    
def check_irc_opt_content(dir_path, direction = 'forward'):
    reactant_path = os.path.join(dir_path, 'reactant.xyz')
    irc_path = path.join(dir_path, "IRC")
    xyzname = '{}.xyz'.format(direction)
    output = path.join(irc_path, xyzname)
    
    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])
    
    outputname = path.join(irc_path, '{}_opt.out'.format(direction))
    
    with open(outputname, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 4096, 0), 0)  # Read last  4kB of file
        lines = f.readlines()
        
    if lines[-5] == '        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n':
        for idx, i in enumerate(lines):
            if i.startswith('                       Coordinates (Angstroms)\n'):
                break

        geo = []
        for i in lines[idx + 2 : idx + 2 + atom_number]:
            atom = i.split()[1:]
            geo.append('  '.join(atom))

        with open(output, 'w') as f2:
            f2.write(str(len(geo)))
            f2.write('\n\n')
            f2.write('\n'.join(geo))
            
        return 'job_success'
    else:
        return 'job_fail'
        
    
def check_irc_opt_job():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_irc_opt_target(direction = 'forward')
    irc_opt_jobid = 'irc_{}_opt_jobid'.format(str('forward'))
    qm_collection = db['qm_calculate_center']

    # 2. check the job pbs_status
    for target in targets:
        job_id = target[irc_opt_jobid]
        new_status = check_irc_opt_job_status(job_id)
        if new_status == "off_queue":
                # 3. check job content
                new_status = check_irc_opt_content(target['path'], direction = 'forward')

                # 4. check with original status which
                # should be job_launched or job_running
                # if any difference update status
                irc_status = 'irc_{}_status'.format(str('forward'))
                irc_opt_status = 'opt_{}_status'.format(str('forward'))
                orig_status = target[irc_opt_status]
                if orig_status != new_status:
                    if new_status == 'job_success':
                        update_field = {
                                        irc_status: 'opt_success', irc_opt_status: new_status, 'irc_equal':'waiting for check'
                                    }
                        qm_collection.update_one(target, {"$set": update_field}, True)
                    else:
                        update_field = {
                                        irc_status: 'opt_fail', irc_opt_status: new_status
                                    }
                        qm_collection.update_one(target, {"$set": update_field}, True)
                               
    # 1. select jobs to check
    targets = select_irc_opt_target(direction = 'reverse')
    irc_opt_jobid = 'irc_{}_opt_jobid'.format(str('reverse'))
    
    # 2. check the job pbs_status
    for target in targets:
        job_id = target[irc_opt_jobid]
        new_status = check_irc_opt_job_status(job_id)
        if new_status == "off_queue":
                # 3. check job content
                new_status = check_irc_opt_content(target['path'], direction = 'reverse')

                # 4. check with original status which
                # should be job_launched or job_running
                # if any difference update status
                irc_status = 'irc_{}_status'.format(str('reverse'))
                irc_opt_status = 'opt_{}_status'.format(str('reverse'))
                orig_status = target[irc_opt_status]
                if orig_status != new_status:
                    if new_status == 'job_success':
                        update_field = {
                                        irc_status: 'opt_success', irc_opt_status: new_status, 'irc_equal':'waiting for check'
                                    }
                        qm_collection.update_one(target, {"$set": update_field}, True)
                    else:
                        update_field = {
                                        irc_status: 'opt_fail', irc_opt_status: new_status
                                    }
                        qm_collection.update_one(target, {"$set": update_field}, True)

"""
OPT check
"""
def select_opt_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    reg_query = {'opt_status':
                    {"$in":
                        ["job_launched", "job_running", "job_queueing"]
                    }
                }
    targets = list(qm_collection.find(reg_query))
    return targets

def check_opt_job_status(job_id):
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
    
def check_opt_content(dir_path):
    reactant_path = os.path.join(dir_path, 'reactant.xyz')
    OPT_path = path.join(dir_path, "OPT")
    
    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])
    
    outputname = path.join(OPT_path, 'opt.out')
    
    with open(outputname, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 4096, 0), 0)  # Read last  4kB of file
        lines = f.readlines()
        
    if lines[-5] == '        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n':
        for idx, i in enumerate(lines):
            if i.startswith('                       Coordinates (Angstroms)\n'):
                break

        geo = []
        for i in lines[idx + 2 : idx + 2 + atom_number]:
            atom = i.split()[1:]
            geo.append('  '.join(atom))

        with open(reactant_path, 'w') as f2:
            f2.write(str(len(geo)))
            f2.write('\n\n')
            f2.write('\n'.join(geo))
            
        return 'job_success'
    else:
        return 'job_fail'

def check_opt_job():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_opt_target()
    qm_collection = db['qm_calculate_center']

    # 2. check the job pbs_status
    for target in targets:
        job_id = target['opt_jobid']
        new_status = check_opt_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_opt_content(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['opt_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                                'opt_status': new_status, 'ssm_status': 'job_unrun'
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                'opt_status': new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)

"""
Low level OPT check
"""
def select_low_opt_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    reg_query = {'low_opt_status':
                    {"$in":
                        ["job_launched", "job_running", "job_queueing"]
                    }
                }
    targets = list(qm_collection.find(reg_query))
    return targets

def check_low_opt_job_status(job_id):
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
    
def check_low_opt_content(dir_path):
    reactant_path = os.path.join(dir_path, 'reactant.xyz')
    OPT_path = path.join(dir_path, "OPT")
    
    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])
    
    outputname = path.join(OPT_path, 'low_opt.out')
    
    with open(outputname, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 4096, 0), 0)  # Read last  4kB of file
        lines = f.readlines()
        
    if lines[-5] == '        *  Thank you very much for using Q-Chem.  Have a nice day.  *\n':
        for idx, i in enumerate(lines):
            if i.startswith('                       Coordinates (Angstroms)\n'):
                break

        geo = []
        for i in lines[idx + 2 : idx + 2 + atom_number]:
            atom = i.split()[1:]
            geo.append('  '.join(atom))

        with open(reactant_path, 'w') as f2:
            f2.write(str(len(geo)))
            f2.write('\n\n')
            f2.write('\n'.join(geo))

        for idx2, i in enumerate(lines):
            if i.startswith(' **  OPTIMIZATION CONVERGED  **\n'):
                break
        
        low_opt_energy = lines[idx2-4].split()[-1]
        return 'job_success', float(low_opt_energy)
    else:
        return 'job_fail', float(0.0)

def check_low_opt_job():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_low_opt_target()
    qm_collection = db['qm_calculate_center']

    # 2. check the job pbs_status
    for target in targets:
        job_id = target['low_opt_jobid']
        new_status = check_opt_job_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status, energy = check_low_opt_content(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['low_opt_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                                'low_opt_status': new_status, 'low_energy': energy, 'check_binding_status': 'need check'
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                'low_opt_status': new_status, 'low_energy': energy
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)

"""
Barrier check.
"""
def select_barrier_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    query = {'barrier':
                {"$in":
                    ["None"]
                }
            }
    targets = list(qm_collection.find(query))
    return targets

def update_barrier_information():
    qm_collection = db['qm_calculate_center']
    targets = select_barrier_target()
    for target in targets:
        reactant_energy = float(target['reactant_scf_energy'])
        ts_energy = float(target['ts_energy'])
        barrier = (ts_energy - reactant_energy) * 627.5095
        update_field = {'barrier':barrier, 'insert reaction':"need insert"}
        qm_collection.update_one(target, {"$set": update_field}, True)

"""
After irc check, insert the reaction into reaction collection.
Here we select the lowest activation energy one.  #TO DO
If we insert to reaction collection, next generation is ready to run. (In the other words, create ard_status: job_unrun)  <-- the lowest activation energy one
BTW, when insert reaction, we may consider the same reaction had been generated by early generation.
# To make sure the lowest barrier add to reaction collection.
# So we need to know all of the job are finished.
# The basic thought is that there is not any unrun job except ard job.
# def create ard job  #TO DO
"""
def select_insert_reaction_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    query = {'insert reaction':
                {"$in":
                    ["need insert"]
                }
            }
    targets = list(qm_collection.find(query))
    return targets

def insert_reaction():
    reactions_collection = db['reactions']
    qm_collection = db['qm_calculate_center']
    targets = select_insert_reaction_target()
    # new one not mean the lowest barrier (so the lowest may in duplicate)
    for target in targets:
        reactant_inchi_key = target['reactant_inchi_key']
        product_inchi_key = target['product_inchi_key']
        reactant_smi = target['Reactant SMILES']
        product_smi = target['Product SMILES']
        path = target['path']
        generations = target['generations']
        barrier = target['barrier']
        check1 = {'reaction':[reactant_inchi_key, product_inchi_key]}
        checker1 = list(reactions_collection.find(check1))
        if len(checker1) == 0:
            reactions_collection.insert_one({
                                'reaction':[reactant_inchi_key, product_inchi_key],
                                'reactant_inchi_key': reactant_inchi_key,
                                'product_inchi_key':product_inchi_key,
                                'reactant_smi':reactant_smi,
                                'product_smi':product_smi,
                                'path':path,
                                'generations':generations,
                                'unique': 'new one',
                                'barrier_energy':barrier})
            qm_collection.update_one(target, {"$set": {'insert reaction':"Already insert"}}, True)
        else:
            reactions_collection.insert_one({
                                'reaction':[reactant_inchi_key, product_inchi_key],
                                'reactant_inchi_key': reactant_inchi_key,
                                'product_inchi_key':product_inchi_key,
                                'reactant_smi':reactant_smi,
                                'product_smi':product_smi,
                                'path':path,
                                'generations':generations,
                                'unique': 'duplicate',
                                'barrier_energy':barrier})
            qm_collection.update_one(target, {"$set": {'insert reaction':"Already insert"}}, True)

"""
ARD check unrun
"""
def insert_ard():
    qm_collection = db['qm_calculate_center']
    statistics_collection = db['statistics']
    reactions_collection = db['reactions']
    energy_query = {"energy_status":
                    {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing"]
                    }
                }
    ssm_query = {"ssm_status":
                    {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing"]
                    }
                }
    ts_query = {"ts_status":
                    {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing"]
                    }
                }
    irc_query_1 = {"irc_forward_status":
                    {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing", "need opt"]
                    }
                }
    irc_query_2 = {"irc_reverse_status":
                    {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing", "need opt"]
                    }
                }
    opt_query_1 = {"opt_forward_status":
                    {"$in":
                        ["job_unrun", "opt_job_launched", "opt_job_running", "opt_job_queueing"]
                    }
                }
    opt_query_2 = {"opt_reverse_status":
                    {"$in":
                        ["job_unrun", "opt_job_launched", "opt_job_running", "opt_job_queueing"]
                    }
                }
    not_finished_number = len(list(qm_collection.find(energy_query))) + len(list(qm_collection.find(ssm_query))) + len(list(qm_collection.find(ts_query))) + len(list(qm_collection.find(irc_query_1))) + len(list(qm_collection.find(irc_query_2))) + len(list(qm_collection.find(opt_query_1))) + len(list(qm_collection.find(opt_query_2)))
    ard_had_add_number = qm_collection.count_documents({})
    ard_should_add_number = sum(statistics_collection.distinct("add how many products"))
    make_sure_not_check_again = reactions_collection.distinct("ard_status")

    if int(not_finished_number) == 0 and int(ard_had_add_number) == int(ard_should_add_number) and make_sure_not_check_again == []:
        targets = list(reactions_collection.find({'unique':'new one'}))
        for target in targets:
            dirpath = target['path']
            ard_qm_target = list(qm_collection.find({'path':dirpath}))[0]
            update_field_for_qm_target = {'ard_status':'job_unrun'}
            update_field_for_reaction_target = {'ard_status':'aleady insert to qm'}
            qm_collection.update_one(ard_qm_target, {"$set": update_field_for_qm_target}, True)
            reactions_collection.update_one(target, {"$set": update_field_for_reaction_target}, True)


"""
Binding energy cutoff filter
"""
def check_bindind_cutoff():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    running_query = {"low_opt_status":
                    {"$in":
                        ["job_unrun", "job_launched", "job_running", "job_queueing"]
                    }
                }
    targets = list(qm_collection.find(running_query))
    if len(targets) == 0:
        query = {'check_binding_status': 'need check'}
        targets = list(qm_collection.find(query))
        cutoff_target = list(qm_collection.find({'Reactant':'initial reactant'}))
        cutoff_energy = float(cutoff_target[0]['binding_mode_energy_cutoff'])

        if cutoff_target[0]['binding_cutoff_select'] == 'starting':
            qchem_energy = cutoff_target[0]['low_energy']
        elif cutoff_target[0]['binding_cutoff_select'] == 'lowest':
            query = [{'$match':{'reactant_inchi_key':cutoff_target[0]['reactant_inchi_key']}},
                    {'$group':{'_id':'$reactant_inchi_key', 'low_energy':{'$min':'$low_energy'}}}]
            qchem_energy = list(qm_collection.aggregate(query))[0]['low_energy']
        else:
            print('Please select the lowest or starting as reference')
            print('Note that "Should manually add into database"')
            raise

        for target in targets:
            delta = (float(target['low_energy']) - float(qchem_energy)) * 627.5095
            if delta > cutoff_energy:
                qm_collection.update_one(target, {"$set": {'check_binding_status': 'greater than cutoff'}}, True)
            else:
                qm_collection.update_one(target, {"$set": {'check_binding_status': 'smaller than cutoff', 'opt_status':'job_unrun'}}, True)
    
    return targets



check_energy_jobs()
check_ssm_jobs()
check_low_opt_job()
check_opt_job()
check_ts_jobs()
check_irc_jobs()
check_irc_equal()
check_irc_opt_job()
update_barrier_information()
insert_reaction()
insert_ard()