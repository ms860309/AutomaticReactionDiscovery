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
from qchem import QChem

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
            q = QChem(outputfile=energy_path)
            energy = q.get_energy()
            zpe = q.get_zpe()
            energy += zpe
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
            if new_status == 'job_success':
                update_field = {
                                'energy_status': new_status,
                                'reactant_energy':energy
                            }
            else:
                update_field = {
                                'energy_status': new_status
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
            f.seek(max(fsize - 4096, 0), 0)  # Read last 4 kB of file
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

def ard_prod_and_ssm_prod_checker(rxn_dir, refine = False):
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
                if refine:
                    update_field = {'product_inchi_key':pyMol_1.write('inchiKey').strip(), 
                                    'initial_dir_name':rxn_dir, 
                                    'path':new_path, 
                                    'ssm_status': 'job_success',                                # Though the reactant equal to product but TS maybe fine. The SSM somehow do not use constrained optimize in product
                                    'ard_ssm_equal':'ssm reactant equal to product',            # This is a special case but sometimes it will happen.
                                    'Product SMILES': prod_smi,
                                    "ts_refine_status":"job_unrun"}                                    # If check ts and get the success maybe need to check irc but QChem's irc is not robust
                else:
                    update_field = {'product_inchi_key':pyMol_1.write('inchiKey').strip(), 
                                    'initial_dir_name':rxn_dir, 
                                    'path':new_path, 
                                    'ssm_status': 'job_success',                               
                                    'ard_ssm_equal':'ssm reactant equal to product',            
                                    'Product SMILES': prod_smi,
                                    "ts_status":"job_unrun"}                                    
                qm_collection.update_one(i, {"$set": update_field}, True)
            else:
                dirname = dir_check(path.dirname(i['path']), pyMol_1.write('inchiKey').strip(), num + 1)
                new_path = path.join(path.dirname(i['path']), dirname)
                os.rename(rxn_dir, new_path)
                prod_smi = pyMol_1.write('can').split()[0]
                if refine:
                    update_field = {'product_inchi_key':pyMol_1.write('inchiKey').strip(), 
                                    'initial_dir_name':rxn_dir, 
                                    'path':new_path, 
                                    'ssm_status': 'job_success', 
                                    "ts_refine_status":"job_unrun", 
                                    'ard_ssm_equal':'not_equal',
                                    'Product SMILES': prod_smi}
                else:
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

def check_ssm_jobs(refine = False):
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
                equal = ard_prod_and_ssm_prod_checker(target['path'], refine = refine)
                if equal == 'equal':
                    if refine:
                        update_field = {
                                        'ssm_status': new_status, "ts_refine_status":"job_unrun", 'ard_ssm_equal':equal
                                    }
                    else:
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
TS refine check
Only update whether job is finished. (Do not check in converge or fail ...)
"""

def select_ts_refine_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    reg_query = {"ts_refine_status":
                    {"$in": 
                        ["job_launched", "job_running", "job_queueing"] 
                    }
                }
    targets = list(qm_collection.find(reg_query))

    return targets

def check_ts_refine_job_status(job_id):
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


def check_ts_refine_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_ts_refine_target()

    qm_collection = db['qm_calculate_center']
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ts_refine_jobid']
        # 2. check the job pbs status
        new_status = check_ts_refine_job_status(job_id)
        if new_status == "off_queue":
            update_field = {
                            'ts_status': 'job_unrun', 'ts_refine_status': 'job_success'
                            }
            qm_collection.update_one(target, {"$set": update_field}, True)
        else:
            orig_status = target['ts_refine_status']
            if orig_status != new_status:
                update_field = {
                                'ts_refine_status': new_status
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

    try:
        q = QChem(outputfile=ts_out_path)
        freqs = q.get_frequencies()
        nnegfreq = sum(1 for freq in freqs if freq < 0.0)
        if nnegfreq > 1:
            return "Have more than one imaginary frequency", 0.0
        elif nnegfreq == 0:
            return "All positive frequency", 0.0
        else:
            energy = q.get_energy()
            zpe = q.get_zpe()
            ts_energy = energy + zpe
            q.create_geo_file(ts_geo_path)
            return "job_success", float(ts_energy)
    except:
        return "job_fail", 0.0

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
                try:
                    reactant_energy = float(target['reactant_energy'])
                    barrier = (ts_energy - reactant_energy) * 627.5095
                except:
                    barrier = 'do not have reactant_energy'
                if target['use_irc'] == '0':
                    if barrier == 'do not have reactant_energy':
                        update_field = {
                                        'ts_status': new_status, 'ts_energy':ts_energy, 'energy_stauts':'job_unrun'
                                        }
                    else:
                        update_field = {
                                        'ts_status': new_status, 'ts_energy':ts_energy, 'barrier':barrier, 'insert reaction':'need insert'
                                        }    
                else:
                    if barrier == 'do not have reactant_energy':
                        update_field = {
                                        'ts_status': new_status, 'ts_energy':ts_energy, 'irc_status':'job_unrun', 'energy_stauts':'job_unrun'
                                        }
                    else:
                        update_field = {
                                        'ts_status': new_status, 'ts_energy':ts_energy, 'irc_status':'job_unrun', 'barrier':barrier
                                        }
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
    opt_lot = path.join(base_dir_path, 'opt_freq.lot')
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
        f.seek(max(fsize - 12280, 0), 0)  # Read last 12 kB of file
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
    output_path = path.join(irc_path, '{}_opt.out'.format(direction))

    try:
        q = QChem(outputfile = output_path)
        q.create_geo_file(output)
        return 'job_success'
    except:
        return 'job_fail'
        
    
def check_irc_opt_jobs():
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
    output_path = path.join(OPT_path, 'opt.out')

    try:
        q = QChem(outputfile=output_path)
        freqs = q.get_frequencies()
        nnegfreq = sum(1 for freq in freqs if freq < 0.0)
        if nnegfreq > 0:
            return 'Have negative frequency', 0, 0.0
        else:
            opt_cycle = q.get_opt_cycle()
            q.create_geo_file(reactant_path)
            energy = q.get_energy()
            zpe = q.get_zpe()
            energy += zpe
            return 'job_success', opt_cycle, energy
    except:
        return 'job_fail', 0, 0.0

def check_opt_jobs():
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
            new_status, opt_cycle, energy = check_opt_content(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['opt_status']
        if orig_status != new_status:
            if new_status == 'job_success' and target['Reactant SMILES'] != 'initial reactant':
                update_field = {
                                'opt_status': new_status, 'ssm_status': 'job_unrun', 'opt_iter':opt_cycle, 'reactant_energy':energy
                            }
            elif new_status == 'job_success' and target['Reactant SMILES'] == 'initial reactant':
                update_field = {
                                'opt_status': new_status, 'opt_iter':opt_cycle, 'reactant_energy':energy
                            }
            elif new_status == "job_running" or new_status == "job_queueing":
                update_field = {
                                'opt_status': new_status
                            }
            else:
                update_field = {
                                'opt_status': new_status, 'opt_iter':opt_cycle
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
    output_path = path.join(OPT_path, 'low_opt.out')

    try:
        q = QChem(outputfile=output_path)
        opt_cycle = q.get_opt_cycle()
        q.create_geo_file(reactant_path)
        energy = q.get_energy()
        return 'job_success', opt_cycle, energy
    except:
        return 'job_fail', 0, 0.0

def check_low_opt_jobs():
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
            new_status, low_opt_cycle, energy = check_low_opt_content(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['low_opt_status']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                                'low_opt_status': new_status, 'low_energy': energy, 'low_opt_iter':low_opt_cycle, 'check_binding_status': 'need check'
                            }
            elif new_status == "job_running" or new_status == "job_queueing":
                update_field = {
                                'low_opt_status': new_status
                            }
            else:
                update_field = {
                                'low_opt_status': new_status, 'low_energy': energy, 'low_opt_iter':low_opt_cycle
                            }
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
    use_irc_query = {'Reactant SMILES':'initial reactant'}
    use_irc = list(qm_collection.find(use_irc_query))[0]['use_irc']
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
    if use_irc == '0':
        not_finished_number = len(list(qm_collection.find(energy_query))) + len(list(qm_collection.find(ssm_query))) + len(list(qm_collection.find(ts_query)))
    else:
        not_finished_number = len(list(qm_collection.find(energy_query))) + len(list(qm_collection.find(ssm_query))) + len(list(qm_collection.find(ts_query))) + len(list(qm_collection.find(irc_query_1))) + len(list(qm_collection.find(irc_query_2))) + len(list(qm_collection.find(opt_query_1))) + len(list(qm_collection.find(opt_query_2)))
    ard_had_add_number = qm_collection.count_documents({})  # Should -1 because the initial reactant
    ard_should_add_number = sum(statistics_collection.distinct("add how many products"))

    if int(not_finished_number) == 0 and int(ard_had_add_number) -1 == int(ard_should_add_number):
        targets = list(reactions_collection.find({'unique':'new one'}))
        for target in targets:
            try:
                ard_stauts = target['ard_status']
                if ard_stauts == 'aleady insert to qm':
                    continue
            except:
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
        cutoff_target = list(qm_collection.find({'Reactant SMILES':'initial reactant'}))
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
                qm_collection.update_one(target, {"$set": {'check_binding_status': 'greater than cutoff', 'deltaH':delta}}, True)
            elif delta <= cutoff_energy:
                qm_collection.update_one(target, {"$set": {'check_binding_status': 'smaller than cutoff', 'opt_status':'job_unrun', 'deltaH':delta}}, True)
            elif target['Reactant SMILES'] == "initial reactant":
                qm_collection.update_one(target, {"$set": {'check_binding_status': 'had checked', 'deltaH':delta}}, True)
            else:
                print('have unknow error')
                raise

    return targets

"""
Check barrier which do not have reactant energy
"""

def select_barrier_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    query = {'$and': 
                    [
                    { "energy_stauts":
                        {"$in":
                        ['job_success']}},
                    {'barrier':
                        {'$in':
                            ['do not have reactant_energy']}}
                    ]
                }
    targets = list(qm_collection.find(query))
    return targets

def check_barrier():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_barrier_target()

    qm_collection = db['qm_calculate_center']
    # 2. check the job pbs status
    for target in targets:
        reactant_energy = float(target['reactant_energy'])
        ts_energy = float(target['ts_energy'])
        barrier = (ts_energy - reactant_energy) * 627.5095
        update_field = {
                        'barrier':barrier
                        }
        qm_collection.update_one(target, {"$set": update_field}, True)


check_energy_jobs()
check_ssm_jobs(refine = True)  # If the ssm perform by orca with xtb GFN2-xtb, then refine the TS is a good choice.  Get a better initial guess
#check_low_opt_jobs()
#check_bindind_cutoff()
#check_opt_jobs()
check_ts_refine_jobs()
check_ts_jobs()
#check_irc_jobs()
#check_irc_equal()
#check_irc_opt_jobs()
check_barrier()
insert_reaction()
insert_ard()
