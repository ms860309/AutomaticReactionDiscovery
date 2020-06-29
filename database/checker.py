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
    elif check_ssm(ssm_path) == 'Exiting early':
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
    pyMol_1 = xyz_to_pyMol(ssm_prod_path)
    rmg_mol_1 = toRMGmol(pyMol_1)
    pyMol_2 = xyz_to_pyMol(ard_prod_path)
    rmg_mol_2 = toRMGmol(pyMol_2)
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
                                'ssm_status': 'job_success', 
                                'ard_ssm_equal':'not_equal but reactant equal to product',
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
                                "energy_status":"job_unrun", 
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

def toRMGmol(pyMol):
    rmg_mol = from_ob_mol(rmgpy.molecule.molecule.Molecule(), pyMol.OBMol)
    return rmg_mol

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
    query = {'reaction':[reactant_inchi_key, product_inchi_key]}
    targets = list(reactions_collection.find(query))
    if len(targets) == 0:
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
                #pool_collection.insert_one({'reactant_inchi_key':target['product_inchi_key']})
            else:
                update_field = {
                                'ts_status': new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)


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
    
def check_irc_content_status(target_path, direction = 'forward'):
    reactant_path = os.path.join(target_path, 'reactant.xyz')
    irc_path = os.path.join(target_path, 'IRC/')
    opt_name = '{}_opt.in'.format(direction)
    opt_in = os.path.join(irc_path, opt_name)
    if direction == 'forward':
        irc_output = os.path.join(irc_path, 'irc_forward.out')
    else:
        irc_output = os.path.join(irc_path, 'irc_reverse.out')

    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])
    with open(irc_output, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 10240, 0), 0)  # Read last 10 kB of file
        lines = f.readlines()

    if lines[-2] == ' IRC backup failure\n' or lines[-2] == ' IRC failed final bisector step\n':
    
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
            f.write('$rem\n')
            f.write('JOBTYPE OPT\n')
            f.write('METHOD B97-D3\n')
            f.write('DFT_D D3_BJ\n')
            f.write('BASIS def2-mSVP\n')
            f.write('SCF_ALGORITHM DIIS\n')
            f.write('MAX_SCF_CYCLES 150\n')
            f.write('SCF_CONVERGENCE 8\n')
            f.write('SYM_IGNORE TRUE\n')
            f.write('SYMMETRY FALSE\n')
            f.write('GEOM_OPT_MAX_CYCLES 150\n')
            f.write('GEOM_OPT_TOL_GRADIENT 100\n')
            f.write('GEOM_OPT_TOL_DISPLACEMENT 400\n')
            f.write('GEOM_OPT_TOL_ENERGY 33\n')
            f.write('WAVEFUNCTION_ANALYSIS FALSE\n')
            f.write('$end\n')
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
            new_status = check_irc_content_status(target['path'], direction='forward')

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        irc_status = 'irc_{}_status'.format('forward')
        orig_status = target[irc_status]
        if orig_status != new_status:
            if new_status == 'job_success':
                generate_xyz(target, direction='forward')
                update_field = {
                    irc_status:new_status, 'irc_equal':'waiting for check'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            elif new_status == 'need opt':
                opt_status = 'opt_{}_status'.format('forward')
                update_field = {
                    irc_status:new_status, 'opt_status':'job_unrun'
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
            new_status = check_irc_content_status(target['path'], direction='reverse')

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        irc_status = 'irc_{}_status'.format('reverse')
        orig_status = target[irc_status]
        if orig_status != new_status:
            if new_status == 'job_success':
                generate_xyz(target, direction='reverse')
                update_field = {
                    irc_status:new_status, 'irc_equal':'waiting for check'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            elif new_status == 'need opt':
                opt_status = 'opt_{}_status'.format('forward')
                update_field = {
                    irc_status:new_status, 'opt_status':'job_unrun'
                                }
                qm_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                irc_status: new_status
                            }
                qm_collection.update_one(target, {"$set": update_field}, True)

def generate_xyz(target, direction='forward'):
    irc_path = os.path.join(target, 'IRC')
    reactant_path = os.path.join(target, 'reactant.xyz')
    
    with open(reactant_path, 'r') as f1:
        lines = f1.readlines()
    atom_number = int(lines[0])
    output_name = 'irc_{}.out'.format(direction)
    output = os.path.join(irc_path, output_name)
    with open(output, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 6144, 0), 0)  # Read last  kB of file
        lines = f.readlines()
    for idx, i in enumerate(lines):
        if i.startswith('             Standard Nuclear Orientation (Angstroms)\n'):
            break
    geo = []
    for i in lines[idx + 3 : idx + 3 + atom_number]:
        atom = i.split()[1:]
        geo.append('  '.join(atom))
    name = '{}.xyz'.format(direction)
    with open(name, 'w') as f:
        f.write(str(atom_number))
        f.write('\n\n')
        f.write('\n'.join(geo))
        
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
                            ['job_success', 'job_fail']}}
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
        min_energy = list(qm_collection.aggregate(
            [
                {"$match": {'reactant_inchi_key':target['reactant_inchi_key']}},
                {"$group":{"_id":{}, 'reactant_scf_energy':{'$min':"$reactant_scf_energy"}}}
            ]
        ))
        reactant_energy = float(min_energy[0]['reactant_scf_energy'])
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
                }, {'ard_ssm_equal':
                    {'$nin':
                        ['not_equal but reactant equal to product']}
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
            query_2 = {'$and':
                            [{"ard_status":
                            {"$in":
                                ['job_unrun']
                            }
                        },
                        {
                            'product_inchi_key':target[0]['product_inchi_key']
                        }]}
            
            if len(list(qm_collection.find(query_2))) == 0:
                qm_collection.update_one(target[0], {"$set": {'ard_status':'job_unrun'}}, True)
                
            insert_exact_rxn(target[0]['reactant_inchi_key'],
                            target[0]['product_inchi_key'],
                            target[0]['Reactant SMILES'],
                            target[0]['Product SMILES'],
                            target[0]['path'],
                            target[0]['generations'],
                            barrier)

def select_same_ssm_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    reactions_collection = db['reactions']
    reg_query = {"ssm":
                    {"$in": 
                        ["job_launched", "job_running", "job_queueing"] 
                    }
                }
    targets = list(reactions_collection.find(reg_query))

    return targets

def check_same_ssm_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_ssm_target()

    reactions_collection = db['reactions']
    
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
                equal = same_ard_prod_and_ssm_prod_checker(target['path'])
                if equal == 'equal':
                    update_field = {
                                    'ssm': new_status, "ts":"job_unrun", 'ard_ssm_equal':equal
                                }
                    reactions_collection.update_one(target, {"$set": update_field}, True)
                else:
                    update_field = {
                                    'ssm': new_status, 'ard_ssm_equal':equal, 'unique': 'fail'
                                }
                    reactions_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                'ssm': new_status, 'unique': 'fail'
                            }
                reactions_collection.update_one(target, {"$set": update_field}, True)

def same_ard_prod_and_ssm_prod_checker(rxn_dir):
    # Use ssm product xyz to check whether ssm prod equal to ard product
    # If equal, insert the inchi key into products pool
    # If not equal, use ssm product as the product and insert inchi key into products pool
    # Next generation use the ssm product to generate
    
    ard_prod_path = path.join(rxn_dir, 'product.xyz')
    ssm_prod_path = path.join(rxn_dir, 'ssm_product.xyz')
    pyMol_1 = xyz_to_pyMol(ssm_prod_path)
    rmg_mol_1 = toRMGmol(pyMol_1)
    pyMol_2 = xyz_to_pyMol(ard_prod_path)
    rmg_mol_2 = toRMGmol(pyMol_2)
    if pyMol_1.write('inchiKey').strip() != pyMol_2.write('inchiKey').strip():
        return 'not_equal'
    else:
        return 'equal'


def select_same_ts_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    reactions_collection = db['reactions']
    reg_query = {"ts":
                    {"$in": 
                        ["job_launched", "job_running", "job_queueing"] 
                    }
                }
    targets = list(reactions_collection.find(reg_query))

    return targets

def check_same_ts_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_same_ts_target()

    reactions_collection = db['reactions']
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
        orig_status = target['ts']
        if orig_status != new_status:
            if new_status == 'job_success':
                update_field = {
                                'ts': new_status, 'ts_energy':ts_energy, 'irc':'job_unrun'
                                }
                reactions_collection.update_one(target, {"$set": update_field}, True)
            else:
                update_field = {
                                'ts': new_status, 'unique': 'fail'
                            }
                reactions_collection.update_one(target, {"$set": update_field}, True)

def check_same_barrier():
    reactions_collection = db['reactions']
    qm_collection = db['qm_calculate_center']
    
    query = {'$and': 
                    [
                    { "ssm":
                        {"$in":
                        ['job_success']}},
                    {'ts':
                        {'$in':
                            ['job_success']}}
                    ]
                }
    
    targets = list(reactions_collection.find(query))
    
    for target in targets:
        reactant_inchi_key = target['reaction'][0]
        
        min_energy = list(qm_collection.aggregate(
            [
                {"$match": {'reactant_inchi_key':reactant_inchi_key}},
                {"$group":{"_id":{}, 'reactant_scf_energy':{'$min':"$reactant_scf_energy"}}}
            ]
        ))
        reactant_scf_energy = float(min_energy[0]['reactant_scf_energy'])
        barrier = (float(target['ts_energy']) - reactant_scf_energy)*627.5095
        update_field = {'barrier_energy':barrier}
        reactions_collection.update_one(target, {"$set": update_field}, True)
        
    clean_same(targets)

def clean_same(targets):
    reactions_collection = db['reactions']
    for target in targets:
        min_barrier = list(reactions_collection.aggregate(
            [
                {"$match": {'unique': 'waiting for check'}},
                {"$group":{"_id":'$reaction', 'barrier_energy':{'$min':"$barrier_energy"}}}
            ]
        ))
        for i in min_barrier:
            query = {'$and': 
                            [
                            { "reaction":i[_id]},
                            {'unique':'waiting for check'}
                            ]
                        }
            for j in list(reactions_collection.find(query)):
                if j['barrier_energy'] == i['barrier_energy']:
                    update_field = {'unique':'new one'}
                    reactions_collection.update_one(target, {"$set": update_field}, True)
                else:
                    update_field = {'unique':'have another lower barrier reaction'}
                    reactions_collection.update_one(target, {"$set": update_field}, True)


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
                        ['job_success']}},
                    {'irc_reverse_status':
                        {'$in':
                            ['job_success']}}
                    ]
                }
    targets = list(qm_collection.find(reg_query))

    return targets

def check_irc_equal():

    targets = select_irc_equal_target()
    qm_collection = db['qm_calculate_center']
    
    for target in targets:
        new_status = check_irc_equal_status(target)
        orig_status = target['irc_equal']
        if orig_status != new_status:
            update_field = {
                                'irc_equal': new_status
                            }
            qm_collection.update_one(target, {"$set": update_field}, True)
            
def check_irc_equal_status(target):
    
    irc_path = os.path.join(target, 'IRC/')
    reactant_path = os.path.join(target, 'reactant.xyz')
    product_path = os.path.join(target, 'ssm_product.xyz')
    forward_output = os.path.join(irc_path, 'forward.xyz')
    reverse_output = os.path.join(irc_path, 'reverse.xyz')
    
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
    elif pyMol_1.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip():
        return 'forward equal to reactant'
    elif pyMol_1.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'reverse equal to reactant'
    elif pyMol_2.write('inchiKey').strip() == pyMol_3.write('inchiKey').strip():
        return 'forward equal to product'
    elif pyMol_2.write('inchiKey').strip() == pyMol_4.write('inchiKey').strip():
        return 'reverse equal to product'
    else:
        return 'unknown'
    
def select_irc_opt_target(direction = 'forward'):
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    qm_collection = db['qm_calculate_center']
    irc_status = 'irc_{}'.format(str(direction))
    targets = list(qm_collection.find(reg_query))
    reg_query = {irc_status:
                    {"$in":
                        ["job_launched", "job_running", "job_queueing"]
                    }
                }
    targets = list(qm_collection.find(reg_query))

    return targets

def check_irc_opt_status(job_id):
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
    
def check_irc_opt_content_status(dir_path, direction = 'forward'):
    
    irc_path = path.join(dir_path, "IRC")
    outputname = '{}_opt.out'.format(direction)
    with open(outputname, 'w') as f:
        # TO DO

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
        new_status = check_irc_opt_status(job_id)
        if new_status == "off_queue":
                # 3. check job content
                new_status = check_irc_opt_content_status(target['path'])

                # 4. check with original status which
                # should be job_launched or job_running
                # if any difference update status
                irc_status = 'irc_{}'.format(str('forward'))
                orig_status = target[irc_status]
                if orig_status != new_status:
                    update_field = {
                                       irc_status: new_status,
                                   }
                    qm_collection.update_one(target, {"$set": update_field}, True)
    # 1. select jobs to check
    targets = select_irc_opt_target(direction = 'reverse')
    irc_opt_jobid = 'irc_{}_opt_jobid'.format(str('reverse'))
    qm_collection = db['qm_calculate_center']
    
    # 2. check the job pbs_status
    for target in targets:
        job_id = target[irc_opt_jobid]
        new_status = check_irc_opt_status(job_id)
        if new_status == "off_queue":
                # 3. check job content
                new_status = check_irc_opt_content_status(target['path'])

                # 4. check with original status which
                # should be job_launched or job_running
                # if any difference update status
                irc_status = 'irc_{}'.format(str('reverse'))
                orig_status = target[irc_status]
                if orig_status != new_status:
                    update_field = {
                                       irc_status: new_status,
                                   }
                    qm_collection.update_one(target, {"$set": update_field}, True)
                    
check_energy_jobs()
check_ssm_jobs()
check_ts_jobs()
check_irc_jobs()
check_irc_equal()
check_irc_opt_job()
#check_same_ssm_jobs()
#check_same_ts_jobs()
#check_same_barrier()


qm_collection = db['qm_calculate_center']
max_gen = qm_collection.find_one(sort=[("generations", -1)])
max_gen = max_gen['generations']

check_ard_barrier_jobs(max_gen)
