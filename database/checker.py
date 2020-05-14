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

def select_energy_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """

    collect = db['molecules']
    reg_query = {"energy_status":
                    {"$in":
                        ["job_launched", "job_running"]
                    }
                }
    targets = list(collect.find(reg_query))

    return targets

def check_energy_status(job_id):
    """
    This method checks slurm status of a job given job_id
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
    else:
        return "job_launched"
    
def check_energy_content_status(path):

    energy_path = path.join(path, "energy.out")
    if not path.exists(energy_path):
        return "job_aborted"
    else:
        with open(energy_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.strip().startswith('SCF   energy in the final basis set'):
                    SCF_Energy = line.split()[-1]
                    energy = float(SCF_Energy)*627.5095
        return "job_success", energy

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
    collect = db['molecules']
    
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
                                       'SCF_energy':energy
                                   }
                    collect.update_one(target, {"$set": update_field}, True)

def select_ssm_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    collect = db['molecules']
    reg_query = {"ssm_status":
                    {"$in": 
                        ["job_launched", "job_running"] 
                    }
                }
    targets = list(collect.find(reg_query))

    return targets

def check_ssm_status(job_id):
    """
    This method checks slurm status of a job given job_id
    Returns off_queue or job_launched or job_running
    """
    commands = ['qstat', '-f', job_id]
    process = subprocess.Popen(commands,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if "Unknown Job Id" in stderr.decode():
        return "off_queue"
    assert "JobId={0}".format(job_id) in stdout, 'PBS cannot show details for job_id {0}'.format(job_id)

    # in pbs stdout is byte, so we need to decode it at first.
    stdout = stdout.decode().strip().split()
    idx = stdout.index('job_state')
    if stdout[idx+2] == 'R':
        return "job_running"
    else:
        return "job_launched"
    
def check_ssm_content_status(target_path):

    ssm_path = path.join(target_path, 'SSM')
    ssm_pic_path = path.join(ssm_path, '0000_string.png')
    if not path.exists(ssm_pic_path):
        return "job_fail"
    else:
        if check_ssm(ssm_path):
            return "Exiting early"
        else:
            generate_ssm_product_xyz(ssm_path)
            return "job_success"

def check_ssm(ssm_path):
    status_path = path.join(target_path, 'status.log')
    
    with open(status_path, 'r') as f:
        f.seek(0, 2)
        fsize = f.tell()
        f.seek(max(fsize - 1024, 0), 0)  # Read last 1 kB of file
        lines = f.readlines()
        
    for idx, i in enumerate(lines):
        if i.startswith('Finished GSM!'):
            break
    
    if lines[idx-1] == 'Exiting early\n'
        return 1
    else:
        return 0
    
    
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

    collect = db['molecules']
    
    # 2. check the job slurm-status
    for target in targets:
        job_id = target['ssm_jobid']
        # 2. check the job slurm_status
        new_status = check_ssm_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = check_ssm_content_status(target['path'])

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ssm_status']
        if orig_status != new_status:
            update_field = {
                            'ssm_status': new_status
                           }

            collect.update_one(target, {"$set": update_field}, True)
            # update reactions collection information
            collect1 = db['reactions']
            dir_name = target['dir']
            reg_query = {"for_ssm_check":dir_name}
            tt = collect1.find_one(reg_query)
            collect1.update_one(tt, {"$set": update_field}, True)
	

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
    parent_ssm_product_path  = path.join(path.abspath(os.path.join(target_path, '../')), 'ssm_product.xyz')
    with open(parent_ssm_product_path,'w') as q:
        q.write('{}\n{}'.format(lines[idx-1], ''.join(lines[idx+1:])))

#check_energy_jobs()
check_ssm_jobs()