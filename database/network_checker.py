from connect import db
import subprocess
import os
from os import path
import sys
import pybel

def select_ard_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    collect = db['qm_calculate_center']
    reg_query = {"ard_status":
                    {"$in": 
                        ["job_launched", "job_running"] 
                    }
                }
    targets = list(collect.find(reg_query))

    return targets

def check_ard_status(job_id):
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
    else:
        return "job_launched"

    
def check_ard_jobs():
    """
    This method checks job with following steps:
    1. select jobs to check
    2. check the job pbs-status, e.g., qstat -f "job_id"
    3. check job content
    4. update with new status
    """
    # 1. select jobs to check
    targets = select_ard_target()

    collect = db['qm_calculate_center']
    
    # 2. check the job pbs status
    for target in targets:
        job_id = target['ard_jobid']
        # 2. check the job pbs status
        new_status = check_ard_status(job_id)
        if new_status == "off_queue":
            # 3. check job content
            new_status = 'job_finished'

        # 4. check with original status which
        # should be job_launched or job_running
        # if any difference update status
        orig_status = target['ard_status']
        if orig_status != new_status:

            if new_status == 'job_finished':
                update_field = {
                                'ard_status': new_status
                            }
            else:
                update_field = {
                                'ard_status': new_status
                            }

            collect.update_one(target, {"$set": update_field}, True)

def print_information(generations):
    """
    For a given generations (int) print the information in database
    """
    collect = db['qm_calculate_center']
    
    gen_query = {"generations":
                    {"$in": 
                        [generations] 
                    }
                }
    ssm_query_1 = {'$and': 
                    [
                    { "ssm_status":
                        {"$in":
                        ["job_running"]}
                        },
                    {'generations':generations}
                    ]
                }
            
    ssm_query_2 = {'$and': 
                    [
                    { "ssm_status":
                        {"$in":
                        ["job_success"]}
                        },
                    {'generations':generations}
                    ]
                }
    ssm_query_3 = {'$and': 
                    [
                    { "ssm_status":
                        {"$in":
                        ['job_fail', 'total dissociation', 'Exiting early']}
                        },
                    {'generations':generations}
                    ]
                }
    ts_query_1 = {'$and': 
                    [
                    { "ts_status":
                        {"$in":
                        ['job_running']}
                        },
                    {'generations':generations}
                    ]
                }
    ts_query_2 = {'$and': 
                    [
                    { "ts_status":
                        {"$in":
                        ['job_success']}
                        },
                    {'generations':generations}
                    ]
                }
    ts_query_3 = {'$and': 
                    [
                    { "ts_status":
                        {"$in":
                        ['job_fail']}
                        },
                    {'generations':generations}
                    ]
                }
    gen_targets = list(collect.find(gen_query))
    ssm_targets_1 = list(collect.find(ssm_query_1))
    ssm_targets_2 = list(collect.find(ssm_query_2))
    ssm_targets_3 = list(collect.find(ssm_query_3))
    ts_targets_1 = list(collect.find(ts_query_1))
    ts_targets_2 = list(collect.find(ts_query_2))
    ts_targets_3 = list(collect.find(ts_query_3))
    reactant_target = list(collect.find({'generations':generations}))
    smi = []
    for target in reactant_target:
        reactant_smi = target['Reactant SMILES']
        smi.append(reactant_smi)
    smi = set(smi)
    print('Extracting information from database....')
    print('-----------------------------------------')
    print('Generations : {}'.format(generations))
    print('Reactant SMILES : {}'.format(smi))
    print('Nodes : {}'.format(len(gen_targets)))
    print('{} nodes running SSM'.format(len(ssm_targets_1)))
    print('{} nodes success in SSM'.format(len(ssm_targets_2)))
    print('{} nodes fail in SSM'.format(len(ssm_targets_3)))
    print('{} nodes running TS'.format(len(ts_targets_1)))
    print('{} nodes success in TS'.format(len(ts_targets_2)))
    print('{} nodes fail in TS'.format(len(ts_targets_3)))
    print('-----------------------------------------')


def update_network_status():
    
    status = db['status']
    qm_cal_center = db['qm_calculate_center']
    statistics = db['statistics']
    total_nodes = list(statistics.find({}, {'add how many products':1}))
    qm_nodes = list(qm_cal_center.find({}))
    count = 0
    
    for i in total_nodes:
        count += i['add how many products']
        
    query = {'$or': 
                    [
                    {'ts_status':
                        {'$in':
                        ['job_fail', 'job_success']}
                        },
                    {'ssm_status':
                        {'$in':
                            ['Exiting early', 'total dissociation', 'job_fail']
                        }}
                    ]
                }
    targets = list(qm_cal_center.find(query))
    
    if len(targets) == count:
        print('Network converged')
        
        target = list(status.find({}))
        update_field = {
                        'status':'Network converged'
                    }
        status.update_one(target[0], {"$set": update_field}, True)
    
    
    
    
check_ard_jobs()

qm_cal_center = db['qm_calculate_center']
max_gen = qm_cal_center.find_one(sort=[("generations", -1)])
max_gen = max_gen['generations']
for i in range(max_gen):
    print_information(i+1)
    
update_network_status()
