from connect import db
import subprocess
import os
from os import path
import sys
from openbabel import pybel

def select_ard_target():
    """
    This method is to inform job checker which targets 
    to check, which need meet one requirement:
    1. status is job_launched or job_running
    Returns a list of targe
    """
    qm_collection = db['qm_calculate_center']
    reg_query = {"ard_status":
                    {"$in": 
                        ["job_launched", "job_running", "job_queueing"] 
                    }
                }
    targets = list(qm_collection.find(reg_query))

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
    elif stdout[idx+2] == 'Q':
        return "job_queueing"
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

    qm_collection = db['qm_calculate_center']
    
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

            qm_collection.update_one(target, {"$set": update_field}, True)


def print_information(generations):
    """
    For a given generations (int) print the information in database
    """
    qm_collection = db['qm_calculate_center']
    
    gen_query = {"generations":
                    {"$in": 
                        [generations] 
                    }
                }
    opt_query_1 = {'$and': 
                    [
                    { "opt_status":
                        {"$in":
                        ["job_running", "job_queueing"]}
                        },
                    {'generations':generations}
                    ]
                }
            
    opt_query_2 = {'$and': 
                    [
                    { "opt_status":
                        {"$in":
                        ["job_success"]}
                        },
                    {'generations':generations}
                    ]
                }
    opt_query_3 = {'$and': 
                    [
                    { "opt_status":
                        {"$in":['job_fail']}
                        },
                    {'generations':generations}
                    ]
                }
    ssm_query_1 = {'$and': 
                    [
                    { "ssm_status":
                        {"$in":
                        ["job_running", "job_queueing"]}
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
                        ['job_running', "job_queueing"]}
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
    irc_query_1 = {'$and': 
                    [
                    { "irc_forward_status":
                        {"$in":
                        ['job_running', "job_queueing"]}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_2 = {'$and': 
                    [
                    { "irc_reverse_status":
                        {"$in":
                        ['job_running', "job_queueing"]}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_3 = {'$and': 
                    [
                    { "irc_forward_status":
                        {"$in":
                        ['job_success', 'opt_success']}
                        },
                    { "irc_reverse_status":
                        {"$in":
                        ['job_success', 'opt_success']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_4 = {'$and': 
                    [
                    { "irc_forward_status":
                        {"$in":
                        ['Bad initial gradient', "Failed line search", 'Error in gen_scfman', 'unknown fail information']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_5 = {'$and': 
                    [
                    { "irc_reverse_status":
                        {"$in":
                        ['Bad initial gradient', "Failed line search", 'Error in gen_scfman', 'unknown fail information']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_6 = {'$and': 
                    [
                    { "irc_forward_status":
                        {"$in":
                        ['need opt']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_7 = {'$and': 
                    [
                    { "irc_reverse_status":
                        {"$in":
                        ['need opt']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_8 = {'$and': 
                    [
                    { "irc_equal":
                        {"$in":
                        ['forward equal to reactant and reverse equal to product', 'reverse equal to reactant and forward equal to product']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_9 = {'$and': 
                    [
                    { "irc_equal":
                        {"$in":
                        ['reverse equal to reactant but forward does not equal to product', 'reverse does not equal to reactant but forward equal to product', 
                        'forward equal to reactant but reverse does not equal to product', 'forward does not equal to reactant but reverse equal to product']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_query_10 = {'$and': 
                    [
                    { "irc_equal":
                        {"$in":
                        ['waiting for check']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_opt_query_1 = {'$and': 
                    [
                    { "opt_forward_status":
                        {"$in":
                        ["opt_job_running", "opt_job_queueing"]}
                        },
                    {'generations':generations}
                    ]
                }
    irc_opt_query_2 = {'$and': 
                    [
                    { "opt_reverse_status":
                        {"$in":
                        ["opt_job_running", "opt_job_queueing"]}
                        },
                    {'generations':generations}
                    ]
                }
    irc_opt_query_3 = {'$and': 
                    [
                    { "opt_forward_status":
                        {"$in":
                        ['job_success']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_opt_query_4 = {'$and': 
                    [
                    { "opt_reverse_status":
                        {"$in":
                        ['job_success']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_opt_query_5 = {'$and': 
                    [
                    { "opt_forward_status":
                        {"$in":
                        ['job_fail']}
                        },
                    {'generations':generations}
                    ]
                }
    irc_opt_query_6 = {'$and': 
                    [
                    { "opt_reverse_status":
                        {"$in":
                        ['job_fail']}
                        },
                    {'generations':generations}
                    ]
                }
    gen_targets = list(qm_collection.find(gen_query))
    opt_targets_1 = list(qm_collection.find(opt_query_1))
    opt_targets_2 = list(qm_collection.find(opt_query_2))
    opt_targets_3 = list(qm_collection.find(opt_query_3))
    ssm_targets_1 = list(qm_collection.find(ssm_query_1))
    ssm_targets_2 = list(qm_collection.find(ssm_query_2))
    ssm_targets_3 = list(qm_collection.find(ssm_query_3))
    ts_targets_1 = list(qm_collection.find(ts_query_1))
    ts_targets_2 = list(qm_collection.find(ts_query_2))
    ts_targets_3 = list(qm_collection.find(ts_query_3))
    irc_targets_1 = list(qm_collection.find(irc_query_1))
    irc_targets_2 = list(qm_collection.find(irc_query_2))
    irc_targets_3 = list(qm_collection.find(irc_query_3))
    irc_targets_4 = list(qm_collection.find(irc_query_4))
    irc_targets_5 = list(qm_collection.find(irc_query_5))
    irc_targets_6 = list(qm_collection.find(irc_query_6))
    irc_targets_7 = list(qm_collection.find(irc_query_7))
    irc_targets_8 = list(qm_collection.find(irc_query_8))
    irc_targets_9 = list(qm_collection.find(irc_query_9))
    irc_targets_10 = list(qm_collection.find(irc_query_10))
    irc_opt_targets_1 = list(qm_collection.find(irc_opt_query_1))
    irc_opt_targets_2 = list(qm_collection.find(irc_opt_query_2))
    irc_opt_targets_3 = list(qm_collection.find(irc_opt_query_3))
    irc_opt_targets_4 = list(qm_collection.find(irc_opt_query_4))
    irc_opt_targets_5 = list(qm_collection.find(irc_opt_query_5))
    irc_opt_targets_6 = list(qm_collection.find(irc_opt_query_6))
    reactant_target = list(qm_collection.find({'generations':generations}))
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
    print('{} nodes running or queueing pre-OPT'.format(len(opt_targets_1)))
    print('{} nodes success in pre-OPT'.format(len(opt_targets_2)))
    print('{} nodes fail in pre-OPT'.format(len(opt_targets_3)))
    print('{} nodes running or queueing SSM'.format(len(ssm_targets_1)))
    print('{} nodes success in SSM'.format(len(ssm_targets_2)))
    print('{} nodes fail in SSM'.format(len(ssm_targets_3)))
    print('{} nodes running or queueing TS'.format(len(ts_targets_1)))
    print('{} nodes success in TS'.format(len(ts_targets_2)))
    print('{} nodes fail in TS'.format(len(ts_targets_3)))
    print('{} nodes running or queueing IRC'.format(len(irc_targets_1) + len(irc_targets_2)))
    print('{} nodes success in IRC (including after opt)'.format(len(irc_targets_3) * 2))
    print('{} nodes fail in OPT'.format(len(irc_targets_4) + len(irc_targets_5)))
    print('{} nodes need to OPT (IRC)'.format(len(irc_targets_6) + len(irc_targets_7)))
    print('{} nodes are waiting for checking irc equal'.format(len(irc_targets_10)))
    print('{} nodes are intended'.format(len(irc_targets_8)))
    print('{} nodes are unintended'.format(len(irc_targets_9)))
    print('{} nodes are running or queueing irc opt job'.format(len(irc_opt_targets_1) + len(irc_opt_targets_2)))
    print('{} nodes are success irc opt job'.format(len(irc_opt_targets_3) + len(irc_opt_targets_4)))
    print('{} nodes are fail irc opt job'.format(len(irc_opt_targets_5) + len(irc_opt_targets_6)))
    print('-----------------------------------------')


def update_network_status():
    status_collection = db['status']
    qm_collection = db['qm_calculate_center']
    statistics_collection = db['statistics']
    ard_had_add_number = qm_collection.count_documents({})
    ard_should_add_number = sum(statistics_collection.distinct("add how many products"))
        
    running_query = {"ard_status":
                    {"$in": 
                        ["job_unrun", "job_launched", "job_running", "job_queueing"] 
                    }
                }
    ard_nodes = list(qm_collection.find(running_query))

    energy_query = {"energy_status":
                    {"$in":
                        ["job_launched", "job_running", "job_queueing", 'job_unrun']
                    }
                }
    ssm_query = {"ssm_status":
                    {"$in":
                        ["job_launched", "job_running", "job_queueing", 'job_unrun']
                    }
                }
    opt_query = {"opt_status":
                    {"$in":
                        ["job_launched", "job_running", "job_queueing", 'job_unrun']
                    }
                }
    ts_query = {"ts_status":
                    {"$in":
                        ["job_launched", "job_running", "job_queueing", 'job_unrun']
                    }
                }
    irc_query_1 = {"irc_forward_status":
                    {"$in":
                        ["job_launched", "job_running", "job_queueing", "need opt", 'job_unrun']
                    }
                }
    irc_query_2 = {"irc_reverse_status":
                    {"$in":
                        ["job_launched", "job_running", "job_queueing", "need opt", 'job_unrun']
                    }
                }
    opt_query_1 = {"opt_forward_status":
                    {"$in":
                        ["opt_job_launched", "opt_job_running", "opt_job_queueing", 'job_unrun']
                    }
                }
    opt_query_2 = {"opt_reverse_status":
                    {"$in":
                        ["opt_job_launched", "opt_job_running", "opt_job_queueing", 'job_unrun']
                    }
                }
    not_finished_number = len(list(qm_collection.find(energy_query))) + len(list(qm_collection.find(ssm_query))) + len(list(qm_collection.find(ts_query))) + len(list(qm_collection.find(irc_query_1))) + len(list(qm_collection.find(irc_query_2))) + len(list(qm_collection.find(opt_query_1))) + len(list(qm_collection.find(opt_query_2))) + len(list(qm_collection.find(opt_query)))


    if ard_had_add_number == ard_should_add_number and len(ard_nodes) == 0 and not_finished_number == 0:
        print('Network converged')
        
        target = list(status_collection.find({}))
        update_field = {
                        'status':'Network converged'
                    }
        status_collection.update_one(target[0], {"$set": update_field}, True)




check_ard_jobs()

qm_collection = db['qm_calculate_center']
max_gen = qm_collection.find_one(sort=[("generations", -1)])
max_gen = max_gen['generations']
for i in range(max_gen):
    print_information(i+1)
    
update_network_status()