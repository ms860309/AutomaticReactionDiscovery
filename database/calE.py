from connect import Connector
from subprocess import Popen

def select_calE_target():
    remote_client = Connector().client
    db = remote_client['network']
    collect = db['molecules']
    reg_query = {"energy_status":"job_unrun"}
    targets = list(collect.find(reg_query))

    selected_targets = []
    for target in targets:
        dir_path = target['path']
        selected_targets.append(dir_path)

    return selected_targets

def push_jobs():
    targets = select_calE_target()
    
    for target in targets:
        if os.path.exists(target):
            os.chdir(target)
        subfile = create_sub_file(target)
        try:
            cmd = 'qsub {}'.format(subfile)
            p = Popen([cmd], shell = True)
            p.wait()
        except:
            continue
        
        update_status(target)
        
def create_sub_file(path, ncpus = 1, mpiprocs = 1):
    subfile = os.path.join(path, 'calculate_energy.job')
    
    shell = '#!/usr/bin/bash'
    pbs_setting = '#PBS -l select=1:ncpus={}:mpiprocs={}'.format(ncpus, mpiprocs)
    target_path = 'cd {}'.format(path)
    nes1 = 'source ~/.bashrc_qchem'
    nes2 = 'export QCSCRATCH=/tmp/ypli/$PBS_JOBID'
    nes3 = 'mkdir -p $QCSCRATCH'
    nes4 = 'qchem -nt 2 energy.in energy.out'
    nes5 = 'rm -r $QCSCRATCH'
    with open(subfile, 'w') as f:
        f.write('{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}').format(shell, pbs_setting, target_path, nes1, nes2, nes3, nes4, nes5)
    return subfile
    
def update_status(target):
    remote_client = Connector().client
    db = remote_client['network']
    collect = db['molecules']
    reg_query = {"path":target}
    collect.update(reg_query,{"energy_status":"finished"})
    
