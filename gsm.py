import subprocess
import os

abspath = os.getcwd()
#run ard job
def run_ard():
    ard_path = abspath +'/'+ 'ard.py' 
    input_path = abspath +'/'+ 'input.txt' 
    p = subprocess.Popen(['python', ard_path, input_path])
    p.wait()
#run gsm job
def run_GSM():
    reactions_path = abspath + '/' + 'reactions'
    dirs = os.listdir(reactions_path)
    gsm_path = abspath +'/pyGSM/examples/QChem/DE_GSM/'

    for i in dirs:
        gsm_input_path = reactions_path + '/' + i + '/'
        p =subprocess.Popen(['gsm -xyzfile ' + gsm_input_path +'input.xyz ' +
                                    '-mode DE_GSM ' +
                                   '-package QChem ' + 
                                   '-lot_inp_file ' + gsm_path + 'qstart ' + 
                                   '-ID 1'
                                   ])
        p.wait()

#run ts search job
#run irc job
p = subprocess.Popen('ls -l', shell=True)
p.wait()

run_ard()
run_GSM()