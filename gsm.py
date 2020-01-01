import subprocess
import os

#get the previous dir path
abspath_pardir = os.path.abspath(os.pardir)
abspath_current = os.getcwd()
#run ard job
def run_ard():
    ard_path = abspath_pardir +'/'+ 'ard.py' 
    input_path = abspath_pardir +'/'+ 'input.txt' 
    p = subprocess.Popen(['python', ard_path, input_path])
    p.wait()
#run gsm job
def run_GSM():
    reactions_path = abspath_pardir + '/' + 'reactions'
    dirs = os.listdir(reactions_path)
    gsm_path = abspath_pardir +'/pyGSM/examples/QChem/DE_GSM/'

    for i in dirs:
        gsm_input_path = reactions_path + '/' + i + '/'
	run_cmd = 'gsm -xyzfile '+gsm_input_path+'input.xyz '+'-mode DE_GSM '+'-package QChem '+'-lot_inp_file '+gsm_path+'qstart '+'-ID 1'
        p =subprocess.Popen([run_cmd],shell=True)
        p.wait()
	#clear scratch
	clear_cmd='rm -rf scratch'
	p=subprocess.Popen([clear_cmd],shell=True)
	p.wait()
        mv_cmd = 'mv '+ abspath_current +'/'+'* '+gsm_input_path
        p =subprocess.Popen([mv_cmd],shell=True)
        p.wait()
#run ts search job
#run irc job

run_ard()
#run_GSM()
