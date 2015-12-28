import math
import numpy as np

def template2str(filename, mapping):
    """
    read in a template file and substitute variables using mapping, return a str
    @param filename str
           absolute path and filename
    @param mapping dict
           for example. {'username','junyi'}
           key must be str; value can be str or type that can be converted to str
    @return subsitituted_file str
    """
    from string import Template
    filein  = open(filename)
    filestr = filein.read()

    for ikey in mapping:
        mapping[ikey] = str(mapping[ikey])

    file_temp = Template(filestr)
    subsitituted_file = file_temp.substitute(mapping)
    return subsitituted_file

def write2file(filename, str2write):
    """
    write str in `str2write` to file `filename`
    """
    file_handle = file(filename,'w')
    file_handle.write(str2write)
    file_handle.close()
    return

template_dir = '../run_env/'
f1_dir      = '/lustre/f1/Junyi.Chai/'
exe_file    = '/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/bld/qg_run.x'

#--------------- experiment parameters -------------------
exp_name    = 'Dec12_kf16_drag8e-4_m16'
num_procs   = 64
walltime    = '02:30:00'

#resolution: default kmax=511
kmax = 511 
Rh       = 1./8. # Rhine scale
F        = 0.0 
beta     = 16.0*np.pi 
kf_min   = 14 
kf_max   = 18 

bot_drag_base = 8.*1e-4 
epsilon = bot_drag_base* (beta)**2 * Rh**4
bot_drag = bot_drag_base*16
e_o      = 1e-4
rmf_norm_min = 1.e-7

#---------------------------------------------------------

exp_dir = f1_dir + exp_name + '/'

print "Bottom friction = ", bot_drag 
print "Beta = ", beta 
print "Rhine scale", Rh
print "Epsilon = ", epsilon
print "Expriment will be ran at ", exp_dir

exp_dir = f1_dir + exp_name + '/'

runscript_template = template_dir + 'runscript.template'
transfer_template  = template_dir + 'transferScript.template'
input_template     = template_dir + 'input_1layer_wtracer.nml.template'

runscript_params = {'num_procs':num_procs, 'exp_name':exp_name, 'walltime':walltime}
transfer_params  = {'exp_name':exp_name}
input_params     = {'F':F, 'beta':beta, 'bot_drag':bot_drag, 'kmax':kmax, 
                    'epsilon':epsilon, 'kf_min':kf_min, 'kf_max':kf_max,
                    'e_o':e_o, 'rmf_norm_min':rmf_norm_min}

runscript_str = template2str(runscript_template, runscript_params)
transfer_str  = template2str(transfer_template,  transfer_params)
input_str     = template2str(input_template,     input_params)

import os
import stat
os.mkdir(exp_dir)
READ_WRITE_EXE = stat.S_IREAD|stat.S_IEXEC|stat.S_IWRITE
write2file(exp_dir + 'runscript', runscript_str)
os.chmod(  exp_dir + 'runscript', READ_WRITE_EXE)
write2file(exp_dir + 'transferScript', transfer_str)
os.chmod(  exp_dir + 'transferScript', READ_WRITE_EXE)
write2file(exp_dir + 'input.nml', input_str)
os.chmod(  exp_dir + 'input.nml', READ_WRITE_EXE)

import shutil
shutil.copyfile(template_dir + 'restartNum', exp_dir + 'restartNum')
shutil.copyfile(exe_file, exp_dir + 'qg_run.x')
os.chmod(template_dir + 'restartNum', READ_WRITE_EXE)
os.chmod(exp_dir + 'qg_run.x', READ_WRITE_EXE)
os.mkdir(exp_dir + 'INPUT/')
os.mkdir(exp_dir + 'RESTART/')

print 'Run environment preparation success'
