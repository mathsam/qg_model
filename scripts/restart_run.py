old_exp_name    = 'Nov5_Sc1.8_drag5e-3'
new_exp_name    = 'Nov5_Sc1.8_drag5e-4' # this is the new experiment name

#match the line with regex `regex2match` and replace the entire line with
#`line2replacew`
regex2match   = '\s+BOT_DRAG\s+='
line2replacew = 'BOT_DRAG        =  3.97887357730000E-2,'

num_procs   = 64 
walltime    = '02:00:00'

f1_dir    = '/lustre/f1/Junyi.Chai/'
f1u_dir   = '/lustre/f1/unswept/Junyi.Chai/RestartFiles/QG_model/2015/'
template_dir = '../run_env/'
exe_file    = '/lustre/f1/unswept/Junyi.Chai/qg_model_mpi/bld/qg_run.x'


import os
import stat
import shutil
import re
import math
old_exp_dir = os.path.join(f1_dir, old_exp_name) + '/'
bk_dir  = os.path.join(f1u_dir,old_exp_name) + '/'

#back up current experiment
if not os.path.exists(bk_dir):
    os.mkdir(bk_dir)
shutil.copyfile(old_exp_dir + 'INPUT/restart.nc', bk_dir + 'restart.nc')
if os.path.isfile(old_exp_dir + 'input.nml'):
    shutil.copyfile(old_exp_dir + 'input.nml',   bk_dir + 'input.nml') 
shutil.copyfile(old_exp_dir + 'restart.nml', bk_dir + 'restart.nml')
shutil.copyfile(old_exp_dir + 'runscript',   bk_dir + 'runscript')
shutil.copyfile(old_exp_dir + 'transferScript', bk_dir + 'transferScript')
print 'backed up last environments in %s' %bk_dir

def match_and_replace(filename, regex2match, line2replacew):
    """find line in file `filename` that matches the regular expression
    `regex2match` and replace this line with `line2replacew`

    Args:
        filename: string for absolute path for file
        regex2match: string for regular expression to match. For example:
                     '\s+BOT_DRAG\s+='
        line2replace:string to replace the line that matches the above regex

    Return:
        newfile: string for the new file
    """
    f = open(filename, 'r')
    newfile = ''
    for line in f:
        if re.match(regex2match, line):
            newfile += line2replacew + '\n'
        else:
            newfile += line

    return newfile


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

#---------------------------------------------------------
new_exp_dir = f1_dir + new_exp_name + '/'

runscript_template = template_dir + 'runscript.template'
transfer_template  = template_dir + 'transferScript.template'
input_template     = template_dir + 'input.nml.template'

runscript_params = {'num_procs':num_procs, 'exp_name':new_exp_name, 'walltime':walltime}
transfer_params  = {'exp_name':new_exp_name}

runscript_str = template2str(runscript_template, runscript_params)
transfer_str  = template2str(transfer_template,  transfer_params)
new_restart_nml = match_and_replace(old_exp_dir+'/restart.nml', regex2match, line2replacew)

os.mkdir(new_exp_dir)
READ_WRITE_EXE = stat.S_IREAD|stat.S_IEXEC|stat.S_IWRITE
READ_WRITE     = stat.S_IREAD|stat.S_IWRITE
write2file(new_exp_dir + 'runscript', runscript_str)
os.chmod(  new_exp_dir + 'runscript', READ_WRITE_EXE)
write2file(new_exp_dir + 'transferScript', transfer_str)
os.chmod(  new_exp_dir + 'transferScript', READ_WRITE_EXE)
write2file(new_exp_dir + 'restart.nml', new_restart_nml)
os.chmod(  new_exp_dir + 'restart.nml', READ_WRITE)

shutil.copyfile(template_dir + 'restartNum', new_exp_dir + 'restartNum')
shutil.copyfile(exe_file, new_exp_dir + 'qg_run.x')
os.chmod(template_dir + 'restartNum', READ_WRITE_EXE)
os.chmod(new_exp_dir + 'qg_run.x', READ_WRITE_EXE)
os.mkdir(new_exp_dir + 'INPUT/')
os.mkdir(new_exp_dir + 'RESTART/')
shutil.copyfile(old_exp_dir + 'INPUT/restart.nc', new_exp_dir + 'INPUT/restart.nc')

print 'New run environment preparation success at %s' %new_exp_dir
