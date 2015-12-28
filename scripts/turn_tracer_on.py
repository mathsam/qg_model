"""
Turn on tracer from a restart run
"""
exp_name    = 'Dec12_kf64_drag1e-4'

#a dictionary with key to be the regex to match and value to be the line to
#replace with
match_and_replace_dict = {'\s+USE_TRACER_Y\s+=': 'USE_TRACER_Y=T,',
    '\s+RESET_TRACER\s+='     : 'RESET_TRACER=T,',
    '\s+USE_MEAN_GRAD_T\s+='  : 'USE_MEAN_GRAD_T=T,',
    '\s+NZT\s+='              : 'NZT=1,',
    '\s+MAXMODE\s+='          : 'MAXMODE=0,',
    '\s+FILTER_TYPE_T\s+='    : "FILTER_TYPE_T='hyperviscous',",
    '\s+FILTER_EXP_T\s+='     : 'FILTER_EXP_T=4.00000000000000,'}

num_procs   = 64 
walltime    = '02:30:00'

f1_dir    = '/lustre/f1/Junyi.Chai/'
template_dir = '../run_env/'


import os
import stat
import shutil
import re
import math
exp_dir = os.path.join(f1_dir, exp_name) + '/'

#back up current restart.nml 
shutil.copyfile(exp_dir + 'restart.nml', exp_dir + '_old_restart.nml')

def match_and_replace(filename, match_and_replace_dict):
    """find line in file `filename` that matches the regular expression
    `regex2match` and replace this line with `line2replacew`

    Args:
        filename: string for absolute path for file
        match_and_replace_dict: a dict with key to be a string for regular 
            expression to match and value for the string to be replaced with.
            For example:
                {'\s+USE_TRACER_Y\s+=': 'USE_TRACER_Y=T,'}

    Return:
        newfile: string for the new file
    """
    f = open(filename, 'r')
    newfile = ''
    for line in f:
        regex_found = False
        for regex2match in match_and_replace_dict:
            if not regex_found and re.match(regex2match, line):
            	newfile += match_and_replace_dict[regex2match] + '\n'
                regex_found = True
        if not regex_found:
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
runscript_template = template_dir + 'runscript.template'
runscript_params = {'num_procs':num_procs, 'exp_name':exp_name, 'walltime':walltime}

runscript_str = template2str(runscript_template, runscript_params)
new_restart_nml = match_and_replace(exp_dir+'/restart.nml', match_and_replace_dict)

READ_WRITE_EXE = stat.S_IREAD|stat.S_IEXEC|stat.S_IWRITE
READ_WRITE     = stat.S_IREAD|stat.S_IWRITE
write2file(exp_dir + 'runscript', runscript_str)
os.chmod(  exp_dir + 'runscript', READ_WRITE_EXE)
write2file(exp_dir + 'restart.nml', new_restart_nml)
os.chmod(  exp_dir + 'restart.nml', READ_WRITE)

print 'Prepared new restart.nml succeed %s' %exp_name
