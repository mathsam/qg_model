old_exp_name    = 'Dec12_kf64_drag64e-4'

#match the line with regex `regex2match
f1_dir    = '/lustre/f1/Junyi.Chai/'
f1u_dir   = '/lustre/f1/unswept/Junyi.Chai/RestartFiles/QG_model/2015/'

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
    print "created bk dir %s" %bk_dir
shutil.copyfile(old_exp_dir + 'INPUT/restart.nc', bk_dir + 'restart.nc')
if os.path.isfile(old_exp_dir + 'input.nml'):
    shutil.copyfile(old_exp_dir + 'input.nml',  bk_dir + 'input.nml') 
shutil.copyfile(old_exp_dir + 'restart.nml',    bk_dir + 'restart.nml')
shutil.copyfile(old_exp_dir + 'runscript',      bk_dir + 'runscript')
shutil.copyfile(old_exp_dir + 'transferScript', bk_dir + 'transferScript')
shutil.copyfile(old_exp_dir + 'restartNum',     bk_dir + 'restartNum')
print 'backed up last environments in %s' %bk_dir
