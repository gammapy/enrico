"""Helper functions to call command line tools from python

@todo:
- Make it flexible so that a dict(a=42) can be given and
'a=42' or '-a 42' or '--a 42' is generated depending on an option!
- Make it a class, so that options like dryrun or options formattting
don't have to be handed around all the time!!!"""

import sys
import logging
import subprocess
import os
#import stat
from os.path import join, dirname
import commands
import time
import datetime
import tempfile
from enrico import environ

def _cmd_to_str(cmd):
    """ Convert entries to strings and join with spaces """
    cmd_string = map(str, cmd)
    return ' '.join(cmd_string)

def _options_to_str(options):
    """ Convert a dictionary of options to a string """
    string = ''
    for key, value in options.items():
        string += ' {0}={1}'.format(key, value)
    return string

def jobs_in_queue():
    """ Returns the number of jobs this user has in the queue """
    import time
    from subprocess import Popen, PIPE
    user = os.environ['USER']
    for trial in xrange(60):
        try:
            fh = Popen("qstat -u {user}".format(user=user), \
                    stdout=PIPE, shell=True)
        except OSError:
            # Wait 5 minutes and try again
            time.sleep(300)
            if trial==59:
                raise
        else:
            break

    njobs = len(fh.stdout.readlines())
    # If there are no jobs we will get 0 lines.
    # If there are jobs there will be two extra header lines.
    # So this works for both cases:
    return max(0, njobs - 2)

def wait_for_slot(max_jobs):
    """Wait until you have less that max_jobs in the queue"""
    njobs = jobs_in_queue()
    while not (njobs < max_jobs):
        time_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        logging.info('{0}, njobs = {1}, max_jobs = {2}'
                     ''.format(time_str, njobs, max_jobs))
        time.sleep(10) # 10 seconds
        njobs = jobs_in_queue()

##Function to chose the Farm commands
def GetSubCmd():
  queuetext = ""
  if [ environ.QUEUE != "" ]:
      queueoptions = "%s " %(environ.TORQUE_RESOURCES)
      queuetext = "-q %s" %(environ.QUEUE)
  cmd = {'LAPP' :    ['qsub -V','-l mem=4096mb'],
         'MPIK' :    ['qsub'],
         'DESY' :    ['qsub','-R y -V -terse -l h_rss=20G -V %s %s'%(queueoptions,queuetext)],
         'DESY_quick' : ['qsub','-js 10 -R y -V -terse -l h_rss=4G -l s_cpu=01:00:00 -l h_cpu=02:00:00 -V %s %s'%(queueoptions,queuetext)],
         'LOCAL' :   ['qsub','-V %s %s'%(queueoptions,queuetext)],
         'CCIN2P3' : ['qsub','-l ct=24:00:00 -l vmem=4G -l fsize=20G -l sps=1 -l os=sl6 -P P_hess']}
  return cmd[environ.FARM]

def GetSubOutput(qsub_log):
  cmd = {'LAPP' :    ['-o', qsub_log, '-j', 'oe'],
         'MPIK' :    ['-o', qsub_log, '-j', 'y'],
         'LOCAL' :   ['-o', qsub_log, '-j', 'oe'],
         'DESY' :    ['-o', qsub_log, '-j', 'y'],
         'DESY_quick' :  ['-o', qsub_log, '-j', 'y'],
         'CCIN2P3' : ['-o', qsub_log, '-e', qsub_log, '-j', 'yes']}
  return cmd[environ.FARM]
###


def call(cmd,
         enricoDir,
         fermiDir,
         scriptfile=None,
         qsub_log=None,
         jobname=None,
	 submit=True,
	 max_jobs=50,
         #logfile=None,
         check_present=None,
         clobber=False,
         exec_dir=None,
         dry=False,
         options=None):
    """Run a command line tool either directly
    or submit to the queue"""
    if check_present and not clobber:
        if os.path.exists(check_present):
            logging.info('{0} exists. Skipping.'
                         ''.format(check_present))
            return

 #   if logfile:
  #      cmd += '>'+ logfile+ '2>&1'

    if not isinstance(cmd, str):
        cmd = _cmd_to_str(cmd)
    if options:
        cmd += _options_to_str(options)
    logging.info(cmd)

    #Number of Max jobs in the queue
    max_jobs = 50
    if environ.FARM=="LAPP":
        max_jobs = 1000
    elif environ.FARM in ["DESY", "DESY_quick"]:
        max_jobs = 50000
    elif environ.FARM=="LOCAL":
        max_jobs = 100
    elif environ.FARM=="CCIN2P3":
        max_jobs = 3500

    # The following steps are different if you submit or not
    if submit:
        wait_for_slot(max_jobs)

        # Note that qsub needs a shell script which sets
        # up the environment and then executes cmd.
        template = join(dirname(__file__),
                        'qsub_'+environ.FARM+'.sh')
        fh = file(template)
        text = fh.read()
        fh.close()

        # Changes to home dir by default, which happens
        # anyway in a new shell.
        if exec_dir:
            text += '\ncd {0}\n\n'.format(exec_dir)

        text +='conda activate fermi'+'\n'
        #text +='export FERMI_DIR='+fermiDir+'\n'
        #text +='export HEADAS_DIR='+fermiDir+'\n'
        text +='export ENRICO_DIR='+enricoDir+'\n'
        #text +='source $HEADAS_DIR/headas-init.sh\n'
        #text +='source $FERMI_DIR/fermi-init.sh\n'
        text +='source $ENRICO_DIR/enrico-init.sh\n'
        #text +='export PYTHONPATH=/usr/local/lib/python2.7/dist-packages/:$PYTHONPATH\n'
        text +='export LATEXDIR=/tmp/aux\n'
        #text +='env\n'
        text +='#PBS -o '+qsub_log+'\n'

        ''' # Trying to move it to the qsub job template
        if environ.FARM in ['DESY','DESY_quick']:
            text += 'export TMPDIR=/tmp/pfiles/$(date +"%s"); mkdir -p $TMPDIR\n'

        if environ.FARM in ['LAPP','DESY','DESY_quick']:
            text += 'cp -R $FERMI_DIR/syspfiles/* $TMPDIR/.\n'
            text += 'export PFILES=$TMPDIR\n'
            text += 'cd $TMPDIR\n'

        text += cmd

        if environ.FARM in ['LAPP','DESY','DESY_quick']:
            text += '\ncd ..\n'
            text += 'rm -fr $TMPDIR\n'
        '''
        
        if environ.FARM in ['LAPP','DESY','DESY_quick']:
            text += 'if [ "$PFILES" != "" ]; then cp -R $FERMI_DIR/syspfiles/* ${PFILES}/; fi\n'
        
        text += cmd

        # Now reset cmd to be the qsub command
        cmd = GetSubCmd()
        if jobname:
            if environ.FARM in ["CCIN2P3","DESY","DESY_quick"]:
                if jobname[0].isdigit():
                    jobname='_'+jobname
            cmd += ['-N', jobname]

        if scriptfile == None:
            # Note that mkstemp() returns an int,
            # which represents an open file handle,
            # which we have to close explicitly to
            # avoid running out of file handles foer
            # > 100s of jobs.
            (outfd,scriptfile)=tempfile.mkstemp()
            outsock=os.fdopen(outfd,'w')
            outsock.close()
            del outfd

        if qsub_log == None:
            (outfd,qsub_log)=tempfile.mkstemp()
            outsock=os.fdopen(outfd,'w')
            outsock.close()
            del outfd

        cmd += GetSubOutput(qsub_log)
        #if environ.FARM in ["DESY"]:
        #    cmd += ['/usr/bin/singularity','exec','/project/singularity/images/SL6.img','sh']
        cmd += [scriptfile]

        cmd = _cmd_to_str(cmd)
        logging.info(cmd)
    else:
        if exec_dir:
            os.chdir(exec_dir)

        text = cmd

    # Now the following steps are again identical
    # for submitting or not
    if scriptfile:
        logging.debug('Saving command in file: {0}'
                      ''.format(scriptfile))
        fh = file(scriptfile, 'w')
        fh.write(text + '\n')
        fh.close()
        #os.chmod(scriptfile, stat.S_IRWXU)

    if not dry:
        print("Running: %s" %cmd)
        os.system(cmd)
