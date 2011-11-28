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
from os.path import join, dirname
import commands
import time
import datetime
import tempfile

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
    user = os.environ['USER']
    #cmd = 'qstat -u {0} | wc -l'.format(user)
    #njobs = int(commands.getoutput(cmd))

    # Dirty hack, because commands.getoutput() seems to
    # open more and more file handles, causing a crash at some point.
    fn = '/tmp/_temp_qsub.txt'
    os.system('qstat -u {user} | wc -l > {fn}'
              ''.format(user=user, fn=fn))
    njobs = int(file(fn).read())
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

def call(cmd,
         scriptfile=None,
         logfile=None,
         check_present=None,
         clobber=False,
         exec_dir=None,
         dry=False,
         options=None,
         submit=False,
         jobname=None,
         max_jobs=100,
         qsub_log=None):
    """Run a command line tool either directly
    or submit to the queue"""
    if check_present and not clobber:
        if os.path.exists(check_present):
            logging.info('{0} exists. Skipping.'
                         ''.format(check_present))
            return

    if logfile:
        cmd += ['>', logfile, '2>&1']

    if not isinstance(cmd, str):
        cmd = _cmd_to_str(cmd)
    if options:
        cmd += _options_to_str(options)
    logging.info(cmd)

    # The following steps are different if you submit or not
    if submit:
        wait_for_slot(max_jobs)

        # Note that qsub needs a shell script which sets
        # up the environment and then executes cmd.
        template = join(dirname(__file__), 
                        'qsub_template.sh')
        fh = file(template)
        text = fh.read()
        fh.close()

        # Changes to home dir by default, which happens
        # anyway in a new shell.
        if exec_dir:
            text += '\ncd {0}\n\n'.format(exec_dir)

        text += cmd

        # Now reset cmd to be the qsub command
        cmd = ['qsub']
        if jobname:
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
            
        cmd += ['-o', qsub_log, '-j', 'y']
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

    if not dry:
        os.system(cmd)

def ask_permission_to_proceed(steps, available_steps):
    """ For windows users ... :-) """
    print
    print '---> Available steps:', available_steps
    print '---> Selected  steps:', steps

    for step in steps:
        if not step in available_steps:
            print('---> USAGE ERROR: "{0}" is not an available step.'
                  ' Exiting.'.format(step))
            sys.exit(-1)

    print '---> Hit Enter to proceed or Ctrl+C to exit.'
    raw_input()

def root(options):
    cmd = ['root -l -q']
    return call(cmd, options=options)
