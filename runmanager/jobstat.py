#!/usr/bin/env python3

import getpass
import logging
import logging.config
import os
import subprocess
import sys
import time
import yaml

top_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(top_dir, 'module'))

from utility import pycolor as cl

logger = logging.getLogger(__name__)
#log_format = 'Last update: %(asctime)s [%(levelname)s] %(message)s'
log_format = 'Last update: %(asctime)s %(message)s'
logging.basicConfig(level=logging.INFO,
                    format=log_format)

interval = 10 # [sec]
target_jobs = ['s', 'l', 'h', 'sx', 'lx', 'hx', 'a']
n_jobs = dict()

#______________________________________________________________________________
def bjobs():
  ''' bjobs '''
  ret = str()
  # ret += ('=== bjobs ' + '='*80)[:80] + '\n'
  buf = process('bjobs')
  global n_jobs
  n_jobs = dict()
  buf_keys = None
  for line in buf.splitlines():
    buf_info = line.split()
    if buf_keys is None:
      buf_keys = buf_info
      continue
    buf_id = buf_info[buf_keys.index('JOBID')]
    user = buf_info[buf_keys.index('USER')]
    stat = buf_info[buf_keys.index('STAT')]
    queue = buf_info[buf_keys.index('QUEUE')]
    from_host = buf_info[buf_keys.index('FROM_HOST')]
    exec_host = buf_info[buf_keys.index('EXEC_HOST')]
    job_name = buf_info[buf_keys.index('JOB_NAME')]
    submit_time = buf_info[buf_keys.index('SUBMIT_TIME')]
    if queue not in n_jobs:
      n_jobs[queue] = dict()
    if stat not in n_jobs[queue]:
      n_jobs[queue][stat] = 0
    n_jobs[queue][stat] += 1
  ret += ('NJOBS = ' + str(n_jobs) + '\n' +
          cl.end)
  return ret

#______________________________________________________________________________
def bqueues():
  ''' bqueues '''
  ret = str()
  # ret += ('=== bqueues '+'='*80)[:80] + '\n'
  buf = process('bqueues')
  njobs = dict()
  buf_keys = None
  for line in buf.splitlines():
    buf_info = line.split()
    if buf_keys is None:
      buf_keys = buf_info
      key_line = line
      continue
    queue = buf_info[buf_keys.index('QUEUE_NAME')]
    prio = buf_info[buf_keys.index('PRIO')]
    status = buf_info[buf_keys.index('STATUS')]
    maxjob = buf_info[buf_keys.index('MAX')]
    jl_u = buf_info[buf_keys.index('JL/U')]
    jl_p = buf_info[buf_keys.index('JL/P')]
    jl_h = buf_info[buf_keys.index('JL/H')]
    n = buf_info[buf_keys.index('NJOBS')]
    pend = buf_info[buf_keys.index('PEND')]
    run = buf_info[buf_keys.index('RUN')]
    susp = buf_info[buf_keys.index('SUSP')]
    if queue in target_jobs and queue not in njobs:
      njobs[queue] = line
    # if pend not in njobs[queue]:
    #   njobs[queue][stat] = 0
    # njobs[queue][stat] += 1
  ret += (cl.blue + key_line + '\n' + cl.end)
  for q in target_jobs:
    if q in n_jobs:
      ret += njobs[q] + '\n'
  buf = process('bqueues -l s')
  if buf is None:
    return ret
  buf_keys = None
  myline = None
  myprio = None
  nusers = 0
  for line in buf.splitlines():
    buf_info = line.split()
    if len(buf_info) == 0:
      continue
    if buf_info[0] == 'USER/GROUP':
      buf_keys = buf_info
      key_line = line
      continue
    if buf_keys is None:
      continue
    if len(buf_keys) == len(buf_info):
      nusers += 1
    if buf_info[buf_keys.index('USER/GROUP')] == getpass.getuser():
      myline = line
      myprio = nusers
  ret += (cl.blue + key_line[:80] + '\n' + cl.end)
  if myline is not None:
    prio = f'({myprio}/{nusers})'
    ret += myline[:80].replace(getpass.getuser() + len(prio)*' ',
                               (getpass.getuser() + f'({myprio}/{nusers})'))
    ret += '\n'
  return ret

#______________________________________________________________________________
def chk_runjob():
  ''' chk_runjob '''
  ret = str()
  # ret += ('=== chk_runjob '+'='*80)[:80] + '\n'
  buf = process('chk_runjob')
  buf_keys = None
  for line in buf.splitlines():
    buf_info = line.split()
    if buf_keys is None:
      buf_keys = buf_info
      # logger.info(buf_keys)
      continue
    jobs = int(buf_info[buf_keys.index('JOBS')])
    slot = int(buf_info[buf_keys.index('SLOT')])
    cputime = int(buf_info[buf_keys.index('CPUTIME')])
    runtime = int(buf_info[buf_keys.index('RUNTIME')])
    cpurun = float(buf_info[buf_keys.index('CPU/RUN')])
  ret += buf
  return ret

#______________________________________________________________________________
def monitor():
  ''' monitor '''
  while True:
    buf = '\n'
    buf += bjobs()
    buf += bqueues()
    buf += chk_runjob()
    while buf[-1] == '\n':
      buf = buf[:-1]
    os.system('clear')
    logger.info(buf)
    time.sleep(3)

#______________________________________________________________________________
def process(command):
  ''' process '''
  try:
    proc = subprocess.run(command,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          shell=True, check=True)
  except subprocess.CalledProcessError as e:
    logger.error(f'command "{e.cmd}" returned error code ({e.returncode})')
    return None
  return proc.stdout.decode()

#______________________________________________________________________________
if __name__ == '__main__':
  try:
    monitor()
  except KeyboardInterrupt as e:
    logger.info('Terminated')
