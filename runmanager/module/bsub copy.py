#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import logging
import psutil
import resource
import shlex
import subprocess
import time

import bjob
import runmanager

logger = logging.getLogger('__main__').getChild(__name__)
rsrc = resource.RLIMIT_NPROC
soft, hard = resource.getrlimit(rsrc)
resource.setrlimit(rsrc, (hard, hard))
MAX_NPROC = 100 #hard
rsrc = resource.RLIMIT_NOFILE
soft, hard = resource.getrlimit(rsrc)
resource.setrlimit(rsrc, (hard, hard))
MAX_NOFILE = hard - 10

#______________________________________________________________________________
class BSub(object):
  ''' BSub class throwing subprocess of "bsub". '''

  #____________________________________________________________________________
  def __init__(self, run, tag, conf, out, log):
    self.__main_process = psutil.Process()
    self.__run = run
    self.__tag = tag
    self.__conf = conf
    self.__out  = out
    self.__log  = log
    self.__pid  = None
    self.__proc = None
    self.__jid  = None
    self.__bjob = None
    self.__stime = None
    self.__rtime = None
    # True: success, False: failure,
    # 0: process thrown, 1: process return and bsub running, 2: killed, -1: unknown
    self.__status   = 'INIT'
    self.__process_status = 'INIT' # UNIX process.
    self.__bjob_status = None # bsub process.
    self.__is_dst = False
    self.__dstin_list = None

  #____________________________________________________________________________
  def check_limits(self):
    ''' Check user limits of fds/proc. '''
    nfds = self.__main_process.num_fds()
    nproc = len(self.__main_process.children(recursive=True))
    logger.debug(f'{self.__out} nfds={nfds}, nproc={nproc}')
    return (nfds < MAX_NOFILE and nproc < MAX_NPROC)

  #____________________________________________________________________________
  def execute(self):
    ''' Execute bsub job. '''
    self.__update_status()
    if self.__process_status != 'INIT':
      return
    while not self.check_limits():
      logger.debug(f'Releasing fds/proc ...')
      time.sleep(1)
    # pf = self.__run.getPreFetchPath()
    cmd = shlex.split('bsub' + ' '
                      + '-q' + ' ' + self.__run.get_queue() + ' '
                      + self.__run.get_option() + ' '
                      + '-o' + ' ' + self.__log + ' '
                      # + '-a \"prefetch (' + pf + ')\"'
                      )
    if self.__is_dst:
      cmd.extend([self.__run.get_bin_path(),
                  self.__conf])
      cmd.extend(self.__dstin_list)
      cmd.append(self.__out)
    else:
      cmd.extend([self.__run.get_bin_path(),
                  self.__conf,
                  self.__run.get_data_path(),
                  self.__out])
    self.__proc = subprocess.Popen(cmd,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    self.__pid = self.__proc.pid
    logger.debug(f'submit <pid={self.__pid}>')
    self.__process_status = 'RUNNING'
    # self.__status   = 0

  #____________________________________________________________________________
  def get_job_id(self):
    ''' Get job id. '''
    self.__update_process_status()
    return self.__jid

  #____________________________________________________________________________
  def get_process_id(self):
    ''' Get process id. '''
    return self.__pid

  #____________________________________________________________________________
  def get_run_time(self):
    ''' Get elapsed running time. '''
    self.__update_job_status()
    if self.__rtime is None:
      return 0
    else:
      return self.__rtime

  #____________________________________________________________________________
  def get_status(self):
    ''' Get status. '''
    if (self.__status is True or
        self.__status is False or
        self.__status is 2):
      return self.__status
    self.__update_status()
    return self.__status

  #____________________________________________________________________________
  def get_tag(self):
    ''' Get tag. '''
    return self.__tag

  #____________________________________________________________________________
  def kill(self):
    ''' Kill job. '''
    self.__update_process_status()
    if (self.__process_status == 'RUNNING' or
        self.__process_status == 'UNKNOWN'):
      logger.info(f'Killing process [pid: {self.__pid}]')
      self.__proc.kill()
      self.__process_status = 'TERMINATED'
    self.__update_job_status()
    if (self.__bjob_status == 0 or self.__bjob_status == -1):
      logger.info(f'Killing job [jid: {self.__pid}]')
      self.__bjob.kill()
      self.__bjob_status = 1

  #____________________________________________________________________________
  def set_dstin_list(self, dstin_list):
    ''' Set list of input files for dst analysis '''
    self.__is_dst = True
    self.__dstin_list = dstin_list

  #____________________________________________________________________________
  def __register_job(self):
    ''' Register job. '''
    self.__update_status()
    if self.__jid is None and self.__process_status == 'DONE':
      buff  = self.__proc.communicate()[0]
      self.__jid = bjob.BJob.read_job_id(buff.decode())
      self.__bjob = bjob.BJob(self.__jid)
      self.__bjob_status = 0
      # self.__status   = 1

  #____________________________________________________________________________
  def __update_job_status(self):
    ''' Update job status. '''
    if self.__bjob_status != 0:
      return
    status = self.__bjob.get_status()
    if status == 'INIT':
      pass
    elif status == 'PEND':      # PEND
      pass
    elif status == 'RUN':      # RUN
      self.__rtime = time.time() - self.__stime
    elif status == 'DONE':      # DONE
      self.__bjob_status = True
      # self.__status   = True
    elif status == 'EXIT':      # EXIT
      self.__bjob_status = False
      # self.__status   = False
    elif status == 'TERMINATED':      # killed
      pass
    else:                  # unknow case
      self.__bjob_status = -1
      # self.__status   = -1

  #____________________________________________________________________________
  def __update_process_status(self):
    ''' Update process status. '''
    if self.__process_status != 'RUNNING':
      return
    if self.__proc.poll() is None:
      return
    elif self.__proc.poll() == 0:
      self.__process_status = 'DONE'
      if self.__bjob is None:
        self.__register_job()
        self.__stime = time.time()
        self.__rtime = self.__stime
    elif self.__proc.poll() == 1:
      logger.error(f'bsub command failed at {self.__tag}')
      outs, errs = proc.communicate()
      self.__process_status = 'FAILED'
    else:
      self.__process_status = 'UNKNOWN'

  #____________________________________________________________________________
  def __update_status(self):
    ''' Update status. '''
    if (self.__status is True or
        self.__status is False or
        self.__status is 2):
      return
    self.__update_process_status()
    self.__update_job_status()
    # proc
    if self.__process_status == 'INIT':
      pass
    elif self.__process_status == 'DONE':
      self.__status = 0
    elif self.__process_status == 'FAILED':
      self.__status = False
    elif self.__process_status == 'TERMINATED':
      self.__status = 2
    else:
      self.__status = -1
    # bjob
    if self.__bjob_status == 'INIT':
      pass
    elif self.__bjob_status is True:
      self.__status = True
    elif self.__bjob_status is False:
      self.__status = False
    elif self.__bjob_status is 0:
      self.__status = 1
    elif self.__bjob_status is 1:
      self.__status = 2
    else:
      self.__status = -1
