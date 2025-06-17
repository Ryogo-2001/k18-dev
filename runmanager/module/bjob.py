#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import bjobmanager
import logging
import shlex
import subprocess
import time

BSUB_RESPONCE = 'Job <JobID> is submitted to queue <queue>'
logger = logging.getLogger('__main__').getChild(__name__)

#______________________________________________________________________________
class BJob(object):
  ''' BJob class throwing subprocess of "bjobs" and "bkill". '''

  #____________________________________________________________________________
  def __init__(self, job_id):
    self.__job_id = job_id
    self.__status = 'INIT'
    self.get_status()

  #______________________________________________________________________________
  def get_job_id(self):
    ''' Get job id. '''
    return self.__job_id

  #______________________________________________________________________________
  def get_status(self):
    ''' Get job status. '''
    if (self.__status == 'DONE' or
        self.__status == 'EXIT' or
        self.__status == 'TERMINATED'):
      return self.__status
    status = bjobmanager.BJobManager().get_job_status(self.__job_id)
    logger.debug(f'{self.__job_id} is {status}')
    # cmd = f'bjobs {self.__job_id}'
    # proc = None
    # try:
    #   proc = subprocess.run(shlex.split(cmd),
    #                         stdout=subprocess.PIPE,
    #                         stderr=subprocess.PIPE,
    #                         check=True)
    # except subprocess.CalledProcessError as e:
    #   logger.error(e)
    #   self.__status = 'FAILED'
    #   return self.__status
    # buff = proc.stdout.splitlines()[1].decode().split()
    # status = 'UNKNOWN'
    # if int(buff[0]) == self.__job_id:
    #   if buff[2] == 'PEND':
    #     status = 'PEND'
    #   elif buff[2] == 'RUN':
    #     status = 'RUN'
    #   elif buff[2] == 'DONE':
    #     status = 'DONE'
    #   elif buff[2] == 'EXIT':
    #     status = 'EXIT'
    self.__status = status
    return status

  #____________________________________________________________________________
  def kill(self):
    ''' Kill job. '''
    self.get_status()
    if (self.__status == 'PEND' or
        self.__status == 'RUN'):
      cmd = f'bkill {self.__job_id}'
      proc = subprocess.run(shlex.split(cmd),
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
      self.__status = 'TERMINATED'

  #____________________________________________________________________________
  @staticmethod
  def read_job_id(buff):
    ''' Read job id. '''
    job_id = None
    flag = True
    words = buff.split()
    for i, item in enumerate(BSUB_RESPONCE.split()):
        if i == 1 or i == 6:
          continue
        if item != words[i]:
          flag = False
    if flag is True:
      job_id = int(words[1][1:-1])
    return job_id
