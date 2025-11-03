#!/usr/bin/env python3

#______________________________________________________________________________
# RAYRAW用に改訂akao/2025/6/23
__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import logging
import shlex
import subprocess
import time
import bjob
import runmanager
import resource
import psutil
logger = logging.getLogger('__main__').getChild(__name__)

#______________________________________________________________________________
class BSub(object):
  ''' Wrapper class of bsub. '''

  #____________________________________________________________________________
  def __init__(self, parent, tag, conf, root, log):
    self.__parent = parent
    self.__tag = tag
    self.__conf = conf
    self.__root = root
    self.__log = log
    self.__proc = None
    self.__job_id = None
    self.__status = 0 # 0: running, 1: submitted, true: done, false: fail

  #____________________________________________________________________________
  def execute(self):
    ''' Execute bsub. '''
    if self.__proc is not None:
      return
    bin_path = self.__parent.get_bin_path()
    data_path = self.__parent.get_data_path()
    queue = self.__parent.get_queue()
    option = self.__parent.get_option()

    # ======================================================================== #
    # vv これが修正されたコマンド組み立て部分です vv
    #
    # Rayraw の正しい引数の順番（conf, data, root）に修正し、
    # 不要な -c と -o を削除しました。
    #
    run_cmd = (f'{bin_path} '
               f'{self.__conf} '
               f'{data_path} '
               f'{self.__root} '
               f'{option}')
    # ======================================================================== #
    
    cmd = shlex.split(f'bsub -q {queue} -o {self.__log} "{run_cmd}"')
    
    try:
      self.__proc = subprocess.Popen(cmd,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE
                                     )
    except subprocess.CalledProcessError as e:
      logger.error(f'command "{e.cmd}" returned error code ({e.returncode})')
      logger.error(proc.stderr)
      self.__status = False
    
  #____________________________________________________________________________
  def get_job_id(self):
    ''' Get job id. '''
    return self.__job_id

  #____________________________________________________________________________
  def get_process_id(self):
    ''' Get job id. '''
    if self.__proc is not None:
      return self.__proc.pid
    else:
      return None

  #____________________________________________________________________________
  def get_status(self):
    ''' Get status. '''
    if self.__status is True or self.__status is False:
      return self.__status
    if self.__proc is not None:
      if self.__proc.poll() is not None:
        buff = self.__proc.stdout.read().decode()
        self.__job_id = self.read_job_id(buff)
        if self.__job_id is None:
          self.__status = False
        else:
          self.__status = 1
    return self.__status

  #____________________________________________________________________________
  def get_tag(self):
    ''' Get tag. '''
    return self.__tag

  #____________________________________________________________________________
  def kill(self):
    ''' Kill job. '''
    if self.__job_id is None:
      return
    subprocess.run(['bkill', str(self.__job_id)])
    
  #____________________________________________________________________________
  @staticmethod
  def read_job_id(buff):
    ''' Read job id from bsub output. '''
    words = buff.split()
    cand = [word for word in words if word.isdigit()]
    return int(cand[0]) if len(cand) == 1 else None