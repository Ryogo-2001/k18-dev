#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import copy
import datetime
import fcntl
import json
import logging
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

import bjob
import bjobmanager
import bsub
import classimpl
import singlerun
import threading
import utility

logger = logging.getLogger('__main__').getChild(__name__)

HSTAGE_PATH = '/ghi/fs02/hstage/requests'
HSTAGE_NMIN = 20
HSTAGE_NMAX = 10000

RUN_PERIOD = 5 # second

#______________________________________________________________________________
class RunManager(metaclass=classimpl.Singleton):
  ''' Manager class for handling run jobs. '''

  #____________________________________________________________________________
  def __init__(self):
    self.start_time = None
    self.__runlist_manager = None
    self.__stat_path = None
    self.__tag = None
    self.__nruns = None
    self.__stage_list_path = None
    self.__is_ready = False
    self.__is_done = list()
    self.__run_job_list = list()
    self.__stage_list = list()
    self.__updater_thread = threading.Thread(target=self.__updater)
    self.__updater_thread.name = 'updater'
    self.__updater_thread.daemon = True
    self.__updater_status = 'IDLE'

  #____________________________________________________________________________
  def dump_status(self):
    ''' Dump status. '''
    if self.__is_ready is False:
      return
    info = dict()
    for run_job in self.__run_job_list:
      info.update(run_job.get_info())
    buff = json.dumps(info, indent=4)
    # with open(self.__stat_path, 'w+') as f:
    with open(self.__stat_path, 'r+') as f:
      try:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
      except IOError as e:
        logger.error(e)
      else:
        f.write(buff)
        f.truncate()
        f.flush()
      finally:
        fcntl.flock(f.fileno(), fcntl.LOCK_UN)

  #____________________________________________________________________________
  def finalize(self, flag=False):
    ''' Finalize run jobs. '''
    if self.__is_ready is False:
      return
    for run_job in self.__run_job_list:
      run_job.finalize(flag)
    self.__updater_status = 'END'
    while self.__updater_thread.is_alive():
      self.__updater_thread.join()

  #____________________________________________________________________________
  def initialize(self):
    ''' Initialize. '''
    if (self.__runlist_manager is None or
        self.__stat_path is None):
      return
    self.start_time = time.time()
    self.__stage_list_path = os.path.join(
      HSTAGE_PATH,
      self.__tag+'.lst.'+datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    for i in range(self.__nruns):
      run = self.__runlist_manager.get_run_info(i)
      run_job = singlerun.SingleRun(self.__tag, run)
      self.__run_job_list.append(run_job)
    self.__stage_list = list()
    for run_job in self.__run_job_list:
      self.__is_done.append(0)
      if not run_job.is_staged():
        self.__stage_list.append(run_job.get_data_path())
    self.__is_ready = True

  #____________________________________________________________________________
  def kill(self):
    ''' Kill all jobs. '''
    if self.__is_ready is False:
      return
    for run_job in self.__run_job_list:
      run_job.update_status()
      status = run_job.get_status()
      if status != 'DONE' or status != 'FAILED':
        run_job.kill_all()

  #____________________________________________________________________________
  def print_info(self):
    logger.info(f'{self.start_time}')
    logger.info(f'{self.__runlist_manager}')
    logger.info(f'{self.__stat_path}')
    logger.info(f'{self.__tag}')
    logger.info(f'{self.__nruns}')
    logger.info(f'{self.__stage_list_path}')
    logger.info(f'{self.__is_ready}')
    logger.info(f'{self.__is_done}')
    logger.info(f'{self.__run_job_list}')
    logger.info(f'{self.__stage_list}')

  #____________________________________________________________________________
  def register_staging(self):
    ''' Register staging. '''
    if self.__is_ready is False:
      return
    self.__stage_list_path = os.path.join(
      HSTAGE_PATH,
      self.__tag+'.lst.'+datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    if len(self.__stage_list) > HSTAGE_NMAX:
      self.__is_ready = False
      logger.error('Number of runs in runlist is too much.')
      exit(1)
    if len(self.__stage_list) > HSTAGE_NMIN:
      logger.info(f'Put {self.__stage_list_path}')
      with open(self.__stage_list_path, 'w') as f:
        for item in self.__stage_list:
          f.write(item + '\n')

  #____________________________________________________________________________
  def run(self):
    ''' Run. '''
    if self.__is_ready is False:
      return
    ptime = time.time()
    self.__updater_status = 'RUNNING'
    self.__updater_thread.start()
    bjob_man = bjobmanager.BJobManager()
    bjob_man.start()
    while sum(self.__is_done) != self.__nruns:
      logger.debug('running')
      self.__single_cycle()
      dtime = RUN_PERIOD - (time.time() - ptime)
      if dtime > 0:
        time.sleep(dtime)
      ptime = time.time()
    bjob_man.stop()
    self.__updater_status = 'END'
    while self.__updater_thread.is_alive():
      self.__updater_thread.join()

  #____________________________________________________________________________
  def set_runlist_manager(self, runlistman):
    ''' Set runlist manager. '''
    self.__runlist_manager = runlistman
    self.__tag = self.__runlist_manager.get_tag()
    self.__nruns = self.__runlist_manager.get_nruns()

  #____________________________________________________________________________
  def set_status_output_path(self, fstat_path):
    ''' Set status output path. '''
    self.__stat_path = fstat_path

  #____________________________________________________________________________
  def __single_cycle(self):
    ''' Single cycle for run. '''
    for i, run_job in enumerate(self.__run_job_list):
      status = run_job.get_status()
      if status == 'INIT':
        pass
      elif status == 'PURGED':
        if len(self.__stage_list) <= HSTAGE_NMIN:
          run_job.access_data_stream()
      elif status == 'STAGED':
        run_job.execute()
      elif status == 'BJOB-RUNNING':
        pass
      elif status == 'BJOB-DONE':
        run_job.merge_root_files()
      elif status == 'MERGING':
        pass
      elif status == 'TERMINATED':
        pass
      elif status == 'DONE':
        if self.__is_done[i] == 0:
          run_job.finalize()
          self.__is_done[i] = 1
      elif status == 'FAILED':
        if self.__is_done[i] == 0:
          run_job.finalize()
          self.__is_done[i] = 1
      else:
        logger.error('Detect unknown status.')

  #____________________________________________________________________________
  def __updater(self):
    logger.debug('updater thread starts')
    while self.__updater_status == 'RUNNING':
      logger.debug(f'updater is running')
      for run in self.__run_job_list:
        run.update_status()
        self.dump_status()
      time.sleep(1)
    logger.debug('updater thread ends')
