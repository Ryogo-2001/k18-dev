#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
import argparse
import json
import logging
import logging.config
import os
import sys
import time
import signal
import yaml

top_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(top_dir, 'module'))

import utility
import runlist
import runmanager

logger = logging.getLogger(__name__)

#______________________________________________________________________________
def main(runlist_path, status_path):
  ''' Main function. '''
  runlist_manager = runlist.RunlistManager()
  runlist_manager.set_run_list(runlist_path)
  run_manager = runmanager.RunManager()
  run_manager.set_runlist_manager(runlist_manager)
  run_manager.set_status_output_path(status_path)
  run_manager.initialize()
  signal.signal(signal.SIGINT, signal_handler)
  logger.info("Press 'Ctrl-C' to terminate processes.")
  run_manager.register_staging()
  run_manager.run()
  logger.info('done')

#______________________________________________________________________________
def signal_handler(num, frame):
  ''' Signal handler for Ctrl-C. '''
  logger.info('KeyboardInterrupt')
  logger.info('Terminating processes...')
  run_manager = runmanager.RunManager()
  run_manager.kill()
  run_manager.dump_status()
  time.sleep(1) # waiting until bsub log files are generated
  keep_log = False
  ret = input('Keep log files? [y/-] >> ')
  if ret == 'y':
    logger.info('Deleting intermediate files...')
    keep_log = True
  else:
    logger.info('Deleting files...')
  run_manager.finalize(keep_log)
  exit(1)

#______________________________________________________________________________
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('run_list', help='run list YAML')
  parsed, unparsed = parser.parse_known_args()
  log_conf = os.path.join(top_dir, 'logging_config.yml')
  with open(log_conf, 'r') as f:
    logging.config.dictConfig(yaml.safe_load(f))
  if not os.path.isfile(parsed.run_list):
    logger.error(f'No such file: {parsed.run_list}')
    exit(1)
  stat_dir = os.path.join(top_dir, 'stat')
  if not os.path.exists(stat_dir):
    os.mkdir(stat_dir)
  stat_file = os.path.join(
    stat_dir,
    os.path.splitext(os.path.basename(parsed.run_list))[0]+'.json')
  if not os.path.exists(stat_file):
    with open(stat_file, 'w') as f:
      pass
  main(parsed.run_list, stat_file)
