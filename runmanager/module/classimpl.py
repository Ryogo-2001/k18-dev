#!/usr/bin/env python3

__author__ = 'Y.Nakada <nakada@ne.phys.sci.osaka-u.ac.jp>'
__version__ = '4.1'
__date__ = '16 Feb. 2021'

#______________________________________________________________________________
class Singleton(type):
  __instances = {}
  def __call__(cls, *args, **kwargs):
    if cls not in cls.__instances:
      cls.__instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
    return cls.__instances[cls]
