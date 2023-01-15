#!/usr/bin/env python

from setuptools import setup

setup(name='sgoal',
      version='1.0',
      # list folders, not files
      packages=['sgoal'],
      scripts=['sgoal/bin/sgoal_script.py'],
      )