#!/usr/bin/env python

from setuptools import setup

setup(name='sgoal',
      version='1.0',
      # list folders, not files
      packages=['sgoal',
                'capitalize.test'],
      scripts=['capitalize/bin/cap_script.py'],
      package_data={'capitalize': ['data/cap_data.txt']},
      )