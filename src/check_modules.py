#!/usr/bin/env python
# -*- coding: utf-8 -*-

from platform import python_version
from funtools import bomb
import subprocess
import sys


def py_version():
    vers_ = python_version().split('.')
    if int(vers_[0]) > 2 and int(vers_[1]) < 6:
        return bomb("Python version not satisfied, install Python V3.6 or later\n")


def pkg_requirements():
    reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
    installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
    required = ['biopython', 'matplotlib']
    if not set(required).issubset(installed_packages):
        bomb(f"install package(s) {required} with: pip install {' '.join(list(required))}\n")
