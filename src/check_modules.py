#!/usr/bin/env python
# -*- coding: utf-8 -*-

from platform import python_version
from funtools import bomb
import subprocess
import sys


def py_version():
    if not float(python_version()[0:3]) >= 3.7:
        float(python_version()[0:3]) >= 3.7
        return bomb("Python version not satisfied, install Python V3.7 or later\n")


def pkg_requirements():
    reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
    installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
    required = ['biopython', 'matplotlib']
    if not set(required).issubset(installed_packages):
        bomb(f"install package(s) {required} with: pip3 install {' '.join(list(required))}\n")
