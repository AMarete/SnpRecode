#!/usr/bin/env python
# -*- coding: utf-8 -*-

from platform import python_version
import subprocess
import sys


def py_version():
    try:
        float(python_version()[0:3]) >= 3.5
    except SystemError:
        print('''Python version not satisfied, install Python V3.5 or later''')
        raise SystemExit


'''
def pkg_requirements():
    required = {'biopython'}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        print(f"install package(s) {list(missing)} with: pip3 install {' '.join(list(missing))}")
        raise SystemExit
'''


def pkg_requirements():
    reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
    installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
    required = 'biopython'
    if required not in installed_packages:
        print(f"install package(s) {required} with: pip3 install {''.join(list(required))}")
        raise SystemExit
