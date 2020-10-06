#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from os.path import isdir, isfile
from check_modules import py_version, pkg_requirements
from funtools import welcome_message, bomb

'''
pkg_requirements()
py_version()
welcome_message()
args0 = ['-h', '--help']
args1 = ['-F', '-O']
args2 = ['-g', '-s', '-o', '-n', '-t', '-a']
'''


def main():
    pkg_requirements()
    py_version()
    welcome_message()
    args1 = ['-F', '-O']
    args2 = ['-g', '-s', '-o', '-n', '-t', '-a']

    if len(sys.argv) == 1:
        from parse_args import msg
        print("\nrun `./snprecode.py -h` for full argument list\nbasic usage: ", msg())
        del msg
        raise SystemExit

    elif sys.argv[1] == '-h' or sys.argv[1] == '-help':
        import geno2fi
        geno2fi
        del geno2fi
        raise SystemExit

    elif any(x in sys.argv[1:] for x in args1):
        gg = set([i for i in sys.argv[1:] if i.startswith('-')])
        if set(args1) - gg:
            bomb(f'Missing argument {list(set(args1) - gg)}\n'
                 'Error: run `./snprecode -h` for complete arguments list\n')
        if sys.argv[1] == "-F" or sys.argv[1] == "--FILEPATH":
            if not isdir(sys.argv[1]):
                bomb(f"Check if directory path (-F, --FILEPATH) exists and non-empty files in dir path\n"
                     f"       Don't forget -O prefix or --OUT prefix\n")

        from garbage import check_dups
        if not check_dups():
            print("File check complete...\n")
        del check_dups

        import geno2fi
        geno2fi
        del geno2fi
        raise SystemExit

    elif any(x in sys.argv[1:] for x in args2):
        gg = set([i for i in sys.argv[1:] if i.startswith('-')])
        if set(args2) - gg:
            bomb(f"Missing argument {list(set(args2) - gg)}\n"
                 "Error: run `./snprecode -h` for complete arguments list\n")

        #if [file_ for file_ in sys.argv[1:] if not isfile(file_) and '-o' not in sys.argv[1:]]:
        if [file_ for file_ in sys.argv[1:] if not isfile(file_) and all(x not in sys.argv[1:] for x in ['-o', '--out'])]:
            bomb(f"Check if all input files exists in file path and are non-empty\n")

        import fi2geno
        fi2geno
        del fi2geno
        raise SystemExit

    else:
        pass


if __name__ == "__main__":
    main()
