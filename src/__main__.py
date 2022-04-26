#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import itertools
from check_modules import py_version, pkg_requirements
from funtools import welcome_message, bomb, MyCols

welcome_message()
print(f"{MyCols.OKCYAN}checking system requirements...")
py_version()
pkg_requirements()
print(f"{MyCols.OKGREEN}system requirements met...")


def main():
    args1 = [('-D', '--DIR'), ('-O', '--OUT')]
    args2 = [('-g', '--geno'), ('-s', '--snps'), ('-o', '--out'),
             ('-n', '--samples'), ('-t', '--type'), ('-a', '--alleles')]

    if len(sys.argv) == 1:
        from parse_args import msg
        print(msg())
        del msg
        raise SystemExit

    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        import geno2fi
        geno2fi
        del geno2fi
        raise SystemExit

    elif any(x in sys.argv[1:] for x in list(itertools.chain(*args1))):
        x = set([item for item in args1 for a in sys.argv[1:] if a in item])
        y = set([item for item in args1 for a in sys.argv[1:] if a not in item])
        z = list(x.symmetric_difference(y))
        if z:
            bomb(f'Missing argument when trying to convert to fimpute:\n'
                 f'       required args: {z}\n'
                 f'       try `./snprecode -h` for complete arguments list\n')

        from check_path import check_path
        if not check_path():
            print("File path OK...")
        else:
            print(check_path())
        del check_path

        from garbage import check_dups, file_list
        if not check_dups():
            print("File check complete...\nFiles to be processed...\n")
            print('\n'.join(map(str, file_list)))
        del check_dups

        import geno2fi
        geno2fi
        del geno2fi
        raise SystemExit

    elif any(x in sys.argv[1:] for x in list(itertools.chain(*args2))):
        x = set([item for item in args2 for a in sys.argv[1:] if a in item])
        y = set([item for item in args2 for a in sys.argv[1:] if a not in item])
        z = list(x.symmetric_difference(y))
        if z:
            bomb(f"Missing argument when trying to convert from fimpute\n"
                 f"       required args: {z}\n"
                 f"       try: `./snprecode -h` for complete arguments list\n")

        import fi2geno
        fi2geno
        del fi2geno
        raise SystemExit

    elif [item for item in sys.argv[1:] if item.endswith(('bim', 'map'))]:
        import snpinfo
        snpinfo
        del snpinfo
        raise SystemExit

    elif [item for item in sys.argv[1:] if item.endswith(('vcf', 'vcf.gz'))]:
        import geno_corr
        geno_corr
        del geno_corr
        raise SystemExit

    else:
        bomb('Unknown argument(s)\n       try ./snprecode -h')


if __name__ == "__main__":
    main()
