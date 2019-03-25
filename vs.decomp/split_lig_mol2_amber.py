#!/usr/bin/env python

'''
./script.py prefix path/to/multimol2/file
'''

import os
import sys

def splitmol2(multiMOL2F, delimiter='@<TRIPOS>MOLECULE'):
    #Generator for splitting multiMOL2 files
    mol2Lines = []
    for mol2Line in multiMOL2F:
        if mol2Line.startswith('###'):
            continue
        if mol2Line.startswith(delimiter) and mol2Lines:
            yield mol2Lines
            mol2Lines = []
        mol2Lines.append(mol2Line)
    if mol2Lines:
        yield mol2Lines

def main(name, ligfile, outdir='lig', size=250):
    ligfile = os.path.abspath(ligfile)
    l = []
    count = 0
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(ligfile) as f:
        for mol2 in splitmol2(f):
            l.append(''.join(mol2))
            if len(l) >= size:
                with open('{}/{}-{:0>4d}.mol2'.format(outdir, name, count), 'w') as fout:
                    fout.write(''.join(l))
                l = []
                count += 1
        else:
            with open('{}/{}-{:0>4d}.mol2'.format(outdir, name, count), 'w') as fout:
                fout.write(''.join(l))
    return

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print __doc__
        raise SystemExit
    main(sys.argv[1], sys.argv[2])
