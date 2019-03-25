#!/usr/bin/env python

import os

#./mmpbsa_decomp.py rec/ ligs.split/5cm8A-I466V-000001.mol2 out/5cm8A-I466V-000001.csv

def chunk(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

for x in chunk(os.listdir('ligs.split'), 16):

    s = '''#PBS -S /bin/bash
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime=6:00:00

CWD="/N/dc2/scratch/yx5/vs.decomp"

cd $CWD

source /N/u/yx5/Karst/amber16/amber.sh
'''

    for i in x:
        s += './mmpbsa_decomp.py rec ligs.split/{0}.mol2 out/{0}.csv &\n'.format(i[:-5])
    s += 'wait\n'

    with open('mmpbsa.sh', 'w') as fout:
        fout.write(s)

    cmd = 'qsub mmpbsa.sh'
    os.system(cmd)
