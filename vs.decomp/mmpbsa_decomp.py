#!/usr/bin/env python

'''
VS Decomposition for Karst amber16
./script.py /path/to/recs/dir/ /path/to/single/multimol2/file [/path/to/output/file]
'''

import os
import sys
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile, mkdtemp
from shutil import rmtree
import numpy as np

sys.path.append(os.path.join(os.getenv('AMBERHOME'), 'bin'))
from MMPBSA_mods import API as MMPBSA_API

#Global working directory
global cDir
cDir = os.getcwd()

def read_mol2(multiMOL2F, delimiter='@<TRIPOS>MOLECULE'):
    '''
    Generator object for splitting multiMOL2 files
    '''
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

def run_lig_prep(mol2Name, mol2Lines):
    '''
    Calculate net charge and run Antechamber
    '''

    #Calculate net charge
    seek1 = mol2Lines.index('@<TRIPOS>ATOM\n')
    seek2 = mol2Lines.index('@<TRIPOS>BOND\n')
    nc = 0.0
    for x in mol2Lines[seek1 + 1 : seek2]:
        nc += float(x.strip().split()[8])
    ncr = int(round(nc))

    #Write the single mol2 file
    with NamedTemporaryFile(suffix = '.mol2',
                            delete = False) as ligtmpFile:
        ligtmpFile.write(''.join(mol2Lines))

    #antechamber
    args = [
        'antechamber',
        '-i',
        ligtmpFile.name,
        '-o',
        os.path.join('lig_amber.mol2'),
        '-fi',
        'mol2',
        '-fo',
        'mol2',
        '-at',
        'gaff',
        '-c',
        'gas', #Gasteiger charge instead of AM1-BCC
        #'bcc'
        '-nc',
        str(ncr),
        ]
    p = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    
    if not os.path.exists('lig_amber.mol2'):
        print out, err
        return 1, [], []
    
    #parmchk
    args = [
        'parmchk',
        '-i',
        'lig_amber.mol2',
        '-o',
        'frcmod',
        '-f',
        'mol2'
        ]
    p = Popen(args)
    out, err = p.communicate()

    if not os.path.exists('frcmod'):
        print out, err
        return 1, [], []

def run_tleap(recPath):
    '''
    Run tleap to generate the topology files and single frame instance
    '''

    tleapLines = '''source leaprc.protein.ff14SB
source leaprc.gaff2
ligforce = loadamberparams frcmod
rec = loadpdb {0}
lig = loadmol2 lig_amber.mol2
com = combine {{rec lig}}
set default PBradii mbondi2
saveamberparm com com.top com.crd
saveamberparm rec rec.top rec.crd
saveamberparm lig lig.top lig.crd
quit
'''.format(recPath)
    
    with open('tleap.in', 'w') as tleapFile:
        tleapFile.write(''.join(tleapLines))

    args = [ 'tleap', '-f', 'tleap.in' ]
    p = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    outFilesPrefix = ( 'com', 'lig', 'rec' )
    outFilesSuffix = ( 'top', 'crd' )
    for prefix in outFilesPrefix:
        for suffix in outFilesSuffix:
            outFile = '.'.join((prefix, suffix))
            if not os.stat(outFile).st_size:
                print out, err
                return 2, [], []

def run_mmpbsa(decompConfigFile='/N/dc2/scratch/yx5/vs.decomp/in/decomp.in'):
    '''
    Run MMPBSA.py and decomposition
    '''
    args = [
        'MMPBSA.py',
        '-O',
        '-i',
        decompConfigFile,
        '-cp',
        'com.top',
        '-rp',
        'rec.top',
        '-lp',
        'lig.top',
        '-y',
        'com.crd'
        ]

    p = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()

    if not os.path.exists('FINAL_DECOMP_MMPBSA.dat'):
        print out, err
        return 3, [], []
    else:
        with open('FINAL_DECOMP_MMPBSA.dat') as f:
            lines = [ x for x in f ]
            if len(lines) < 10:
                print out, err
                return 3, [], []

def read_mmpbsa_dict(data, term, descriptor='gb'):
    try:
        dataCom = data[descriptor]['complex'][term]
        dataRec = data[descriptor]['receptor'][term]
        dataLig = data[descriptor]['ligand'][term]
        return np.mean(dataCom - dataRec - dataLig)
    except:
        return 9999.9

def read_mmpbsa():
    '''
    Read the MMPBSA.py output using MMPBSA.MPI
    '''
    try:
        data = MMPBSA_API.load_mmpbsa_info('_MMPBSA_info')
    except:
        return 4, [], []
    dETotal = read_mmpbsa_dict(data, 'TOTAL') #Total energy
    dEVDW = read_mmpbsa_dict(data, 'VDWAALS') #van der Waals
    dEEEL = read_mmpbsa_dict(data, 'EEL') #Electrostatic
    dEEGB = read_mmpbsa_dict(data, 'EGB') #Polar solvation
    dEESURF = read_mmpbsa_dict(data, 'ESURF') #Non-polar solvation

    dESummary = [ dEVDW, dEEEL, dEEGB, dEESURF, dETotal ]

    lResiduesDecomp = []
    for residue in data['decomp']['gb']['receptor']['TDC']:
        try:
            dEComRes = data['decomp']['gb']['complex']['TDC'][residue]['tot']
            dERecRes = data['decomp']['gb']['receptor']['TDC'][residue]['tot']
            dERes = np.mean(dEComRes - dERecRes)
        except:
            dERes = 0.0
        lResiduesDecomp.append(dERes)
        #print residue, np.mean(dERes)
    if residue != len(lResiduesDecomp):
        return 5, [], []
    else:
        return 0, dESummary, lResiduesDecomp

def decomp_wrapper(recName, recPath, mol2Name, mol2Lines):
    '''
    Decomposition workflow for a single receptor and ligand
    '''
    #Run in a temporary working directory
    tmpDir = mkdtemp()
    os.chdir(tmpDir)

    #Run antechamber and parmchk on single mol2
    err = run_lig_prep(mol2Name, mol2Lines)
    if err:
        os.chdir(cDir)
        rmtree(tmpDir)
        return err

    #tleap
    err = run_tleap(recPath)
    if err:
        os.chdir(cDir)
        rmtree(tmpDir)
        return err

    #MMPBSA.py
    err = run_mmpbsa()
    if err:
        os.chdir(cDir)
        rmtree(tmpDir)
        return err

    #Read _MMPBSA_info using API
    err, ldeltaE, lResiduesDecomp = read_mmpbsa()

    os.chdir(cDir)
    rmtree(tmpDir)

    if err:
        return err
    else:
        return 0, ldeltaE, lResiduesDecomp

def main(recDir, ligPath, outFilePath=None):
        ligFile = os.path.basename(ligPath)
        recName = ligFile.split('-')[0]
        recPath = os.path.join(recDir,
                               recName,
                               '{}_amber.pdb'.format(recName))

        if outFilePath:
            #If log file exists, read finished and append new
            completedMol2s = []
            if os.path.exists(outFilePath):
                with open(outFilePath) as outFile:
                    for outFileLine in outFile:
                        outFileLine = outFileLine.strip().split(',')
                        if len(outFileLine) > 8:
                            completedMol2s.append(outFileLine[1])
                outFile = open(outFilePath, 'a', 1)

            else:
                outFile = open(outFilePath, 'w', 1)
            completedMol2s = set(completedMol2s)

        with open(ligPath) as ligFileObject:
            for mol2Lines in read_mol2(ligFileObject):
                mol2Name = mol2Lines[1].strip()

                #Skip if finished in output file
                if outFilePath and os.path.exists(outFilePath) and mol2Name in completedMol2s:
                    continue

                err, ldeltaE, lResiduesDecomp = decomp_wrapper(recName, recPath, mol2Name, mol2Lines)
                outputLine = [ recName, mol2Name ] + ldeltaE + lResiduesDecomp

                print '\t'.join(map(str, outputLine))

                if outFilePath:
                    outFile.write(','.join(map(str, outputLine)) + '\n')
        if outFilePath:
            outFile.close()

if __name__ == '__main__':
    args = sys.argv
    if len(args) < 3:
        print __doc__
        raise SystemExit
    recDir = os.path.abspath(args[1])
    ligPath = os.path.abspath(args[2])
    if len(args) != 4:
        main(recDir, ligPath)
    else:
        outFilePath = os.path.abspath(args[3])
        main(recDir, ligPath, outFilePath)
