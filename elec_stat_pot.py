#!/bin/env python
from subprocess import Popen, PIPE, STDOUT
from ase import Atom
from ase.io import read
from math import ceil

def grab_esp(poscar,outcar):
    '''Finds electrostatic potentials printed in a VASP OUTCAR (outcar) corresponding to the system described by the POSCAR (poscar).
Electrostatic potentials are indexed with atom indexes defined with ASE.'''
    slab = read(poscar,format='vasp')
    atoms = float(len(slab))
    rows = int(ceil(atoms/5))
    add = 3
    elec = []
    with open(outcar, 'r') as out:
        lines = out.readlines()
    for i,l in enumerate(lines):
        s = l.split()
        for j,c in enumerate(s):
            if c == '(electrostatic)':
                d = i
    for k in range(d+add,d+add+rows):
        for e in lines[k].split()[1::2]:
            elec.append(float(e))
    return elec

def espdiff(elec1,elec2,pair):
    '''Calculating differences in electrostatic potential per atom between systems 1 and 2.
Diff = esp_2 - esp_1
Obtain electrostatic potential array with grab_esp().

First and second arguments are arrays of electrostatic potentials for systems 1 and 2. 
Third is index pairs between elec1 and elec2
'''
    diffs = []
    for p in pair:
        diff = elec2[p[1]]-elec1[p[0]]
        diffs.append(diff)
    return diffs

def remove_duplicate_espdiffs(dexlist,espdiffs):
    '''Removes indexes of atoms in transmute list with repeated electrostatic potential differences'''
    #Pairing espdiff values with corresponding atom index in dexlist
    tol = 0.01
    espd_dict = {}
    for i in dexlist:
        espd_dict[str(i)] = espdiffs[i]
    #Grabbing unique espdiff values
    unique_diffs = {}
    for dex,diff in espd_dict.items():
        if diff not in unique_diffs.values():
            unique_diffs[dex] = diff
    print('Unique Electrostatic Potential Differences:')
    print(unique_diffs.values())
    print('')

    copy_unique_diffs = unique_diffs.copy()
    for dex,diff in copy_unique_diffs.items():
        c = 0
        for dex2,diff2 in unique_diffs.items():
            if abs(diff2 - diff) < tol and abs(diff2 - diff) != 0:
                c += 1
        if c > 0:
            del unique_diffs[dex]
    
    print('Unique Electrostatic Potential Differences WITHIN A TOLERANCE VALUE:')
    unique_vals = [u for u in unique_diffs.values()]
    print(unique_vals)
    print(' ')
    unique_dexes = []
    for dex in unique_diffs.keys():
        unique_dexes.append(int(dex))
    
    return unique_dexes

def heatmap(poscar,dexlist,espdiffs):
    '''Assigns electrostatic potential differences to each atom in system as an initial charge in ASE GUI to visualize.'''
    system = read(poscar,format='vasp')
    system.set_initial_charges(espdiffs)
    return system
        