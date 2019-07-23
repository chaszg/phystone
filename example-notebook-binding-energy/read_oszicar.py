#!/bin/env python

def grab_energy(oszicar):
    '''Grabs energy (E0) from OSZICAR file.'''
    with open(oszicar, 'r') as os:
        energy = float(os.readlines()[-1].split()[4])
    return energy