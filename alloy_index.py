#!/usr/bin/env python

import os
from ase import Atom

def index_transmuted(slab,transmute_atom_sym,counter_atom_sym,transmute_num,counter_num,symmetric):
    '''Identifies in slab the indexes of atoms to be transmuted near the surface by knowing chemical symbol (transmute_atom_sym) and finds indexes of atoms far from the surface to be counter transmuted by knowing the symbol (counter_atom_sym). 
    The number of indexes found controlled by inputs (transmute_num and counter_num).
    Works for symmetric (symmetric = TRUE) and non-symmetric (symmetric = FALSE) slab.'''
    #Dictionary of all metal atoms in the surface
    transmute_atom = {}
    counter_atom = {}
    #This currently works if you loop through either ads or slab
    for s in slab:
        #Identifying metal atoms
        if s.symbol == transmute_atom_sym:
            transmute_atom[str(s.index)] = s.position[2]
        if s.symbol == counter_atom_sym:
            counter_atom[str(s.index)] = s.position[2]
            
    #List of atom indexes transmuted at surface
    #List of atom indexes counter transmuted from from surface
    transmute = []
    counter = []
       
    if symmetric:
        for i in range(0,transmute_num):
            topmax = max(transmute_atom, key = transmute_atom.get)
            transmute.append(int(topmax))
            del transmute_atom[topmax]
        
        
        for j in range(0,transmute_num):
            botmin = min(transmute_atom, key = transmute_atom.get)
            transmute.append(int(botmin))
            del transmute_atom[botmin]
            
        for dex in transmute_atom.keys():
            counter.append(int(dex))
            
    else:
        for k in range(0,transmute_num):
            topmax = max(transmute_atom, key = transmute_atom.get)
            transmute.append(int(topmax))
            del transmute_atom[topmax]
        
        
        for l in range(0,counter_num):
            botmin = min(counter_atom, key = counter_atom.get)
            counter.append(int(botmin))
            del counter_atom[botmin]
            
    return transmute, counter

def transmuter(slab,atomdex,trans):
	'''Takes an atoms object of a surface, and transmutes atoms specified by an array of indexes and an array of atom objects.
    len(atomdex) must equal len(trans)'''
	slabcopy = slab.copy()
	for i in range(0,len(atomdex)):
		slabcopy[atomdex[i]].symbol = trans[i].symbol
	return slabcopy

def transmuted_directory_names(bdex,tdex,dexes,atoms_array):
	ad = ''
	for i in range(0,len(dexes)):
		ad = ad + '.' + atoms_array[i].symbol + str(dexes[i])
	dir_name = str(bdex) + '.' + str(tdex) + ad
	return dir_name