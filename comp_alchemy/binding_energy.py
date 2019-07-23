#!/bin/env python
from numpy import dot, zeros

def alc_be(transmute,counter,espdiffs,charge):
    '''Calculates delta binding energy with alchemy.'''
    transmute_charges = [charge for i in range(len(transmute))]
    counter_charges = [-charge for j in range(len(counter))]
    
    dn = zeros(len(espdiffs))
    
    for k,t in enumerate(transmute):
        dn[t] = transmute_charges[k]
        
    for l,c in enumerate(counter):
        dn[c] = counter_charges[l]
    
    be = dot(espdiffs, dn)
        
    return [dn, be]
        