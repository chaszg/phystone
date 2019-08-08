#!/usr/bin/env python

'''
Utility modules in this file.
1. Find sites/atoms in  with unique alchemical derivatives. 

'''

def esp_layer(esp,natoms,layer):
    st = natoms*(layer-1)
    end = natoms*layer
    return [float(e) for e in esp[st:end]]

def uniquefy(layer):
    print layer
    ind = [] 
    for i in range(len(layer)):
        for j in range(i+1,len(layer)):
            diff = layer[i] - layer[j]
            if abs(round(diff,4)) < 0.001:
                ind.append(j)
    return [i for i in range(len(layer)) if i not in ind]

def unique_esp(esp,n,natoms,facet,layer):

    if facet == 211:
        inert_layer = esp_layer(esp, natoms, 4) 
        active_layer = esp_layer(esp, natoms, layer) 
    else:
        inert_layer = esp_layer(esp, natoms, 1) 
        active_layer = esp_layer(esp, natoms, layer)

    unique_inert_sites = uniquefy(inert_layer)    
    unique_active_sites = uniquefy(active_layer)    
    return unique_inert_sites, unique_active_sites

