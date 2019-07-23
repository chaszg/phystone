#!/usr/bin/env python
import os
from ase import Atom
from ase.io import read as r
from ase.visualize import view
from math import sqrt


def pairs(slab, ads):
	"""Using atomic distance in unit cell to match the indices of atoms in ads with atoms in slab"""
	sp = []
	for s in slab:
		#Calculating the distance of each atom in slab from reference points
		d = sqrt(s.position[0]**2 + s.position[1]**2 + (s.position[2])**2)
		d1 = sqrt((s.position[0]-10)**2 + (s.position[1]-10)**2 + (s.position[2])**2)
		d2 = sqrt((s.position[0]-0)**2 + (s.position[1]-10)**2 + (s.position[2])**2)
		d3 = sqrt((s.position[0]-10)**2 + (s.position[1]-0)**2 + (s.position[2])**2)
		#List of slab atom index and distances
		sp.append([s.index, d, d1, d2, d3])

	ap = []
	for a in ads:
		#Calculating the distance of each atom in ads from reference points
		d = sqrt(a.position[0]**2 + a.position[1]**2 + (a.position[2])**2)
		d1 = sqrt((a.position[0]-10)**2 + (a.position[1]-10)**2 + (a.position[2])**2)
		d2 = sqrt((a.position[0]-0)**2 + (a.position[1]-10)**2 + (a.position[2])**2)
		d3 = sqrt((a.position[0]-10)**2 + (a.position[1]-0)**2 + (a.position[2])**2)
		#List of ads atom index and distances
		ap.append([a.index, d, d1, d2, d3])
	
	#Tolerance used when calculating difference in distances to obtain matches
	d_tol = 0.3
	pair = []
	for s in sp:
		mult = []
		for a in ap:
			if abs(s[1]-a[1]) <= d_tol and abs(s[2]-a[2]) <= d_tol and abs(s[3]-a[3]) <= d_tol and abs(s[4]-a[4]) <= d_tol:
				#List of matching atom index in ads with two differences in distance
				mult.append([a[0], (s[1]-a[1]), s[2]-a[2]])
		
		#Placeholder in mult if no matches are found
		if len(mult) == 0:
			mult.append(['n'])
		#List of index pairs found (slab, ads)
		pair.append([s[0], mult[0][0]])

	return pair