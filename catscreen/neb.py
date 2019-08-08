#!/usr/bin/env python

import os
from jasp import *
from ase import Atom, Atoms
import ase.io as io
from ase.lattice.surface import fcc111, fcc100, fcc211, fcc111_root, add_adsorbate
from ase.constraints import FixAtoms
import itertools as it
import random
from subprocess import Popen, PIPE, STDOUT
from dict_definitions import *

import dft 

class alloy:
    def __init__(self, **kwargs):
        """Initialization of class for alloy system
        calcdir: the directory to run calculation in
        **kwargs: all the required keywords

        """
        self.cwd = os.getcwd()
        self.wdir = kwargs['wdir']

    def setup(self, deltaz, sites, layer, conc_values,
              q, clus, **args):
        refatoms = ase.io.read('{}/CONTCAR'.format(self.wdir),format='vasp')
        host = args['host']
        cov = args['cov']
        facet = args['fac']
        #natoms = n_layer['{}-{}x{}'.format(111,cov[0],cov[1])]
        natoms = args['natoms'] 
        plusz_atom = atnum_dict[(ele_dict[host] + deltaz)]
        minusz_atom = atnum_dict[(ele_dict[host] - deltaz)]
        trans_atoms = [[plusz_atom, minusz_atom],[minusz_atom,plusz_atom]]
        all_conf = dft.gen_sites(conc_values, inert_sites=sites[0],
                                 active_sites=sites[1])

        ss = {1:6, 2:6, 3: 4, 4:1}
        #ss = {1:2, 2:2, 3: 3, 4:1}

        for nc in conc_values:
            conf = all_conf[nc]

            for num, trans in enumerate(trans_atoms):
                inert_sol = trans[0]
                active_sol = trans[1]
                all_sites = conf[0]
                if num == 1:
                    all_sites = conf[1]
                
                random.seed(500)
                rand_sites = random.sample(all_sites,ss[nc])

                #for sites in all_sites:
                for sites in rand_sites:
                    #print nc, trans, sites
                    inert_str = '.'.join(str(a) for a in sites[:nc])
                    active_str = '.'.join(str(a) for a in sites[nc:])
                    alloy_str = [inert_str, inert_sol,active_str,active_sol]
                    alloy = '{}_{}-{}_{}'.format(*alloy_str)
                        
                    if facet == 211:
                        act_site = map(lambda x: x + (natoms * (layer-1)),
                                                     sites[nc:])
                        in_site = map(lambda x: x + (natoms * 3),
                                                     sites[:nc])
                    else:
                        in_site = sites[:nc]
                        act_site = map(lambda x: x + (natoms * (layer-1)),
                                                      sites[nc:])

                    print alloy, in_site, act_site
 
                    atoms = refatoms.copy()
                    fn = '{0}/{1}/'.format(self.wdir, alloy)
                    for i in in_site:
                        atoms[i].symbol = inert_sol
                    for j in act_site:
                        atoms[j].symbol = active_sol

                    #io.write('test/{}.vasp'.format(alloy),atoms,format='vasp',
                    #          vasp5=True)
                    if not os.path.isdir(fn):
                        dft.call_jasp(atoms = atoms,
                                      kp = kp_dict[cov[1]],
                                      calcdir = fn,
                                      spin = 1,
                                      nsteps = 0,
                                      qe = q, cluster=clus)

    #def read_energy(self, calcdir):
    #    with jasp(calcdir) as calc:
    #        try:
    #            atoms = calc.get_atoms()
    #            energy = atoms.get_potential_energy()
    #            return energy
    #        except:
    #            print slabdir, 'Calc not completed'

    def __exit__(self, etype, evalue, traceback):
        """On exit, change back to the original directory.

        """
        os.chdir(self.cwd)
