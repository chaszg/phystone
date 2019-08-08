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

def add_sbr(atoms,ads,lat_cov):
    at_index = lat_cov * 4 - 1
    xy_atom = atoms[at_index]
    x = xy_atom.position[0]
    y = xy_atom.position[1]
    h = 1.0
    return x,y,h

def add_br100(atoms,ads,lat_cov):
    y = atoms.cell[1][1]/(lat_cov*2)
    x = atoms.cell[0][0] - 0.6
    h = 1.0
    return x,y,h

def add_tt(atoms,ads,lat_cov):
    apl = lat_cov * 1
    x = (atoms[2*apl - 1].position[0] + atoms[3*apl - 1].position[0]) / 2
    y = atoms[2*apl - 1].position[1]
    h = 0.5
    return x,y,h

def add_hollow(atoms,ads,lat_cov):
    apl = lat_cov * 1
    x = atoms[apl].position[0] + 0.45
    y = atoms.cell[1][1]/(lat_cov*2)
    h = 0.5
    return x,y,h

def root3_111(atoms,ads,site):
    if site == 'fcc':
        x = atoms[7].position[0]
        y = atoms[7].position[1]
        h = 1.5
    elif site == 'bridge':
        x = atoms[3].position[0] - 0.25
        y = atoms[3].position[1] + 0.55
        h = 1.5
    elif site == 'ontop':
        x = atoms[9].position[0]
        y = atoms[9].position[0]
        h = 1.5
    return x,y,h

def call_jasp(atoms, kp, calcdir, spin, nsteps, qe, cluster):

    JASPRC['queue.nodes'] = q_dict[qe][0]
    JASPRC['queue.ppn'] = q_dict[qe][1]
    JASPRC['queue.cluster'] = q_dict[qe][2]

    with jasp(calcdir,
        xc = 'PBE',
        kpts = [kp, kp, 1],
        encut = 350,
        nsw = nsteps,
        ibrion = 2,
        lwave = False,
        lcharg = False,
        ispin = spin,
        idipol = 3,
        ediff = 1e-05,
        ediffg = -0.05,
        atoms = atoms) as calc:
	    try:
	    	return atoms.get_potential_energy()
	    except (VaspSubmitted, VaspQueued):
	    	print "Submitted to queue"

def ref_calc_is_ok(slabdir,adsdir):
    slabenergy = 0
    with jasp(slabdir) as calc:
        try:
            atoms = calc.get_atoms()
            slabenergy = atoms.get_potential_energy()
        except:
            print slabdir, 'Calc not completed'

    adsenergy = 0
    with jasp(adsdir) as calc:
        try:
            atoms = calc.get_atoms()
            adsenergy = atoms.get_potential_energy()
        except:
            print adsdir, 'Calc not completed'

    if (slabenergy == 0) or (adsenergy == 0):
        return False
    else:
        return True

def facet_setup(host,fac,c,lat,root=None):
    if fac == 111:
        if root:
            atoms = fcc111_root(host, root=root, a=lat, size=(c[0],c[1],4),
                                vacuum = 12.0)
        else:
            atoms = fcc111(host, a=lat, size=(c[0],c[1],4),
                           vacuum=12.0)
    elif fac == 100:
        atoms = fcc100(host, a=lat, size=(c[0],c[1],4),
                       vacuum=12.0)

    elif fac == 211:
        atoms = fcc211(host, a=lat, size=(c[0],c[1],4),
                       vacuum=12.0)

    return atoms

    
def gen_sites(conc_values, inert_sites, active_sites):
    #all_combos = []
    all_conf = {}
    for conc in conc_values:
        all_combos = []
        for is_comb in it.combinations(inert_sites, conc):
            for as_comb in it.combinations(active_sites, conc):
                all_combos.append(list(is_comb)+list(as_comb))

        if len(all_combos) <= 15:
            sample_size =  len(all_combos)
        else:
            sample_size = 15
        random.seed(500)

        # Choose random alloy configurations of given sample size from all
        # possible configs. Making sure all are unique combinations.
        # combos and inv_combos are mutually exclusive site selections
        combos = random.sample(all_combos,sample_size)
        ind = [i for i, a in enumerate(all_combos) if a in combos]
        for id in sorted(ind,reverse=True):
            del all_combos[id]
        if len(all_combos) < sample_size:
            inv_combos = all_combos + random.sample(combos,
                                  sample_size-len(all_combos))
        else:
            inv_combos = random.sample(all_combos, sample_size)
        combos = [c for c in combos if len(c) == conc*2]
        inv_combos = [c for c in inv_combos if len(c) == conc*2]
        all_conf[conc] = [combos, inv_combos]
    return all_conf


class ref(object):
    def __init__(self, **kwargs):
        """Initialization of class for reference system
        calcdir: the directory to run calculation in
        **kwargs: all the required keywords

        """
        self.cwd = os.getcwd()
        #self.calcdir = calcdir

    def slab_setup(self, calcdir, q, clus, **kwargs):
        c = kwargs['cov']
        atoms = facet_setup(host=kwargs['host'],fac=kwargs['fac'],
                            c=kwargs['cov'],lat=kwargs['lat'],
                            root=kwargs['root'])

        if kwargs['fac'] == 211:
            constraint = FixAtoms(indices=[atom.index for atom in atoms
                                           if atom.z < 14.0])
        else:
            constraint = FixAtoms(indices=[atom.index for atom in atoms
                                           if atom.z < 12.2])
                                           #if atom.z < 12.2])
        atoms.set_constraint(constraint)
        call_jasp(atoms = atoms,
                  kp = kp_dict[c[1]],
                  calcdir = calcdir,
                  spin = 1,
                  nsteps = 500,
                  qe = q, cluster=clus)

    def ads_setup(self, calcdir, q, clus, **kwargs):
        c = kwargs['cov']
        ads = kwargs['ads']

        atoms = facet_setup(host=kwargs['host'],fac=kwargs['fac'],
                            c=kwargs['cov'],lat=kwargs['lat'],
                            root=kwargs['root'])

        if kwargs['fac'] == 211:
            mol = Atoms(ads, coord[ads_dict[ads][0]])
            ads_call = {'sbr':add_sbr, 'br':add_br100,
                        'tt':add_tt, 'hollow':add_hollow}
            x,y,h = ads_call[kwargs['site']](atoms,mol,c[1])
            add_adsorbate(slab=atoms, adsorbate=mol, height=h,
                          position=(x,y))
            constraint = FixAtoms(indices=[atom.index for atom in atoms
                                           if atom.z < 14])

        else: #111 or 100
            mol = Atoms(ads, coord[ads_dict[ads][0]])
            if kwargs['root'] is None:
                print kwargs['site']
                add_adsorbate(atoms, mol, height=2.0, position=kwargs['site'])
            else:
                x,y,h = root3_111(atoms,mol,kwargs['site'])
                add_adsorbate(slab=atoms, adsorbate=mol, height=h,
                              position=(x,y))

            constraint = FixAtoms(mask=[True for atom in atoms
                                        if atom.position[2] < 12.2])
                                        #if atom.position[2] < 12.2])

        atoms.set_constraint(constraint)
        call_jasp(atoms = atoms,
                  kp = kp_dict[c[1]],
                  calcdir = calcdir,
                  spin = 1,
                  nsteps = 500,
                  qe = q, cluster=clus)

    def read_energy(self, calcdir):
        with jasp(calcdir) as calc:
            try:
                atoms = calc.get_atoms()
                energy = atoms.get_potential_energy()
                return energy
            except:
                print slabdir, 'Calc not completed'

    def ret_esp(self, calcdir, cov, ads=None):
        '''
        Returns the electrostatic potential of reference calculations.
        Current code works only for single host reference calccs.
        This will not work for Pt3Ni skin as the reference.

        '''
        nl_dict = {1:6, 2:8, 3:10}
        nlines = nl_dict[cov]
        n = 0
        if (ads):
            n = ads_dict[ads][1]
        cmd = 'grep -A {} electrostatic {}/OUTCAR | tail -{}'.format(nlines,
                                                                     calcdir,
                                                                     nlines-2)
        s = Popen(cmd,shell=True,
                  stdin=PIPE,stdout=PIPE,
                  stderr=STDOUT,close_fds=True)
        lis = s.communicate()[0]
        esp = []
        for lines in lis.split('\n'):
            for num, pot in enumerate(lines.split(),1):
                if (num%2) == 0:
                    esp.append(pot)
        last = len(esp)-n
        return esp[:last]

    def __exit__(self, etype, evalue, traceback):
        """On exit, change back to the original directory.

        """
        os.chdir(self.cwd)

class alloy:
    def __init__(self, **kwargs):
        """Initialization of class for alloy system
        calcdir: the directory to run calculation in
        **kwargs: all the required keywords

        """
        self.cwd = os.getcwd()
        self.slabdir = kwargs['slabdir']
        self.adsdir = kwargs['adsdir']

        if not ref_calc_is_ok(self.slabdir, self.adsdir):
            print 'Fix ref calcs'
            sys.exit()
        # Object created without any issues

    def setup(self, calctype, deltaz, sites, layer, conc_values,
              q, clus, **args):
        with jasp(self.slabdir) as calc:
            slabatoms = calc.get_atoms()

        with jasp(self.adsdir) as calc:
            adsatoms = calc.get_atoms()

        host = args['host']
        cov = args['cov']
        facet = args['fac']
        #natoms = n_layer['{}-{}x{}'.format(111,cov[0],cov[1])]
        natoms = args['natoms'] 
        plusz_atom = atnum_dict[(ele_dict[host] + deltaz)]
        minusz_atom = atnum_dict[(ele_dict[host] - deltaz)]
        trans_atoms = [[plusz_atom, minusz_atom],[minusz_atom,plusz_atom]]
        all_conf = gen_sites(conc_values, inert_sites=sites[0],
                             active_sites=sites[1])

        for nc in conc_values:
            conf = all_conf[nc]

            for num, trans in enumerate(trans_atoms):
                inert_sol = trans[0]
                active_sol = trans[1]
                all_sites = conf[0]
                if num == 1:
                    all_sites = conf[1]
                for sites in all_sites:
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
 
                    if calctype == 'slab':
                        atoms = slabatoms.copy()
                        fn = 'slab/{}x{}/{}'.format(cov[0],cov[1],alloy)
                        for i in in_site:
                            atoms[i].symbol = inert_sol
                        for j in act_site:
                            atoms[j].symbol = active_sol

                    elif calctype == 'ads':
                        atoms = adsatoms.copy()
                        ads = args['ads']
                        site = args['site']
                        n = ads_dict[ads][1] # no. of atoms in ads
                        alloy_str = [ads,site,cov[0],cov[1],alloy,ads]
                        fn = '{}_BE/{}/{}x{}/{}_{}'.format(*alloy_str)
                        for i in in_site:
                            atoms[i].symbol = inert_sol
                        for j in act_site:
                            atoms[j].symbol = active_sol

                    #io.write('test/{}.vasp'.format(alloy),atoms,format='vasp',
                    #          vasp5=True)
                    if not os.path.isdir(fn):
                        call_jasp(atoms = atoms,
                                  kp = kp_dict[cov[1]],
                                  calcdir = fn,
                                  spin = 1,
                                  nsteps = 0,
                                  qe = q, cluster=clus)

    def read_energy(self, calcdir):
        with jasp(calcdir) as calc:
            try:
                atoms = calc.get_atoms()
                energy = atoms.get_potential_energy()
                return energy
            except:
                print slabdir, 'Calc not completed'

    def __exit__(self, etype, evalue, traceback):
        """On exit, change back to the original directory.

        """
        os.chdir(self.cwd)
