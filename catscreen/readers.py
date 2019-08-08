#!/usr/bin/env python

'''
Reader modules in this file.
1. Read and store electrostatic potentials from reference calculations. These
are the alchemical derivatives.
'''

def ret_esp(calc_type,ads,site,cov,lat_cov,host):
    if calc_type == 'slab':
        fn = '{2}-slab-{0}x{1}'.format(cov,lat_cov,host)
    else:
        fn = '{4}-slab-{0}x{3}-{1}-{2}'.format(cov,ads,site,lat_cov,host)
    with jasp('{0}/{1}'.format(dn,fn)) as calc:
        atoms = calc.get_atoms()
    if int(lat_cov) == 2:
        s = Popen("grep -A 8 electrostatic {0}/{1}/OUTCAR | tail -6".format(dn,fn), shell= True, stdin=PIPE,stdout=PIPE, stderr=STDOUT, close_fds=True)
    elif int(lat_cov) == 3:
        s = Popen("grep -A 10 electrostatic {0}/{1}/OUTCAR | tail -8".format(dn,fn), shell= True, stdin=PIPE,stdout=PIPE, stderr=STDOUT, close_fds=True)
    elif int(lat_cov) == 1:
        s = Popen("grep -A 6 electrostatic {0}/{1}/OUTCAR | tail -4".format(dn,fn), shell= True, stdin=PIPE,stdout=PIPE, stderr=STDOUT, close_fds=True)
    lis = s.communicate()[0]
    esp = []
    i = 0
    data = []
    fs = []
    for lines in lis.split('\n'):
        for num, pot in enumerate(lines.split(),1):
            if (num%2) == 0 and atoms.get_atomic_numbers()[i] > 10:
                atom_pos = (atoms.get_scaled_positions()[i]).tolist()
                data.append([i, atoms.get_atomic_numbers()[i],atom_pos, pot, atoms[i].index])
                esp.append(pot)
                i = i+1
    return esp

