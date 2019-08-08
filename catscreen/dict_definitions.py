#!/usr/bin/env python

q_dict = {'test':[1, 12,'--partition=test'],
          'mpi-opa':[2, 28,'--cluster=mpi\n#SBATCH --partition=opa'],
          'mpi-ib':[2, 20,'--cluster=mpi\n#SBATCH --partition=ib'],
          'smp':[1, 12,'--cluster=smp'],
          'mpi-legacy':[2,16,'--cluster=mpi\n#SBATCH --partition=legacy']}

# number of atoms per layer
n_layer = {'111-3x3':9,'100-3x3':9,'211-3x3':9,
           '111-2x2':4,'100-2x2':4,'211-3x2':6,
           '111_root3-1x1':3}

# kPoint dictionary
kp_dict = {1: 8, 2: 8, 3: 4}

# Coordinates for molecules
coord = {1: [(0.,0.,0.)], # O,N,C,H
         2: [(0.,0.,0.),(0.,0.,1.5)], #OH, CH, NH
         3: [(0.,0.,0.),(0.,-0.75,0.75),(0.,0.75,0.75)], #CH2,NH2,OH2
         4: [(0.,0.,0.),(-0.5,-0.75,0.75),(-0.5,0.75,0.75), # CH3, NH3
             (0.8,0.20,0.75)],
         5: [(0.,0.,0.),(-1.3,0.,0.3),(-1.3,0.,1.2)]} #OOH

# 'Ads':[molid in coord, no. of atoms in ads]
ads_dict = {'C':[1,1], 'N':[1,1], 'O':[1,1], 'H':[1,1],
            'CH':[2,2], 'NH': [2,2], 'OH': [2,2],
            'CH2':[3,3], 'NH2': [3,3], 'OH2': [3,3], 'OOH': [5,3],
            'CH3':[4,4], 'NH3':[4,4], 'CH4':[4,5]}

atnum_dict = {20:'Ca', 21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn',
              26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 39:'Y',
              40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh',
              46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 73:'Ta', 74:'W',
              75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg',
              81:'Tl', 82:'Pb', 83:'Bi'}

ele_dict = {'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79,
            'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Ca':20, 'Sc':21, 'Ti':22,
            'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29,
            'Zn':30, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44,
            'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49}

# Energy of molecules
e_h2 = -6.75846269
e_o2 = -9.84818417
e_c = -18.572761/2  # Value copied from eric_graphene/pby/
e_n = -16.641257/2 # Value copied from eric_graphene/pby/
e_oh = -7.52462863
e_ch = -5.58103363
e_nh = -6.26796846
e_ch2 = -11.34463559
e_nh2 = -12.76256159
e_ch3 = -17.52129266
e_h2o = -14.21910042
e_nh3 = -19.51488381

e_mol_dict = {'H':0.5*e_h2, 'O':e_h2o-2*0.5*e_h2, 'C':e_c, 'N':e_n,
              'OH':e_h2o-1*0.5*e_h2, 'CH':e_ch, 'NH':e_nh, 'CH2':e_ch2, 
              'NH2':e_nh2, 'OH2':e_h2o,'CH3':e_ch3, 'NH3':e_nh3, 
              'OOH': 2*e_h2o-3*0.5*e_h2, 'O2': e_o2}
