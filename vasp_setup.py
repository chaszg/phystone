#!/usr/bin/env python
import os
import sys
potdir = '~/potcars/'
def make_potcar(wdir,potdir):
	with open('{0}/POSCAR'.format(wdir),'r') as poscar:
		poslines = poscar.readlines()
	atomlist = poslines[0].split()
	command = 'cat'
	for a in atomlist:
		command = command + ' ' + '{0}{1}/POTCAR'.format(potdir,a)
	command = command + ' >> {0}/POTCAR'.format(wdir)
	os.system('{0}'.format(command))
    
def make_kpoints(wdir,name,kspacing):
    with open('{0}/KPOINTS'.format(wdir),'w') as kpt:
        kpt.write('''{0}
0
M
{1} {2} {3}
0 0 0
'''.format(name,kspacing[0],kspacing[1],kspacing[2]))

def make_incar(wdir,name,encut,ismear,sigma,ibrion,nsw):
    with open('{0}/INCAR'.format(wdir),'w') as inc:
        inc.write('''SYSTEM = {0}

#-- Start parameters --
ISTART  = 0
LCHARG  = .FALSE.
LWAVE   = .FALSE.
LVTOT   = .FALSE.

#-- Electronic parameters --
NELM   = 1000
NELMIN = 5
ENCUT  = {1}
ISMEAR = {2}
SIGMA  = {3}
EDIFF  = 1.0e-5
LREAL  = Auto
LPLANE = .TRUE.
NPAR   = 4
LSCALU = .FALSE.
NSIM   = 4
IDIPOL = 3

#-- MD
IBRION ={4}
NSW={5}'''.format(name,encut,ismear,sigma,ibrion,nsw))
        
def make_job_script(wdir,name,jobnum,nodes,cores,cluster,partition,hours):
    with open('{0}/job_sub.slurm'.format(wdir),'w') as js:
        js.write('''#!/bin/bash

#SBATCH --job-name="{0}-{1}"
#SBATCH --nodes={2}
#SBATCH --ntasks={3}
#SBATCH --cluster={4}
#SBATCH --partition={5}
#SBATCH --error=VASP-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cdg36@pitt.edu
#SBATCH --time={6}:00:00

set -v
ulimit -s unlimited

module purge
module load intel/2017.1.132
module load intel-mpi/2017.1.132
module load mkl
module load fftw
module load vasp/5.4.4


# BEFORE running section
echo "JOB_ID: $SLURM_JOB_ID JOB_NAME: $SLURM_JOB_NAME" >> runstats.out
before =$(date +%s)
echo "The JOB started on : $(date)" >> runstats.out

# RUN secion
srun --mpi=pmi2 vasp_std  >& stdout.prod

# AFTER running section
after=$(date +%s)
elapsed_seconds=$(expr $after - $before)
echo "The JOB ended on: $(date)" >> runstats.out
echo "The JOB ran for: $elapsed_seconds seconds" >> runstats.out'''.format(name,jobnum,nodes,cores,cluster,partition,hours))