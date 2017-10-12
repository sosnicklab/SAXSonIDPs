import numpy as np
import tables as tb
from matplotlib.pyplot import *
import subprocess as sp
import collections
from glob import glob
import re
import sys,os
import cPickle as cp
import IPython.display as disp
import uuid, json
import time
import argparse


upside_dir = '/home/fillin/upside/'
params_dir = '/home/fillin/upside-parameters/'
run_dir    = '/scratch/fillin/runupside/'

if upside_dir + 'src' not in sys.path: sys.path = [upside_dir+'src'] + sys.path
import PDB_to_initial_structure as ps

deg = np.pi/180.

from IPython.display import clear_output


def upside_config(fasta, output, n_system, dimer=False, backbone=True, rotamer=False, 
                  sidechain=None,sidechain_exclude=[], hbond=None, sheet_mix_energy=None, environment=None,
                  init=None, rama_pot='rama_libraries.h5',rama_pot_product=False, reference_rama=None, restraint_groups=[], restraint_spring=None,
                  hbond_coverage_radius=None,target=[]):
    
    args = [upside_dir + 'src/upside_config.py', '--fasta=%s'%fasta, '--output=%s'%output, '--n-system=%i'%n_system]

    if init:
        args.append('--initial-structures=%s'%init)
    for targ in target:
        args.append('--target-structures=%s'%targ)
    if rama_pot!=None:
        args.append('--rama-library=%s'%(params_dir+rama_pot))
    if rama_pot_product:
        args.append('--rama-library-combining-rule=product')
    if sheet_mix_energy is not None:
        args.append('--rama-sheet-mixing-energy=%f'%sheet_mix_energy)
    if dimer:
        args.append('--dimer-basin-library=%s'%(params_dir+'TCB_count_matrices.pkl'))
    if not backbone:
        args.append('--no-backbone')
    if hbond:
        args.append('--hbond-energy=%f'%hbond)
    if reference_rama:
        args.append('--reference-state-rama=%s'%reference_rama)
    for rg in restraint_groups:
        args.append('--restraint-group=%s'%rg)
    if restraint_spring is not None:
        args.append('--restraint-spring-constant=%f'%restraint_spring)
        
    if rotamer or environment:
        args.append('--rotamer-placement=%s'%(params_dir+'rotamer-MJ-1996-with-direc.h5'))

    if rotamer:
        args.append('--rotamer-interaction=%s'%rotamer_interaction_param)
    if environment:
        args.append('--environment=%s'%environment)
    
    if sidechain!=None:
        args.append('--sidechain-radial=%s'%(params_dir+sidechain))
        #args.append('--backbone-dependent-point=%s'%(params_dir+'backbone_dependent_com.h5'))
    if len(sidechain_exclude)>0:
        args.append('--sidechain-radial-exclude-residues=%s' % sidechain_exclude[0])
        
    # print ' '.join(args)

    return ' '.join(args) + '\n' + sp.check_output(args)

def compile():
    return sp.check_output(['/bin/bash', '-c', 'cd %s; make -j4'%(upside_dir+'obj')])


UpsideJob = collections.namedtuple('UpsideJob', 'job config output'.split())


def run_upside(queue, config, duration, frame_interval, n_threads=1, hours=36, temperature=1., seed=None,
               replica_interval=None, anneal_factor=1., anneal_duration=-1., pivot_interval=None, 
               time_step = None, swap_sets = None,
               log_level='basic', account=None, estimate_rama_bandwidth=None):
    if isinstance(config,str): config = [config]
    
    upside_args = [upside_dir+'obj/upside', '--duration', '%f'%duration, '--frame-interval', '%f'%frame_interval] + config

    try:
        upside_args.extend(['--temperature', ','.join(map(str,temperature))])
    except TypeError:  # not iterable
        upside_args.extend(['--temperature', str(temperature)])

    if replica_interval is not None:
        upside_args.extend(['--replica-interval', '%f'%replica_interval])
        for s in swap_sets:
            upside_args.extend(['--swap-set', s])
    if pivot_interval is not None:
        upside_args.extend(['--monte-carlo-interval', '%f'%pivot_interval])
    if anneal_factor != 1.:
        upside_args.extend(['--anneal-factor', '%f'%anneal_factor])
    if anneal_duration != -1.:
        upside_args.extend(['--anneal-duration', '%f'%anneal_duration])
    upside_args.extend(['--log-level', log_level])
    
    if time_step is not None:
        upside_args.extend(['--time-step', str(time_step)])

    upside_args.extend(['--seed','%li'%(seed if seed is not None else np.random.randint(1<<31))])
    
    output_path = config[0]+'.output'
    
    if estimate_rama_bandwidth:
        for configi in config:
            fn = configi + '.rama.pkl'
            if os.path.exists(fn):
                os.remove(fn)

            upside_args.extend(['&&', 
                                upside_dir + '/src/estimate_rama_distributions.py', 
                                '--system=0',
                                '--bandwidth=%f'%estimate_rama_bandwidth,
                                configi,
                                fn])

    if queue == '': 
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(n_threads)
        output_file = open(output_path,'w')
        job = sp.Popen(upside_args, stdout=output_file, stderr=output_file)
    elif queue == 'srun':
        # set num threads carefully so that we don't overwrite the rest of the environment
        # setting --export on srun will blow away the rest of the environment
        # afterward, we will undo the change

        old_omp_num_threads = os.environ.get('OMP_NUM_THREADS', None)

        try:
            os.environ['OMP_NUM_THREADS'] = str(n_threads)
            args = ['srun', '--ntasks=1', '--nodes=1', '--cpus-per-task=%i'%n_threads, 
                    '--slurmd-debug=0', '--output=%s'%output_path] + upside_args
            job = sp.Popen(args, close_fds=True)
        finally:
            if old_omp_num_threads is None:
                del os.environ['OMP_NUM_THREADS']
            else:
                os.environ['OMP_NUM_THREADS'] = old_omp_num_threads
    else:
        #print upside_args
        args = ['sbatch', '-p', queue, '--time=0-%i'%hours, '--ntasks=1', 
                '--cpus-per-task=%i'%n_threads, '--export=OMP_NUM_THREADS=%i'%n_threads,
                '--output=%s'%output_path, '--parsable', '--wrap', ' '.join(upside_args)]
        if account is not None:
            args.append('--account=%s'%account)
        
        print args
        job = sp.check_output(args).strip()
    print job,config,output_path
    return UpsideJob(job,config,output_path)


def status(job):
    try:
        job_state = sp.check_output(['/usr/bin/env', 'squeue', '-j', job.job, '-h', '-o', '%t']).strip()
    except sp.CalledProcessError:
        job_state = 'FN'
        
    if job_state == 'PD':
        status = ''
    else:
        status = sp.check_output(['/usr/bin/env','tail','-n','%i'%1, job.output])[:-1]
    return '%s %s' % (job_state, status)

def statusi(job):
    try:
        job_state = sp.check_output(['/usr/bin/env', 'squeue', '-j', job, '-h', '-o', '%t']).strip()
    except sp.CalledProcessError:
        job_state = 'FN'
        
    if job_state == 'PD':
        status = ''
    
    return '%s' % (job_state)

def swap_table2d(nx,ny):
    idx = lambda xy: xy[0]*ny + xy[1]
    good = lambda xy: (0<=xy[0]<nx and 0<=xy[1]<ny)
    swap = lambda i,j: '%i-%i'%(idx(i),idx(j)) if good(i) and good(j) else None
    horiz0 = [swap((a,b),(a+1,b)) for a in range(0,nx,2) for b in range(0,ny)]
    horiz1 = [swap((a,b),(a+1,b)) for a in range(1,nx,2) for b in range(0,ny)]
    vert0  = [swap((a,b),(a,b+1)) for a in range(0,nx)   for b in range(0,ny,2)]
    vert1  = [swap((a,b),(a,b+1)) for a in range(0,nx)   for b in range(1,ny,2)]
    sets = (horiz0,horiz1,vert0,vert1)
    sets = [[y for y in x if y is not None] for x in sets]
    return [','.join(x) for x in sets if x]


jobs = dict()
compile()
