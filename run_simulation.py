#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from openmm import app
import openmm as mm
from openmm import unit
import numpy as np
import pandas as pd
import json
import sys

from force_field import nonbonded_terms
from utils import gene_psf,CustomLangevinIntegrator

_kcal_to_kj = 4.1840
_A_to_nm = 0.1
print('Loading simulation parameters')
## Parameters about simulaiton
simu_para = json.load(open('simulation_parameters.json')) 
T = simu_para['Temperature']
friction = simu_para['Friction']
timestep = simu_para['Timestep']
equi_steps = simu_para['Equilibrium_time']
dcdperiod = simu_para['DCD_period']
logperiod = simu_para['log_period']
tot_simu_steps = simu_para['total_simu_steps']
platform_type = simu_para['platform_type']
dimension = simu_para['dimension']
initial_pdb_path = simu_para['initial_pdb_path']
output_folder_name = simu_para['out_folder_name']
output_file_name = simu_para['out_file_name']
cutoffdis = simu_para['cutoffdistance'] * _A_to_nm
box_vec = simu_para['box_size'] * _A_to_nm

print('Reading pdb')
pdb = app.PDBFile(initial_pdb_path)
top = pdb.topology
atoms = list(top.atoms())
# get Exclusion
exclusions = [[ai.index,aj.index] for i,ai in enumerate(atoms) 
                                  for j,aj in enumerate(atoms[i+1:]) 
                                  if ai.index-aj.index <=4 
                                  and ai.residue.chain.index == ai.residue.chain.index]

print('Initial system')
box_vec = np.array([[box_vec[3],0,0], [0,box_vec[4],0], [0,0,box_vec[5]]])*unit.nanometer
top.setPeriodicBoxVectors(box_vec)
forcefield = app.ForceField('force_field/force.xml')
system = forcefield.createSystem(top,nonbondedMethod=app.CutoffPeriodic,
                                    nonbondedCutoff=cutoffdis*unit.nanometer,
                                    removeCMMotion=False)

use_pbc = system.usesPeriodicBoundaryConditions()

## Create the excluded force and lj-12-6 potential
# excluded force
atom_type = list(set([ai.name for ai in atoms]))
num_atom_type = len(atom_type)
ex_epsilon_map = pd.read_csv('params/ex_epsilon.csv').to_numpy()[:,1:] * _kcal_to_kj
ex_sigma_map = pd.read_csv('params/ex_sigma.csv').to_numpy()[:,1:] * _A_to_nm
ex_atomtype = [] 
ex_chain_index = []
for ai in atoms:
    if ai.name == "DN":
        ex_atomtype.append(0)
        # there are excluded force between all DNA atoms whether in inter-chain or intra-chain,
        # so the chain index of DNA atoms are identical. 
        ex_chain_index.append(-1)
    elif ai.name == "FN":
        ex_atomtype.append(1)
        ex_chain_index.append(ai.residue.chain.index)
    elif ai.name == "FM":
        ex_atomtype.append(2)
        ex_chain_index.append(ai.residue.chain.index)
    elif ai.name == "FC":
        ex_atomtype.append(3)
        ex_chain_index.append(ai.residue.chain.index)
excluded_force = nonbonded_terms.excluded_term(ex_atomtype,ex_epsilon_map,ex_sigma_map,exclusions,
                                               chain_index=ex_chain_index,use_pbc=use_pbc,
                                               cutoff=2.0,force_group=1)

lj_epsilon_map = pd.read_csv('params/lj_epsilon.csv').to_numpy()[:,1:] * _kcal_to_kj
lj_sigma_map = pd.read_csv('params/lj_sigma.csv').to_numpy()[:,1:] * _A_to_nm
lj_atomtype = [] 
lj_chain_index = []
for ai in atoms:
    lj_chain_index.append(ai.residue.chain.index)
    if ai.name == "DN":
        lj_atomtype.append(0)
    elif ai.name == "FN":
        lj_atomtype.append(1)
    elif ai.name == "FM":
        lj_atomtype.append(2)
    elif ai.name == "FC":
        lj_atomtype.append(3)
lj_force = nonbonded_terms.LJ_12_6_term(lj_atomtype,lj_epsilon_map,lj_sigma_map,exclusions,
                                        chain_index=lj_chain_index,use_pbc=use_pbc,
                                        cutoff=2.0,force_group=2)

# Add force
print('Add excluded force')
system.addForce(excluded_force)
print('Add lj-12-6 force')
system.addForce(lj_force)

print('Create simulation')
integrator = CustomLangevinIntegrator(T*unit.kelvin,friction/unit.picosecond,timestep*unit.femtosecond,dimension)
if platform_type == "CPU":
   platform = mm.Platform.getPlatformByName('CPU')
   simulation = app.Simulation(top,system,integrator,platform)
elif platform_type == "CUDA":
  platform = mm.Platform.getPlatformByName('CUDA')
  properties = {'DeviceIndex':'0','Precision':'mixed'}
  simulation = app.Simulation(top,system,integrator,platform,properties)

# Set the parameter for exculded force between same type atoms
# Set the initial positions and velocities
positions = pdb.positions
# move conter of system geometry to box center 
box_center = 0.5*(box_vec[3] + box_vec[4] + box_vec[5])
center_of_mass = np.average(positions, axis=0)
positions = positions - center_of_mass + box_center
simulation.context.setPositions(positions*unit.nanometer)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(T*unit.kelvin)
simulation.minimizeEnergy()

# Equilibrating
simulation.step(equi_steps)
# saving the trajectory
positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
gene_psf('%s/%s'%(output_folder_name,output_file_name),top)
app.PDBFile.writeFile(top,positions,open('%s/%s.pdb'%(output_folder_name,output_file_name),'w'))
simulation.reporters.append(app.DCDReporter('%s/%s.dcd'%(output_folder_name,output_file_name),dcdperiod))
simulation.reporters.append(app.StateDataReporter('%s/%s.log'%(output_folder_name,output_file_name),logperiod,step=True,
                            potentialEnergy=True, temperature=True, progress=True,remainingTime=True,    
                            speed=True, totalSteps=tot_simu_steps,separator='\t'))

print('Running simulation!')
simulation.step(tot_simu_steps)
