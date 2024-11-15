#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 09:28:54 2021
Descripations:
    simulation 1d: dimension = 1
    simulation 2d: dimension = 2

@author: windyer
"""
import openmm as mm
from openmm import unit
import numpy as np

def CustomLangevinIntegrator(temperature=298.0*unit.kelvin, collision_rate=91.0/
                             unit.picoseconds, timestep=1.0*unit.femtoseconds,dimensions=2):
    
    kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
    # Compute constants.
    kT = kB * temperature

    gamma = collision_rate
    # Create a new custom integrator
    integrator = mm.CustomIntegrator(timestep) 
    
    if dimensions == 2:
        integrator.addPerDofVariable("dumv",1.0)
        integrator.setPerDofVariableByName("dumv",[mm.Vec3(x=1.0,y=1.0,z=0.0)])
    # Integrator initialization
    integrator.addGlobalVariable("kT",kT) # thermal energy
    integrator.addGlobalVariable("T",temperature) # temperature 
    integrator.addGlobalVariable("b", np.exp(-gamma*timestep)) # velocity mixing parameter
    integrator.addPerDofVariable("sigma_i",0)
    integrator.addPerDofVariable("x1",0) # position before application of constraints

    # Allow context updating here
    integrator.addUpdateContextState();
    
    # Pre-computation.
    # This only needs to be done once, but it needs to be done for each degree of freedom.
    # Could move this to initialization?
    #
    integrator.addComputePerDof("sigma_i", "sqrt(kT/m)")
    
    # Veclocity perturbation
    integrator.addComputePerDof("v","sqrt(b)*v + sqrt(1-b)*sigma_i*gaussian")
    integrator.addConstrainVelocities();
    
    # Metropolized symplectic step
    integrator.addComputePerDof("v", "v + 0.5*dt*f/m")
    if dimensions == 2: # To get a 1D system, make y and z velocities zero when moving x
       integrator.addComputePerDof("x","x + v*dumv*dt")
    else:
        integrator.addComputePerDof("x","x+v*dt")
    integrator.addComputePerDof("x1","x")
    integrator.addConstrainPositions()
    integrator.addComputePerDof("v","v + 0.5*dt*f/m + (x-x1)/dt")
    integrator.addConstrainVelocities();
    
    # Velocity randomization
    integrator.addComputePerDof("v", "sqrt(b)*v + sqrt(1-b)*sigma_i*gaussian")
    integrator.addConstrainVelocities();
    
    if dimensions == 2: # Remove the resulting y and z velocities to get the correct Kinectic Energy
       integrator.addComputePerDof("v","v*dumv")
       
    return integrator
