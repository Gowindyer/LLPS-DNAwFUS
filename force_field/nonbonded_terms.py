import openmm as mm
from openmm import unit
import pandas as pd
import numpy as np

def excluded_term(atom_types,epsilon_map,sigma_map,exclusions,extra_exclusions=None,chain_index=None,use_pbc=False,cutoff=2.5,force_group=1):
    """
    An excluded interaction between atoms is implemented to prevent overlap.

    Parameter
    ---------
    atom_type: list, int 
        Atom type.
    
    epsilon_map: array
        the strength of specific residue pair

    sigma_map: array
        the rest length of specific residue pair
    
    exclusions: list
        particular pairs of particles whose interactions should be omitted from force and energy calculations.
    
    extra_exclusions: list
        This represents the extra particle pairs that are neglected in the calculation of the potential.
    
    chain_index: list
        the chain index of any atom.

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.

    cutoff: float
        cutoff = truncated distance / rest length. If cutoff is None, the potential will not be trancated.
    
    force_group: int
        Force group.
    
    Return
    ------
    excluded_force: Force
        OpenMM Force object
    """    
    energyfunc = '''eps*(sig/r)^12;
                    eps=epsilon(atom_type1,atom_type2);
                    sig=sigma(atom_type1,atom_type2);'''
    num_type_atom = np.shape(epsilon_map)[0]
    if chain_index != None:
        energyfunc = '(delta(chidx1-chidx2)+step(-(chidx1*chidx2)))*' +energyfunc
    excluded_force = mm.CustomNonbondedForce(energyfunc)#;ex_resi=(step(abs(ex_resi1-ex_resi2)-3))') 
    excluded_force.addTabulatedFunction('epsilon', mm.Discrete2DFunction(num_type_atom,num_type_atom,epsilon_map.flatten()))
    excluded_force.addTabulatedFunction('sigma', mm.Discrete2DFunction(num_type_atom,num_type_atom,sigma_map.flatten()))
    excluded_force.addPerParticleParameter('atom_type')
    if chain_index != None:
        excluded_force.addPerParticleParameter('chidx')
        for i,ai in enumerate(atom_types):
            excluded_force.addParticle([ai,chain_index[i]])
    else:
        for ai in atom_types:
            excluded_force.addParticle([ai])
    for i_exclu in exclusions:
        excluded_force.addExclusion(i_exclu[0],i_exclu[1])
    if extra_exclusions != None:
        for i_exclu in extra_exclusions:
            if ([i_exclu[0],i_exclu[1]] == np.array(exclusions)).all(1).any() or ([i_exclu[1],i_exclu[0]] == np.array(exclusions)).all(1).any():
                continue
            excluded_force.addExclusion(i_exclu[0],i_exclu[1])
    if use_pbc:
        excluded_force.setNonbondedMethod(excluded_force.CutoffPeriodic)
        excluded_force.setCutoffDistance(cutoff*np.max(sigma_map))
    elif cutoff == None:
        excluded_force.setNonbondedMethod(excluded_force.NoCutoff)
    else:
        excluded_force.setNonbondedMethod(excluded_force.CutoffNonPeriodic)
        excluded_force.setCutoffDistance(cutoff*np.max(sigma_map))
    excluded_force.setForceGroup(force_group)
    return excluded_force

def LJ_12_6_term(atom_types,epsilon_map,sigma_map,exclusions,extra_exclusions=None,chain_index=None,use_pbc=False,cutoff=2.5,force_group=1):
    """
    The Lennard-Jones 12-6 potential.

    Parameter
    --------
    atom_type: list, int 
        Atom type.
    
    epsilon_map: array
        the strength of specific residue pair

    sigma_map: array
        the rest length of specific residue pair

    chain_index: list
        the chain index of any atom.
    
    exclusions: list
        particular pairs of particles whose interactions should be omitted from force and energy calculations.
    
    extra_exclusions: list
        This represents the extra particle pairs that are neglected in the calculation of the potential.

    use_pbc: bool, optional
        Whether use periodic boundary conditions.  If False (default),then pbc would not apply to force.

    cutoff: float
        cutoff = truncated distance / rest length. If cutoff is None, the potential will not be trancated.
    
    force_group: int
        Force group.
    
    Return
    ------
    excluded_force: Force
        OpenMM Force object    
    """

    num_type_atom = np.shape(epsilon_map)[0]
    energyfunc = '''4*eps*((sig/r)^12-(sig/6)^6);
                    eps=epsilon(atomtype1,atomtype2);
                    sig=sigma(atomtype1,atomtype2);'''

    if chain_index != None:
        energyfunc = '(1-delta(chidx1-chidx2))*' +energyfunc      
    lj12_6_force = mm.CustomNonbondedForce(energyfunc)
    lj12_6_force.addTabulatedFunction('epsilon', mm.Discrete2DFunction(num_type_atom,num_type_atom,epsilon_map.flatten()))
    lj12_6_force.addTabulatedFunction('sigma', mm.Discrete2DFunction(num_type_atom,num_type_atom,sigma_map.flatten()))
    lj12_6_force.addPerParticleParameter('atomtype')
    if chain_index != None:
        lj12_6_force.addPerParticleParameter('chidx')
        for i,ai in enumerate(atom_types):
            lj12_6_force.addParticle([ai,chain_index[i]])
    else:      
        for ai in atom_types:
            lj12_6_force.addParticle([ai])
    for i_exclu in exclusions:
        lj12_6_force.addExclusion(i_exclu[0],i_exclu[1])
    if extra_exclusions != None:
        for i_exclu in extra_exclusions:
            if ([i_exclu[0],i_exclu[1]] == np.array(exclusions)).all(1).any() or ([i_exclu[1],i_exclu[0]] == np.array(exclusions)).all(1).any():
                continue
            lj12_6_force.addExclusion(i_exclu[0],i_exclu[1])
    if use_pbc:
        lj12_6_force.setNonbondedMethod(lj12_6_force.CutoffPeriodic)
        lj12_6_force.setCutoffDistance(cutoff*np.max(sigma_map))
    elif cutoff == None:
        lj12_6_force.setNonbondedMethod(lj12_6_force.NoCutoff)
    else:
        lj12_6_force.setNonbondedMethod(lj12_6_force.CutoffNonPeriodic)
        lj12_6_force.setCutoffDistance(cutoff*np.max(sigma_map))
    lj12_6_force.setForceGroup(force_group)
    return lj12_6_force
    
    
