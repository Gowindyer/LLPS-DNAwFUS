#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Descripation:
    Construct pdb file for fus segment
"""

import openmm as mm
import openmm.app as app
from openmm import unit
import mdtraj as md
import numpy as np
import json

configdic = json.load(open('../../simulation_parameters.json'))

# bond length(angstrom)
l0 = configdic["bond_length_fus"]
# the number of bead in N terminal and mediate segment and C terminal.
num_bead_N_medi_C = configdic["num_bead_N_mediate_C_fus"]
name_residue = configdic["residue_name_of_fus"]
name_atom = configdic["atom_name_of_fus"]
tot_num_bead = sum(num_bead_N_medi_C) 
## filename of PDB
pdb_file_name = "protein_fus"
## positions
positions = np.array([[0,l0*i,4] for i in range(tot_num_bead)])
## the number of bond in topology
bondinfo = np.array([[i,i+1] for i in range(tot_num_bead-1)])

## CONSTRUCT PDB 
topology = mm.app.topology.Topology()
topology.addChain('FUS')
fus_chain = list(topology.chains())[0]

for i,num_bead_i in enumerate(num_bead_N_medi_C):
    for _ in range(num_bead_i):
        topology.addResidue(name_residue[i],fus_chain)
        fus_chain_resi = list(topology.residues())[-1]
        topology.addAtom(name_atom[i],app.Element.getBySymbol('Ti'),fus_chain_resi)

atoms = list(topology.atoms())
for i in range(len(bondinfo)):
    topology.addBond(atoms[bondinfo[i,0]],atoms[bondinfo[i,1]],type='Single',order=1)
    
positions_dimension = np.array(positions)*0.1*unit.nanometer
## OUTPUT PDB
app.PDBFile.writeFile(topology,positions_dimension,open('%s.pdb'%pdb_file_name,'w'))


res = list(topology.residues())
bond = list(topology.bonds())
space = ' '
with open('%s.psf'%pdb_file_name,'w') as psf:
    psf.write('PSF\n')
    # ATOM
    #psf.write(space*6)
    psf.write('%8d !NATOM\n'%len(atoms))
    for i in range(len(atoms)):
        psf.write('%8d%7d%7s%s%4s%5s%15d%14d%8d\n'%
                  (atoms[i].index+1,1,atoms[i].residue.name,
                   space*2,atoms[i].name,atoms[i].name,
                   0,100,0))
    psf.write('\n')
    psf.write('%8d !NBOND: bonds\n'%len(bond))
    for i in range(len(bond)):
        psf.write('%8d%8d'%(bond[i][0].index+1,bond[i][1].index+1))
        if (i+1)%4 == 0:
           psf.write('\n')
        elif i == int(len(bond)-1):
           psf.write('\n')
