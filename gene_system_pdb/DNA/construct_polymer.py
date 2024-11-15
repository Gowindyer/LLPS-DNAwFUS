#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Descripation:
    Construct the position and pdb file of single chain
"""

import openmm as mm
import openmm.app as app
from openmm import unit
import mdtraj as md
import numpy as np
import json

configdic = json.load(open('../../simulation_parameters.json'))

# bond length(angstrom)
l0 = configdic["bond_length_DNA"]
# number of bead in chain.
num_bead = configdic["num_bead_in_DNA"]
# filename of PDB
pdb_file_name = "DNA"
# name of residue
resi_name = configdic["residue_name_of_DNA"]
# name of atom 
atom_name = configdic["atom_name_of_DNA"]

# postions
positions = np.array([[0,l0*i,4]for i in range(num_bead)])
## bond
bondinfo = np.array([[i,i+1] for i in range(num_bead-1)])

# create PDB 
topology = mm.app.topology.Topology()
topology.addChain('LINEAR')
sin_chain = list(topology.chains())[0]

for i in range(0,num_bead):
    topology.addResidue(resi_name[0],sin_chain)
    sin_chain_resi = list(topology.residues())[-1]
    topology.addAtom(atom_name[0],app.Element.getBySymbol('Ti'),sin_chain_resi)

atoms = list(topology.atoms())
for i in range(len(bondinfo)):
    topology.addBond(atoms[bondinfo[i,0]],atoms[bondinfo[i,1]],type='Single',order=1)
    
positions_dimension = np.array(positions)*0.1*unit.nanometer
## output PDB
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
        if (i+1)%4==0:
           psf.write('\n')
        elif i == int(len(bond)-1):
           psf.write('\n')
