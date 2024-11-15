#!/bin/bash

path="$PWD"
# the json file that contains all Parameters about constructing the system
cd gene_system_pdb/DNA 
# Get the pdb and psf of single chain 
python construct_polymer.py
cd $path
cd gene_system_pdb/protein_fus 
# Get the pdb and psf of star-shaped chain 
python construct_polymer.py

cd $path
cd gene_system_pdb/pack_mixed_system/ 
# Mix the two type molecular to gain the system for simulation
python pack_system.py 
