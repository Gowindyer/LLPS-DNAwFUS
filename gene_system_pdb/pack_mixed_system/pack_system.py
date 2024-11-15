import json
import sys
import os
import numpy as np
import pandas as pd

# This JSON file loads up the user input parameters about packing system.
configdic = json.load(open('../../simulation_parameters.json'))
# This paramters is the executed file of packmol
path_to_packmol = configdic['execute_packmol_path']#'/home/windyer/software/packmol/packmol'
# parameters about the name of output packmol.inp file 
name_pack_inp = configdic['packmol_inp']
ex_force_para = pd.read_csv('../../params/ex_sigma.csv').to_numpy()[:,1:]
lj_para =pd.read_csv('../../params/lj_sigma.csv').to_numpy()[:,1:]
# The seprated tolerance between intermelecular atoms
all_sigma = np.stack((ex_force_para[:],lj_para[:]))
tolerance_pack = 2.5*np.max(all_sigma)
# The pdb name of mixed system
output_pdb_name = configdic['output_name_sys']
# Path of type 1 structure 
path_type_1_pdb = configdic['path_type1_pdb']
# Path of type 2 structure 
path_type_2_pdb = configdic['path_type2_pdb']
# Parameters that determine the number of type 1 and type 2, box size
num_type_1 = configdic['num_type_1']
num_type_2 = configdic['num_type_2']
box_size = configdic['box_size']
with open(name_pack_inp,'w+') as pack_sys:
    pack_sys.write('tolerance %.2f\n'%tolerance_pack)
    pack_sys.write('add_box_sides = 1.0\n')
    pack_sys.write('filetype pdb\n')
    pack_sys.write('output %s.pdb\n'%output_pdb_name)
    
    pack_sys.write('\n')
    pack_sys.write('# Here is Randomly place a certain number of type1 molecules into a given region\n')
    pack_sys.write('structure %s\n'%path_type_1_pdb)
    pack_sys.write('   number %d\n'%num_type_1)
    pack_sys.write('   changechains\n')
    pack_sys.write('   inside box %.4f %.4f %.4f %.4f %.4f %.4f\n'%\
            (box_size[0],box_size[1],box_size[2],box_size[3],\
                box_size[4],box_size[5]))
    pack_sys.write('end structure\n')


    pack_sys.write('\n')
    pack_sys.write('# Here is Randomly place a certain number of type2 molecules into a given region\n')
    pack_sys.write('structure %s\n'%path_type_2_pdb)
    pack_sys.write('   changechains\n')
    pack_sys.write('   number %d\n'%num_type_2)
    pack_sys.write('   inside box %.4f %.4f %.4f %.4f %.4f %.4f\n'%\
            (box_size[0],box_size[1],box_size[2],box_size[3],\
                box_size[4],box_size[5]))
    pack_sys.write('end structure\n')

# generate the pdb of mixed system
os.system('%s < %s'%(path_to_packmol,name_pack_inp))
