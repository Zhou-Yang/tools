import numpy as np
import mdtraj as md
from itertools import combinations
from itertools import product

NATIVE_CUTOFF = 0.45  # nanometers

ipdb = md.load_pdb('traj_comp_nowater_nojump.pdb')

ipdb_heavy = ipdb.topology.select_atom_indices('heavy')

ipdb_pro = np.intersect1d(ipdb.topology.select('protein'), ipdb_heavy)
ipdb_lig = ipdb.topology.select("(resname LIG) and not type H")



lig_pro_pairs = np.array([(i,j) for (i,j) in product(ipdb_pro,ipdb_lig)])

lig_pro_pairs_distance = md.compute_distances(ipdb[0], lig_pro_pairs)[0]

native_lig_pro_contact = lig_pro_pairs[lig_pro_pairs_distance < NATIVE_CUTOFF]


num=1
for i in native_lig_pro_contact:
    dis = float(md.compute_distances(ipdb[0], np.array([i]))[0])
    print "ATOMS%d=%d,%d SWITCH%d={RATIONAL R_0=1.2 D_0=0.0 NN=6 MM=10} WEIGHT%d=1" %(num,i[0],i[1],num,num)
    num += 1
