#!/usr/bin/env python
from test import * 

# Processes count
cores = 12 

# Load trj to MDTraj
traj = md.load("npt_final-28.dcd", top="npt_final-28.pdb") 

# residue and atom index list for relevant groups
atom_dict = get_atoms(traj, ["OCT"])
cation_residues, cation_atoms = atom_dict["OCT"][0], atom_dict["OCT"][1]

# index list from 0
res_index = [i for i in range(len(cation_residues))]

frames = len(traj)

# np array of cation atoms
cation_xyz = get_xyz(traj, cation_atoms)

# check for out-of-box coords and remove
cation_xyz = check_coords(cation_xyz, traj.unitcell_lengths[0,0])

# divide box into cells 
cell_map = make_map(10)

num_frames = len(traj)

edges = Parallel(n_jobs=cores, backend="multiprocessing")\
            (delayed(get_network)\
            (cation_xyz[frame], traj.unitcell_lengths[0,0], cell_map, 10, 0.4, frame)\
            for frame in range(num_frames))

#With the edges obtained from get_network, add to a NetworkX graph object
graphs = []
for edge in edges:
    graph = make_graph(edge, res_index)
    graphs.append(graph)

#Uses the formed graphs to compute graph properties
diams = Parallel(n_jobs=cores,backend="multiprocessing")(
        delayed(
            ###Change the "compute_metric" function here if you want to compute a different property
            compute_metric
            )(graphs[frame]) for frame in range(len(graphs)))

#Change this function depending on which metric you want to plot or analyze
plt_metric(diams)


