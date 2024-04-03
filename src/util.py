import mdtraj as md
import sys
import numpy as np

def initialize(dcd, pdb, residue_names):
    """Initializes parameters for networkx functions.
    Parameters
    ----------
    dcd : 
        .dcd trajectory file
    pdb : 
        .pdb topology file
    residue_names : 
        species name to match in pdb file for inter-species analysis
    
    Returns   
    -------
    xyz : 
        ...
    traj : 
        ...
    cell_map : 
        ...
    residues : 
        ...
    atoms : 
        ...
    num_frames : 
        ...
    """
    traj = md.load(pdb, dcd)
    residues, atoms = get_atoms(traj, residue_names)
    xyz = traj.xyz[:, atoms, :]
    #print(len(xyz[0]))
    #print(len(xyz[0])/len(residues))
    cell_map = None
    num_frames = None
    return xyz, traj, cell_map, residues, atoms, num_frames
def get_atoms(traj, residue_name):
    """Filter out list of indices for residue and respective atoms
    Parameters
    ----------
    traj : mdtraj.core.trajectory.Trajectory 
        MDTraj trajectory object
    residue_name : str
        name of residue of interest matching pdb topology file 
    Returns
    -------
    residues : list  
        list of residue indices for selected residue name
    atoms : list 
        list of atom indices for selected atom name 
    
    """
    residues, atoms = [], []
    for res in traj.top.residues: 
        if res.name == residue_name: 
            residues.append(res.index)
            for atom in res.atoms: 
                atoms.append(atom.index)
    #print(f"({len(atoms)/len(residues)}) atoms in residue")
    return residues, atoms
def make_head():
    pass 
def set_cell():
    pass
def make_map():
    pass 

def get_distance():
    pass
def get_angle():
    pass
def hbond_criteria():
    pass
def check_xyz():
    pass

