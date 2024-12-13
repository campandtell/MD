#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import sys
import glob
from pathlib import Path

def get_com_velocities(atoms_obj, chunk_start, chunk_end):
    """Calculate center of mass velocities for given atom selection."""
    masses = atoms_obj.atoms.masses
    com_vels = []
    for _ in atoms_obj.universe.trajectory[chunk_start:chunk_end]:
        vel = atoms_obj.velocities
        com_vels.append(np.matmul(np.transpose(vel), masses) / np.sum(masses))
    return com_vels

def get_adjusted_vels_per_atom(atom_obj, com_vels, chunk_start, chunk_end):
    """Get COM-adjusted velocities for a single atom."""
    adj_vels = []
    for _ in atom_obj.universe.trajectory[chunk_start:chunk_end]:
        adj_vels.append(atom_obj.velocities[0])
    adj_vels = np.array(adj_vels)
    adj_vels = adj_vels - com_vels
    return adj_vels

def get_correlation(Nf, Nw, dN, atom_obj, com_vels, chunk_start, chunk_end, cov_save):
    """Calculate velocity autocorrelation for a given atom."""
    adj_vels = get_adjusted_vels_per_atom(atom_obj, com_vels, chunk_start, chunk_end)
    cov = []
    chunk_size = chunk_end - chunk_start
    
    a = 0
    while a <= (chunk_size - Nw):
        init_vel = adj_vels[a]
        vel_chunk = adj_vels[a:a+Nw]
        if len(vel_chunk) != Nw:
            print(a)
        cov.append(np.dot(vel_chunk, init_vel) / np.square(np.linalg.norm(init_vel)))
        a += dN
    
    avg_cov = np.mean(cov, axis=0)
    np.save(file=f"{cov_save}_chunk_{chunk_start}_{chunk_end}.npy", arr=avg_cov)
    return

def make_chunks(Nf, n_chunks):
    """Divide trajectory into chunks."""
    chunk_i = list(range(0, Nf, int(Nf/n_chunks)))
    chunk_i.append(Nf)
    return chunk_i

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python vel_autocorr.py <pdb_file> <trajectory>")
        sys.exit(1)

    pdb_struct = sys.argv[1]
    traj = sys.argv[2]
    
    # Load with PDB format
    md1 = mda.Universe(pdb_struct, traj)
    
    # Parameters
    Nw = 500  # window size
    dN = 1    # stride
    Nf = len(md1.trajectory)
    
    # Select alpha carbons
    ca_obj = md1.select_atoms("name CA")
    resi_list = ca_obj.residues.resids
    N_res = len(resi_list)
    
    # Divide trajectory into 25 chunks
    n_chunks = 25
    chunk_i = make_chunks(Nf, n_chunks)
    
    # Check existing calculations
    get_finished = glob.glob("*.npy")
    
    # Process each chunk
    for i in range(0, len(chunk_i)-1):
        chunk_start = chunk_i[i]
        chunk_end = chunk_i[i+1]
        
        # Check which residues need processing
        resi_to_run = []
        for i in resi_list:
            cov_save = f"velcorr_{i}"
            name = f"{cov_save}_chunk_{chunk_start}_{chunk_end}.npy"
            if name in get_finished:
                continue
            else:
                resi_to_run.append(i)
        
        if len(resi_to_run) == 0:
            continue
            
        # Calculate COM velocities once per chunk
        com_vels = get_com_velocities(ca_obj, chunk_start, chunk_end)
        
        # Process each residue
        for i in resi_to_run:
            cov_save = f"velcorr_{i}"
            name = f"{cov_save}_chunk_{chunk_start}_{chunk_end}.npy"
            atom_obj = md1.select_atoms(f"name CA and resid {i}")
            get_correlation(Nf, Nw, dN, atom_obj, com_vels, chunk_start, chunk_end, cov_save)
