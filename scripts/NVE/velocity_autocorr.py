#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import sys
import glob
from pathlib import Path
import argparse

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

def get_correlation(Nw, dN, atom_obj, com_vels, chunk_start, chunk_end, cov_save):
    """Calculate velocity autocorrelation for a given atom."""
    adj_vels = get_adjusted_vels_per_atom(atom_obj, com_vels, chunk_start, chunk_end)
    cov = []
    chunk_size = chunk_end - chunk_start
    
    a = 0
    while a <= (chunk_size - Nw):
        init_vel = adj_vels[a]
        vel_chunk = adj_vels[a:a+Nw]
        if len(vel_chunk) != Nw:
            print(f"Warning: Inconsistent chunk size at index {a}")
            continue
        cov.append(np.dot(vel_chunk, init_vel) / np.square(np.linalg.norm(init_vel)))
        a += dN
    
    avg_cov = np.mean(cov, axis=0)
    np.save(f"{cov_save}_chunk_{chunk_start}_{chunk_end}.npy", avg_cov)

def make_chunks(Nf, n_chunks):
    """Divide trajectory into chunks."""
    chunk_i = list(range(0, Nf, int(Nf/n_chunks)))
    chunk_i.append(Nf)
    return chunk_i

def main():
    parser = argparse.ArgumentParser(description='Calculate velocity autocorrelation functions')
    parser.add_argument('topology', help='Topology file (parm7)')
    parser.add_argument('trajectory', help='Trajectory file (nc)')
    parser.add_argument('--window', type=int, default=500, help='Window size for correlation')
    parser.add_argument('--stride', type=int, default=1, help='Stride for correlation calculation')
    parser.add_argument('--chunks', type=int, default=25, help='Number of chunks to divide trajectory')
    args = parser.parse_args()

    # Load trajectory
    u = mda.Universe(args.topology, args.trajectory)
    
    # Setup calculation parameters
    Nw = args.window
    dN = args.stride
    Nf = len(u.trajectory)
    
    # Select alpha carbons
    ca_obj = u.select_atoms("name CA")
    resi_list = ca_obj.residues.resids
    
    # Create output directory
    output_dir = Path('velcorr_output')
    output_dir.mkdir(exist_ok=True)
    
    # Get chunks
    chunk_i = make_chunks(Nf, args.chunks)
    
    # Check existing calculations
    get_finished = list(output_dir.glob("*.npy"))
    
    # Process each chunk
    for i in range(len(chunk_i)-1):
        chunk_start = chunk_i[i]
        chunk_end = chunk_i[i+1]
        
        # Check which residues need processing
        resi_to_run = []
        for resid in resi_list:
            cov_save = output_dir / f"velcorr_{resid}"
            name = f"{cov_save}_chunk_{chunk_start}_{chunk_end}.npy"
            if not Path(name).exists():
                resi_to_run.append(resid)
        
        if not resi_to_run:
            continue
            
        # Calculate COM velocities once per chunk
        com_vels = get_com_velocities(ca_obj, chunk_start, chunk_end)
        
        # Process each residue
        for resid in resi_to_run:
            cov_save = output_dir / f"velcorr_{resid}"
            atom_obj = u.select_atoms(f"name CA and resid {resid}")
            get_correlation(Nw, dN, atom_obj, com_vels, chunk_start, chunk_end, cov_save)
            print(f"Processed residue {resid} for chunk {chunk_start}-{chunk_end}")

if __name__ == "__main__":
    main()
