#!/usr/bin/env python
import mdanalysis as mda
from mdanalysis.analysis import align
from mdanalysis.analysis.rms import RMSD
import numpy as np
import argparse
import os

"""
Sample Usage:
./cluster_traj.py --traj 1btl_production_npt_gpu.nc --top 1btl.parm7 --rmsd 2.5 --frames 1000
"""

def cluster_trajectory(traj_file, top_file, rmsd_cutoff, n_frames='all', output_dir='clusters'):
    """
    Cluster trajectory frames based on RMSD.
    
    Parameters:
    -----------
    traj_file : str
        Path to trajectory file (.nc, .dcd, etc.)
    top_file : str
        Path to topology file (.parm7, .psf, etc.)
    rmsd_cutoff : float
        RMSD cutoff in Angstroms for clustering
    n_frames : str or int
        Number of frames from the end to analyze ('all' for entire trajectory)
    output_dir : str
        Directory to save output files
    """
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load trajectory
    u = mda.Universe(top_file, traj_file)
    
    # Select frames to analyze
    if n_frames == 'all':
        frames = range(len(u.trajectory))
    else:
        frames = range(max(0, len(u.trajectory) - int(n_frames)), len(u.trajectory))
    
    # Prepare protein selection (backbone atoms)
    protein = u.select_atoms('protein and backbone')
    
    # Initialize arrays for clustering
    frame_coords = []
    frame_indices = []
    
    # Align trajectory to first frame
    ref_coords = protein.positions
    for ts in u.trajectory[frames]:
        # Align current frame to reference
        _, rmsd = align.rotation_matrix(protein.positions, ref_coords)
        frame_coords.append(protein.positions.copy())
        frame_indices.append(ts.frame)
    
    # Convert to numpy array
    frame_coords = np.array(frame_coords)
    
    # Perform clustering
    clusters = []
    used_frames = set()
    
    # Stats file
    with open(f'{output_dir}/cluster_statistics.txt', 'w') as stats:
        stats.write(f"RMSD cutoff: {rmsd_cutoff} Å\n\n")
        
        # Clustering algorithm
        cluster_id = 0
        while len(used_frames) < len(frame_coords):
            current_cluster = []
            
            # Find first unused frame
            for i in range(len(frame_coords)):
                if i not in used_frames:
                    reference = frame_coords[i]
                    current_cluster.append(i)
                    used_frames.add(i)
                    break
            
            # Find all frames within RMSD cutoff
            for i in range(len(frame_coords)):
                if i not in used_frames:
                    rmsd = np.sqrt(np.mean(np.sum((frame_coords[i] - reference)**2, axis=1)))
                    if rmsd < rmsd_cutoff:
                        current_cluster.append(i)
                        used_frames.add(i)
            
            clusters.append(current_cluster)
            
            # Write cluster statistics
            stats.write(f"Cluster {cluster_id}:\n")
            stats.write(f"  Size: {len(current_cluster)} frames\n")
            stats.write(f"  Representative frame: {frame_indices[current_cluster[0]]}\n")
            stats.write(f"  Frames: {[frame_indices[i] for i in current_cluster]}\n\n")
            
            # Save representative structure
            u.trajectory[frame_indices[current_cluster[0]]]
            protein.write(f'{output_dir}/cluster_{cluster_id}_rep.pdb')
            
            cluster_id += 1
        
        # Write summary
        stats.write(f"\nTotal number of clusters: {len(clusters)}\n")
        stats.write(f"Largest cluster size: {max(len(c) for c in clusters)} frames\n")
        stats.write(f"Smallest cluster size: {min(len(c) for c in clusters)} frames\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cluster trajectory frames based on RMSD')
    parser.add_argument('--traj', required=True, help='Trajectory file')
    parser.add_argument('--top', required=True, help='Topology file')
    parser.add_argument('--rmsd', type=float, required=True, help='RMSD cutoff (Å)')
    parser.add_argument('--frames', default='all', help='Number of frames from end to analyze')
    parser.add_argument('--out', default='clusters', help='Output directory')
    
    args = parser.parse_args()
    cluster_trajectory(args.traj, args.top, args.rmsd, args.frames, args.out)
