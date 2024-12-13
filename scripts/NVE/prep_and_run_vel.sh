#!/bin/bash
#SBATCH -n 8
#SBATCH -p general
#SBATCH -t 12:00:00
#SBATCH -o velcorr_%j.out
#SBATCH -e velcorr_%j.err

module purge all
module load amber/22v3
module load anaconda3/5.3.0

# Get command line arguments
parm7=$1  # Path to parm7 file
traj=$2   # Path to trajectory

# Create working directory if it doesn't exist
mkdir -p vel_analysis
cd vel_analysis

# Extract first frame and convert to PDB
cpptraj $parm7 << EOF
    trajin $traj 0 1
    trajout struct_model.pdb pdb
EOF

# Run velocity autocorrelation calculation
python vel_autocorr.py struct_model.pdb $traj
