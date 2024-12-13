Velocity Autocorrelation Analysis Setup
Setup Instructions

Create analysis directory and copy scripts

bashCopy# Create directory for analysis
mkdir -p vel_analysis
cd vel_analysis

# Copy scripts into directory
# Save prepare_and_run.sh and vel_autocorr.py here
# Make the submission script executable
chmod +x prepare_and_run.sh
Running the Analysis
From your main simulation directory:
bashCopysbatch vel_analysis/prepare_and_run.sh 1efa.parm7 outputs/gpuProduction/1efa_production_nve_gpu.nc
Workflow
The script performs the following steps:

Creates a working directory for the analysis
Uses cpptraj to extract the first frame as PDB structure
Runs velocity autocorrelation calculations using:

Generated PDB as structure file
Original trajectory for velocity data


Outputs correlation data as .npy files for each residue and chunk

Expected Outputs
Copyvel_analysis/
├── struct_model.pdb           # Structure file from cpptraj
├── velcorr_output/           # Directory containing results
│   ├── velcorr_1_chunk_0_*.npy
│   ├── velcorr_2_chunk_0_*.npy
│   └── ...                   # One file per residue per chunk
├── velcorr_*.out             # SLURM output file
└── velcorr_*.err             # SLURM error file
