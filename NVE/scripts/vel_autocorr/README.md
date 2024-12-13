# Velocity Autocorrelation Analysis Setup

## Setup Instructions

1. Create analysis directory and copy scripts
```bash
# Create directory for analysis
mkdir -p vel_analysis
cd vel_analysis

# Copy scripts into directory
# Save prepare_and_run.sh and vel_autocorr.py here
# Make the submission script executable
chmod +x prepare_and_run.sh
```

## Running the Analysis

From your main simulation directory:
```bash
sbatch vel_analysis/prepare_and_run.sh 1efa.parm7 outputs/gpuProduction/1efa_production_nve_gpu.nc
```

## Workflow
The script performs the following steps:
1. Creates a working directory for the analysis
2. Uses cpptraj to extract the first frame as PDB structure
3. Runs velocity autocorrelation calculations using:
   - Generated PDB as structure file
   - Original trajectory for velocity data
4. Outputs correlation data as .npy files for each residue and chunk

## Expected Outputs
```
vel_analysis/
├── struct_model.pdb           # Structure file from cpptraj
├── velcorr_output/           # Directory containing results
│   ├── velcorr_1_chunk_0_*.npy
│   ├── velcorr_2_chunk_0_*.npy
│   └── ...                   # One file per residue per chunk
├── velcorr_*.out             # SLURM output file
└── velcorr_*.err             # SLURM error file
```

## Understanding Chunk Outputs

### File Naming Convention
Files are named using the pattern:
```
velcorr_[RESIDUE_NUMBER]_chunk_[START_FRAME]_[END_FRAME].npy
```

### Example Output Structure
For a trajectory with:
- 50,000 frames
- 25 chunks
- 100 residues

You would get files like:
```
velcorr_1_chunk_0_2000.npy      # Residue 1, frames 0-2000
velcorr_1_chunk_2000_4000.npy   # Residue 1, frames 2000-4000
velcorr_2_chunk_0_2000.npy      # Residue 2, frames 0-2000
velcorr_2_chunk_2000_4000.npy   # Residue 2, frames 2000-4000
...
```

### Content of Output Files
Each .npy file contains:
- A numpy array of length 500 (the window size)
- Velocity autocorrelation function values for that residue in that chunk
- Data showing correlation of motion over time within the chunk

### Purpose of Chunking
The trajectory is divided into chunks for:
1. Memory management - processes smaller pieces at a time
2. Restart capability - can resume if calculation fails
3. Parallel processing potential
4. Statistical sampling across different time periods

## Understanding .npy File Contents

### File Structure
Each .npy file contains a 1-dimensional numpy array with 500 elements (default window size). You can load it using:
```python
import numpy as np
data = np.load('velcorr_1_chunk_0_2000.npy')
```

### Data Interpretation
- Array index represents time lag (τ)
- Array values are normalized autocorrelation values in range [-1, 1]
  - 1.0: Perfect correlation
  - 0.0: No correlation
  - -1.0: Perfect anti-correlation

### Example Data Structure
```python
array([1.000000,    # τ = 0, always 1.0 (perfect self-correlation)
       0.982340,    # τ = 1, correlation after 1 timestep
       0.937219,    # τ = 2
       ...,
       0.023145,    # τ = 498
       0.021678])   # τ = 499, correlation at maximum lag
```

### Physical Interpretation
- Fast decay to zero: Rapid loss of velocity correlation (more random motion)
- Slow decay: Velocity remains correlated for longer (more directed motion)
- Negative values: Motion reversal
- Oscillations: Periodic motion

### Reading and Processing Example
```python
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.load('velcorr_1_chunk_0_2000.npy')

# Create time axis (assuming 2fs timestep)
time = np.arange(len(data)) * 2  # time in fs

# Plot
plt.plot(time, data)
plt.xlabel('Time (fs)')
plt.ylabel('Velocity Autocorrelation')
plt.title('Residue 1 Velocity Autocorrelation')
plt.show()
```
