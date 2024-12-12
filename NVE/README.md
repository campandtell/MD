# Molecular Dynamics NVE Simulation Protocol
A comprehensive guide for running NVE (microcanonical ensemble) molecular dynamics simulations using AMBER on HPC clusters.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Initial Setup](#initial-setup)
- [System Minimization](#system-minimization)
- [System Heating](#system-heating)
- [NVE Equilibration](#nve-equilibration)
- [Production Runs](#production-runs)
- [Analysis Tools](#analysis-tools)
- [Troubleshooting](#troubleshooting)

## Prerequisites
- AMBER (version 22 or later)
- Access to HPC cluster with GPU capabilities
- Initial PDB structure
- Python with MDAnalysis (for analysis)

## Initial Setup
Create your working directory and copy your PDB file:
```bash
cd /scratch/YOUR_USERNAME
mkdir -pv nve_simulation
cd nve_simulation
cp /path/to/your/protein.pdb ./1efa.pdb
```

### System Preparation (tleap)
Create `tleap.in`:
```bash
source leaprc.protein.ff14SB
source leaprc.water.tip3p
pdb_file = loadpdb 1efa.pdb
solvateBox pdb_file TIP3PBOX 16.0
charge pdb_file
addions pdb_file Na+ 0
addions pdb_file Cl- 0
saveAmberParm pdb_file 1efa.parm7 1efa.crd
quit
```

Run tleap (make sure to load your HPC's version of amber22 or higher):
```bash
module load amber/22v3
tleap -f tleap.in
```

## System Minimization

### Step 1: Solvent Minimization
Create `minimization_solvent.in`:
```bash
# Energy minimization solvent
&cntrl
imin=1,
ntpr=10000,
ntr=1,
restraintmask='!(:WAT,Na+,Cl-)',
restraint_wt=10.0,
maxcyc=50000,
ncyc=25000,
ntc=1,
ntf=1,
ntb=1,
cut=12.0,
&end
/
```

Submit solvent minimization:
```bash
sbatch minimization_solvent_sbatch
```

### Step 2: Full System Minimization
Create `minimization_solution.in`:
```bash
# Energy minimization solution
&cntrl
imin=1,
ntpr=10000,
maxcyc=100000,
ncyc=50000,
ntc=1,
ntf=1,
ntb=1,
cut=12.0,
&end
/
```

Submit full minimization:
```bash
sbatch minimization_solution_sbatch
```

## System Heating
Create `heatup_nve.in`:
```bash
# NVE-preparatory heating
&cntrl
imin=0,
ntx=1,
irest=0,
ntpr=500,
ntwr=5000,
iwrap=1,
ntwx=500,
ntwv=500,
ntr=1,
restraintmask='!(:WAT,Na+,Cl-)',
restraint_wt=10.0,
nstlim=100000,
dt=0.002,
ntt=3,
temp0=300.0,
tempi=0.0,
ig=-1,
gamma_ln=1.0,
ntc=2,
ntf=2,
ntb=1,
cut=12.0,
ioutfm=1,
&end
/
```

Key differences for NVE preparation:
- Constant volume heating (ntb=1)
- Longer duration (200ps)
- Velocity output enabled (ntwv=500)
- Reduced Langevin coupling (gamma_ln=1.0)
- More frequent output for monitoring

## NVE Equilibration
Create `equilibration_nve.in`:
```bash
# Initial NVE equilibration
&cntrl
imin=0,
ntx=5,
irest=1,
ntpr=500,
ntwr=5000,
iwrap=1,
ntwx=500,
ntwv=500,
nstlim=50000,
dt=0.002,
ntt=0,
temp0=300.0,
ig=-1,
ntc=2,
ntf=2,
ntb=1,
cut=12.0,
ioutfm=1,
&end
/
```

Important considerations:
- No temperature coupling (ntt=0)
- Constant volume (ntb=1)
- Frequent output for monitoring stability
- Monitor total energy conservation

## Production Runs

### CPU Production
Create `production_nve_cpu.in`:
```bash
# CPU NVE production
&cntrl
imin=0,
ntx=5,
irest=1,
ntpr=5000,
ntwr=5000,
iwrap=1,
ntwx=5000,
ntwv=5000,
nstlim=50000,
dt=0.002,
ntt=0,
temp0=300.0,
ig=-1,
ntc=2,
ntf=2,
ntb=1,
cut=12.0,
ioutfm=1,
&end
/
```

### GPU Production
Create `production_nve_gpu.in`:
```bash
# GPU NVE production - 20ns
&cntrl
imin=0,
ntx=5,
irest=1,
ntpr=5000,
ntwr=5000,
iwrap=1,
ntwx=5000,
ntwv=5000,
nstlim=10000000,
dt=0.002,
ntt=0,
temp0=300.0,
ig=-1,
ntc=2,
ntf=2,
ntb=1,
cut=12.0,
ioutfm=1,
&end
/
```

## Analysis Tools

### Energy Conservation Check
Create an energy analysis script `check_energy.cpptraj`:
```bash
parm 1efa.parm7
trajin 1efa_production_nve_gpu.nc
energy output energy.dat
run
```

Run analysis:
```bash
cpptraj -i check_energy.cpptraj
```

### Trajectory Clustering
Save the provided clustering script as `cluster_traj.py` and run:
```bash
./cluster_traj.py --traj 1efa_production_nve_gpu.nc --top 1efa.parm7 --rmsd 2.5 --frames 1000
```

## Troubleshooting

### Common Issues and Solutions

1. **Energy Conservation Problems**
   - Check for appropriate timestep (dt)
   - Verify minimization was complete
   - Ensure heating was stable
   - Consider longer equilibration

2. **Temperature Drift**
   - Extend heating phase
   - Reduce Langevin coupling during heating
   - Check for bad contacts in structure

3. **Simulation Crashes**
   - Verify box size is appropriate
   - Check for bad contacts
   - Ensure proper neutralization
   - Review restart files for anomalies

### Quality Control Metrics

Monitor these values throughout your simulation:
- Total energy conservation (should be < 0.01% drift)
- Temperature stability (should be ±5K around target)
- System density (should be ~1.0 g/cm³)
- RMSD from starting structure
- Velocity distribution (should be Maxwell-Boltzmann)

## Notes
- Always check energy conservation before proceeding to production
- Save velocities for complete phase space analysis
- Consider running multiple replicates
- Monitor system properties regularly during production

For additional help, consult the AMBER manual or contact your HPC support team.
