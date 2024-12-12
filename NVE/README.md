# Molecular Dynamics NVE Simulation Protocol
A comprehensive guide for running NVE (microcanonical ensemble) molecular dynamics simulations using AMBER on HPC clusters.

## Table of Contents
1. [Directory Setup](#directory-setup)
2. [Initial System Setup](#step-1-initial-system-setup)
   - [Create tleap input](#create-tleap-input)
   - [Run tleap](#run-tleap)
3. [Solvent Minimization](#step-2-solvent-minimization)
   - [Create minimization input](#create-minimization-input)
   - [Create submission script](#create-submission-script)
4. [Full System Minimization](#step-3-full-system-minimization)
   - [Create minimization input](#create-minimization-input-1)
   - [Create submission script](#create-submission-script-1)
5. [System Heating](#step-4-system-heating)
   - [Create heating input](#create-heating-input)
   - [Create submission script](#create-submission-script-2)
6. [NVE Equilibration](#step-5-nve-equilibration)
   - [Create equilibration input](#create-equilibration-input)
   - [Create submission script](#create-submission-script-3)
7. [CPU Production](#step-6-cpu-production)
   - [Create production input](#create-production-input)
   - [Create submission script](#create-submission-script-4)
8. [GPU Production](#step-7-gpu-production)
   - [Create production input](#create-production-input-1)
   - [Create submission script](#create-submission-script-5)
9. [Monitoring Jobs](#monitoring-jobs)
10. [Directory Structure](#directory-structure)

## Directory Setup
```bash
mkdir -p nve_simulation/{scripts,inputs,outputs}
cd nve_simulation
```

## Step 1: Initial System Setup

### Create tleap input
```bash
cat > inputs/tleap.in << 'EOF'
source leaprc.protein.ff14SB
source leaprc.water.tip3p
pdb_file = loadpdb 1efa.pdb
solvateBox pdb_file TIP3PBOX 16.0
charge pdb_file
addions pdb_file Na+ 0
addions pdb_file Cl- 0
saveAmberParm pdb_file 1efa.parm7 1efa.crd
quit
EOF
```

Run tleap:
```bash
module load amber/22v3
tleap -f inputs/tleap.in
```

## Step 2: Solvent Minimization

### Create minimization input
```bash
cat > inputs/minimization_solvent.in << 'EOF'
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
EOF
```

### Create submission script
```bash
cat > scripts/minimization_solvent_sbatch << 'EOF'
#!/usr/bin/env bash
#SBATCH -n 8
#SBATCH -p general
#SBATCH -t 7-00:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load amber/22v3

mpiexec.hydra -n 8 pmemd.MPI -O \
-i inputs/minimization_solvent.in \
-o outputs/minimization_solvent.out \
-p 1efa.parm7 \
-c 1efa.crd \
-r outputs/1efa_minimization_solvent.rst7 \
-x outputs/1efa_minimization_solvent.crd \
-ref 1efa.crd
EOF

chmod +x scripts/minimization_solvent_sbatch
```

Submit job:
```bash
sbatch scripts/minimization_solvent_sbatch
```

## Step 3: Full System Minimization

### Create minimization input
```bash
cat > inputs/minimization_solution.in << 'EOF'
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
EOF
```

### Create submission script
```bash
cat > scripts/minimization_solution_sbatch << 'EOF'
#!/usr/bin/env bash
#SBATCH -n 8
#SBATCH -p general
#SBATCH -t 7-00:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load amber/22v3

mpiexec.hydra -n 8 pmemd.MPI -O \
-i inputs/minimization_solution.in \
-o outputs/minimization_solution.out \
-p 1efa.parm7 \
-c outputs/1efa_minimization_solvent.rst7 \
-r outputs/1efa_minimization_solution.rst7 \
-x outputs/1efa_minimization_solution.crd
EOF

chmod +x scripts/minimization_solution_sbatch
```

Submit job:
```bash
sbatch scripts/minimization_solution_sbatch
```

## Step 4: System Heating

### Create heating input
```bash
cat > inputs/heatup_nve.in << 'EOF'
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
EOF
```

### Create submission script
```bash
cat > scripts/heatup_nve_sbatch << 'EOF'
#!/usr/bin/env bash
#SBATCH -n 8
#SBATCH -p general
#SBATCH -t 7-00:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load amber/22v3

mpiexec.hydra -n 8 pmemd.MPI -O \
-i inputs/heatup_nve.in \
-o outputs/heatup_nve.out \
-p 1efa.parm7 \
-c outputs/1efa_minimization_solution.rst7 \
-r outputs/1efa_heatup_nve.rst7 \
-x outputs/1efa_heatup_nve.nc \
-v outputs/1efa_heatup_nve.vel \
-ref outputs/1efa_minimization_solution.rst7
EOF

chmod +x scripts/heatup_nve_sbatch
```

Submit job:
```bash
sbatch scripts/heatup_nve_sbatch
```

## Step 5: NVE Equilibration

### Create equilibration input
```bash
cat > inputs/equilibration_nve.in << 'EOF'
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
EOF
```

### Create submission script
```bash
cat > scripts/equilibration_nve_sbatch << 'EOF'
#!/usr/bin/env bash
#SBATCH -n 8
#SBATCH -p general
#SBATCH -t 7-00:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load amber/22v3

mpiexec.hydra -n 8 pmemd.MPI -O \
-i inputs/equilibration_nve.in \
-o outputs/equilibration_nve.out \
-p 1efa.parm7 \
-c outputs/1efa_heatup_nve.rst7 \
-r outputs/1efa_equilibration_nve.rst7 \
-x outputs/1efa_equilibration_nve.nc \
-v outputs/1efa_equilibration_nve.vel
EOF

chmod +x scripts/equilibration_nve_sbatch
```

Submit job:
```bash
sbatch scripts/equilibration_nve_sbatch
```

## Step 6: CPU Production

### Create production input
```bash
cat > inputs/production_nve_cpu.in << 'EOF'
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
EOF
```

### Create submission script
```bash
cat > scripts/production_nve_cpu_sbatch << 'EOF'
#!/usr/bin/env bash
#SBATCH -n 8
#SBATCH -p general
#SBATCH -t 7-00:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load amber/22v3

mpiexec.hydra -n 8 pmemd.MPI -O \
-i inputs/production_nve_cpu.in \
-o outputs/production_nve_cpu.out \
-p 1efa.parm7 \
-c outputs/1efa_equilibration_nve.rst7 \
-r outputs/1efa_production_nve_cpu.rst7 \
-x outputs/1efa_production_nve_cpu.nc \
-v outputs/1efa_production_nve_cpu.vel
EOF

chmod +x scripts/production_nve_cpu_sbatch
```

Submit job:
```bash
sbatch scripts/production_nve_cpu_sbatch
```

## Step 7: GPU Production

### Create production input
```bash
cat > inputs/production_nve_gpu.in << 'EOF'
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
EOF
```

### Create submission script
```bash
cat > scripts/production_nve_gpu_sbatch << 'EOF'
#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -G 1
#SBATCH -p general
#SBATCH -t 1-00:00:00
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

module load amber/22v3

pmemd.cuda -O \
-i inputs/production_nve_gpu.in \
-o outputs/production_nve_gpu.out \
-p 1efa.parm7 \
-c outputs/1efa_production_nve_cpu.rst7 \
-r outputs/1efa_production_nve_gpu.rst7 \
-x outputs/1efa_production_nve_gpu.nc \
-v outputs/1efa_production_nve_gpu.vel
EOF

chmod +x scripts/production_nve_gpu_sbatch
```

Submit job:
```bash
sbatch scripts/production_nve_gpu_sbatch
```

## Monitoring Jobs
Check job status:
```bash
squeue -u $USER
```

## Directory Structure
After setup, your directory should look like:
```
nve_simulation/
├── 1efa.pdb
├── 1efa.parm7
├── 1efa.crd
├── inputs/
│   ├── tleap.in
│   ├── minimization_solvent.in
│   ├── minimization_solution.in
│   ├── heatup_nve.in
│   ├── equilibration_nve.in
│   ├── production_nve_cpu.in
│   └── production_nve_gpu.in
├── outputs/
│   └── (simulation outputs will appear here)
└── scripts/
    ├── minimization_solvent_sbatch
    ├── minimization_solution_sbatch
    ├── heatup_nve_sbatch
    ├── equilibration_nve_sbatch
    ├── production_nve_cpu_sbatch
    └── production_nve_gpu_sbatch
```
