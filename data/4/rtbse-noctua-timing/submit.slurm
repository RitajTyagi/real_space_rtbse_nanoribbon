#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH --exclusive
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH -N 1
#SBATCH -J "rirs"
#SBATCH -p largemem
#SBATCH --account=hpc-prf-rqsued

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=true

ulimit -s unlimited

module reset
module load chem/CP2K/2026.2-foss-2025b-cpu
export CP2K_DATA_DIR=/pc2/groups/hpc-prf-metdyn/metdyn07_Ritaj/implementation/github/cp2k/data
srun /pc2/groups/hpc-prf-metdyn/metdyn07_Ritaj/implementation/github/cp2k-noctua/cp2k/build/bin/cp2k.psmp input.inp &> output.out

