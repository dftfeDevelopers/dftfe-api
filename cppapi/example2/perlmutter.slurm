#!/bin/bash
#SBATCH -A m3916_g
#SBATCH -C gpu
#SBATCH -q early_science
#SBATCH --job-name dftfeapitest-example2
#SBATCH -t 00:05:00
#SBATCH -n 64
#SBATCH --ntasks-per-node=4
#SBATCH -c 32
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=map_gpu:0*1,1*1,2*1,3*1

export SLURM_CPU_BIND="cores"
module load PrgEnv-gnu
module load cudatoolkit
module unload cray-libsci/21.08.1.2
module load nccl/2.11.4
export OMP_NUM_THREADS=1

srun ./dftfeapitest > output
