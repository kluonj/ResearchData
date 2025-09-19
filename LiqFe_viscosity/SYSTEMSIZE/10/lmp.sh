#!/bin/bash
#SBATCH --job-name=lammps-sys
#SBATCH --nodes=1
#SBATCH --ntasks=36 #8 #12 #24
#SBATCH --cpus-per-task=1
#SBATCH --time=10-01:00:00
#SBATCH --partition=medium #small #medium
#SBATCH --output=%j.log

#source ~kluo/.bashrc

module load mpi/intelmpi/2017.4.239 
module load compiler/intel/intel-compiler-2017.5.239

ulimit -s unlimited
export OMP_NUM_THREADS=1

#export LAMMPS_POTENTIALS=/path/to/lammps/potentials
#EXE=~kluo/miniconda3/envs/deepmd/bin/lmp
EXE=~/.conda/envs/deepmd/bin/lmp

mpirun -n $SLURM_NTASKS $EXE -in in.lammps
