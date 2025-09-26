#!/bin/bash -l
#SBATCH --parsable
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --gres=gpu:0
#SBATCH --partition amd_a8_384

source ~/.bashrc
conda activate ~/soft/deepmd-kit/


export OMP_NUM_THREADS=4
#export TF_INTRA_OP_PARALLELISM_THREADS=8
#export TF_INTER_OP_PARALLELISM_THREADS=2


dp train input.json 

bash freeze.sh
