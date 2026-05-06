#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --partition=gpumd
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#module load gcc/8.3  cuda/11.0  
module load  gcc/9.3.0  cuda/12.2
source ~/.bashrc 
conda activate calorine

mkdir dump thermo movie mc_energyLog mc_summaryLog

# loop the graphene growth simulation
for i in {1..10}; do
  for j in {1..10}; do
    /data/home/yfzhao/GPUMD-5.0-edmkc/src/gpumd
    rm model.xyz
    python add-carbon-gpumd.py 
    cp dump.xyz ./dump/dump.$(( (i-1)*10+j )).xyz
    cp thermo.out ./thermo/thermo.$(( (i-1)*10+j )).out
    cp movie.xyz ./movie/movie.$(( (i-1)*10+j )).xyz
    rm dump.xyz thermo.out movie.xyz 
  done  
  
  python edkmc.py
  cp mc_energy.log ./mc_energyLog/mc_energy.$i.log
  cp mc_summary.log ./mc_summaryLog/mc_summary.$i.log
  rm mc_energy.log mc_summary.log
done


