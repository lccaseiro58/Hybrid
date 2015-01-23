##!/bin/bash 
#PBS -l walltime=03:00 
#PBS -j oe 
#PBS -N mpiP 
#PBS -l nodes=4:ppn=4
#PBS -q mei
 
cat $PBS_NODEFILE 
cd ~/PCP/Hybrid/MPIP 

module add gnu/4.9.0
module add gnu/openmpi_eth/1.8.4
for proc in 4 9 16
do
mpirun -np $proc -machinefile $PBS_NODEFILE -mca bt1 ^openib ./par_mpiP

done