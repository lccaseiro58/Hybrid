#!/bin/bash
#PBS -l walltime=03:00
#PBS -j oe
#PBS -N omp
#PBS -q mei

cat $PBS_NODEFILE
cd ~/Hybrid/

module add gnu/openmpi_eth/1.8.4
module add gnu/4.9.0


for proc in 4 9 16
do
	for th in 2 4 8 9 12 16
	do
	echo "PROC=$proc TH=$th"
	export OMP_NUM_THREADS=$th
	mpirun -np $proc -machinefile $PBS_NODEFILE -mca bt1 ^openib ./par_hybrid input_hybrid.txt
	mpirun -np $proc -machinefile $PBS_NODEFILE -mca bt1 ^openib ./par_hybrid big_input_hybrid.txt
	done
done
