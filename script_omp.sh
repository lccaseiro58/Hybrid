##!/bin/bash 
#PBS -l walltime=03:00 
#PBS -j oe 
#PBS -N omp 
#PBS -l nodes=4:ppn=2 
#PBS -q mei
 
cat $PBS_NODEFILE 
cd ~/PCP/Hybrid/OMP

module add gnu/4.9.0

for th in 2 4 8 9 12 16
do
echo export OMP_NUM_THREADS = $th
export OMP_NUM_THREADS=$th

echo ./par_ompN input_omp.txt
./par_ompN input_omp.txt

echo ./par_ompN big_input_omp.txt
./par_ompN big_input_omp.txt
done
