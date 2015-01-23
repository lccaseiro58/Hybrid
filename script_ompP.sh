##!/bin/bash 
#PBS -l walltime=03:00 
#PBS -j oe 
#PBS -N ompP 
#PBS -l nodes=4:ppn=2 
#PBS -q mei
 
cat $PBS_NODEFILE 
cd ~/PCP/Hybrid/OMPP

module add gnu/4.9.0

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/papi/5.3.2/lib/

for th in 8
do
echo export OMP_NUM_THREADS = $th
export OMP_NUM_THREADS=$th

export OMPP_APPNAME=par_omp_input
echo ./par_omp input_omp.txt
./par_omp input_omp.txt

export OMPP_APPNAME=par_omp_biginput
echo ./par_omp big_input_omp.txt
./par_omp big_input_omp.txt
done