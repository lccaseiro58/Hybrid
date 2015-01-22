##!/bin/bash 
#PBS -l walltime=03:00 
#PBS -j oe 
#PBS -N a 
#PBS -l nodes=4:ppn=2 
#PBS -q mei
 
cat $PBS_NODEFILE 
cd ~/PCP/Hybrid 

#valgrind --toold=callgrind --cache-sim=yes --separate-threads=yes omp_mat_vect_rand_split
module add gnu/4.9.0

g++ -g -Wall -O3 -fopenmp -lm jacobi_par_omp.c -o par_omp

for th in 2 4 8 9 12 16
do
echo export OMP_NUM_THREADS = $th
export OMP_NUM_THREADS=$th

export OMPP_APPNAME=par_omp_input
echo ./kinst-par_omp input_omp.txt
./kinst-par_omp input_omp.txt

export OMPP_APPNAME=par_omp_biginput
echo ./kinst-par_omp big_input_omp.txt
./kinst-par_omp big_input_omp.txt

