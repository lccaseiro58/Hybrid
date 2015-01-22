##!/bin/bash 
#PBS -l walltime=03:00 
#PBS -j oe 
#PBS -N a 
#PBS -l nodes=4:ppn=2 
#PBS -q mei
 
cat $PBS_NODEFILE 
cd ~/PCP 

#valgrind --toold=callgrind --cache-sim=yes --separate-threads=yes omp_mat_vect_rand_split
module add gnu/4.9.0

export OMPP_CTR1 = PAPI_L2_DCW
export OMPP_CTR4 = PAPI_L3_DCW
export OMPP_CTR2 = PAPI_L2_DCR
export OMPP_CTR3 = PAPI_L1_DCM

for th in 1 2 4
do
echo export OMP_NUM_THREADS = $th
export OMP_NUM_THREADS=$th

export OMPP_APPNAME=omp_mat_vect_rand_split-8000000x8
echo ./kinst-omp_mat_vect $th 8000000 8
./kinst-omp_mat_vect $th 8000000 8

export OMPP_APPNAME=omp_mat_vect_rand_split-8000x8000
echo ./kinst-omp_mat_vect $th 8000 8000
./kinst-omp_mat_vect $th 8000 8000

export OMPP_APPNAME=omp_mat_vect_rand_split-8x8000000
echo ./kinst-omp_mat_vect $th 8 8000000
./kinst-omp_mat_vect $th 8 8000000
