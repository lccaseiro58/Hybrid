!/bin/bash
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
mpirun -np $proc -machinefile $PBS_NODEFILE -mca bt1 ^openib ./par_mpi input_mpi.txt
mpirun -np $proc -machinefile $PBS_NODEFILE -mca bt1 ^openib ./par_mpi big_input_mpi.txt
done
