#!/bin/bash --login
#SBATCH --account=pawsey0420
#SBATCH --partition=work
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=0:20:00

#####  benchmarking bookkeeping  #####

echo $SLURM_JOB_NODELIST  # print node id to the beginning output file for benchmarking / identifying troublesome nodes
echo "CLUSTERID=SETONIX-CPU" # print which cluster was used, so if a system is run on multiple systems I can separate the log files for benchmarking
echo "BENCHMARKING RUN: [RUN DETAILS]" # text field with parameters being examined in benchmarking (system size, build, node configuration, mpi configuration, etc)

##### load the gromacs version to be used in test #####

module load gromacs/2023
module list


# Temporal workaround for avoiding Slingshot issues on shared nodes:
export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)
export MPICH_OFI_STARTUP_CONNECT=1
export MPICH_OFI_SKIP_NIC_SYMMETRY_TEST=1  # failing due to uneven NICs presenting a performance risk.  For EQ, I am not concerned about performance.  Lets try.

export OMP_NUM_THREADS=8

export NRANK=$(($SLURM_NTASKS/$OMP_NUM_THREADS))

#####  set up gromacs variables  #####

GMX="gmx_mpi"  # used to call gmx for anything that isn't mdrun

# GMXMDRUN is used to call MDRUN.  This is where you should set parameters about mpi configuration, load balancing, and mdrun optimisation (see https://manual.gromacs.org/current/user-guide/mdrun-performance.html ) 
# maxh specifies how long to run the simulation, in this case 15 minutes.
# this uses my default parameters for a 256 core run

GMXMDRUN="srun -N $SLURM_JOB_NUM_NODES -n $NRANK -c $OMP_NUM_THREADS -m block:block:block $GMX mdrun  -maxh 0.25 " 
echo "'GMXMDRUN COMMAND IS "$GMXMDRUN"'"


#####  START MD WITH A GROMPP  #####

$GMX grompp -f waterprod.mdp -c waterbox15_EQ.gro -o waterbox15.tpr -p waterbox15.top  -maxwarn 2 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt

#####  RUN MD  #####

$GMXMDRUN -v -deffnm waterbox15

#####  END  #####


