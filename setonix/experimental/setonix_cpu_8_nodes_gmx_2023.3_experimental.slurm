#!/bin/bash --login
#SBATCH --account=pawsey0420
#SBATCH --partition=work
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
##SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --time=10:00:00
#SBATCH --exclude=nid00[2024-2055],nid00[2792-2823]
#### nodelist is to only use nodes with a higher number of NIC. This is a mismatched hardware fault of Setonix where some nodes were installed with 1 NIC and 2 NICs and they don't talk to eachother with MPI.

#### run with multiple configurations where 8 was replaced with some multiple of 128 betwen 128 and 1024 inclusive

#####  benchmarking bookkeeping  #####


echo $SLURM_JOB_NODELIST  # print node id to the beginning output file for benchmarking / identifying troublesome nodes
echo "CLUSTERID=SETONIX-CPU" # print which cluster was used, so if a system is run on multiple systems I can separate the log files for benchmarking
echo "BENCHMARKING RUN: purempi" # text field with parameters being examined in benchmarking (system size, build, node configuration, mpi configuration, etc)

##### load the gromacs version to be used in test #####

module load gromacs/2023.3-mx33m5h

module list


export OMP_NUM_THREADS=1

export NRANK=$(($SLURM_NTASKS/$OMP_NUM_THREADS)) # 128

#####  set up gromacs variables  #####

GMX="gmx_mpi"  # used to call gmx for anything that isn't mdrun

# GMXMDRUN is used to call MDRUN.  This is where you should set parameters about mpi configuration, load balancing, and mdrun optimisation (see https://manual.gromacs.org/current/user-guide/mdrun-performance.html ) 
# maxh specifies how long to run the simulation, in this case 15 minutes.
# this uses my default parameters for a 256 core run

GMXMDRUN="srun $GMX mdrun -maxh 9.95 " 
echo "'GMXMDRUN COMMAND IS "$GMXMDRUN"'"


errexit ()

{
    errstat=$?
    if [ $errstat != 0 ]; then
        # A brief nap so slurm kills us in normal termination
        # Prefer to be killed by slurm if slurm detected some resource excess
        sleep 5
        echo "Job returned error status $errstat - stopping job sequence $SLURM_JOB_NAME at job $SLURM_JOB_ID"
        exit $errstat
    fi
}

#####  START MD WITH A GROMPP  #####

## GROMPP if there's no TPR file (eg, this is the first submission)
if [ ! -f ${SLURM_JOB_NAME}.tpr ]; then
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top  -n ${SLURM_JOB_NAME}.ndx -maxwarn 2 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
fi

#####  RUN MD FOR 9.95 HOURS  #####

$GMXMDRUN -v -deffnm ${SLURM_JOB_NAME} -cpi ${SLURM_JOB_NAME}.cpt || errexit

#####  CHECK IF JOB IS DONE; IF NOT DONE RESUBMIT THIS SCRIPT  #####


# Check the log file for the number of steps completed
steps_done=`perl -n -e'/Statistics over (\d+) steps using (\d+) frames/ && print $1' ${SLURM_JOB_NAME}.log`
# Check the mdp file for the number of steps we want
steps_wanted=`perl -n -e'/nsteps\s*=\s*(\d+)/ && print $1' ${SLURM_JOB_NAME}.mdp`
# Resubmit if we need to
if (( steps_done < steps_wanted )); then
    echo "Job ${SLURM_JOB_NAME} terminated with ${SLURM_JOB_NAME}/${SLURM_JOB_NAME} steps finished." 
    echo "Submitting next job in sequence ${SLURM_JOB_NAME}."
    sbatch ${SLURM_JOB_NAME}
fi

# ## backup logfile for live debugging

# cp ${SLURM_JOB_NAME}.log lastlog.txt

### append with some info onthe run for debugging and benchmarking

echo $SLURM_JOB_NODELIST  # print node id to the end of output file for benchmarking / identifying troublesome nodes
echo "CLUSTERID=SETONIX GPU" # print which cluster was used, so if a system is run on multiple systems I can separate the log files for benchmarking

#####  END  #####
