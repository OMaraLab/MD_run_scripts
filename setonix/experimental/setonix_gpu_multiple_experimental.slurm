#!/bin/bash --login
#SBATCH --account=pawsey0420-gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1               # Total number of nodes
##SBATCH --ntasks-per-node=8     # 7 MPI ranks per node
#SBATCH --gres=gpu:1       # Allocate one gpu per MPI rank
##SBATCH --exclusive
#SBATCH --time=12:00:00

#### IMPORTANT:  THIS SCRIPT IS EXPERIMENTAL #####

# This script runs MD for 12 hours, then checks whether the job is complete by 
# comparing the number of steps in the log file with nsteps in the mdp.  
# If the job is not done, it resubmits automatically.
#
# Do not specify a job name, instead rename this file with the intended job name / system name
#
# Requires the following files
#
#   This script saved with file name :      <SYSTEMNAME>              NB: note no extenstion
#   MDP file with filename :                <SYSTEMNAME>.mdp
#
# and then either:
#
#   Input coordinates with filename :       <SYSTEMNAME>_start.gro
#   Index file with filename :              <SYSTEMNAME>.ndx
#   Topology file with filename :           <SYSTEMNAME>.top
#   All force field include files as specified in the topology file
#
# or:
#
#   A .tpr file with filename :             <SYSTEMNAME>.top


echo $SLURM_JOB_NODELIST  # print node id to the end of output file for benchmarking / identifying troublesome nodes
echo "CLUSTERID=SETONIX-GPU" # print which cluster was used, so if a system is run on multiple systems I can separate the log files for benchmarking

module load craype-accel-amd-gfx90a
#module load rocm/5.2.3
module load gromacs-amd-gfx90a/2023 #amdgromacs/2023-deolboo
module list
export MPICH_GPU_SUPPORT_ENABLED=1

export GMX_GPU_DD_COMMS=true
export GMX_GPU_PME_PP_COMMS=true
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1
export OMP_NUM_THREADS=8


GMX="gmx_mpi"
GMXMDRUN="srun --export=ALL -n 7 gmx_mpi mdrun -pin on  -maxh 11.95  -ntomp 8 -noconfout -nb gpu -bonded gpu -pme gpu -npme 1 -v -nstlist 400 -gpu_id 0123456"



# Define error function so we can see the error code given when something
# important crashes

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


#####  END  #####


