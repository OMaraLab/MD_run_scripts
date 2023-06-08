#!/bin/bash --login
#SBATCH --partition=gpu
#SBATCH --nodes=1              #1 nodes in this example
#SBATCH --ntasks-per-node=1    #1 tasks for the 1 GPUs in this job
#SBATCH --gpus-per-node=1      #1 GPUs in this job
#SBATCH --sockets-per-node=1   #Use the 1 slurm-sockets in this job
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --account=pawsey0420-gpu

# This script runs MD for 10 hours, then checks whether the job is complete by 
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
#
#
#
#
#   !!! IMPORTANT !!!           
#                                       
#   Check your mdp options.  
#
#   I am finding that my UA systems will only run with constraints=h-bonds and constraint-algorithm=lincs
# 
#   constraint-algorithm=shake fails at mdrun with an error that shake is not supported
#
#   constraints=all-bonds and constraint-algorithm=lincs fails at mdrun with an error that there are 
#   more coupled constaints that supported by the GPU LINCS code
#
#   constraints=none causes my test systems to blow up and hard crash ; looks like it's unstable even at a 2 fs timestep.
#
#   That's all I got.  I would love if there was more detailed documentation on this amd gromacs fork
#
#   Using PME with a 220k system, I get 2-30 ns/day in early testing
#   Using RF with a 400k particle system, Martin and Sharif were getting ~65 ns day.
#
#   Good luck!  
#       -Ada
#

#####  SETUP ENVIRONMENT  #####

# load modules

module load PrgEnv-cray
module load rocm craype-accel-amd-gfx90a

module load gcc/12.1.0
module load gromacs-amd-gfx90a/2022.3.amd1_174

module list
export GMX_MAXBACKUP=-1

unset OMP_NUM_THREADS

GMX='srun -l -u -c 8 gmx'
GMXMDRUN='srun -l -u -c 8 gmx mdrun -nb gpu -bonded gpu -pin on -update gpu -ntomp 8 -ntmpi 1 -maxh 9.95'


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

echo $SLURM_JOB_NODELIST  # print node id to the end of output file for benchmarking / identifying troublesome nodes
echo "CLUSTERID=SETONIX GPU" # print which cluster was used, so if a system is run on multiple systems I can separate the log files for benchmarking

#####  END  #####
