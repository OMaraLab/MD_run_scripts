#!/bin/bash
#PBS -P bd36
#PBS -q gpuvolta
#PBS -l walltime=04:00:00
#PBS -l mem=32GB
#PBS -l jobfs=16000MB
#PBS -l ngpus=1
#PBS -l ncpus=12
#PBS -l other=mpi:hyperthread
#PBS -l wd
#PBS -r y
#PBS -l storage=scratch/bd36+gdata/q95`
# THIS WAS ADA TEST 01       
# GOAL:  SEE IF IT WORKS

## General-purpose resubmit script for GROMACS jobs on Bunya

## this runs jobs in four hour blocks with checkpinting and resubmission, 
## edit $GMXMDRUN to vary the runtime
## nsteps is set in the mdp
## Starting structure, .top, .ndx and .mdp should have the same name as the
## script, and all be in the same folder. Output will also have 
## the same name.
## Eg:  if script name is GlyT2_POPC_CHOL_r1, mdp is GlyT2_POPC_CHOL_r1.mdp
##
## IMPORTANT: SCRIPT NAME SHOULD NOT END WITH .sh
##
## IMPORTANT: If SCRIPTNAME.mdp does not exist, automatic resubmission will fail
##
## IMPORTANT:  YOU NEED TO HAVE ACCESS TO A GROMACS SINGULARITY CONTAINER.  

 
# Load module, always specify version number.
module load singularity
 
# Must include `#PBS -l storage=scratch/ab12+gdata/yz98` if the job
# needs access to `/scratch/ab12/` and `/g/data/yz98/`. Details on:
# https://opus.nci.org.au/display/Help/PBS+Directives+Explained
 
# Run program
singularity exec mycontainer.sif /path/to/my/program

#####  SETUP ENVIRONMENT  #####

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

# srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l 

export SINGULARITY_TMPDIR=/scratch/bd23/aq8103/tmp/ 



GMX='singularity run --nv /g/data/q95/SHARED/gromacs_2022.3.sif gmx'
GMXMDRUN='singularity run --nv /g/data/q95/SHARED/gromacs_2022.3.sif gmx mdrun -ntmpi 1 -ntomp 12 -pin on -dlb yes -maxh 3.95'


# GMXMDRUN='singularity run --nv /scratch/user/uqhaitke/containers/gromacs_2022.3.sif gmx mdrun -ntmpi 4 -ntomp 8 -pin on -nb gpu -dlb yes -maxh 9.95'

#####  START MD WITH A GROMPP  #####

## GROMPP if there's no TPR file (eg, this is the first submission)
if [ ! -f ${SLURM_JOB_NAME}.tpr ]; then
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top  -n ${SLURM_JOB_NAME}.ndx -maxwarn 2 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
fi

#####  RUN MD FOR 3.95 HOURS  #####

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
    qsub ${SLURM_JOB_NAME}
fi

#####  END  #####