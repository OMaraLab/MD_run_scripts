#!/bin/bash
#PBS -P q95
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=380GB
#PBS -l jobfs=16000MB
#PBS -l ncpus=96
#PBS -l other=mpi:hyperthread
#PBS -l wd
#PBS -r y

## General-purpose resubmit script for GROMACS jobs on Wiener

## this runs jobs in ten hour blocks with checkpinting and resubmission, 
## edit $GMXMDRUN to vary the runtime
## nsteps is set in the mdp
## Starting structure, .top, .ndx and .mdp should have the same name as the
## script, and all be in the same folder. Output will also have 
## the same name.
## Eg:  if script name is GlyT2_POPC_CHOL_r1, mdp is GlyT2_POPC_CHOL_r1.mdp
##
## IMPORTANT: SCRIPT NAME SHOULD NOT END WITH .sh


## The script automatically chooses the multithreading 
## mode based on the number of CPUs it finds. 


# Define error function so we can see the error code given when something
# important crashes
errexit ()
{
    errstat=$?
    if [ $errstat != 0 ]; then
        # A brief nap so PBS kills us in normal termination
        # Prefer to be killed by PBS if PBS detected some resource
        # excess
        sleep 5
        echo "Job returned error status $errstat - stopping job sequence $PBS_JOBNAME at job $PBS_JOBID"
        exit $errstat
    fi
}

# Guarantee GPUs are visible
# export CUDA_VISIBLE_DEVICES=$(seq 0 $(( $PBS_NGPUS-1 )) | tr  '\r\n' ',')

# Change to working directory - not necessary with #PBS -l wd
cd $PBS_O_WORKDIR

# Terminate the job sequence if the file STOP_SEQUENCE is found in pwd
if [ -f STOP_SEQUENCE ]; then
    echo "STOP_SEQUENCE file found - terminating job sequence $PBS_JOBNAME at job $PBS_JOBID"
    exit 0
fi

####  GROMACS time!  ####

# Load the gromacs module

module load gromacs/2021.4

# Let GROMACS choose how to parallelise everything, unless we specify something later:
unset OMP_NUM_THREADS

GMX='gmx'
# GMXMDRUN='gmx mdrun -maxh 9.95'
GMXMDRUN='gmx mdrun -maxh 9.95 -ntmpi 8 -ntomp -12'
#GMXMDRUN='mpirun gmx mdrun -maxh 9.95'



## GROMPP if there's no TPR file (eg, this is the first submission)
if [ ! -f ${PBS_JOBNAME}.tpr ]; then
    $GMX grompp -f ${PBS_JOBNAME}.mdp -c ${PBS_JOBNAME}_start.gro -o ${PBS_JOBNAME}.tpr -p ${PBS_JOBNAME}.top  -n ${PBS_JOBNAME}.ndx -maxwarn 2 || errexit
fi

## run MD!
$GMXMDRUN -v -deffnm ${PBS_JOBNAME} -cpi ${PBS_JOBNAME}.cpt || errexit
 # Notes:
 # -cpi: Continue from checkpoint if available, otherwise start new simulation
# -maxh n: Write a checkpoint and terminate after 0.99 * n hours
# -nb gpu: Die if GPU not usable (unfortunately, won't die if GPU isn't found)

# Check the log file for the number of steps completed
steps_done=`perl -n -e'/Statistics over (\d+) steps using (\d+) frames/ && print $1' ${PBS_JOBNAME}.log`
# Check the mdp file for the number of steps we want
steps_wanted=`perl -n -e'/nsteps\s*=\s*(\d+)/ && print $1' ${PBS_JOBNAME}.mdp`
# Resubmit if we need to
if (( steps_done < steps_wanted )); then
    echo "Job ${PBS_JOBID} terminated with ${steps_done}/${steps_wanted} steps finished."
    echo "Submitting next job in sequence $PBS_JOBNAME."
    qsub $PBS_JOBNAME
fi