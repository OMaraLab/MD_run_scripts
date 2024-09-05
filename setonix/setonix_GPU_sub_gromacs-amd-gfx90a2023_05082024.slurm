#!/bin/bash --login
#SBATCH --partition=gpu
#SBATCH --nodes=1              #1 nodes in this example
##SBATCH --ntasks-per-node=1    #1 tasks for the 1 GPUs in this job
#SBATCH --gpus-per-node=1      #1 GPUs in this job
##SBATCH --sockets-per-node=1   #Use the 1 slurm-sockets in this job
#SBATCH --time=10:00:00
#SBATCH --account=pawsey0420-gpu

#   !!! IMPORTANT !!!           
#  
#   We are still working on improving the GPU Gromacs submission script and gromacs module as a part of the Pawsey Uptake Project. 
#                                       
#

#####  SETUP ENVIRONMENT  #####

# load modules (adapt this for your own needs):

module load gromacs-amd-gfx90a/2023
module load rocm/5.7.3 # to report GPU
module load fftw/3.3.10
module list

# export GMX_GPU_DD_COMMS=true # only included when there are multiple GPUs.
export GMX_GPU_PME_PP_COMMS=true
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GMX_FORCE_GPU_AWARE_MPI=1
export GMX_FORCE_UPDATE_DEFAULT_GPU=true
export OMP_NUM_THREADS=8
#

### Report GPU utilisation ###

cat << EOF > select_gpu
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

### Create GMX variables ###


#GMX='singularity exec -B <moduleDir>/99-shpc.sh:/.singularity.d/env/99-shpc.sh <container> gmx "$@"'
GMX='gmx_mpi'  # srun -l -u -c 8 
GMXMDRUN=" srun -N1 -n 1 -c 8 --gres=gpu:1 --gpus-per-task=1 --gpu-bind=closest $GMX -u -c 8 mdrun -nb gpu -pme gpu -pin on -maxh 9.95" # -ntmpi 1 (include in mdrun with multiple GPUs)


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
    $GMX grompp -f ${SLURM_JOB_NAME}.mdp -r ${SLURM_JOB_NAME}_start.gro -c ${SLURM_JOB_NAME}_start.gro -o ${SLURM_JOB_NAME}.tpr -p ${SLURM_JOB_NAME}.top  -n ${SLURM_JOB_NAME}.ndx -maxwarn 2 &> ${SLURM_JOB_NAME}_grompp_${SLURM_JOB_ID}.txt
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


### Report GPU utilisation ###

nsamples=100
sleeptime=5
for i in $(seq 1 $nsamples); do 
  rocm-smi -u >> smi_output_partial.txt
  sleep $sleeptime
done

wait

#####  END  #####
