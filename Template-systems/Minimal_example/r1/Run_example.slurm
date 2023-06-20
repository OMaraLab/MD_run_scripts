#!/bin/bash

#SBATCH --account=a_omara
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH --partition=general


# This script runs MD for 4 hours
#
# Do not specify a job name, instead rename this file with the intended job name / system name
#
# Requires the following files
#
#   This script, saved with the filename = <SYSTEMNAME>              NB: note no extenstion
#   MDP file with filename =               <SYSTEMNAME>.mdp          This is a control file, it contains instructions for gromacs
#
# and then either:
#
#   Input coordinates with filename :       <SYSTEMNAME>_start.gro
#   Index file with filename :              <SYSTEMNAME>.ndx (make with gmx make_ndx)
#   Topology file with filename :           <SYSTEMNAME>.top
#   All force field include files as specified in the topology file
#
# or:
#
#   A .tpr file with filename :             <SYSTEMNAME>.tpr (this is the entire simulation in a single binary, made with gmx grompp)
#
#
# change the accountat the top of the file to your accounting group on bunya
#
# to run this job, go to the folder that contains your script and type sbatch <SYSTEMNAME>
# so ith the example system, 
#
#    $ cd path/to/Minimal_example/r1
#    $ sbatch Run_example.slurm
#

#####  SETUP ENVIRONMENT  #####

module load gromacs/2021.3-foss-2021a  # this choses which gromacs to use.  
# 2021.3-foss-2021a is the default module installed on bunya.  I have no idea how efficient it is; I use a custom install.
module list # print a list to the logfile

export OMP_NUM_THREADS=1

#####  START MD WITH A GROMPP  #####

## GROMPP if there's no TPR file (eg, this is the first submission)
if [ ! -f Example.tpr ]; then
    gmx grompp -f Example.mdp -c Example_start.gro -o Example.tpr -p Example.top  -n Example.ndx -maxwarn 2 &> Example_grompp_${SLURM_JOB_ID}.txt
fi

#####  RUN MD  #####

GMXOPTIONS="-ntmpi 48 -ntomp 1 -pin on -pme cpu -dlb yes " # this is a bunch of configuration settings telling gromacs how to use the CPU and GPU resources.  I don't use 2021.3-foss-2021a so I don't know what's optimal.

srun gmx mdrun -v -deffnm Example -cpi Example.cpt $GMXOPTIONS # run MD of the system Example.tpr with the options from GMXOPTIONS

#####  END  #####
