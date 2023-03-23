#!/bin/bash

# This is a script to run equilibritation on an md system, slowly relating the position restriction over a number of runs.

# load intended gromacs version gromacs

module unload gromacs
module load gromacs/2021.4

# EDIT THESE

STARTFILE="path/to/coordinates_after_EM2.gro"    # eg: "05_hOCT1_THII_p2_POPC_CLR_water_ions_EM2.gro"
SYSNAME="sysname"                                # eg: "hOCT1_THII_p2_POPC_CLR"

# YOU SHOULD ADJUST YOUR STEP NUMBERS TO BE SEQUENTIAL FOLLOWING ON FROM $STARTFILE

gmx grompp -f equil_memb_1000.mdp -c ${STARTFILE} -r ${STARTFILE} -n ${SYSNAME}.ndx -p ${SYSNAME}.top -o 06_${SYSNAME}_equil_1000.tpr  -maxwarn 2  2>&1 | tee 06_grompp.txt
gmx mdrun -v -deffnm 06_${SYSNAME}_equil_1000

gmx grompp -f equil_memb_500.mdp -c 06_${SYSNAME}_equil_1000.gro -r 06_${SYSNAME}_equil_1000.gro -n ${SYSNAME}.ndx -p ${SYSNAME}.top -o 07_${SYSNAME}_equil_500.tpr  -maxwarn 2  2>&1 | tee 07_grompp.txt
gmx mdrun -v -deffnm 07_${SYSNAME}_equil_500

gmx grompp -f equil_memb_100.mdp -c 07_${SYSNAME}_equil_500.gro -r 07_${SYSNAME}_equil_500.gro -n ${SYSNAME}.ndx -p ${SYSNAME}.top -o 08_${SYSNAME}_equil_100.tpr  -maxwarn 2  2>&1 | tee 08_grompp.txt
gmx mdrun -v -deffnm 08_${SYSNAME}_equil_100

gmx grompp -f equil_memb_50.mdp -c 08_${SYSNAME}_equil_100.gro -r 08_${SYSNAME}_equil_100.gro -n ${SYSNAME}.ndx -p ${SYSNAME}.top -o 09_${SYSNAME}_equil_50.tpr  -maxwarn 2  2>&1 | tee 09_grompp.txt
gmx mdrun -v -deffnm 09_${SYSNAME}_equil_50

gmx grompp -f equil_memb_10.mdp -c 09_${SYSNAME}_equil_50.gro -r 09_${SYSNAME}_equil_50.gro -n ${SYSNAME}.ndx -p ${SYSNAME}.top -o 10_${SYSNAME}_equil_10.tpr  -maxwarn 2  2>&1 | tee 10_grompp.txt
gmx mdrun -v -deffnm 10_${SYSNAME}_equil_10