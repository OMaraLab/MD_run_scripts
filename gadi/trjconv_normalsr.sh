#!/bin/bash
#PBS -P q95
#PBS -q normalsr
#PBS -l walltime=10:00:00
#PBS -l mem=50GB
#PBS -l jobfs=50GB
#PBS -l ncpus=13
#PBS -l other=mpi:hyperthread
#PBS -l wd
#PBS -r y
#PBS -l storage=scratch/q95

# example script to run trajectory processing on gadi
# you will need to adapt all commands and scripts to match your file structure and system design
module load gromacs/2021.4

ligs=( "DAMP" "MF1A" "MF1B" "PFMA" "PFMB")
sites=( "D474" "E386")
apos=(  "apo_HE357")

GMX="gmx"
path="path/to/prod_data/"
length="300" # desired final frame in ns
mkdir clean
mkdir clean/ff
mkdir clean/ff/top
mkdir clean/ndx
mkdir clean/mdp
mkdir clean/tpr
mkdir trjconv
cp -r $path/gromos54a7_atb.ff clean/

GMX_MAXBACKUP=-1

function run_trjconv() {

    i=1
    cp  $path/${sys}/r$i/hOCT1_${sys}*_r$i.mdp clean/mdp/hOCT1_${sys}.mdp
    cp  $path/${sys}/r$i/hOCT1_${sys}*_r$i.ndx clean/ndx/hOCT1_${sys}.ndx
    cp  $path/${sys}/r$i/hOCT1_${sys}*_r$i.top clean/ff/top/hOCT1_${sys}.top
    $GMX grompp -f  clean/mdp/hOCT1_${sys}.mdp -p clean/ff/top/hOCT1_${sys}.top -n clean/ndx/hOCT1_${sys}.ndx -c $path/${sys}/r$i/hOCT1_${sys}*_r${i}_start.gro -maxwarn 2 -o clean/tpr/hOCT1_${sys}.tpr
    # regrompp with the version of gromacs you intend to use for your analysis
    # often you might wait until you have data locally before running this step
    
    while [ $i -le 3 ]; do
        printf "1\n0\n" | GMX_MAXBACKUP=-1 $GMX trjconv -f  $path/${sys}/r${i}/hOCT1_${sys}*_r${i}.xtc -o trjconv/nojump.xtc -s clean/tpr/hOCT1_${sys}.tpr -dt 1 -tu ns -pbc nojump -center -n clean/ndx/hOCT1_${sys}.ndx -b 0 -e ${length}
        printf "1\n0\n" | GMX_MAXBACKUP=-1 $GMX trjconv -f  trjconv/nojump.xtc -o trjconv/atom.xtc  -s clean/tpr/hOCT1_${sys}.tpr -pbc atom -center -n clean/ndx/hOCT1_${sys}.ndx
        printf "1\n0\n" | GMX_MAXBACKUP=-1 $GMX trjconv -f  trjconv/atom.xtc -o clean/hOCT1_${sys}_r${i}_1ns.xtc  -s clean/tpr/hOCT1_${sys}.tpr -pbc mol -center -n clean/ndx/hOCT1_${sys}.ndx
        printf "1\n0\n" | GMX_MAXBACKUP=-1 $GMX trjconv -f  trjconv/atom.xtc -o clean/hOCT1_${sys}_r${i}_start.gro  -s clean/tpr/hOCT1_${sys}.tpr -tu ns -pbc mol -center -n clean/ndx/hOCT1_${sys}.ndx -dump 0
        printf "1\n0\n" | GMX_MAXBACKUP=-1 $GMX trjconv -f  trjconv/atom.xtc -o clean/hOCT1_${sys}_r${i}_end.gro  -s clean/tpr/hOCT1_${sys}.tpr -tu ns -pbc mol -center -n clean/ndx/hOCT1_${sys}.ndx -dump ${length}
        i=$(( i + 1 ))
    done

    GMX_MAXBACKUP=-1 $GMX trjcat -f  clean/hOCT1_${sys}_r1_1ns.xtc clean/hOCT1_${sys}_r2_1ns.xtc clean/hOCT1_${sys}_r3_1ns.xtc -o clean/hOCT1_${sys}_cat_1ns.xtc -cat

}


for l in "${ligs[@]}"; do
  for s in "${sites[@]}"; do
    sys="${l}_${s}"
    run_trjconv
  done
done
