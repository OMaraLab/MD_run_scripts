#!/usr/bin/env bash

# makes gromos54a7 waterboxs boxes from 5 nm^3 to 20 nm^3
# to be used for benchmarking

module purge
module load gromacs/2021.4

$GMX="gmx"
$GMXMDRUN="gmx mdrun"

sizes=("5" "7.5" "10" "12.5" "15" "17.5" "20")


for i in ${sizes[@]}; do

# make water box

mkdir 54A7_water_${i}

cp water.mdp 54A7_water_${i}/
cp waterprod.mdp 54A7_water_${i}/54A7_water_${i}.mdp
cp template.top 54A7_water_${i}/54A7_water_${i}.top

cd 54A7_water_${i}/

$GMX solvate -cs spc216.gro -o 00_54A7_water_${i}.gro -box 5 5 5 -p 54A7_water_${i}.top    2>&1 | tee 00_solvate.txt

$GMX grompp -f ../minimise.mdp -c 00_54A7_water_${i}.gro -p 54A7_water_${i}.top -o 01_54A7_water_${i}_min.tpr  2>&1 | tee 01_grompp.txt
$GMXMDRUN -v -deffnm 01_54A7_water_${i}_min

$GMX grompp -f ../EQ.mdp -c  01_54A7_water_${i}_min.gro -p 54A7_water_${i}.top -o 02_54A7_water_${i}_EQ.tpr  2>&1 | tee 02_grompp.txt
$GMXMDRUN -v -deffnm  02_54A7_water_${i}_EQ

$GMX grompp -f 54A7_water_${i}.mdp -c 02_54A7_water_${i}_EQ.gro  -o 03_54A7_water_${i}_prod.tpr -p 54A7_water_${i}.top -maxwarn 1   2>&1 | tee 03_grompp.txt

cd ..

done





