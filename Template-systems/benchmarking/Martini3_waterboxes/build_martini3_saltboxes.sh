#!/usr/bin/env bash

# Make martini 3.0 water boxes with new martini3 water and 150 mM NaCl
# For these large boxes I'm using an series of EQ stages with increasing timesteps to even out problems in packing
# makes boxes from 10 nm^3 to 55 nm^3, intended to be used for benchmarking

#afterwards you can use grompp_Martini3.sh to turn them into tprs so you only have one file per system for benchmarking

i="10" # start with a 10 nm cubic box

while [ $i -le 55 ] # iterate up to a 40 nm cubic box  # total system range goes up to 55 but I already have the 40 nmto 55 nm systems
do

# make M3 water box with 150 mM NaCl

mkdir M3_W_ions_$i
./insane_M3_01_August_2022.py -x $i -y $i -z $i -d 0 -pbc cubic -sol W -salt 0.15 -excl -0.5 -o M3_W_ions_$i/M3_W_ions_$i.gro -p M3_W_ions_$i/M3_W_ions_$i.top
cp -r martini_v300.ff M3_W_ions_$i/

#  no M3 parameter files yet

cp M3.mdp  M3_W_ions_$i/M3_W_ions_$i.mdp
cp emin.mdp  M3_W_ions_$i/
cp EQ_05fs.mdp  M3_W_ions_$i/
cp EQ_10fs.mdp  M3_W_ions_$i/
cp EQ_15fs.mdp  M3_W_ions_$i/
cp EQ_20fs.mdp  M3_W_ions_$i/

cd  M3_W_ions_$i/

# fix ion names in gro file
# I've hacked insane to generate the ions with martini 3 residue names, but it seems to pull the residue name for the bead name
# easier to let it to it wrong and fix the gro file for now

if [ -f /opt/homebrew/bin/gsed ] # use gsed if on a mac
then 

    gsed -i "s?Insanely solvated protein.?Martini 3.0 $i nm cubic waterbox with W water and 150 mM NaCl?"  M3_W_ions_$i.top
    gsed -i "s?TCL    TCL?TCL     CL?"  M3_W_ions_$i.gro
    gsed -i "s?TNA    TNA?TNA     NA?"  M3_W_ions_$i.gro

else

    sed -i "s?Insanely solvated protein.?Martini 3.0 $i nm cubic waterbox with W water and 150 mM NaCl?"  M3_W_ions_$i.top
    sed -i "s?TCL    TCL?TCL     CL?"  M3_W_ions_$i.gro
    sed -i "s?TNA    TNA?TNA     NA?"  M3_W_ions_$i.gro

fi


printf "q\n" | gmx make_ndx -f  M3_W_ions_$i.gro -o  M3_W_ions_$i.ndx

gmx grompp -f emin.mdp -c  M3_W_ions_$i.gro -p  M3_W_ions_$i.top -o  M3_W_ions_$i\_min.tpr
gmx mdrun -v -deffnm  M3_W_ions_$i\_min

gmx grompp -f EQ_05fs.mdp -c  M3_W_ions_$i\_min.gro -p  M3_W_ions_$i.top -n  M3_W_ions_$i.ndx -o  M3_W_ions_$i\_EQ_05.tpr
gmx mdrun -v -deffnm   M3_W_ions_$i\_EQ_05

gmx grompp -f EQ_10fs.mdp -c  M3_W_ions_$i\_EQ_05.gro -p  M3_W_ions_$i.top -n  M3_W_ions_$i.ndx -o  M3_W_ions_$i\_EQ_10.tpr
gmx mdrun -v -deffnm   M3_W_ions_$i\_EQ_10

gmx grompp -f EQ_15fs.mdp -c  M3_W_ions_$i\_EQ_10.gro -p  M3_W_ions_$i.top -n  M3_W_ions_$i.ndx -o  M3_W_ions_$i\_EQ_15.tpr
gmx mdrun -v -deffnm   M3_W_ions_$i\_EQ_15

gmx grompp -f EQ_20fs.mdp -c  M3_W_ions_$i\_EQ_15.gro -p  M3_W_ions_$i.top -n  M3_W_ions_$i.ndx -o  M3_W_ions_$i\_EQ_20.tpr
gmx mdrun -v -deffnm   M3_W_ions_$i\_EQ_20

mv M3_W_ions_$i\_EQ_20.gro   M3_W_ions_$i\_start.gro 

gmx grompp -f M3_W_ions_$i.mdp -c M3_W_ions_$i\_start.gro  -o M3_W_ions_$i.tpr -p M3_W_ions_$i.top -n M3_W_ions_$i.ndx -maxwarn 1  2>&1 | tee M3_W_ions_$i\_prod_grompp.txt


cd ..

# increment by 5 nm

i=$[$i+5]

done





