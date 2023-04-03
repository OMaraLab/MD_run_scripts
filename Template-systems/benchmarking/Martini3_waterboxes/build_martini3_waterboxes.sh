#!/usr/bin/env bash

# Make martini 3.0 water boxes with new martini3 water 
# For these large boxes I'm using an series of EQ stages with increasing timesteps to even out problems in packing
# makes boxes from 10 nm^3 to 55 nm^3, intended to be used for benchmarking

i="10" # start with a 10 nm cubic box

while [ $i -le 55 ] # iterate up to a 55 nm cubic box  # testing only do tiny box
do

# make pure water box

mkdir M3_W_$i
./insane_M3_01_August_2022.py -x $i -y $i -z $i -d 0 -pbc cubic -sol W -excl -0.5 -o M3_W_$i/M3_W_$i.gro -p M3_W_$i/M3_W_$i.top
cp -r martini_v300.ff M3_W_$i/

#  no M3 parameter files ye

cp M3.mdp  M3_W_$i/M3_W_$i.mdp
cp emin.mdp  M3_W_$i/
cp EQ.mdp  M3_W_$i/

cd  M3_W_$i/



if [ -f /opt/homebrew/bin/gsed ] # use gsed if on a mac
then 

    gsed -i "s?Insanely solvated protein.?Martini 3.0 $i nm cubic waterbox with W water?"  M3_W_$i.top

else

    sed -i "s?Insanely solvated protein.?Martini 3.0 $i nm cubic waterbox with W water?"  M3_W_$i.top

fi


printf "q\n" | gmx make_ndx -f  M3_W_$i.gro -o  M3_W_$i.ndx

gmx grompp -f emin.mdp -c  M3_W_$i.gro -p  M3_W_$i.top -o  M3_W_$i\_min.tpr
gmx mdrun -v -deffnm  M3_W_$i\_min

gmx grompp -f EQ.mdp -c  M3_W_$i\_min.gro -p  M3_W_$i.top -n  M3_W_$i.ndx -o  M3_W_$i\_EQ.tpr
gmx mdrun -v -deffnm   M3_W_$i\_EQ

mv M3_W_$i\_EQ.gro   M3_W_$i\_start.gro 

gmx grompp -f M3_W_$i.mdp -c M3_W_$i\_start.gro  -o M3_W_$i.tpr -p M3_W_$i.top -n M3_W_$i.ndx -maxwarn 1   2>&1 | tee M3_W_$i\_prod_grompp.txt


cd ..

# increment by 5 nm

i=$[$i+5]

done


