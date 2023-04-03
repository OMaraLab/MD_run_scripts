#!/usr/bin/env bash

# makes martini 2.2 and 2.2P waterboxs boxes from 10 nm^3 to 55 nm^3
# to be used for benchmarking

#afterwards you can use grompp_Martini2.sh to turn them into tprs so you only have one file per system for benchmarking

# Part one: Make martini 2.2 water boxes with 10% antifreeze water or 100 % polarizable water

i="10" # start with a 10 nm cubic box

while [ $i -le 55 ] # iterate up to a 55 nm cubic box
do

# make antifreeze water box

mkdir M2.2_WF_$i
./insane_13July2020.py -x $i -y $i -z $i -d 0 -pbc cubic -sol W:0.9 -sol WF:0.1 -excl -0.5 -o M2.2_WF_$i/M2.2_WF_$i.gro -p M2.2_WF_$i/M2.2_WF_$i.top
cp martini_v2.2.itp M2.2_WF_$i/
cp Martini2.2.mdp M2.2_WF_$i/M2.2_WF_$i.mdp
cp emin_2.2.mdp M2.2_WF_$i/
cp EQ_2.2.mdp M2.2_WF_$i/


cd M2.2_WF_$i/
sed -i 's?#include "martini.ff/martini_v2.2P.itp"?#include "martini_v2.2.itp"?' M2.2_WF_$i.top
sed -i 's?#include "martini.ff/martini_v2.0_ions.itp"?   ?' M2.2_WF_$i.top
sed -i 's?;#include "martini.ff/martini_v2.0_lipids_all_201506.itp"?   ?' M2.2_WF_$i.top
sed -i 's?;#include "martini.ff/martini_v2.0_lipids_all_201506.itp"?   ?' M2.2_WF_$i.top
sed -i 's?#include "martini.ff/martini_v2.0_brain_complex_b5-GMs_old.itp"?   ?' M2.2_WF_$i.top
sed -i 's?#include "martini.ff/martini_v2_PI-lipids.itp"?   ?' M2.2_WF_$i.top
sed -i "s?Insanely solvated protein.?Martini 2.2 $i nm cubic waterbox with antifreeze water?" M2.2_WF_$i.top
printf "q\n" | gmx make_ndx -f M2.2_WF_$i.gro -o M2.2_WF_$i.ndx

gmx grompp -f emin_2.2.mdp -c M2.2_WF_$i.gro -p M2.2_WF_$i.top -o M2.2_WF_$i\_min.tpr
gmx mdrun -v -deffnm M2.2_WF_$i\_min

gmx grompp -f EQ_2.2.mdp -c M2.2_WF_$i\_min.gro -p M2.2_WF_$i.top -n M2.2_WF_$i.ndx -o M2.2_WF_$i\_EQ.tpr
gmx mdrun -v -deffnm  M2.2_WF_$i\_EQ
mv  M2.2_WF_$i\_EQ.gro  M2.2_WF_$i\_start.gro 

gmx grompp -f M2.2_WF_$i.mdp -c M2.2_WF_$i\_start.gro  -o M2.2_WF_$i.tpr -p M2.2_WF_$i.top -n M2.2_WF_$i.ndx -maxwarn 1   2>&1 | tee M2.2_WF_$i\_prod_grompp.txt

cd ..

# Part two: make polarizable water box

mkdir M2.2P_PW_$i
./insane_13July2020.py -x $i -y $i -z $i -d 0 -pbc cubic -sol PW -excl -0.5 -o M2.2P_PW_$i/M2.2P_PW_$i.gro -p M2.2P_PW_$i/M2.2P_PW_$i.top
cp martini_v2.2P.itp M2.2P_PW_$i/
cp Martini2.2P.mdp M2.2P_PW_$i/M2.2P_PW_$i.mdp
cp emin_2.2P.mdp M2.2P_PW_$i/
cp EQ_2.2P.mdp M2.2P_PW_$i/

cd M2.2P_PW_$i/
sed -i 's?#include "martini.ff/martini_v2.2P.itp"?#include "martini_v2.2P.itp"?' M2.2P_PW_$i.top
sed -i 's?#include "martini.ff/martini_v2.0_ions.itp"?   ?' M2.2P_PW_$i.top
sed -i 's?;#include "martini.ff/martini_v2.0_lipids_all_201506.itp"?   ?' M2.2P_PW_$i.top
sed -i 's?;#include "martini.ff/martini_v2.0_lipids_all_201506.itp"?   ?' M2.2P_PW_$i.top
sed -i 's?#include "martini.ff/martini_v2.0_brain_complex_b5-GMs_old.itp"?   ?' M2.2P_PW_$i.top
sed -i 's?#include "martini.ff/martini_v2_PI-lipids.itp"?   ?' M2.2P_PW_$i.top
sed -i "s?Insanely solvated protein.?Martini 2.2P $i nm cubic waterbox with polarizable water?" M2.2P_PW_$i.top
printf "q\n" | gmx make_ndx -f M2.2P_PW_$i.gro -o M2.2P_PW_$i.ndx


gmx grompp -f emin_2.2P.mdp -c M2.2P_PW_$i.gro -p M2.2P_PW_$i.top -o M2.2P_PW_$i\_min.tpr
gmx mdrun -v -deffnm M2.2P_PW_$i\_min

gmx grompp -f EQ_2.2P.mdp -c M2.2P_PW_$i\_min.gro -n M2.2P_PW_$i.ndx -p M2.2P_PW_$i.top -o M2.2P_PW_$i\_EQ.tpr
gmx mdrun -v -deffnm  M2.2P_PW_$i\_EQ
mv  M2.2P_PW_$i\_EQ.gro  M2.2P_PW_$i\_start.gro 

gmx grompp -f M2.2P_PW_$i.mdp -c M2.2P_PW_$i\_start.gro  -o M2.2P_PW_$i.tpr -p M2.2P_PW_$i.top -n M2.2P_PW_$i.ndx -maxwarn 1  2>&1 | tee M2.2P_PW_$i\_prod_grompp.txt


cd ..

# increment by 5 nm

i=$[$i+5]

done





