i=10

while [ $i -le 55 ] # iterate up to a 55 nm cubic box; test build with first box
do


cd M2.2_WF_$i


gmx grompp -f ../Martini2.2.mdp -c M2.2_WF_$i\_start.gro  -o M2.2_WF_$i.tpr -p M2.2_WF_$i.top -n M2.2_WF_$i.ndx -maxwarn 1   2>&1 | tee M2.2_WF_$i\_prod_grompp.txt

cd ..


i=$[$i+5]

done

i=10

while [ $i -le 55 ] # iterate up to a 55 nm cubic box; test build with first box
do


cd M2.2P_PW_$i


gmx grompp -f ../Martini2.2.mdp -c M2.2P_PW_$i\_start.gro  -o M2.2P_PW_$i.tpr -p M2.2P_PW_$i.top -n M2.2P_PW_$i.ndx -maxwarn 1  2>&1 | tee M2.2P_PW_$i\_prod_grompp.txt

cd ..


i=$[$i+5]

done
