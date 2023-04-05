i=10

while [ $i -le 55 ] # iterate up to a 55 nm cubic box; test build with first box
do


cd M3_W_$i


gmx grompp -f ../M3.mdp -c M3_W_$i\_start.gro  -o M3_W_$i.tpr -p M3_W_$i.top -n M3_W_$i.ndx -maxwarn 1   2>&1 | tee M3_W_$i\_prod_grompp.txt

cd ..


i=$[$i+5]

done

i=10

while [ $i -le 55 ] # iterate up to a 55 nm cubic box; test build with first box
do


cd M3_W_ions_$i


gmx grompp -f ../M3.mdp -c M3_W_ions_$i\_start.gro  -o M3_W_ions_$i.tpr -p M3_W_ions_$i.top -n M3_W_ions_$i.ndx -maxwarn 1  2>&1 | tee M3_W_ions_$i\_prod_grompp.txt

cd ..


i=$[$i+5]

done
