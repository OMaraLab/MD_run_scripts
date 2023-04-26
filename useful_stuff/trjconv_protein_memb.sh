
#!/bin/bash

cd /path/to/system/prod/r1 # change to working direcotry of your production replicate

mkdir trjconv  # make a directory where your trjconv steps will go

# step one, turn off the entire pbc, and only save every tenth frame.  Figure out what the right "skip" number is so it saves a frame every 1 ns
gmx trjconv -f prod.xtc -o trjconv/01_your_system_name_trjconv_1ns_nobox.xtc -skip 10  -pbc nojump -center 2>&1 | tee trjconv/01.txt 

# step two, turn the pbc back on for atoms.  
gmx trjconv -f trjconv/01_your_system_name_trjconv_1ns_nobox.xtc -o trjconv/02_your_system_name_trjconv_1ns_pbc_atom.xtc  -center -s prod.tpr 2>&1 | tee trjconv/02.txt 
# chose protein for centering, entire system for output

# step three, turn the pbc back on for entire molecules.  
gmx trjconv -f trjconv/02_your_system_name_trjconv_1ns_pbc_atom.xtc  -o trjconv/03_your_system_name_trjconv_1ns_pbc_mol.xtc  -pbc atom -center -s prod.tpr 2>&1 | tee trjconv/03.txt 
# again, choose protein for centering and entire system for output

# i don't really know why you need step two and step three, but if you don't do them both it breaks

# step four:  look at each stage in vmd.  Does step three look normal?

# step five: repeat for rep 2 and rep 3

# step six:  use trjcat to combine all three replicates into a single long .xtc 

cd /path/to/system/prod/

# make your analysis folder if you haven't yet.  put a folder inside named clean, this is where your clean trajectories will go

gmx trjcat -f r1/trjconv/03*.xtc r2/trjconv/03*.xtc r3/trjconv/03*.xtc -o /path/to/analysis_folder/clean/your_system_name_1ns_cat.xtc -cat 2>&1 | tee cat.txt # concatenate trajectories

gmx trjconv -f /path/to/analysis_folder/clean/your_system_name_1ns_cat.xtc -o /path/to/analysis_folder/clean/your_system_name.pdb -s r1/prod.tpr -dump 0 # save the first frame as a pdb

