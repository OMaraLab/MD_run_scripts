#! usr/bin/env bash

# calculate contact fractions for gromacs trajectory

###     This isn't done yet
###     What I need:  calculate contact fractions between a drug and ligand
###     for each replicate of a simulation system
###     currently it does it manually, requiring me to make a contacts.sh
###     script for every single drug and pose, and then run them all 
###     in sequence with another bash script.  
###
###     That's a great way to generate errors with typos; it needs automation
###

###
###     OVERVIEW
###
###     The whole process is spread across multiple scripts
###
###      * run_contacts.sh runs all the contact.sh scripts
###
###      * contacts.sh runs the contacts for a particular system,
###             for each replicate in turn, specified manually
###
###      * contactFreq.tcl is the program to calculate arbitrary contact 
###             fractions in vmd
###
###      * contactfreq.sh is a little script called by vmd, which specifies 
###             which specific contact fractions are to be calculated,
###             and calls contactFreq.tcl on those atomgroups
###
###      * process_contacts.py converts the contacts for any particular system 
###              to a dataframe in a .csv file, and copies it to one place for
###              further processing
###
###      * concat.py concatenates all the per system .csv files into one
###              big csv file that's easy to work with
###
###     I want to automate all of this.

###     SHORT TERM GOAL:  Automate contacts.sh to run with a command line argument
###     so I don't need to make a contacts.sh for every single system.  
###     It will need to correctly identify the .gro and .xtc file for that system, 
###     and then save the output and run process_contacts.py
###     
###     Ideally this would be able to run both individual replicates and 
###     concatenated trajectories
###     
###     Secondly, this process only generates a contact fraction if there is any
###     contact between drug and residue in a given replicate
###
###     It does NOT generate a zero freuqncy contact if there is no contact.
###
###     To see why this is a problem, consider the following case case:
###
###     A drug has 100% contact with some arbitrary residue in replicates one and
###     two, and no contact with that arbitrary residue in replicate three.  
###     The expected behaviour is that .mean() of contacts with that residue should
###     be df[foo].mean() = mean( [100, 100, 0] ) = 66% , but the observed behaviour
###     is usually df[foo].mean() = mean ( [100, 100]) = 100 % 
###
###     To avoid this issue I currently need to do one of three things:
###
###       * calculate per replicate contacts AND full length contacts (wastes time)
###       * manually add zero contact rows to the final dataframe, but only when
###         a contact is not present in a particular replicate (pain in the ass)
###       * not calculate per replicate contacts at all (only works if trajectories
###         converge, which they should but sometimes they don't and in those cases 
###         I often want to look at them in isolation)   
###
###     I don't have a fix in mind yet but it's something to keep in mind.


###     POSSIBLE LONG TERM GOAL:  
###     Build the entire process into a single python script.
###
###     Eliminates the need to run individual sub scripts.  
###
###     One drawback of this proposed new approach is that contacts take a 
###     long time to run.  As it currently stands, if I need to rerun contacts 
###     for one system,  I can just rerun that system's contacts.sh script 
###     and then rerun concat.py, which doesn't require me waste a bunch of time
###     regenerating contact fractions for systems that are fine.

###     POSSIBLE EXTREMLY LONG TERM GOAL:  Can I switch to calculating contact
###     fractions in python using MDAnalysis?

# contacts template



drugs=(  "BB4R" "BB4S" "C06R" "C06S" "E25R" "E25S" "G07R" "G07S"  "H03R"  "H03S"  ) 

for drug in ${drugs[@]}; do

echo "#### #### #### #### ####"
echo ${drug}
echo "#### #### #### #### ####"

mkdir ${drug}
mkdir ${drug}/upper
cd ${drug}/upper
cp ../../contactFreq.tcl ./
cp ../../contactfreq.sh ./
cp ../../process_contacts.py ./
printf "source contactFreq.tcl \n source contactfreq.sh \n exit \n"|vmd ../../../clean/${drug}_upper*frame0.gro ../../../clean/${drug}_upper*short_1ns_cat.xtc -dispdev text
./process_contacts.py  GlyT2 POPC_CHOL ${drug} LAS upper -1

mkdir ../lower
cd ../lower

cp ../../contactFreq.tcl ./
cp ../../contactfreq.sh ./
cp ../../process_contacts.py ./
printf "source contactFreq.tcl \n source contactfreq.sh \n exit \n"|vmd ../../../clean/${drug}_lower*frame0.gro ../../../clean/${drug}_lower*1ns_cat.xtc -dispdev text
./process_contacts.py  GlyT2 POPC_CHOL ${drug} LAS lower 0

cd ../..

done

cd ../../../contacts/pr
./concat.py
