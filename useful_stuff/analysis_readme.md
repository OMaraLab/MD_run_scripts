# QacA Analysis - Dr. Heather Aitken 
Feb 6 2023
* The following uses GROMACS, CPPTRAJ (Amber Tools), VMD and some web servers
* Each analysis type can be found in its own directory and the concatonated trajectories (+ tpr and gro) can be found in the "trajectories" directory
* * How to concatonate (join) a trajectory: 

``` gmx trajcat -f run1.xtc run2.xtc run3.xtc -o run_all.xtc -s run1.tpr -skip 10 ```

* I made a super small trajectory (skip 10 of a previous skip 5) for speedy analysis, but I wouldn't do this as standard practice:

``` for i in `cat list `; do printf '0\n' | gmx trjconv -f $i\_sk5.xtc -skip 10 -o $i\_sk5\_sk10.xtc -s $i.tpr -pbc whole ; done ```

* The printf allows me to specify my index file option - in this case the 'System' which is index 0 in the default index file. If you want to do anything bespoke you will likely want to create an index file: 

``` gmx make_ndx -f run1.tpr -o whatever.ndx ``` and follow the prompts

### You may also consider this readme a brief introduction to the concept of for loops and why they are so powerful ... and why you should make your file names consistent... ;-)

- "list" is a list of all traj names 
- "drugs" is a list of all of the sims that have substrates 
- I make a list by doing the following ls *gro > list 
- - This will search for everything with the extension gro and create a list of these in a file called "list" 

## RMSD 

``` for i in `cat ../trajectories/list `; do printf '4\n 1\n' | gmx rms -f ../trajectories/$i\_sk5.xtc -s ../trajectories/$i.tpr -xvg none -o $i-p.xvg  ; done ```

* This will do the protein RMSD relative to the backbone. -xvg none gets rid of all the header stuff to make easier to load into python - including to then run the inhouse script stat.py (gives mean and SD)

``` for i in `cat ../trajectories/drugs `; do printf '4\n 13\n' | gmx rms -f ../trajectories/$i\_sk5.xtc -s ../trajectories/$i.tpr -xvg none -o $i-d.xvg  ; done ```

* This will do the drug RMSD relative to the backbone of the protein. I.e. how much does the drug move over the course of the trajectory relative to the backbone of the protein

## Cluster Analysis 

* Using RMSD of 2 nm 

``` for i in `cat ../trajectories/list `; do gmx cluster -f ../trajectories/$i\_sk5.xtc -s ../trajectories/$i.tpr -g $i.log -cl $i-clusters.pdb -n ../trajectories/$i.ndx -cutoff 0.2 ; done ```



## Solvent Accessible Surface Area (SASA)

``` for i in `cat ../trajectories/drugs `; do printf '13\n' | gmx sasa -f ../trajectories/$i\_sk5\_sk10.xtc -s ../trajectories/$i.tpr -xvg none -o $i.xvg -n ../trajectories/$i.ndx ; done ```

## RMSF 

``` for i in `cat ../trajectories/list `; do printf '1\n' | gmx rmsf -f ../trajectories/$i\_sk5\_sk10.xtc -s ../trajectories/$i.tpr -xvg none -o $i.xvg -oq $i-bf.pdb -res ; done ```

* The -oq pdb is a heatmap of RMSF that you can visual in VMD as a surface on the protein. Make sure to change the range to something sensible. 

* The -res averages over the residue rather than calculating per atom (messy) 

## Contacts in VMD (in this case drug-protein contacts)

``` for i in `cat ../trajectories/list `; do printf "\n source contactFreq.tcl \n contactFreq {protein} {resname ETHI CHLX CLHX} 4 0 contact-d-"$i".out \n exit \n" | vmd ../trajectories/$i.gro ../trajectories/$i\_sk5_sk10.xtc -dispdev text ; done ```

- ETHI CHLX and CLHX are the names of my drug residues in this case

- change 4 0 to change the bond length and percentage persistance assessed 

## HBONDS in AmberTools

- make a pdb from the gro so that Amber can read it 

``` for i in `cat list `; do printf '0\n' | gmx trjconv -f $i.gro -o $i.pdb -s $i.tpr ; done ```

- iterate through with cpptraj 

``` for i in `cat ../trajectories/list `; do printf "parm ../trajectories/"$i".pdb \n trajin ../trajectories/"$i"\_sk5_sk10.xtc \n autoimage \n hbond "$i"-prot out hbond-"$i".out angle 120 dist 3.4 donormask :1-514@N*,O* acceptormask :1-514@N*,O* avgout hbond-"$i".dat \n go \n quit \n" | cpptraj ; done ```

## Salt Bridges 

- same as for Hbonds but just specify the residues of interest 

``` for i in `cat ../trajectories/list `; do printf "parm ../trajectories/"$i".pdb \n trajin ../trajectories/"$i"\_sk5_sk10.xtc \n autoimage \n hbond "$i"-prot-salt out salt-"$i".out angle 120 dist 3.4 donormask :ARG,LYS,HIS acceptormask :ASP,GLU avgout salt-"$i".dat \n go \n quit \n" | cpptraj ; done ```

## Electrostatic Potential visualisation 

- We use an online tool for this one https://server.poissonboltzmann.org/ 
- First you need to make a pqr with the PDB2PQR tool and then you submit for APBS calculation 
- Download the pqr and dx files and rename them appropriately 
- You can visualise these in VMD 
* * Load the pqr (essentially a pdb) and then load the dx into that 
* * Colour by volume and view as volume slice
* * Chimera does this better than VMD (a tutorial in and of itself)

## Cavity Size

- A web server that can do this statically (a single pdb) is http://sts.bioe.uic.edu/castp/index.html
- However to do it over the coarse of a sim do it we need to use trj_cavity which means you need to install a OLD version of gromacs (< gromacs 5) 
* * Fun... 