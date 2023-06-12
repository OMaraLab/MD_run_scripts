# MD_run_scripts

Standard scripts, MDP files and other documentation for running MD simulations

The most up to date version of this will always be hosted on github at [https://github.com/OMaraLab/MD_run_scripts](https://github.com/OMaraLab/MD_run_scripts).  This repository is private, so you need to be logged in to github, and be a member of the omara organisation.

## Other resources

* [Gadi user guide](https://opus.nci.org.au/display/Help/Gadi+User+Guide)
* [Pawsey user guide](https://support.pawsey.org.au/documentation/display/US/User+Support+Documentation)
* [UQ RCC user guide](https://github.com/UQ-RCC/hpc-docs)
* [CTCMS guide to working with linux and supercomputers](https://ctcms-uq.github.io/)
  * If these instructions do not make sense to you, read the CTCMS guide first.

## Submission Scripts

`bunya/`, `weiner/`, `setonix/`, `topaz/` and `gadi/` contain standard submission scripts for running jobs on the super computers

Chose a script appropriate for your system, considering the size of your system, the gromacs version you wish to use, and how much compute resources are appropriate.  Most of the submission scripts have been tested and have shown adequate performance, but you should always check the performance of your own system.

These scripts will submit a job to run for a specific block of time, usually either four hours or ten hours.  After every submission, they will check if the job is done by comparing the `.log` and `.mdp` files.  If the job is not complete, they will automatically resubmit for another go.

Automatic resubmission is great, but don't assume everything is fine.  Always keep an eye on your jobs, and check at least once a day to check whether they've crashed.

Choose the appropriate script for your job, and set up a production folder for each replicate.

## MDP files

MDP files for typical system setups will be stored in `/mdp`

Do not blindly choose an `.mdp` file.  Always check if it is appropriate to your project.

## Useful stuff

`useful_stuff/` contains several documents that are generally useful, including:

* a barebones guide to setting up a workstation
* a `.vmdrc` file you can use to set up a good default vmd configuration.
* a tutorial for generic protein analysis
* `.bash_aliases` files that contain settings that are useful when working with pbs or slurm clusters.
* an example of a script for processing production for analysis, using `gmx trjconv`
* an example of the scripts you might use to look at contact fractions in vmd

## The most important thing

**Q:** Hey I deleted a file on the supercomputer by mistake, how do I get it back?

**A:** Unless you are absurdly lucky, you probably can't.

Be careful!  Only delete things when you are 100% sure you are done with them, and that you will never need them again.  Always double check your deletion commands before you use them and make sure they are going to do what you expect.  **Always remember: If you delete something and you don't have a backup, there is a good chance it is *gone!***

## How to get a job running on a computing cluster

### Step one:  Set up a folder for your production run

For any gereric MD system, your production folder should look something like `Template-systems/GlyT2_apo_prod`

If your production run is named GlyT2_apo_POPC_CHOL, you MUST have five things in your production folder

* A submission script named `GlyT2_apo_POPC_CHOL`
  * (**NB:** not `GlyT2_apo_POPC_CHOL.sh`, `GlyT2_apo_POPC_CHOL.slurm` or `GlyT2_apo_POPC_CHOL.pbs`)
* An mdp file named `GlyT2_apo_POPC_CHOL.mdp`
* An index file named `GlyT2_apo_POPC_CHOL.ndx`
* A topology file named `GlyT2_apo_POPC_CHOL.top`
* Starting coordinates named `GlyT2_apo_POPC_CHOL_start.gro`  
  * (**NB:** Note the different filename)

You must also have your forcefield folder and all `.itp` files in the locations referenced by your `.top` file.

If any of these files are missing, our generic submission scripts will not function correctly.

### Step two:  Copy the system to the supercomputer / cluster you wish to use

You can use `scp` or `rsync`

example scp syntax to copy from a local machine to a remote machine:

`$ scp -r path/to/thing_to_copy username@remote.computer.address:/path/to/place-to-put-it`

eg:  if I wanted to copy the system GlyT2_apo_POPC_CHOL from my local machine to gadi, I might use this command:

`$ scp -r /store/users/ada/GlyT2/buildprod/GlyT2_apo_POPC_CHOL aq8103@gadi.nci.org.au:/scratch/q95/aq8103/GlyT2/`

Where you want to put it depends on which cluster you are using.  You will want to put it somewhere inside your user folder.  Make sure your files are structured logically, and that other people in our group can access them.

|Cluster     | Cluster address | Typical path to an appropriate user folder   |
|Bunya (RCC)      |`bunya.rcc.uq.edu.au`     |`/scratch/user/your_username/`   |
|Wiener(RCC)      |`wiener.hpc.dc.uq.edu.au` |`/scratch/aibn/your_username/` |
|Gadi (NCI)       |`gadi.nci.org.au`         |`/scratch/q95/your_username/`          |
|Setonix (Pawsey) | `setonix.pawsey.org.au`  |`/scratch/pawsey0420/your_username/`   |
|Topaz (Pawsey)   | `topaz.pawsey.org.au`    |`/group/pawsey0420/your_username/`   |

**Make sure you read the user guide for the cluster you are using.  Some folders are intended only for short term use!  Do not assume files on the clusters will remain there forever.**

### Step three:  Log into the remote system andque your submission script

Log in to the cluster with ssh

|Bunya (RCC)      |`$ ssh username@bunya.rcc.uq.edu.au`   |
|Wiener(RCC)      |`$ ssh username@wiener.hpc.dc.uq.edu.au` |
|Gadi (NCI)       |`$ ssh username@gadi.nci.org.au`         |
|Setonix (Pawsey) |`$ ssh username@setonix.pawsey.org.au`   |
|Topaz (Pawsey)   |`$ ssh username@topaz.pawsey.org.au`   |

Go to the folder for your first replicate.  

`$ cd /path/to/system/r1'`

Jobs on computing lcusters usually use queuing systems.  There are two common systems, PBS (used by gadi) and slurm (used everywhere else) 

Different clusters use 

submit your job by queueing the submission script

|PBS queing system (gadi) | `$ qsub the_name_of_your_submission_script` |
|Slurm queing system (pawsey and UQ clusters) | `$ sbatch the_name_of_your_submission_script` |

If you don't get an error, it should be queued.

Our submission scripts should run for a specified number of hours, after which they check if the simulation is complete.  If the simulation is not complete, the job **should** automatically resubmit itself.  You still need to keep an eye on it.  Sometimes resubmission fails.  Sometimes clusters go down.  Sometimes a job crashes because something is wrong.  If you have jobs running on a cluster, you should check once a day to see if they are ok.

You can use these commands to check the status of all queued jobs.

|PBS queing system (gadi) | `$ qstat` |
|Slurm queing system (pawsey and UQ clusters) | `$ squeue --me` |

You can always check the end of your gromacs log file to see how far your simulation has gotten, and how quickly it is running.   The queuing systems also generate error logs for every run.  On pbs systems, these will be called something like `systemname.e72234683`. On slurm systems these will be called something `slurm-83245234.out`.  If your job has stopped and you don't know why, read the most recent error log and see what happened.

You can these log files (and any other text file) through the command line using the program `less`.  Run it with ` $ less filename.log`.  Press `G` to scroll to the end of the file.  When you are done reading, press `q` to quit `less`

It is also usually a good idea to look at your system in VMD, to see if everything looks normal.  To do this, use scp to download your system to your local machine, and look the it in VMD.  

`$ scp -r username@remote.computer.address:path/to/system_that_is_running /path/to/local-place-to-put-it`

### Step four:  Once your job is finished

Download the system to your local machine.  Check it is actually finished, and that it ran the full time you expected.  Get rid of the PBS or Slurm error logs.  Then, make a backup to the RDM before you do **anything** else.

Congratulations, you have data!
