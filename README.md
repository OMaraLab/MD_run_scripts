# MD_run_scripts

Standard scripts, MDP files and other documentation for running MD simulations

## Useful stuff

`useful_stuff/` contains a barebones guide to setting up a workstation, and broadly useful files like the group's `.vmdrc` files.

## MDP files

MDP files for standard system setups will be stored in `/mdp`

Do not blindly choose an `.mdp` file.  Always check if it is appropriate to your project.

## Automatic resubmission scripts

`bunya/`, `weiner/`, `setonix/`, `topaz/` and `gadi/` contain standard submission scripts for running jobs on the super computers

These scripts will submit a job to run for a specific block of time, usually either four hours or ten hours.  After every submission, they will check if the job is done by comparing the `.log` and `.mdp` files.  If the job is not complete, they will automatically resubmit for another go.

Automatic resubmission is great, but don't assume everything is fine.  Always keep an eye on your jobs, and check at least once a day to check whether they've crashed.

Choose the appropriate script for your job, and set up a production folder for each replicate.

## Production folder setup

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

