to run jobs:

edit the submission script r1/Run_example.slurm to use the correct account for your groups

copy the whole directory to bunya

got to path/to/Minimal_example/r1/

run the submission script with sbatch Run_example.slurm (ie, sbatch and then the filename of the submission script)

it should run for 10000 steps and then stop.

Performance is not optimised in this example script, it might be slow or wasteful.