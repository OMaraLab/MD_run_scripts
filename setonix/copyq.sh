#!/bin/bash --login
#SBATCH --partition=copy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --account=pawsey0420

rsync -vhsrl --chmod=Dg+s -e ssh path/to/source user@data.qriscloud.org.au:/path/to/dest/
