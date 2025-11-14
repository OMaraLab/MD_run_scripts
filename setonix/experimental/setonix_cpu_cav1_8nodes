!/bin/bash
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks=484
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --account=m72
#SBATCH --exclusive
#SBATCH --partition=work
#SBATCH --exclude=nid00[2024-2055],nid00[2792-2823]

export OMP_NUM_THREADS=2

module unload gcc-native/14.2
module switch pawseyenv/2025.08 pawseyenv/2024.05
module load gcc/12.2.0
module load libfabric/1.15.2.0
module load gromacs/2023.3-mx33m5h

GMX='gmx_mpi'
#GMXMDRUN='srun -n 242 -m block:block:block gmx_mpi mdrun -ntomp 2 -maxh 23.95 -dd 11 11 2' #535.179ns/day at 39.15% efficiency
GMXMDRUN='srun -n 484 -m block:block:block gmx_mpi mdrun -ntomp 2 -pin on -maxh 0.08 -dd 11 11 4' # estimating 991.921ns/day at 39.25% efficiency 

# expected cost for 15us
# normally system with 242 CPUs production run at ~535 ns/day, so if we multiply by 1.854 (per scaling seen in tests) we should get 991.9 ns/day
#-n 128: 15000/308 = 48.7 days = 1169 hours -> 299,264 SUs (exlusive memory)
#-n 242: 15000/535 = 28.0 days = 673 hours -> 344,576 SUs (exlusive memory) = total 1,033,728 SUs for 3 replicates
#-n 484: 15000/992 = 15.0 days = 360 hours -> 368,640 SUs (exlusive memory) = total 1,105,920 SUs for 3 replicates

errexit ()

{
    errstat=$?
    if [ $errstat != 0 ]; then
        # A brief nap so slurm kills us in normal termination
        # Prefer to be killed by slurm if slurm detected some resource excess
        sleep 5
        echo "Job returned error status $errstat - stopping job sequence $SLURM_JOB_NAME at job $SLURM_JOB_ID"
        exit $errstat
    fi
}

if [ ! -f cav1_epithelial_rN_prod.tpr ]; then
    $GMX grompp -f ../production_pr5.mdp -c 27_cav1_epithelial_rN_eq20fs_pr0_renamedlipids_resnr.gro -o cav1_epithelial_rN_prod.tpr -p cav1_epithelial.top -n cav1_epithelial.ndx -maxwarn 2 &> cav1_epithelial_rN_prod_grompp_${SLURM_JOB_ID}.txt
fi

$GMXMDRUN -v -deffnm cav1_epithelial_rN_prod -cpi cav1_epithelial_rN_prod.cpt || errexit


# Check the log file for the number of steps completed
steps_done=`perl -n -e'/Statistics over (\d+) steps using (\d+) frames/ && print $1' cav1_epithelial_rN_prod.log`
# Check the mdp file for the number of steps we want
steps_wanted=`perl -n -e'/nsteps\s*=\s*(\d+)/ && print $1' ../production.mdp`
# Resubmit if we need to
if (( steps_done < steps_wanted )); then
    echo "Job ${SLURM_JOB_NAME} terminated with ${SLURM_JOB_NAME}/${SLURM_JOB_NAME} steps finished."
    echo "Submitting next job in sequence ${SLURM_JOB_NAME}."
    sbatch ${SLURM_JOB_NAME}
fi

