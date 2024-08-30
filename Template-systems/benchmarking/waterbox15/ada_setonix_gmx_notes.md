# Setonix gromacs notes from Ada

## my custom gromacs builds

I have made mixed precision mpi builds of gromacs 2023.0 and gromacs 2023.2.  
Both are stored in `/software/projects/pawsey0420/manual/gromacs/`

I cobbled together some crude module files for these in `/software/projects/pawsey0420/modulefiles`

Both were assembled using the same process:

Download gromacs source code and regression tests

Start interactive session on one of the work nodes with 16 cores

`module load gromacs/2023` to inherit the default environment

`tar xfz gromacs-source-tarball.tar.gz`
`tar xfz regressiontest-tarball.tar.gz`

``` sh
cd /path/to/gromacs_source
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON \
    -DREGRESSIONTEST_DOWNLOAD=OFF -DREGRESSIONTEST_PATH=path/to/regressiontests/ \
    -DCMAKE_INSTALL_PREFIX=path/to/install_dir \
    -DGMX_MPI=on -DGMX_HWLOC=ON -DGMX_GPU=OFF
make -j 16
make check -j 16
make install
```

## benchmarking system

The benchmarking system is four files

* `submit_minimal.slurm` - example submission script to run the system for 15 minutes
* `waterbox15_EQ.gro` - starting coordinates
* `waterbox15.top` - starting topology
* `waterprod.mdp` - starting control file

To make the system run for fifteen minutes you need to specify -maxh 0.25 in the gmx mdrun command in the submission script.  This is specified in the example.  If you don’t, it will run forever.

There are two ways we run benchmarking.

Method one:  Take a given template system and runb with various node, mpi and gromacs configurations to find the optimal settings
Method two:  Take a given node, mpi and gromacs configuration, and run it with multiple systems of increasing sizse to examine how the configuration scales with larger problems

This test system is for method one.  As we proceed with this project and identify optimal configurations, we may wish to switch to method 2 to evaluate how these configurations scale.
