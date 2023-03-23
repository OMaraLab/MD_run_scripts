# how ada set up her workstation

## gromacs 2022 came preinstalled

check version with

~~~s
which gmx
gmx pdb2gmx --version
~~~

## install vmd

go to https://www.ks.uiuc.edu/Research/vmd/alpha/ , downloaded `VMD 1.9.4a55 for RHEL 7+ Linux, 64-bit Intel x86 (x86_64), RTX RTRT`

navigate to download directory

~~~s
tar -xf vmd-1.9.4a55.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185-RTXRTRT.opengl.tar.gz 
cd vmd-1.9.4a55/
vi configure # check that the install path is correct; only important if you have existing vmd installation
./configure
cd src/
sudo make install
which vmd
~~~

if you have a `.vmdrc`file, put it in your home directory

## install slack

go to https://slack.com/intl/en-au/downloads/linux and download an install .deb

navidate to downloaded file using GUI
double click to install

## install vscode

go to https://code.visualstudio.com/download, download .deb file

instructions are here: https://code.visualstudio.com/docs/setup/linux

using terminal, navigate to downloaded .deb file

~~~s
sudo apt install ./vscode_downloaded_file_name.deb
~~~

within vscode, I recommend installing these extensions:

~~~s
Name: Jupyter
Id: ms-toolsai.jupyter
Description: Jupyter notebook support, interactive programming and computing that supports Intellisense, debugging and more.
Version: 2022.7.1102252217
Publisher: Microsoft
VS Marketplace Link: https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter

Name: Python
Id: ms-python.python
Description: IntelliSense (Pylance), Linting, Debugging (multi-threaded, remote), Jupyter Notebooks, code formatting, refactoring, unit tests, and more.
Version: 2022.12.1
Publisher: Microsoft
VS Marketplace Link: https://marketplace.visualstudio.com/items?itemName=ms-python.python

Name: Pylance
Id: ms-python.vscode-pylance
Description: A performant, feature-rich language server for Python in VS Code
Version: 2022.8.40
Publisher: Microsoft
VS Marketplace Link: https://marketplace.visualstudio.com/items?itemName=ms-python.vscode-pylance

Name: gromacs helper
Id: SupernovaZJX.gmx-helper
Description: Make MD easier!
Version: 1.4.1
Publisher: SupernovaZJX
VS Marketplace Link: https://marketplace.visualstudio.com/items?itemName=SupernovaZJX.gmx-helper

Name: indent-rainbow
Id: oderwat.indent-rainbow
Description: Makes indentation easier to read
Version: 8.3.1
Publisher: oderwat
VS Marketplace Link: https://marketplace.visualstudio.com/items?itemName=oderwat.indent-rainbow
~~~

## install zoom

download the .deb from https://zoom.us/download?os=linux

navigate to file, run 

~~~s
sudo apt install ./zoom_downloaded_file_name.deb
~~~

## set up internal 2TB hard drive

guides I referenced:  

 * https://help.ubuntu.com/community/InstallingANewHardDrive
 * https://help.ubuntu.com/community/Fstab#Editing_fstab
 * https://help.ubuntu.com/community/UsingUUID

### first need to format the drive. 

Install partition manager

~~~s
sudo apt-get install gparted
~~~

use partition manager to format the internal drive.  Default file format was ntfs, I changed it to ext3 based on guides

THIS WILL WIPE THE DRIVE!  MAKE SURE YOU ARE WIPING THE DRIVE YOU INTENDED!

launch the partition manager gui with 

~~~s
sudo gparted
~~~

chose the correct drive. For me, it was `/dev/sda`

select the large data partition, and format as the intended file system (I used ext3)

for me this was `/dev/sda2`
there was also a small partition at `/dev/sda1` ; it looks like some kind of microsoft recovery partition; I just left it

Double check everything is correct, then once you aresatisfied click the big green tick to apply the changes.  Takes a few minutes.

once the new partition is ready, exit gparted

### permanently mount the drive to `/store`

make a mount directory

~~~s
sudo mkdir /store
~~~

bet the UUID of the partition you want to mount

~~~s
sudo blkid
~~~

mine was

~~~s
/dev/sda2: UUID="221ae19f-4e10-4c80-98f3-42f30bfad182" SEC_TYPE="ext2" TYPE="ext3" PARTLABEL="Basic data partition" PARTUUID="92a540e3-79cd-4f4e-a175-4263104c34ef"
~~~

add the drive to your fstab 

~~~s
sudo nano /etc/fstab
~~~

for mine I addded this line:

~~~s
UUID=221ae19f-4e10-4c80-98f3-42f30bfad182 /store          ext3    defaults          0       0
~~~

reboot to check it worked

# install conda

instructions:  https://docs.anaconda.com/anaconda/install/linux/

go here, download installer script : https://www.anaconda.com/products/distribution#linux

check installer script

~~~s
shasum -a 256 Anaconda3-2022.05-Linux-x86_64.sh
~~~

run installer

~~~s
bash ./Anaconda3-2022.05-Linux-x86_64.sh 
~~~

you might need to add conda to your path, and then run `conda init yourshell` 

use gui to set up environments

~~~s
anaconda-navigator
~~~

you'll need to add conda forge to the channels

https://ostechnix.com/enable-conda-forge-channel-for-conda-package-manager/

click channels, add `https://conda.anaconda.org/conda-forge/` save, click update index

make a general working environemnt; mine is called mdanalysis

Add the following packages

* pandas
* numpy
* cython
* seaborn
* matplotlib
* mdanalysis
* ipython
* jupyter
* scipy

build the environment.

edit your `.bashrc` or equivalent to load this environment by default

update mdanalysis to use openMP to enable parallel processing (make it go fast)

~~~s
pip install --upgrade MDAnalysis
pip install --upgrade MDAnalysisTests
~~~

## install pymol open source

with the correct conda environment open

ref: https://pymolwiki.org/index.php/Linux_Install

~~~s
sudo apt-get install pymol
conda install pymol
~~~

have to do both otherwise it doesn't work

## install python2 if you need it for things like martinize.py or insane.py

~~~s
sudo apt install python2  
which python2
~~~

## install tmux

~~~s
sudo apt install tmux  
which tmux
~~~

## setup ssh

make a .ssh/config file

~~~s
touch ~/.ssh/config
~~~

mine looks like this

~~~s
Host *
  AddKeysToAgent yes
#  UseKeychain yes
  IdentityFile ~/.ssh/id_rsa

Host gadi
    Hostname gadi.nci.org.au
    User aq8103
    ForwardAgent yes

Host setonix
    Hostname setonix.pawsey.org.au
    User ada
    ForwardAgent yes
~~~

generate an ssh key

~~~s
cd ~/.ssh/
ssh-keygen -b 1024 -t rsa -f id_rsa
~~~

set up simple gadi and setonix login

~~~s
ssh-copy-id ada@setonix

# then input your password

ssh-copy-id aq8103@gadi

# then input your password

alias gadi="ssh aq8103@gadi"  # i think this is right for bash, you might have to check syntax
alias setonix="ssh ada@setonix"
~~~

this should allow you to connect to either supercomputer by running `gadi` or `setonix` from the command line 

## autodock vina


(download vina) [https://vina.scripps.edu/]

extract vina to /packages/vina

symlink into usr/local/bin to add binaries to system path:

~~~s
cd /usr/local/bin
sudo ln -s /packages/vina/bin/vina_split vina_split   
sudo ln -s /packages/vina/bin/vina vina
~~~

(download mgltools) [https://ccsb.scripps.edu/mgltools/downloads/]

change to download directory

~~~s
tar xzvf mgltools_x86_64Linux2_1.5.7p1.tar.gz
cd mgltools_x86_64Linux2_1.5.7/   
bash install.sh -d /packages/mgltools

# symlink into usr/local/bin to add binaries to system path:

cd /usr/local/bin
sudo ln -s /packages/mgltools/bin/adt adt
sudo ln -s /packages/mgltools/bin/pmv pmv
sudo ln -s /packages/mgltools/bin/vision vision
sudo ln -s /packages/mgltools/bin/pythonsh pythonsh
~~~

# alternate gromacs

If you want/need multiple versions of gromacs, I suggest you use module manager

https://modules.readthedocs.io/en/latest/module.html

I needed to install my own gromacs because the version that came with this computer seems to have some problems with numberical instability; I wonder if it was installed through apt-get instead of compiled here

I installed my own version of gromacs 2021.4

download gromacs 2021.4 to `~/Downloads`
untar it

I install alternate gromacs versions in `/packages/gromacs`
I store my modules files in `/packages/modules`

~~~s
tar xfz gromacs-2021.4.tar.gz 
cd gromacs-2021.4
mkdir build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=/packages/gromacs/2021.4/ -DGMX_MPI=off
make
make check
# you need to make sure that gromacs passes all the different tests
# if it doesn't pass the test, DO NOT PROCEED
sudo make install
~~~

you need the right version of cuda for your version of gromacs. hopefully this should be checked during `cmake`, `make`, and `make check`.

Once everything is installed, set up module files in `/packages/modules/gromacs/`

you'll need one for the default install, and one for the new version

default install was version 2022.  the module file was `/packages/modules/gromacs/2022`

~~~s
#%Module

proc ModulesHelp { } {
   puts stderr "This module adds gromacs 2022 to your path"
}

module-whatis "This module adds gromacs 2022 to your path\n"
conflict gromacs

set basedir "/usr/local/gromacs/"

setenv GMXPREFIX "${basedir}"
setenv GMXBIN  "${basedir}/bin"
setenv GMXLDLIB  "${basedir}/lib"
setenv GMXMAN  "${basedir}/share/man"
setenv GMXDATA  "${basedir}/share/gromacs"
setenv GROMACS_DIR "${basedir}"

prepend-path PATH "${basedir}/bin"
prepend-path LD_LIBRARY_PATH "${basedir}/lib"
prepend-path PKG_CONFIG_PATH "${basedir}/lib/pkgconfig"
prepend-path MANPATH "${basedir}/man"
~~~

honestly I should probably just uninstall this at some point

for 2021.4, module file was

~~~s
#%Module

proc ModulesHelp { } {
   puts stderr "This module adds gromacs 2021.4 to your path"
}

module-whatis "This module adds gromacs 2021.4 to your path\n"
conflict gromacs

set basedir "/packages/gromacs/2021.4/"

setenv GMXPREFIX "${basedir}"
setenv GMXBIN  "${basedir}/bin"
setenv GMXLDLIB  "${basedir}/lib"
setenv GMXMAN  "${basedir}/share/man"
setenv GMXDATA  "${basedir}/share/gromacs"
setenv GROMACS_DIR "${basedir}"

prepend-path PATH "${basedir}/bin"
prepend-path LD_LIBRARY_PATH "${basedir}/lib"
prepend-path PKG_CONFIG_PATH "${basedir}/lib/pkgconfig"
prepend-path MANPATH "${basedir}/man"
~~~

used a similar process to build gmx 5.1.5 for some old libraries, comake command was

~~~s
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DCMAKE_INSTALL_PREFIX=/packages/gromacs/5.1.5/ -DREGRESSIONTEST_PATH=/home/uqadaqu1/Downloads/regressiontests-5.1.5/ -DGMX_GPU=off
~~~

I'm not sure to what degree this version actually works; it's not very compatible with my current system, but hopefully it has functional libraries for the programs that will use it. 