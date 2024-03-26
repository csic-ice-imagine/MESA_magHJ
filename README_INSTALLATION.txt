Installation instructions for MESA.

MESA Website: docs.mesastar.org/en/release-r23.05.1/

1. Download MESA version r15140.
File: mesa-r15140.zip
Webpage: zenodo.org/records/4311514

2. Download the compatible MESA SDK version.
File: mesasdk-x86_64-linux-21.4.1.tar.gz
Webpage: user.astro.wisc.edu/~townsend/static.php?ref=mesasdk#Download

3. Create dedicated MESA installation directories.
mkdir MESA
mkdir MESA/version_r15140

4. Copy downloaded files and extract.
SDK:
mv mesasdk-x86_64-linux-21.4.1.tar.gz ~/MESA/version_r15140/
gunzip mesasdk-x86_64-linux-21.4.1.tar.gz
tar xvf mesasdk-x86_64-linux-21.4.1.tar

MESA:
mv mesa-r15140.zip ~/MESA/version_r15140/
unzip mesa-r15140.zip

5. Install any missing libraries listed in the prerequisites table (cf. the SDK webpage).
sudo apt install binutils
sudo apt install make
sudo apt install perl
sudo apt install libx11-6
sudo apt install libx11-dev
sudo apt install zlib1g
sudo apt install zlib1g-dev
sudo apt install tcsh

You also need numpy:
sudo apt install python3-numpy

6. Export paths.
Edit .bashrc and add the following lines:

# MESA path.
export MESA_DIR=~/MESA/version_r15140/mesa-r15140

# Set OMP_NUM_THREADS to be the number of cores on your machine.
export OMP_NUM_THREADS=2

# SDK.
export MESASDK_ROOT=~/MESA/version_r15140/mesasdk
source $MESASDK_ROOT/bin/mesasdk_init.sh

Initialize .bashrc (in all open terminals):
source .bashrc

7. Install MESA.
cd ~/MESA/version_r15140/mesa-r15140
./install

The installation may take 15 min to 1 hour.
If successful, the following message will appear:

************************************************
************************************************
************************************************

MESA installation was successful

************************************************
************************************************
************************************************
