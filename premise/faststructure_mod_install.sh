#!/bin/bash


#Instructions to install a modified version of faststructure to deal with the issue of nan values during logistic prior runs

#First clone the faststructure repo localle
cd repo
cd faststructure
git remote set-url origin http://github.com/ewmorr/faststructure
#create repo on github.com
git push -u origin master

#The cloned copy of the script faststructure.pyx is modied at L68-70, L74 (moved iter advancement down below new if statement), and L146-148

module purge
module load anaconda/colsa

conda create --name faststructure_mod -c bioconda python=2.7 numpy=1.16.5 cython=0.27.3 scipy=1.2.1 pip git

conda activate faststructure_mod

conda install -c conda-forge gsl

conda deactivate
conda activate faststructure_mod

which python
echo $PATH
#check to make sure that correct version of python is being used and .conda/envs/faststructure_mod/bin is in the path
#export PATH=~/.conda/envs/faststructure_mod/bin:$PATH
#This doesn't appear to be necessary after deactivating and reactivating the env (the main conda bin shows up after instead of before the .local after deactivate and activate)

#pip install git+git://github.com/ewmorr/faststructure@master
#The pip install doesn't work because of no python 2.7 support

#Try the install prescribed in docs

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/home/garnas/ericm/.conda/envs/faststructure_mod/lib/
export CFLAGS="-I/mnt/home/garnas/ericm/.conda/envs/faststructure_mod/include"
export LDFLAGS="-L/mnt/home/garnas/ericm/.conda/envs/faststructure_mod/lib/"

source ~/.bashrc

cd ~/repo
git clone https://github.com/ewmorr/faststructure
cd ~/repo/faststructure/vars/
python setup.py build_ext --inplace

cd ..
python setup.py build_ext -f --inplace

python ~/repo/faststructure/structure.py
