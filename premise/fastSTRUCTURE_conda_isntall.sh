#!/bin/bash

module purge
module load anaconda/colsa

conda create --name faststructure -c bioconda python=2.7 faststructure

conda activate faststructure

#For some reason the path in $PATH is set to ~/.conda/envs/faststructure/local/bin but there is no such directory, instead ~/.conda/envs/faststructure/bin

export PATH=~/.conda/envs/faststructure/bin:$PATH
#This still isn't working to call the main python script ...

#clone the git repo to get the test data
mkdir fast_structure_git
cd fast_structure_git/
git clone https://github.com/rajanil/fastStructure

cd ..
#test
#run the commands at http://rajanil.github.io/fastStructure/ and compare to the test files in the git repo
mkdir test

python ~/.conda/envs/faststructure/bin/structure.py -K 3 --input=fastStructure/test/testdata --output=test/testoutput_simple --full --seed=100
python ~/.conda/envs/faststructure/bin/structure.py -K 3 --input=fastStructure/test/testdata --output=test/testoutput_logistic --full --seed=100 --prior=logistic

