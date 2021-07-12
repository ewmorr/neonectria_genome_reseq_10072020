#!/usr/bin/bash

module purge
module load anaconda/colsa

python --version
pip --version

pip3 install structure_threader --user

export PATH=~/.local/bin:$PATH

