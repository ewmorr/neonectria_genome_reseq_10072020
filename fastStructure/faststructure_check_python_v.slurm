#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="faststructure"
#SBATCH --output=faststructure_version_check.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate faststructure

export PATH=~/.conda/envs/faststructure/bin:$PATH

python --version

python -c   "import pkg_resources; print(pkg_resources.get_distribution('numpy').version)"

python structure.py
