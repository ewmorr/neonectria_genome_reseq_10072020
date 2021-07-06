#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="faststructure"
#SBATCH --output=faststructure_logistic.7.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate faststructure
export PATH=~/.conda/envs/faststructure/bin:$PATH

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/

#serial calls of faststructure because `wait` seemingly doesn't work with for loop in slurm

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 7 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 


