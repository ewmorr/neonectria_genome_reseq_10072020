#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="faststructure"
#SBATCH --output=faststructure_logistic.8.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate faststructure_mod
export PATH=~/.conda/envs/faststructure_mod/bin:$PATH

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/

python ~/repo/faststructure/structure.py -K 8 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/logistic \
--full --seed=100 --prior=logistic --tol=1e-5
