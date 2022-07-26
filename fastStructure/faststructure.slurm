#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="faststructure"
#SBATCH --output=faststructure2.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate faststructure

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/

#serial calls of faststructure because `wait` seemingly doesn't work with for loop in slurm

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 1 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 2 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 3 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 4 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 5 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 6 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 7 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 8 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 9 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 10 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 11 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 12 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 13 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 14 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

~/.conda/envs/faststructure/bin/python ~/.conda/envs/faststructure/bin/structure.py -K 15 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure/simple \
--full --seed=100 &

wait
