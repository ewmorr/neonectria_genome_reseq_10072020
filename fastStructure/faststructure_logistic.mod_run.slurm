#!/usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="faststructure"
#SBATCH --output=faststructure_logistic.mod_run.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

conda activate faststructure_mod

export PATH=~/.conda/envs/faststructure_mod/bin:$PATH

cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/

#serial calls of faststructure because `wait` seemingly doesn't work with for loop in slurm

python ~/repo/faststructure/structure.py -K 1 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 2 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 3 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 4 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 5 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 6 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 7 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 8 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 9 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 10 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 11 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 12 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 13 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 14 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

python ~/repo/faststructure/structure.py -K 15 \
--input=Nf.out.filtered.LD_filtered_0.5_10Kb \
--output=faststructure_mod/logistic \
--full --seed=100 --prior=logistic --tol=1e-4 &

wait

