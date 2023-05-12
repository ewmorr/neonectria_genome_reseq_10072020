#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="admixture"
#SBATCH --output=admixture.out
#SBATCH --partition=shared
#SBATCH --exclude=node117,node118

module purge
module load linuxbrew/colsa

cd ~/Nd_SPANDx_all_seqs_admixture
#prev had up to 15 for Nf

#Note that there is some weird stuff going on with the flags. Seems like rodering may be important and the number of dahses before dif commands has changed since we first ran, but haven't done a lot of troubleshooting.
for K in 1 2 3 4 5 6 7 8 9 10 11 12
do(
    srun admixture ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.ped $K -haploid="*" --cv | tee log$K.out
)
done

grep -h CV log*.out > CV_by_K.text
