#!/usr/bin/bash
#SBATCH --job-name="gatk"
#SBATCH --output=gatk.out
#SBATCH --partition=shared
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --array=0-32%24

module purge
module load anaconda/colsa

conda activate spandx

#mkdir Nd_invariant_sites_GVCF
#mkdir Nd_invariant_sites_GVCF/indv_GVCFs

#change array number appropriately
#find Nd_SPANDx_all_seqs/ -name *.gvcf | grep NG | wc -l

ref_dir=/mnt/home/garnas/ericm/Nd_SPANDx_all_seqs
out_dir=/mnt/home/garnas/ericm/Nd_invariant_sites_GVCF

reference=$ref_dir/ref.fasta

gatk CreateSequenceDictionary -R $reference
samtools faidx $reference

array=($(find $ref_dir/Outputs/bams/ -name *.dedup.bam -printf "%f "))
id_arr=("${array[@]%.dedup.bam}")

#indv GVCFs
gatk HaplotypeCaller -R ${reference} --ploidy 1 -ERC GVCF \
    --I $ref_dir/Outputs/bams/${id_arr[$SLURM_ARRAY_TASK_ID]}.dedup.bam \
    -O $out_dir/indv_GVCFs/${id_arr[$SLURM_ARRAY_TASK_ID]}.raw.gvcf

#gatk HaplotypeCaller -R ${reference} --ploidy 1 -ERC GVCF \
#    --I $ref_dir/Outputs/bams/${id_arr[1]}.dedup.bam \
#    -O $out_dir/indv_GVCFs/${id_arr[1]}.raw.gvcf
