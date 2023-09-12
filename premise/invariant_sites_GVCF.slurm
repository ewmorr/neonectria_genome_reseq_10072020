

  process VariantCallingGVCF {

    label "spandx_gatk"
    tag { "$id" }
    //publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: true, pattern: '*.gvcf'

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") from variantcallingGVCF_ch

    output:
    set id, file("${id}.raw.gvcf")
          file("${id}.raw.gvcf") into gvcf_files


    """
    gatk HaplotypeCaller -R ${reference} --ploidy 1 -ERC GVCF --I ${id}.dedup.bam -O ${id}.raw.gvcf
    """
  }


  process Master_vcf {
    label "master_vcf"
    publishDir "./Outputs/Master_vcf", mode: 'copy', overwrite: true

    input:
    file("*.raw.gvcf") from gvcf_files.collect()
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1

    output:
    set file("out.filtered.vcf"), file("out.vcf") into snp_matrix_ch

    script:
    """
    bash Master_vcf.sh ${reference.baseName}
    gatk VariantFiltration -R ${reference} -O out.filtered.vcf -V out.vcf \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "QD < $params.QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "HaplotypeScoreFilter"
    """

  }


#!/bin/bash

ref=$1

echo "Creating master VCF file"

array=($(find *.gvcf -printf "%f "))
n="${#array[@]}"
array2=("${array[@]/#/-V }")
gatk CombineGVCFs -R ${ref}.fasta ${array2[*]} -O master.vcf
gatk GenotypeGVCFs -R ${ref}.fasta -V master.vcf -O out.vcf
 
exit 0


#Call new GVCFs with haplotype caller
#call new master GVCF (with combineGCVFs -- note there are currently two master.vcf files that must be from separate runs which explains the multiple sets of gvcf files... need to recall to ensure we have the right run (with the correct dedeup.bams)
#genotype the new master using the invariants call here https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html (i.e., with the --all-sites flag)
#perform variant filtration as above using the call following the master_vcf bash script call



#!/usr/bin/bash

#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --array=0-116%10

module purge
module load anaconda/colsa

conda activate spandx

ref_dir=/mnt/home/garnas/ericm/Nf_SPANDx_all_seqs
out_dir=/mnt/home/garnas/ericm/Nf_invariant_sites_GVCF

reference=$ref_dir/ref.fasta

mkdir $out_dir/indv_GVCFs

array=($(find $ref_dir/Outputs/bams/ -name *.dedup.bam -printf "%f "))
id_arr=("${array[@]%.dedup.bam}")

#indv GVCFs
gatk HaplotypeCaller -R ${reference} --ploidy 1 -ERC GVCF \
    --I $ref_dir/Outputs/bams/${id_arr[$SLURM_ARRAY_TASK_ID]}.dedup.bam \
    -O $out_dir/indv_GVCFs/${id_arr[$SLURM_ARRAY_TASK_ID]}.raw.gvcf

$SLURM_ARRAY_TASK_ID
