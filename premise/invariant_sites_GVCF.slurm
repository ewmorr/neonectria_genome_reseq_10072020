

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
#SBATCH --job-name="gatk"
#SBATCH --output=gatk.out
#SBATCH --partition=shared
#SBATCH --ntasks=1 --cpus-per-task=1

module purge
module load anaconda/colsa

conda activate spandx

out_dir=/mnt/home/garnas/ericm/Nf_invariant_sites_GVCF
in_dir=/mnt/home/garnas/ericm/Nf_invariant_sites_GVCF/indv_GVCFs

reference=/mnt/home/garnas/ericm/Nf_SPANDx_all_seqs/ref.fasta

#SPANDx gatk configs
    CLUSTER_SNP=3
    CLUSTER_WINDOW_SNP=10
    MLEAF_SNP=0.95
    QD_SNP=10.0
    MQ_SNP=30.0
    FS_SNP=60.0
    QUAL_SNP=30.0
    LOW_DEPTH=2
    HIGH_DEPTH=3

  MLEAF_INDEL = 0.95
  QD_INDEL = 10.0
  MQ_INDEL = 30.0
  FS_INDEL = 200.0
  QUAL_INDEL = 30.0
  LOW_DEPTH_INDEL = 2
  HIGH_DEPTH_INDEL = 3
#end

array=($(find $in_dir/ -name *.gvcf -printf "%f "))
array2=("${array[@]/#/-V }")

gatk CombineGVCFs -R $reference ${array2[*]} -O $out_dir/master.invariant_sites.vcf
gatk GenotypeGVCFs -R $reference -V $out_dir/master.vcf -all-sites -O $out_dir/out.invariant_sites.vcf

gatk VariantFiltration -R ${reference} -O out.invariant_sites.filtered.vcf -V out.invariant_sites.vcf \
    --cluster-size CLUSTER_SNP -window CLUSTER_WINDOW_SNP \
    -filter "QD < QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > FS_SNP" --filter-name "HaplotypeScoreFilter"


