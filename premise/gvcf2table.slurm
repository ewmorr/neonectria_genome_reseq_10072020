#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="snp2table"
#SBATCH --output=snp2table.out
#SBATCH --partition=shared

module purge
module load linuxbrew/colsa

cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/


echo "Creating VCF tables"
gatk VariantsToTable -V out.filtered.LD_filtered_0.5_10Kb.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -O out.filtered.LD_filtered_0.5_10Kb.table
gatk VariantsToTable -V out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -O out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.table

echo "Clean table"
echo "Creating SNP matrix"
awk '$5 ~/SNP/' out.filtered.LD_filtered_0.5_10Kb.table | awk '$4 !~/\.,\*/' | grep -v '\./\.' > out.filtered.LD_filtered_0.5_10Kb.table.clean
awk '$5 ~/SNP/' out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.table | awk '$4 !~/\.,\*/' | grep -v '\./\.' > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.table.clean
