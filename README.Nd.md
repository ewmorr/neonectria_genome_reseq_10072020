# Repo for analysis of genome resequencing data of Neonectria faginata and Neonectria ditissima
## This workflow picks up after SNP calling and is modified from the main README to be performed on the Nd SNP calls

#### Workflows performed on UNH Premise or locally using R as appropriate


## SNP and sample filtering

### Filter out SNPs that SPANDx FAIL filtering
```
cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf
grep "##\|#\|PASS" out.filtered.vcf > out.filtered.PASS.vcf
grep -v "##\|#" out.filtered.PASS.vcf | wc -l
```
1,112,323 SNPs remaining x 19 samples = 21,134,137
```
grep -o "\s\.:" out.filtered.PASS.vcf | wc -l
```
1182518 NA sites
1182518/21134137 = 5.595%


### Post-SNP calling calculations and filtering

#### Note that while the documentation/paper says that SPANDx performs depth filtering, it is not clear from the scripts that this is done (and there are many SNPs with DP == 1 in the filtered VCF)

### Coverage calculations
calculate raw coverage
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/sample_coverage.slurm
```

concatenate coverage results
```
cd ~/Nd_SPANDx_all_seqs/Outputs/bams
for i in *.coverage_by_sequence.txt
do
    cov="$(grep "genome" $i| cut -f 3)"
    sample=${i%.coverage_by_sequence.txt}
    echo -e "$sample\t$cov" >> ~/Nd_SPANDx_all_seqs/Outputs/coverage_by_sample.dedup_bam.txt
done
```

calculate coverage based on DP after spandx filtering. Note that the R script points to the file paths
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/cov_vcfR.slurm
```

Depth results and plots stored locally at `~/repo/neonectria_genome_reseq_10072020/data/coverage_Nd`



### Next will need to perform coverage based filtering based on these results. 6 read minimum 79 read max (0.25 or 3x of mean 26.3 DP. Note that median DP is 16 )
### The VCF file `out.filtered.PASS.vcf` first needs to be modified to exclude commas within quotes (i.e., description fields) in the ##INFO rows as vcftools will not handle these. Write to `out.filtered.PASS.no_comma.vcf` the easiest/most acurate way seems to be to modify manually
```
cp out.filtered.PASS.vcf out.filtered.PASS.no_comma.vcf
vim out.filtered.PASS.no_comma.vcf
# remove commas by hand
```

```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vcftools_DP_filter_min6_max79.slurm

```
Count variants and NA vals
```
cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/
grep -v "##\|#" out.filtered.PASS.DP_filtered.recode.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.recode.vcf | wc -l
```
#### Note that `vcftools` recodes NA values "." as diploid NA "./."
1112323 variants # (first grep) x 19 samples
21134137 total variants
3891541 NA # second grep
3891541/21134137 = 18.4%

#### Missing data filter, maximum 25%. (i.e., no more than 25% of samples are NA at the site)
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_missing_dat_filter_0.25.slurm
```
The slurm script is giving an illegal instrution error if we end up on node 117 or 118.
Add `#SBATCH --exclude=node[117-118]` to fix illegal instruction
```
cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.lt25missing.vcf | wc -l
```
915968 remaining variants x 19 = 17403392
2074428 NA
2074428/17403392 = 11.9% NA

### Remove polyallelic sites
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vcftools_polyalleles_rm.slurm 
cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.biallele.recode.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.lt25missing.biallele.recode.vcf | wc -l
```
871375 remaining variants x 19 = 16556125
1976185 NA
1976185/16556125 = 11.9% NA

### Filtering sites based on minor allele *count* `--mac` of 3 (at least 3 samples) and then filtering samples based on missing data
Run in interactive session! First start a tmux session
```
module purge
tmux new -s vcf_filter
srun --pty bash -i
module purge
module load linuxbrew/colsa
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.recode.vcf --mac 2 --recode --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.recode.vcf --mac 3 --recode --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac3
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.recode.vcf | wc -l
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.biallele.mac3.recode.vcf | wc -l
```
Don't forget to `exit` interactive session!
228931 sites remaining at --mac 2
150533 sites remaining at --mac 3
using mac 2

### need to look at proportion of missing data per individual. Useful to plot. 
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/NA_from_vcfR.slurm 
```
Download plot
```
cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/
```

Download to Nd_SPANDx_all_seqs. 2 samples with >30% missing data

```
srun --pty bash -i
module purge
module load linuxbrew/colsa

vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.recode.vcf --missing-indv
module purge
awk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
```
2 samples flagged. 
```
module load linuxbrew/colsa
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.recode.vcf --remove lowDP.indv --recode --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf | wc -l

```
17 individuals and 228931 sites remaining = 3891827
309744 NA
309744/13357614 = 7.959% NA

#
x sites recognized by R::PopGenome
x sites recongized including poly allelic


### Perform LD filtering before population structure analyses
Using BCFtools
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_10KB.slurm
grep -v "^##" ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.vcf | wc -l
```
11346 SNPs remaining at 10Kb filter

## Format conversions for downstream processes

### ADMIXTURE
Convert VCF to PED format for input to ADMIXTURE or LEA. We use using plink 1.9 after trying several options (see old readme for more). 
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_BED.slurm
```
Note that for ADMIXTURE PED needs to be recoded using `--recode12` and not `--recode`
- ADMIXTURE also requires that the `.map` file containes only numerics in "chromosome" (contig) names

```
cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf
cp out.filtered.LD_filtered_0.5_10Kb.map out.filtered.LD_filtered_0.5_10Kb.map.original
sed 's/[a-z,_,A-Z,\.]*//g' out.filtered.LD_filtered_0.5_10Kb.map.original > out.filtered.LD_filtered_0.5_10Kb.map
```

### LFMM
REcoding entire SNP set with `--recode01 --missing-genotype 9` for input to R::LFMM (or R::LEA, but that package (older v. of LFMM) seems to be broken). Note, this must be filtered for biallelic SNPs only. The resulting PED will need the first 6 columns removed. We also remove every other line because plink is stupid and doubles every SNP assuming diploid data
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.full_dat_LFMM.slurm

cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf

#rm first 6
cut out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.ped -d " " -f 1-6 --complement | awk  '{for (i=1;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.lfmm
#write sample IDs
cut out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.ped -d " " -f 1 > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.sampleIDs

grep "##contig=<ID=" out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf > scaffold_lengths.txt

sed -e 's/##contig=<ID=//' -e 's/length=//' -e 's/>//' scaffold_lengths.txt > scaffold_lengths.csv
```
### Pcadapt
Also producing BED files of non-LD-filtered VCF for pcadapt
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_BED.full_dat.slurm
```

### STRUCTURE format
Download `out.filtered.LD_filtered_0.5_10Kb.fam` to make indfile for structure `make_structure_indfile.r`
Downloaded .fam, .lfmm file to `repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs/`

upload `Nd_SPANDx_all_seqs/ind_file.structure`

```
module purge
module load linuxbrew/colsa
cd ~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/
```
First run pgdspider with no -spid to generate template
```
pgdspider -inputfile out.filtered.LD_filtered_0.5_10Kb.vcf -inputformat VCF -outputfile Nd.out.filtered.LD_filtered_0.5_10Kb.structure -outputformat STRUCTURE 
```
Edits to template_VCF_STRUCTURE.spid
- VCF_PARSER_POP_FILE_QUESTION=ind_file.structure
- VCF_PARSER_PLOIDY_QUESTION=HAPLOID
- VCF_PARSER_POP_QUESTION=true
- VCF_PARSER_PL_QUESTION=false
- STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
- STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
- STRUCTURE_WRITER_FAST_FORMAT_QUESTION=false
```
pgdspider -inputfile out.filtered.LD_filtered_0.5_10Kb.vcf -inputformat VCF -outputfile out.filtered.LD_filtered_0.5_10Kb.structure -outputformat STRUCTURE -spid template_VCF_STRUCTURE.spid
awk '{print NF}' out.filtered.LD_filtered_0.5_10Kb.structure | sort -nu | tail -n 1
```

PGSpider is leaving 9925 SNPs (above awk -2). Set mainparams loci number to 9925,  LABEL to 1 and POPDATA to 1.

### LD filtered VCF, PED, and BED files are at
```
~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.ped
~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.bed
~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.vcf
```
### Full SNP set (i.e., quality filtered but not LD filtered) is at
```
~/Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf
```

For local pcadapt. sftp:
```
lcd repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs
get Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb*
get Nd_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode*
```

# left off here and moved to structure runs

## LD decay
Trying with plink 1.9 (.bim file) (https://www.biostars.org/p/300381/)
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_LD_decay.slurm
```
File is 5.8 G. Header `CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2 `. Checking to maker sure SNPs were not calculated across chromosomes
```
gunzip -c out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay.ld.gz | awk '$1 != $4' | wc -l
# 1 (the header)
```
Test data
```
gunzip -c out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay.ld.gz | head -n 5000 > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_decay.ld.test
```
Run R plotting script with avgs across different distances (1Kb, 2.5Kb, 5Kb, 10Kb)
```
cd neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/LD_decay_plots.R.slurm 
```

## LD blocks 
Gabriel et al. 2002 https://www.biostars.org/p/300381/
https://www.cog-genomics.org/plink/1.9/ld#blocks
Gabriel method is based partly on comparing to phenotypes, which we don't have of course.

```
cd neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_LD_block.slurm 
```
run report
```
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block.log.
Options in effect:
  --allow-extra-chr
  --bfile out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode
  --blocks no-pheno-req
  --blocks-max-kb 200
  --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block
  --threads 24

128714 MB RAM detected; reserving 64357 MB for main workspace.
130957 variants loaded from .bim file.
102 people (0 males, 0 females, 102 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block.nosex
.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 102 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.890081.
130957 variants and 102 people pass filters and QC.
Note: No phenotypes present.
--blocks: 99%^M--blocks: 5612 haploblocks written to
out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block.blocks
.
Extra block details written to
out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block.blocks.det
.
Longest span: 49.861kb.
```
The .det file contains block position and block length (Kb) and n snps per block. Count number of SNPs categorized
```
LD_block_SNP_count.r
``` 
19914 of 130957 variants in blocks. Mean block len = 650 bp, median = 153 bp, not much help for mappping SNPs to genes


## population structure analyses

### Running structure_threader
First generate the params files
```
cd 
mkdir Nd_SPANDx_all_seqs_structure_th
module purge
module load anaconda/colsa
export PATH=~/.local/bin:$PATH

cd Nd_SPANDx_all_seqs_structure_th
structure_threader params -o ~/Nd_SPANDx_all_seqs_structure_th/

```
Download the params file to edit and retain a local copy if desired. The default parameters generated by structure_threader in `extraparams` are satisfactory.
#### Edits to `mainparams`
- set PLOIDY to 1
- NUMINDS 17 #after filtering out seuqncing reps and high missing data
- NUMLOCI 9925 #SNPs remaining after pgspider
- MISSING -9
- LABEL 1
- POPDATA 1
- MAXPOPS 8 #although this is set by -K on the command line; use 8 bc there are only 17 samples, silly to have more than 2 samples per pop

Run structure_threader (started at 3/20/23 3:45p)
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/structure_threader.slurm
```
The  structure threader run started writing results at 11 days. After 14.75 days there are nine results written (K4 and K3). K8 started writing at 20d 14h (6/17 04:00), K1 has finished by 23d5h (06/19 20:00), K10 started writing at at 23d15h (06/20 0600) finished 24d14h (06/21 0400)
#### Note that the run may be killed early due to Premise downtime starting 06/26. If that is the case can restart any jobs that didn't finish and/or use structureharvester outside of the structure threader wrapper to calculate evanno method

#### Plot evanno
Run locally
```
Rscript plot_evanno.r
```

### Running admixture CV
```
cd 
mkdir Nd_SPANDx_all_seqs_admixture
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/admixture_CV.slurm
```
Ran ADMIXTURE with --haploid flag set 
```
grep "CV error" admixture.out
less ~/Nd_SPANDx_all_seqs_admixture/CV_by_K.text 
```
Look for a dip in CV results. There is only a steady increase


### Run PCADAPT locally `repo/pcadapt/pca_LD_filtered.r`
```
#file
out.filtered.LD_filtered_0.5_10Kb.bed
#script
pcadapt/pca_LD_filtered.r
```

### Run IBD locally (popgen/adegenet methods)
```
#file
out.filtered.LD_filtered_0.5_10Kb.bed
#script
R_scripts/IBD.adegenet_metrics.multiple_subsamples.r
partial_mantel_IBD_dur_inf.r
```

### Calculate phylogeny
Convert GVCF to SNP matrix
```
cd neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/gvcf2table.slurm 
```
There are very few SNPs with no NA values so the .clean files are useless. Convert NA i.e., `./.` to gaps for tree calculation
```
sed 's:\./\.:-:g' out.filtered.LD_filtered_0.5_10Kb.table > out.filtered.LD_filtered_0.5_10Kb.table.na2gap
sed 's:\./\.:-:g' out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.table > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.table.na2gap
```
Once NAs are converted convert to fasta multiple sequence alignment (run locally)
```
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/snp_table2fasta.pl out.filtered.LD_filtered_0.5_10Kb.table.na2gap 5 > out.filtered.LD_filtered_0.5_10Kb.fasta
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/snp_table2fasta.pl out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.table.na2gap 5 > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.fasta
```
Calculate ML tree in R
```
ml_tree.adegenet.r
ml_tree.adegenet.plot.r
```



## After pop structure analyses performing analyses of diversity (e.g. nucleotide diversity and Fst) between sites OR across dataset
### R package PopGenome 
### Goal is to perform whole genome and sliding window analyses

## Genome scans for diversity metrics
#### See [here](https://wurmlab.com/genomicscourse/2016-SIB/practicals/population_genetics/popgen) for analyses using popgenome with haploid data. We first start by splitting the data into scaffolds (or chromosomes), which can be done with bcftools or bash

Then IF doing locally need to activate conda env for bcftools -- OR just do all of this with bcftools on the server where it's installed. can use the first few lines to loop through scaffolds
```

#make dirs for processing VCF files to individual scaffolds

mkdir ~/repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs/scaffolds_split
grep -v "#" ~/repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf | cut -f 1 | uniq \
    > repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs/scaffolds.txt
    
cut -f 1 ~/repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.map | uniq

conda activate bcftools

#compress and index VCF file
cd ~/repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs

bgzip out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf
tabix -p vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf.gz

while read line
do(
    bcftools view out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf.gz $line > scaffolds_split/$line
)
done < scaffolds.txt

conda deactivate
```
Need to set up separate dirs for each scaffold
```
cd ~/repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs

#array of files (there are 20)
scf_files=(scaffolds_split/*)
#loop through by number and copy to dir
for i in {0..19}
do(
    dir_name=$((i+1))
    mkdir scf_$dir_name
    cp ${scf_files[$i]} scf_$dir_name
)
done 
```
Genome scans for Pi theta TajD
```
F_stats.PopGenome.read_ind_scf.scan_Pi_theta.no_pops.r
```
NOTE that four of the total 24 contigs 4 are left out of VCF because there were no SNPs called on these contigs. tig00007952_pilon tig00007948_pilon tig00007947_pilon tig00000405_pilon. These are all less than 100kb. Would be interesting to see what's on these contigs



## LFMM
#### Whole SNP set. DO NOT filter to sites with >= 4 samples per site for analyses to be performed in conjuntion with LFMM tests of SNP-climate correlation/GWAS, because the minimum sample number is not necessary when not performing at the population(site) level and we will have more data points
First running LFMM to identify significant SNPs. We will use the automatically calibrated p-values instead of manually adjusted according to Francois et al. 2016 "recalibrating is always subjective" (both auto and manually recalibrated are retained in the summary tables)
```
LFMM_analyses/initial_LFMM_full_data.r
```
We are not picking up many SNPs with freeze-thaw. It is possible the more important var is tmin. Running that as well for comparison
```
LFMM_analyses/LFMM_full_data.tmin.r
```
Also run growing season HDD4
```
LFMM_analyses/LFMM_full_data.growing_hdd4.r
```
#### Comparing SNP position to gene position (from gff) to identify nearest neighbors to significant SNPs of different variables
run locally
```
LFMM_analyses/nearest_neighbor_genes.hdd4.r
LFMM_analyses/nearest_neighbor_genes.freezeThaw.r
LFMM_analyses/nearest_neighbor_genes.ppt.r
LFMM_analyses/nearest_neighbor_genes.tmin.r
LFMM_analyses/nearest_neighbor_genes.growing_hdd4.r
```
Venn diagrams of genes and SNPs shared between different climate metrics
```
LFMM_analyses/gene_SNP_venn.R
```
PCA of climate vars
```
LFMM_analyses/climate_vars_pca.r
```

### The following routine is no longer used in preference for the nearest neighbor method above
### Identify which genes SNPs occur on
get gene IDS to pull from gff
```
cd neonectria_genome_reseq_10072020/maker2_run
grep ">" makerFINAL.all.maker.transcripts.fasta > transcipt_IDs.txt
```
Download GFF and IDs
```
cd repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs/maker2_ann
sftp
get neonectria_genome_reseq_10072020/maker2_run/makerFINAL.all.gff
get neonectria_genome_reseq_10072020/maker2_run/transcipt_IDs.txt
exit
```
pull relevant lines from GFF
```
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_mRNA_IDs_from_GFF.pl makerFINAL.all.gff > makerFINAL.all.mRNA_ONLY.gff
```
count the number of SNPs on genes. files written to `data/Nd_LFMM_tables` dir
```
sig_SNPs_gene_count.hdd4.r
sig_SNPs_gene_count.freezeThaw.r
sig_SNPs_gene_count.ppt.r
```

## GO mapping and enrichment tests for LFMM partitions
### GO mapping and GO slim mapping performed by Tuan D. using blast2GO. 
Some initial exploration of the mapping file
```
grep "NA---" data/blast2GO/makerFINAL.all.maker.proteins-1.txt | wc -l
#1025
grep "NO-BLAST" data/blast2GO/makerFINAL.all.maker.proteins-1.txt | wc -l
#1025
grep -v "MAPPED" data/blast2GO/makerFINAL.all.maker.proteins-1.txt | wc -l
#3652
grep "MAPPED" data/blast2GO/makerFINAL.all.maker.proteins-1.txt | wc -l
#10413
grep "ANNOTATED" data/blast2GO/makerFINAL.all.maker.proteins-1.txt | wc -l
#9631
```
As indicated in Tuan's email, the number of genes with high-quality annotations is 9631 which corresponds to the 'ANNOTATED' flag in the file
```
grep "NO-BLAST" data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.txt | wc -l
#1025
grep "MAPPED" data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.txt | wc -l
#10413
grep "ANNOTATED" data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.txt | wc -l
#6
grep "GO-SLIM" data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.txt | wc -l
#9625
```
In the slim mapping the flag for high-quality mappings is 'GO-SLIM'. The 'ANNOTATED' flag indicates six mappings that are apparently high-quality, but map to `C: obsolete intrinsic/anchored compnent of membrane.
```
grep "ANNOTATED" data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.txt | wc -l
#1469
grep "GO-SLIM" data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.txt | wc -l
#8162
grep "ANNOTATED" data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.txt | less
```
There are many "hypothetical" and "putative" proteins/functions in the ANNOTATED set for the generic slim, but also many with apparently poor slim mappings. Can presumably continue with just the yeast slim, but will still plot up functions

### Convert the B2GO slim mapping to a gene association style format, with one row per gene-GO association instead of one row per gene with multiple GO terms listed comma-separated in a single column. Only high-quality annotations are considered
```
cd repo/neonectria_genome_reseq_10072020/
perl GO_slim_and_enrichment/B2GOslim2long.pl data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.txt GO-SLIM | less
#first run the granular terms for comparison
perl GO_slim_and_enrichment/B2GOslim2long.pl data/blast2GO/makerFINAL.all.maker.proteins-1.txt ANNOTATED > data/blast2GO/makerFINAL.all.maker.proteins-1.long.txt
#slim sets
perl GO_slim_and_enrichment/B2GOslim2long.pl data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.txt GO-SLIM> data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.long.txt
perl GO_slim_and_enrichment/B2GOslim2long.pl data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.txt GO-SLIM > data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.long.txt
```
Note there are a few genes that hit the same slim term multiple times. Presumably these are hits from different child terms (i.e., the gene mapped to multiple children that then mapped to the same slim term), but we should make sure to only count these once
```
sort data/blast2GO/makerFINAL.all.maker.proteins-1.long.txt | uniq > data/blast2GO/makerFINAL.all.maker.proteins-1.long.uniq.txt
sort data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.long.txt | uniq > data/blast2GO/makerFINAL.all.maker.proteins-1_GOSlim_Yeast.long.uniq.txt
sort data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.long.txt | uniq > data/blast2GO/makerFINAL.all.maker.proteins-1_GOslim_generic.long.uniq.txt
```
basic stats and plots
```
GO_slim_and_enrichment/slim_stats_and_plots.r
```
In general F aspect has more genes mapped (about 2k) but lower number of categories (like half) than P for both yeast and generic slims


### This routine no longer needed because switched to GO slim mapping of whole set and then performing enrichment tests using R instead of blast2GO
Extract sequences of genes associated with different SNP sets for GO enrichment test sets (run locally)
```
cd repo/neonectria_genome_reseq_10072020/
perl perl_scripts/get_seqs_by_list_from_fasta.pl data/Nd_SPANDx_all_seqs/maker2_ann/makerFINAL.all.maker.proteins.faa data/Nd_LFMM_tables/hdd4.geneIDs.nearest_neighbors.txt > data/Nd_LFMM_tables/hdd4.nearest_neighbors.faa

perl perl_scripts/get_seqs_by_list_from_fasta.pl data/Nd_SPANDx_all_seqs/maker2_ann/makerFINAL.all.maker.proteins.faa data/Nd_LFMM_tables/freezeThaw.geneIDs.nearest_neighbors.txt > data/Nd_LFMM_tables/freezeThaw.nearest_neighbors.faa

perl perl_scripts/get_seqs_by_list_from_fasta.pl data/Nd_SPANDx_all_seqs/maker2_ann/makerFINAL.all.maker.proteins.faa data/Nd_LFMM_tables/ppt.geneIDs.nearest_neighbors.txt > data/Nd_LFMM_tables/ppt.nearest_neighbors.faa

grep ">" data/Nd_LFMM_tables/hdd4.nearest_neighbors.faa | wc -l
#49
grep ">" data/Nd_LFMM_tables/freezeThaw.nearest_neighbors.faa | wc -l
#5
grep ">" data/Nd_LFMM_tables/ppt.nearest_neighbors.faa | wc -l
#129

wc -l data/Nd_LFMM_tables/hdd4.geneIDs.nearest_neighbors.txt
#49
wc -l data/Nd_LFMM_tables/freezeThaw.geneIDs.nearest_neighbors.txt
#5
wc -l data/Nd_LFMM_tables/ppt.geneIDs.nearest_neighbors.txt
#129
```







#### Have not run pairwise site comparisons for now. May be interesting, but also hard to compare to significant SNP correlations because pariwise stats are run on a subset of the data
## pop level diversity
### Using R::PopGenome for diversity stats
### For population level stats (e.g., Fst comps between sites) will need to filter low n sites

Issues with calculating population level diversity stats in sites with low sample number inclduing ME.N, MI, and NJ, etc. Have output of list of sample IDs associated with these sites and will filter out the samples.
```
list_low_samples_to_retain_n_ge4.r
```
Then rerun scaffold split
```
conda activate bcftools
cd ~/repo/neonectria_genome_reseq_10072020/data/Nd_SPANDx_all_seqs

bcftools view --samples-file retain_samples.txt out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf.gz > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.rm_low_n_sites.vcf
bgzip out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.rm_low_n_sites.vcf
tabix -p vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.rm_low_n_sites.vcf.gz

mkdir scaffolds_split_rm_low_n

while read line
do(
bcftools view out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.rm_low_n_sites.vcf.gz $line > scaffolds_split_rm_low_n/$line
)
done < scaffolds.txt
```
### Now running using popgenome with new script
```
PopGenome/F_stats.PopGenome.rm_low_n_sites.r
```

