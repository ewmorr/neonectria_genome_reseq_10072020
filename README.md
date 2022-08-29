# Repo for analysis of genome resequencing data of Neonectria faginata and Neonectria ditissima
#### Workflows performed on UNH Premise or locally using R as appropriate

## set up sequence files for SNP calling

### Set up dirs for performing SPANDx SNP calling (premise)
```
mkdir Nf_SPANDx_all_seqs
mkdir Nd_SPANDx_all_seqs
```

local. Find matching sequence files by sample to combine
```
Rscript sort_seq_files_for_SPANDx.r
```

### Cp seqs to new dirs using the lists of IDs created above. Some samples have two sets of sequences and those are concatenated
```
#Nf first set
while IFS= read -r line 
do(
    cp neonectria_genome_reseq_10072020/reads/${line}_*R1*.fastq.gz Nf_SPANDx_all_seqs/${line}_1.fastq.gz
    cp neonectria_genome_reseq_10072020/reads/${line}_*R2*.fastq.gz Nf_SPANDx_all_seqs/${line}_2.fastq.gz
)
done < sample_IDs.Nf.10072020.txt

#Nf second set
while IFS= read -r line 
do(
    cp neonectria_genome_reseq_03312022/reads/${line}_*R1*.fastq.gz Nf_SPANDx_all_seqs/${line}_1.fastq.gz
    cp neonectria_genome_reseq_03312022/reads/${line}_*R2*.fastq.gz Nf_SPANDx_all_seqs/${line}_2.fastq.gz
)
done < sample_IDs.Nf.03312022.txt

#Nf second set with additional reads
while IFS= read -r line 
do(
    cat neonectria_genome_reseq_03312022/reads/${line}_*R1*.fastq.gz neonectria_genome_reseq_03312022_adtl_reads/reads/${line}_*R1*.fastq.gz > Nf_SPANDx_all_seqs/${line}_1.fastq.gz
    cat neonectria_genome_reseq_03312022/reads/${line}_*R2*.fastq.gz neonectria_genome_reseq_03312022_adtl_reads/reads/${line}_*R2*.fastq.gz > Nf_SPANDx_all_seqs/${line}_2.fastq.gz
)
done < sample_IDs.Nf.03312022_adtl_reads.txt

#Nd first set
while IFS= read -r line 
do(
    cp neonectria_genome_reseq_10072020/reads/${line}_*R1*.fastq.gz Nd_SPANDx_all_seqs/${line}_1.fastq.gz
    cp neonectria_genome_reseq_10072020/reads/${line}_*R2*.fastq.gz Nd_SPANDx_all_seqs/${line}_2.fastq.gz
)
done < sample_IDs.Nd.10072020.txt

#Nd second set
while IFS= read -r line 
do(
    cp neonectria_genome_reseq_03312022/reads/${line}_*R1*.fastq.gz Nd_SPANDx_all_seqs/${line}_1.fastq.gz
    cp neonectria_genome_reseq_03312022/reads/${line}_*R2*.fastq.gz Nd_SPANDx_all_seqs/${line}_2.fastq.gz
)
done < sample_IDs.Nd.03312022.txt

#Nd second set with additional reads
while IFS= read -r line 
do(
    cat neonectria_genome_reseq_03312022/reads/${line}_*R1*.fastq.gz neonectria_genome_reseq_03312022_adtl_reads/reads/${line}_*R1*.fastq.gz > Nd_SPANDx_all_seqs/${line}_1.fastq.gz
    cat neonectria_genome_reseq_03312022/reads/${line}_*R2*.fastq.gz neonectria_genome_reseq_03312022_adtl_reads/reads/${line}_*R2*.fastq.gz > Nd_SPANDx_all_seqs/${line}_2.fastq.gz
)
done < sample_IDs.Nd.03312022_adtl_reads.txt

```
Also cat the duplicates
```
cat neonectria_genome_reseq_10072020/reads/NG10_*R1*.fastq.gz neonectria_genome_reseq_10072020/reads/NG76_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG10_1.fastq.gz
cat neonectria_genome_reseq_10072020/reads/NG10_*R2*.fastq.gz neonectria_genome_reseq_10072020/reads/NG76_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG10_2.fastq.gz

cat neonectria_genome_reseq_03312022/reads/NG108_*R1*.fastq.gz neonectria_genome_reseq_03312022/reads/NG149_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG108_1.fastq.gz
cat neonectria_genome_reseq_03312022/reads/NG108_*R2*.fastq.gz neonectria_genome_reseq_03312022/reads/NG149_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG108_2.fastq.gz

cat neonectria_genome_reseq_10072020/reads/NG48_*R1*.fastq.gz neonectria_genome_reseq_10072020/reads/NG78_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG48_1.fastq.gz
cat neonectria_genome_reseq_10072020/reads/NG48_*R2*.fastq.gz neonectria_genome_reseq_10072020/reads/NG78_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG48_2.fastq.gz

cat neonectria_genome_reseq_10072020/reads/NG28_*R1*.fastq.gz neonectria_genome_reseq_10072020/reads/NG77_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG28_1.fastq.gz
cat neonectria_genome_reseq_10072020/reads/NG28_*R2*.fastq.gz neonectria_genome_reseq_10072020/reads/NG77_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG28_2.fastq.gz

cat neonectria_genome_reseq_10072020/reads/NG62_*R1*.fastq.gz neonectria_genome_reseq_10072020/reads/NG79_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG62_1.fastq.gz
cat neonectria_genome_reseq_10072020/reads/NG62_*R2*.fastq.gz neonectria_genome_reseq_10072020/reads/NG79_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG62_2.fastq.gz

cat neonectria_genome_reseq_03312022/reads/NG106_*R1*.fastq.gz neonectria_genome_reseq_03312022/reads/NG161_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG106_1.fastq.gz
cat neonectria_genome_reseq_03312022/reads/NG106*R2*.fastq.gz neonectria_genome_reseq_03312022/reads/NG161_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG106_2.fastq.gz

```

## SPANDx SNP calling
Copy the referece genomes into the SPANDx working dirs
```
cp SPANDx_Nf/ref.fasta Nf_SPANDx_all_seqs/
cp N_ditissima_ref_genome/LDPL01.1.fsa_nt.fasta Nd_SPANDx_all_seqs/ref.fasta


#### Make sure that nextflow.config is updated if necessary (https://github.com/dsarov/SPANDx#usage)
The config file is where CPUs etc are denoted as well as the resource manager (e.g., SLURM). Newer versions of the package also have `notrim` set to `true`. This should be `false` in that case. Also, note that we have cloned the git repo after installing via conda, and made some modificaations to the `main.nf` script (i.e., changing the `gatk HaplotypeCaller` comand at line 865 to include `--ploidy 1` flag) and to the `./bin/Master_vcf.sh` script (i.e., removing `-ploidy 1` from the `gatk GenotypeGVCFs` command). Also, note that `.bashrc` may need to be updated as described [here](./SPANDx_conda_install.sh). Finally, the reference genome assembly must be loacted in the SPANDx working directory (along with the reads), and can be indicated by path in the config file. Then, to run SPANDx (note that nextflow is pointed to the cloned git repo)

```

Running full sample set of Nf SPANDx
```
module purge
module load anaconda/colsa

cd ~/Nf_SPANDx_all_seqs
conda activate spandx
source ~/.bashrc

nextflow config ~/SPANDx_git_clone/

screen

nextflow run ~/SPANDx_git_clone/

Ctrl-a Ctrl-d
screen -list
```
The first run on 05132022 is `30523.pts-47.login01`
```
screen -r 30523.pts-47.login01
```
Running Nd
```
cd ~/Nd_SPANDx_all_seqs/

screen

nextflow run ~/SPANDx_git_clone/

Ctrl-a Ctrl-d
screen -list
```
25526.pts-0.login01

## SNP and sample filtering

### Filter out SNPs that SPANDx FAIL filtering
```
cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf
grep "##\|#\|PASS" out.filtered.vcf > out.filtered.PASS.vcf
grep -v "##\|#" out.filtered.PASS.vcf | wc -l
```
516,326 SNPs remaining x 117 samples = 60410142
```
grep -o "\s\.:" out.filtered.PASS.vcf | wc -l
```
2719040 NA sites
2719040/60410142 = 4.50%


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
cd ~/Nf_SPANDx_all_seqs/Outputs/bams
for i in *.coverage_by_sequence.txt
do
    cov="$(grep "genome" $i| cut -f 3)"
    sample=${i%.coverage_by_sequence.txt}
    echo -e "$sample\t$cov" >> ~/Nf_SPANDx_all_seqs/Outputs/coverage_by_sample.dedup_bam.txt
done
```

calculate coverage based on DP after spandx filtering. Note that the R script points to the file paths
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/cov_vcfR.slurm
```

Depth results and plots stored locally at the home dir `coverage/Nf_all_seqs/`



### Next will need to perform coverage based filtering based on these results. 4 read minimum 48 read max (0.25 or 3x of mean 15.8 DP. Note that median DP is 6 and the mean of the first set of samples was 7.5 so setting DP min = 4 allows to capture more valid SNPs across entire sample set)
### The VCF file `out.filtered.PASS.vcf` first needs to be modified to exclude commas within quotes (i.e., description fields) in the ##INFO rows as vcftools will not handle these. Write to `out.filtered.PASS.no_comma.vcf` the easiest/most acurate way seems to be to modify manually
```
cp out.filtered.PASS.vcf out.filtered.PASS.no_comma.vcf
vim out.filtered.PASS.no_comma.vcf
# remove commas by hand
```

```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vcftools_DP_filter_min4_max48.slurm

```
Count variants and NA vals
```
cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/
grep -v "##\|#" out.filtered.PASS.DP_filtered.recode.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.recode.vcf | wc -l
```
#### Note that `vcftools` recodes NA values "." as diploid NA "./."
516326 variants # (first grep) x 117 samples
60410142 total variants
13,130,656 NA # second grep
13130656/60410142 = 21.7%

#### Missing data filter, maximum 25%. (i.e., no more than 25% of samples are NA at the site)
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_missing_dat_filter_0.25.slurm
```
The slurm script is giving an illegal instrution error (WHY?) This is a small enough job to run on the head node but need to figure this out
```
module purge
module load linuxbrew/colsa
cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/
bcftools view -i 'F_MISSING<0.25' out.filtered.PASS.DP_filtered.recode.vcf -Ov -o out.filtered.PASS.DP_filtered.lt25missing.vcf
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.lt25missing.vcf | wc -l
```
433071 remaining variants x 117 = 50669307
8921908 NA
8921908/50669307 = 17.6% NA

### Remove polyallelic sites
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vcftools_polyalleles_rm.slurm 
cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.biallele.recode.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.lt25missing.biallele.recode.vcf | wc -l
```
416688 remaining variants x 117 = 48752496
8583122 NA
8583122/48752496 = 17.6% NA

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
130957 sites remaining at --mac 2
94489 sites remaining at --mac 3
using mac 2

### need to look at proportion of missing data per individual. Useful to plot. This could be run on the server using vcfR conda env and `NA_from_VCF.r` but the server is currently not running conda due to hardware issues. Trying to run locally and have not set up a slurm script
Download to Nf_SPANDx_all_seqs. 15 samples with >30% missing data

```
srun --pty bash -i
module purge
module load linuxbrew/colsa

vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.recode.vcf --missing-indv
module purge
awk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
```
15 samples flagged. 
```
module load linuxbrew/colsa
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.recode.vcf --remove lowDP.indv --recode --out out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf | wc -l

```
102 individuals and 130957 sites remaining = 13357614
1460689 NA
1460689/13357614 = 10.9% NA

#
83542 sites recognized by R::PopGenome
120803 sites recongized including poly allelic


### Perform LD filtering before population structure analyses
Using BCFtools
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_10KB.slurm
grep -v "^##" ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.vcf | wc -l
```
Tried filter at 50Kb and 10Kb orignially, 50Kb seems excessive esp. for genome size. 45019 SNPs remaining at 10Kb filter

## Format conversions for downstream processes

### ADMIXTURE
Convert VCF to PED format for input to ADMIXTURE or LEA. We use using plink 1.9 after trying several options (see old readme for more). 
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_BED.slurm
```
This appears to work.... Note that for ADMIXTURE PED needs to be recoded using `--recode12` and not `--recode`
- ADMIXTURE also requires that the `.map` file containes only numerics in "chromosome" (contig) names


```
cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf
cp out.filtered.LD_filtered_0.5_10Kb.map out.filtered.LD_filtered_0.5_10Kb.map.original
sed 's/[a-z,_]*//g' out.filtered.LD_filtered_0.5_10Kb.map.original > out.filtered.LD_filtered_0.5_10Kb.map
```

### LFMM
REcoding entire SNP set with `--recode01 --missing-genotype 9` for input to R::LFMM (or R::LEA, but that package (older v. of LFMM) seems to be broken). Note, this must be filtered for biallelic SNPs only. The resulting PED will need the first 6 columns removed. We also remove every other line because plink is stupid and doubles every SNP assuming diploid data
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.full_dat_LFMM.slurm

cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf

cut out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.ped -d " " -f 1-6 --complement | awk  '{for (i=1;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' > out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.lfmm

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
Downloaded .fam, .lfmm file to `repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs/`

upload `Nf_SPANDx_all_seqs/ind_file.structure`

```
module purge
module load linuxbrew/colsa
cd ~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/
```
First run pgdspider with no -spid to generate template
```
pgdspider -inputfile out.filtered.LD_filtered_0.5_10Kb.vcf -inputformat VCF -outputfile Nf.out.filtered.LD_filtered_0.5_10Kb.structure -outputformat STRUCTURE 
```
Edits
- VCF_PARSER_PLOIDY_QUESTION=HAPLOID
- VCF_PARSER_POP_QUESTION=true
- VCF_PARSER_POP_FILE_QUESTION=ind_file.structure
- VCF_PARSER_PL_QUESTION=false
- STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
- STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
- STRUCTURE_WRITER_FAST_FORMAT_QUESTION=false
```
pgdspider -inputfile out.filtered.LD_filtered_0.5_10Kb.vcf -inputformat VCF -outputfile out.filtered.LD_filtered_0.5_10Kb.structure -outputformat STRUCTURE -spid template_VCF_STRUCTURE.spid
awk '{print NF}' out.filtered.LD_filtered_0.5_10Kb.structure | sort -nu | tail -n 1
```

PGSpider is leaving 40857 SNPs (above awk -2)... will move forward with structure run but will need to investigate. Set mainparams loci number to 40857,  LABEL to 1 and POPDATA to 1.

### LD filtered VCF, PED, and BED files are at
```
~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.ped
~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.bed
~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb.vcf
```
### Full SNP set (i.e., quality filtered but not LD filtered) is at
```
~/Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf
```

For local pcadapt. sftp:
```
lcd repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs
get Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.LD_filtered_0.5_10Kb*
lcd ..
get Nf_SPANDx_all_seqs/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode*
```

## LD decay
Trying with plink 1.9 (.bim file)
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/plink1.9_LD_deacy.slurm
```


## population structure analyses

### Running structure_threader
First generate the params files
```
cd 
mkdir Nf_SPANDx_all_seqs_structure_th
module purge
module load anaconda/colsa
export PATH=~/.local/bin:$PATH

cd Nf_SPANDx_all_seqs_structure_th
structure_threader params -o ~/Nf_SPANDx_all_seqs_structure_th/

```
Download the params file to edit and retain a local copy if desired. The default parameters generated by structure_threader in `extraparams` are satisfactory.
#### Edits to `mainparams`
- set PLOIDY to 1
- NUMINDS 102 #after filtering out seuqncing reps and high missing data
- NUMLOCI 40857 #SNPs remaining after pgspider
- MISSING -9
- LABEL 1
- POPDATA 1
- MAXPOPS 15 #although this is set by -K on the command line

Run structure_threader (started at 5/27 2p)
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
mkdir Nf_SPANDx_all_seqs_admixture
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/admixture_CV.slurm
```
Ran ADMIXTURE with --haploid flag set 
```
grep "CV error" admixture.out
less ~/Nf_SPANDx_all_seqs_admixture/CV_by_K.text 
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



## After pop structure analyses performing analyses of pairwise diversity (e.g. nucleotide diversity and Fst) between sites OR across dataset
### R package PopGenome 
### Goal is to perform whole genome and sliding window analyses


## Sliding window analysis and LFMM
#### Whole SNP set. DO NOT filter to sites with >= 4 samples per site for analyses to be performed in conjuntion with LFMM tests of SNP-climate correlation/GWAS, because the minimum sample number is not necessary when not performing at the population(site) level and we will have more data points
First running LFMM to identify significant SNPs
```
LFMM_analyses/initial_LFMM_full_data.r
```

#### See [here](https://wurmlab.com/genomicscourse/2016-SIB/practicals/population_genetics/popgen) for analyses using popgenome with haploid data. We first start by splitting the data into scaffolds (or chromosomes), which can be done with bcftools or bash

Then IF doing locally need to activate conda env for bcftools -- OR just do all of this with bcftools on the server where it's installed. can use the first few lines to loop through scaffolds
```

#make dirs for processing VCF files to individual scaffolds

mkdir ~/repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs/scaffolds_split
grep -v "#" ~/repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode.vcf | cut -f 1 | uniq \
    > repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs/scaffolds.txt
    
cut -f 1 ~/repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.recode01missing9.map | uniq

conda activate bcftools

#compress and index VCF file
cd ~/repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs

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
cd ~/repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs

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
NOTE that four of the total 24 contigs are left out of VCF because there were no SNPs called on these contigs. tig00007952_pilon tig00007948_pilon tig00007947_pilon tig00000405_pilon. These are all less than 100kb. Would be interesting to see what's on these contigs






### Identify which genes SNPs occur on
get gene IDS to pull from gff
```
cd neonectria_genome_reseq_10072020/maker2_run
grep ">" makerFINAL.all.maker.transcripts.fasta > transcipt_IDs.txt
```
Download GFF and IDs
```
cd repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs/maker2_ann
sftp
get neonectria_genome_reseq_10072020/maker2_run/makerFINAL.all.gff
get neonectria_genome_reseq_10072020/maker2_run/transcipt_IDs.txt
exit
```
pull relevant lines from GFF
```
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_mRNA_IDs_from_GFF.pl makerFINAL.all.gff > makerFINAL.all.mRNA_ONLY.gff
```
count the number of SNPs on genes. files written to `data/Nf_LFMM_tables` dir
```
sig_SNPs_gene_count.hdd4.r
sig_SNPs_gene_count.freezeThaw.r
sig_SNPs_gene_count.ppt.r
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
cd ~/repo/neonectria_genome_reseq_10072020/data/Nf_SPANDx_all_seqs

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


########
########
########
########
########
Also running new gene mark annotations for LFMM

```
cd neonectria_genome_reseq_10072020/
sbatch ~/repo/ONS_Nf/genemark.pilon_polished.slurm

```
Get a couple of regions of high SNP correlation
```
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_segments_by_pos_from_fasta.pl ~/neonectria_minion/MAT2_polish/pilon_.fasta tig00000025_pilon 2524273 2544922 > neonectria_genome_reseq_10072020/high_LFMM_snp_tig00000025_pilon_2524273-2544922.fasta

perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_segments_by_pos_from_fasta.pl ~/neonectria_minion/MAT2_polish/pilon_.fasta tig00000025_pilon 3070669 3078240 > neonectria_genome_reseq_10072020/high_LFMM_snp_tig00000025_pilon_3070669-3078240.fasta

cat high_LFMM_snp_tig00000025_pilon_2524273-2544922.fasta high_LFMM_snp_tig00000025_pilon_3070669-3078240.fasta > high_LFMM_SNP.fasta

sbatch ~/repo/neonectria_genome_reseq_10072020/premise/blastx_swissprot_high_LFMM.slurm
```
Also blast maker proteins for comparison
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/blastp_swissprot_maker1.slurm #This got overwritten
```
Also blast genemark proteins for comparison
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/blastp_swissprot_genemark_LFMM.slurm
```
Top hits in order of 4085_g, 3399_g, 4236_g, 4236_g (note the last two don't have great hits based on e-value and this also blasts to a transposon region)
```
grep GIS2_YEAST ~/blast_dbs/swiss-prot_fungi_reviewed.fasta
#Zinc finger protein GIS2; GO regulation of translation
grep RHO3_YEAST ~/blast_dbs/swiss-prot_fungi_reviewed.fasta
grep RHO3_ASHGO ~/blast_dbs/swiss-prot_fungi_reviewed.fasta
#GTP-binding protein RHO3
grep YEG5_SCHPO ~/blast_dbs/swiss-prot_fungi_reviewed.fasta
#Ankyrin and IPT/TIG repeat-containing protein C26H5.05; GO regulation of transcription
grep TIP20_YEAST ~/blast_dbs/swiss-prot_fungi_reviewed.fasta
#Protein transport protein TIP20
```


### Extract PRISM data for new sites
On Premise
```
module purge
module load anaconda/colsa
conda activate PRISM
screen
srun Rscript ~/repo/neonectria_genome_reseq_10072020/R_scripts/PRISM_dailies_ppt_download.r
srun Rscript ~/repo/neonectria_genome_reseq_10072020/R_scripts/PRISM_dailies_download.r
```
Once the dailies are downloaded
```
sbatch repo/neonectria_genome_reseq_10072020/premise/PRISM_dalies_site_extract.slurm
```
Also process 30 yr normals locally using `sites_climate_dat.r` and process the created .rds object from `PRISM_dalies_site_extract.slurm` locally
```
#first download the .rds objects. Requires: sites_daily_tmin_tmax_ppt_20072018.df.rds
cd repo/neonectria_genome_reseq_10072020/
Rscript PRISM_calcs/tmin_tmax_dailys_GDD_calcs.r
Rscript PRISM_calcs/sites_climate_dat.r
```
