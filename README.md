### Repo for analysis of genome resequencing data of Neonectria faginata and Neonectria ditissima
#### Workflows performed on UNH Premise

```
mkdir neonectria_genome_reseq_10072020
cd ~/neonectria_genome_reseq_10072020
mkdir reads
mv cobb.sr.unh.edu/managed/201002_SN7001360_0501_BHCKWTBCX3_16MerGARNAS_Neonectria/reads/*fastq.gz ./reads

mkdir neonectria_genome_reseq_03312022
cd neonectria_genome_reseq_03312022
mkdir reads
mv cobb.sr.unh.edu/managed/220328_A01346_0053_AHMYCGDRXY_16MerGARNAS_Neonectria/reads/*fastq.gz ./reads
```

Sample IDs
```
cd neonectria_genome_reseq_10072020

for i in reads/*R1*.fastq.gz
do

    FILE=${i##*/}
    ID=${FILE%_S*}
    echo $ID >> sample_IDs.txt
done

cd neonectria_genome_reseq_03312022

for i in reads/*R1*.fastq.gz
do

    FILE=${i##*/}
    ID=${FILE%_S*}
    echo $ID >> sample_IDs.txt
done

```

BBDUK adapter trimming and quality trimming
```
cd neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bbduk_trim.slurm

#trim new treads
cd neonectria_genome_reseq_03312022/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bbduk_trim_033122_reads.slurm

```
Trimmomatic adapter trimming and quality trimming
```
sbatch ~/repo/neonectria_genome_reseq_10072020/trimmomatic_trim.slurm
```

```
cd ~/neonectria_genome_reseq_10072020/
grep "Result" bbduk.out > bbduk.results
cut -f 2 bbduk.results | cut -f 1 -d ' ' > bbduk_filtered_read_count.txt
grep "Input Read Pairs" trimmomatic.out > trimmomatic.results

cd ~/neonectria_genome_reseq_03312022/
grep "Result" bbduk.out > bbduk_03312022.results
cut -f 2 bbduk_03312022.results | cut -f 1 -d ' ' > bbduk_filtered_read_count_03312022.txt
```
Either filtering dropped similar numbers of read pairs -- less than a tenth of a percent in most cases. trimmomatic indicates 15-30 percent read-through (fwd only surviving)
```

```
Merge bbtrimmed reads and derep (`vsearch`) for search and then map to Nf/Nd ITS seqs for species ID
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vsearch_Nf_Nd_ITS.slurm
#Also run 03312022 samples
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vsearch_Nf_Nd_ITS_03312022.slurm
```
Eight of the genomes are Nd based on this analysis
#### ALSO comparing to MAT1 and MAT2 seqs
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vsearch_Nf_Nd_MAT_otu_tab.slurm
```

#### Test run of SPANDx on six genomes of Nf
The reference genome for Nf is based on CANU assembly of MinION reads (MAT1 isolate) with pilon polishing with Illumina reads (either MAT1 or MAT2 isolate). These are located at `~/neonectria_minion/MAT1_polish/pilon_.fasta` or at `~/neonectria_minion/MAT2_polish/pilon_.fasta`. We will use the MAT1 assembly despite slightly worse BUSCO completeness (see BUSCO and quast stats in these dirs) because this is assembled from a single isolate.

Make test run dir and copy test files with rename of sequence files to follow the format expected by SPANDx
```
mkdir SPANDx_test_run

for i in NG1 NG2 NG4 NG5 NG9 NG10
do(
    cp neonectria_genome_reseq_10072020/reads/${i}_*R1* ~/SPANDx_test_run/${i}_1.fastq.gz
    cp neonectria_genome_reseq_10072020/reads/${i}_*R2* ~/SPANDx_test_run/${i}_2.fastq.gz
)
done
```

#### Make sure that nextflow.config is updated if necessary (https://github.com/dsarov/SPANDx#usage)
The config file is where CPUs etc are denoted as well as the resource manager (e.g., SLURM). Newer versions of the package also have `notrim` set to `true`. THis should be `false` in that case. Also, note that we have cloned the git repo after installing via conda, and made some modificaations to the `main.nf` script (i.e., changing the `gatk HaplotypeCaller` comand at line 865 to include `--ploidy 1` flag) and to the `./bin/Master_vcf.sh` script (i.e., removing `-ploidy 1` from the `gatk GenotypeGVCFs` command). Also, note that `.bashrc` may need to be updated as described [here](./SPANDx_conda_install.sh). Finally, the reference genome assembly must be loacted in the SPANDx working directory (along with the reads), and can be indicated by path in the config file. Then, to run SPANDx (note that nextflow is pointed to the cloned git repo)

```
cd SPANDx_test_run
conda activate spandx
source ~/.bashrc

nextflow config ~/SPANDx_git_clone/

screen

nextflow run ~/SPANDx_git_clone/

```
Test run completed with no errors after `.bashrc` modifications

#### Running full sample set of Nf
First setting up sampes in dirs with names formatted
```
mkdir SPANDx_Nf

for i in neonectria_genome_reseq_10072020/reads/*R1*
do(

    DIR=${i%/*}
    FILE=${i##*/}
    BEFFILE=${FILE%R1*}
    AFTFILE=${FILE##*R1}
    R2FILE=${BEFFILE}R2${AFTFILE}
    SAMPNAM=${FILE%_S*}
    echo $SAMPNAM

    cp $i ./SPANDx_Nf/${SAMPNAM}_1.fastq.gz
    cp $DIR/$R2FILE ./SPANDx_Nf/${SAMPNAM}_2.fastq.gz

)
done

mkdir SPANDx_Nd

for i in NG3_ NG6_ NG18_ NG20_ NG27_ NG42_ NG45_ NG69_
do(
    mv SPANDx_Nf/$i* ./SPANDx_Nd/
)
done
```
Running SPANDx on Nf
```
module purge
module load anaconda/colsa

cd SPANDx_Nf
conda activate spandx
source ~/.bashrc

nextflow config ~/SPANDx_git_clone/

screen

nextflow run ~/SPANDx_git_clone/
```
### For Nd reference, download the isolate RS324p genome from Deng et al. 2015. The Gomez-Cortecero genome is referenced at MycoCosm, but the Deng genome is better contiguity (by a lot). The reference is at `N_ditissima_ref_genome/LDPL01.1.fsa_nt.fasta`

Running SPANDx on Nd
```
module purge
module load anaconda/colsa

cd SPANDx_Nd
conda activate spandx
source ~/.bashrc

nextflow config ~/SPANDx_git_clone/

screen

nextflow run ~/SPANDx_git_clone/
```
NG6 seems to have hung at trimmomatic step (after > 24 hours not finished). Killed and try restarting with `resume`
```
screen
nextflow run ~/SPANDx_git_clone/ --resume
```

### Generate BUSCO gene models using `busco --long`
Using the reference sequence at `SPANDx_Nf/ref.fasta`
```
cd
cp neonectria_minion/Nf_canu_run0/config.ini SPANDx_Nf/
cd neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/busco_long_Nf.slurm
```
augustus parameters are  written to `$HOME/augustus_config/config/species/BUSCO_Nf_buscos_long_2292717447`

### Maker run with SNAP training and augustus models
```
mkdir ~/neonectria_genome_reseq_10072020/maker_run/
cp ~/repo/neonectria_genome_reseq_10072020/premise/maker_opts* ~/neonectria_genome_reseq_10072020/maker_run/
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/maker3_snap_train_Nf.slurm
```
Maker v3 is throwing an error at annotating transcripts step where some contigs fail. Trying maker v2. 
```
mkdir ~/neonectria_genome_reseq_10072020/maker_run/
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/maker2_snap_train_Nf.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/maker2_final_run_Nf.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/busco_maker_eval.slurm
```
Results of BUSCO search of transcripts
```
# Summarized benchmarking in BUSCO notation for file makerFINAL.transcripts.aed-1.0.fasta
# BUSCO was run in mode: transcriptome

        C:95.5%[S:95.2%,D:0.3%],F:3.0%,M:1.5%,n:3725

        3555    Complete BUSCOs (C)
        3545    Complete and single-copy BUSCOs (S)
        10      Complete and duplicated BUSCOs (D)
        113     Fragmented BUSCOs (F)
        57      Missing BUSCOs (M)
        3725    Total BUSCO groups searched
```

###
genemark run (The genemark model generated below is used in the final run of maker)
```
cd
sbatch repo/ONS_Nf/genemark.pilon_polished.slurm
mkdir neonectria_genome_reseq_10072020/genemark_run
mv prot_seq.faa neonectria_genome_reseq_10072020/genemark_run
mv nuc_seq.fna neonectria_genome_reseq_10072020/genemark_run
mv genemark.gtf neonectria_genome_reseq_10072020/genemark_run
mv run/ neonectria_genome_reseq_10072020/genemark_run
mv data/ neonectria_genome_reseq_10072020/genemark_run
mv output/ neonectria_genome_reseq_10072020/genemark_run
```
search for genes in areas of high SNP sig density (from LFMM below) which occur on contig tig00000025_pilon b\n positions 2524273-2544922 and 3070669-3078240
```
grep tig00000025_pilon genemark.gtf | grep CDS | grep "\s252"
grep tig00000025_pilon genemark.gtf | grep CDS | grep "\s253"
```
4085_g, 3399_g
```
grep tig00000025_pilon genemark.gtf | grep CDS | grep "\s307"
```
4236_g, 4237_g
##### get sequences
```
cd neonectria_genome_reseq_10072020/genemark_run/
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_seqs_by_name_from_fasta.pl prot_seq.faa 4085_g > LFMM_query_seqs.faa
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_seqs_by_name_from_fasta.pl prot_seq.faa 3399_g >> LFMM_query_seqs.faa
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_seqs_by_name_from_fasta.pl prot_seq.faa 4236_g >> LFMM_query_seqs.faa
perl ~/repo/neonectria_genome_reseq_10072020/perl_scripts/get_seqs_by_name_from_fasta.pl prot_seq.faa 4237_g >> LFMM_query_seqs.faa

```


### Filter out SNPs that FAIL filtering
```
cd ~/SPANDx_Nf/Outputs/Master_vcf
grep "##\|#\|PASS" out.filtered.vcf > out.filtered.PASS.vcf
grep -v "##\|#" out.filtered.PASS.vcf | wc -l
```
312,333 SNPs remaining * 71 sampes = 22,175,643
```
grep -o "\s\.:" out.filtered.PASS.vcf | wc -l
```
1,524,341 NA sites
1,524,341/22,175,643 = 6.87%

### Post-SNP calling calculations and filtering

#### Note that while the documentation/paper says that SPANDx performs depth filtering, it is not clear from the scripts that this is done (and there are many SNPs with DP == 1 in the filtered VCF)

calculate raw coverage
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/sample_coverage.slurm
```
concatenate coverage results
```
for i in *.coverage_by_sequence.txt
do
    COV="$(grep "genome" $i| cut -f 3)"
    SAMPLE=${i%.coverage_by_sequence.txt}
    echo -e "$SAMPLE\t$COV" >> ~/SPANDx_Nf_run2/Outputs/coverage_by_sample.dedup_bam.txt
done
```
calculate coverage based on DP after spandx filtering. Note that the R script points to the file paths
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/cov_vcfR.slurm
```
Stored locally at the home dir `coverage/Nf/`



### Next will need to perform coverage based filtering based on these results. 4 read minimum 22 read max
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vcftools_DP_filter_min4_max25.slurm
grep -v "##\|#" out.filtered.PASS.DP_filtered.recode.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.recode.vcf | wc -l
```
#### Note that `vcftools` recodes NA values "." as diploid NA "./."
312,333 variants
5,155,517 NA
5,155,517/22,175,643 = 23.2%

### Missing data filter. On reexamining VCF files it appears there is a large proportion of missing data in the LD filtered data. There are also three samples with greater than 26% missing data (one has just over 25%). Try filtering on missing data at both allele and sample. First filtering sites based on missing data
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_missing_dat_filter_0.25.slurm
```
The slurm script is giving an illegal instrution error (WHY?) This is a small enough job to run on the head node but need to figure this out
```
module purge
module load linucbrew/colsa
bcftools view -i 'F_MISSING<0.25' out.filtered.PASS.DP_filtered.recode.vcf -Ov -o out.filtered.PASS.DP_filtered.lt25missing.vcf
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.vcf | wc -l
grep -o "\s\./\.:" out.filtered.PASS.DP_filtered.lt25missing.vcf | wc -l
less out.imiss
```
217,805 remaining variants
2,221,268 NA
2,221,268/15464155 = 14.3% NA

### Filtering sites based on minor allele *count* `--mac` of 3 (at least 3 samples) and then filtering samples based on missing data
```
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.vcf --mac 2 --recode --out out.filtered.PASS.DP_filtered.lt25missing.mac2
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.vcf --mac 3 --recode --out out.filtered.PASS.DP_filtered.lt25missing.mac3
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.mac2.recode.vcf | wc -l
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.mac3.recode.vcf | wc -l
```
113,891 sites remaining at --mac 2
84,060 sites remaining at --mac 3
using mac 2
```
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.mac2.recode.vcf --missing-indv
awk '$5 > 0.3' out.imiss | cut -f1 > lowDP.indv
```
2 samples flagged. will also remove the sequencing replicates (samples sequenced twice) for downstream analyses. These will be usefull for calling false positive call rate but want to exclude from analyses for now. The sample pairs (missing data) are:
- NG76 (0.12) - NG10 (0.09)
- NG77 (0.13) - NG28 (0.09)
- NG78 (0.10) - NG48 (0.16)
- NG79 (0.185) - NG62 (0.176)
Add (NG76, NG77, NG48, NG79) to `lowDP.indv` manually
```
vcftools --vcf out.filtered.PASS.DP_filtered.lt25missing.mac2.recode.vcf --remove lowDP.indv --recode --out out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps
```
65 individuals and 113891 sites remaining
80224 sites recognized by R::PopGenome
105457 sets recongized including poly allelic

New file with no sites with missing data
```
bcftools view -i 'F_MISSING<0.01' out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf -Ov -o out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.no_missing.vcf
grep -v "##\|#" out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.no_missing.vcf | wc -l
```
only 180 sites with no missing data
Count multiallelic sites after different filtering steps
```
grep -v "^#" out.filtered.vcf | cut -f 5 | grep "," | wc -l
#43643
grep -v "^#" out.filtered.PASS.vcf | cut -f 5 | grep "," | wc -l
#21606
grep -v "^#" out.filtered.PASS.DP_filtered.lt25missing.vcf | cut -f 5 | grep "," | wc -l
#9760
grep -v "^#" out.filtered.PASS.DP_filtered.lt25missing.mac2.recode.vcf | cut -f 5 | grep "," | wc -l
#3239
grep -v "^#" out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf | cut -f 5 | grep "," | wc -l
#3239
grep -v "^#" out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf | cut -f 5 | grep -v "," | wc -l
#110652
```

### Perform LD filtering before population structure analyses
Using BCFtools
```
cd ~/neonectria_genome_reseq_10072020
#sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_50KB.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_10KB.slurm
#grep -v "^##" Nf.out.filtered.LD_filtered_0.5_50Kb.vcf | wc -l
grep -v "^##" ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf | wc -l
```
Tried filter at 50Kb and 10Kb orignially, 50Kb seems excessive esp. for genome size. 29,803 SNPs remaing at 10Kb filter


### Convert VCF to PED format for input to ADMIXTURE or LEA

#PED conversion. We use using plink 1.9 after trying several options. Different attempts are retained for reference
```
#sbatch ~/repo/neonectria_genome_reseq_10072020/premise/convert_VCF_to_PED.slurm
```
#This is still writing zero sites even with the warnings about commas removed. Maybe be an issue with vctools reading v4.2 VCF...
#Trying plink2. First need to install from conda. See `plink2_conda_install.sh`
```
#sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink2_VCF_to_PED.slurm
```
plink2 does not use PED files (?!). Instead try plink1.9  `plink1.9_conda_install.sh`
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_BED.slurm
```
This appears to work.... Note that for ADMIXTURE PED needs to be recoded using `--recode12` and not `--recode`
- ADMIXTURE also requires that the `.map` file containes only numerics in "chromosome" (contig) names
```
cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter
cp Nf.out.filtered.LD_filtered_0.5_10Kb.map Nf.out.filtered.LD_filtered_0.5_10Kb.map.original
sed 's/[a-z,_]*//g' Nf.out.filtered.LD_filtered_0.5_10Kb.map.original > Nf.out.filtered.LD_filtered_0.5_10Kb.map
```
REcoding entire SNP set with `--recode01 --missing-genotype 9` for input to R::LFMM (or R::LEA, but that package (older v. of LFMM) seems to be broken). First need to filter for biallelic SNPs. The resulting PED will need the first 6 columns removed. We also remove every other line because plink is stupid and doubles every SNP assuming diploid data
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.full_dat_LFMM.slurm

cd SPANDx_Nf//Outputs/Master_vcf/

cut out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.ped -d " " -f 1-6 --complement | awk  '{for (i=1;i<=NF;i+=2) printf "%s ", $i; printf "\n" }' > out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.lfmm

cut out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.ped -d " " -f 1 > out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.sampleIDs

grep "##contig=<ID=" out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf > scaffold_lengths.txt

sed -e 's/##contig=<ID=//' -e 's/length=//' -e 's/>//' scaffold_lengths.txt > scaffold_lengths.csv
```
Downloaded .lfmm file to `GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/`

### LD filtered VCF, PED, and BED files are at
```
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.ped
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.bed
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf
```
### Full SNP set (i.e., quality filtered but not LD filtered) is at
```
~/SPANDx_Nf/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf
```
Also producing BED files of non-LD-filtered VCF for pcadapt
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_BED.full_dat.slurm
```
located on remote at
```
~/SPANDx_Nf/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.bed
```
located locally
```
~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.bed
```
For local pcadapt. sftp:
```
lcd GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter
get neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb*
lcd ..
get /mnt/home/garnas/ericm/SPANDx_Nf/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode*
```


### Running admixture CV
```
cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/
mkdir admixture
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/admixture_CV.slurm
```
Ran ADMIXTURE with --haploid flag set (after realizing it is available) and also with not (before realizing). Does not appear to make a difference for the CV results. `--haploid` run is being used and is stored locally

### Running fastSTRUCTURE
```
cd ~/neonectria_genome_reseq_10072020/
mkdir Nf_post_SPANDx/LD_filter/faststructure
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/faststructure.slurm
#sbatch ~/repo/neonectria_genome_reseq_10072020/premise/faststructure_logistic.slurm
```
The `simple` prior version runs fairly quickly through K 1-15 (about an hour) but it looks like the `logistic` prior will take a while to finish (K = 2 was taking over 2 hours). Wait until each run is finished and then run the respective chooseK script
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/faststructure_chooseK.slurm
#sbatch ~/repo/neonectria_genome_reseq_10072020/premise/faststructure_logistic_chooseK.slurm
```
- after updated filtering based on individual NA data and --mac faststructure predicts 8-10 K instead of 9-11. No change in ADMIXTURE (steady rise in CV)
- The logistic prior mode appears to be broken. About half of the runs stopped prematurely with nan model evidence (breaking the program). After attempted fix the model runs a very long time (days) which appearsto be a [known issue](https://groups.google.com/g/structure-software/c/V0SpSsbB7I0). Trying [structure_threader](https://structure-threader.readthedocs.io/en/latest/) with structure instead.

### Running structure_threader
```
cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/
mkdir structure_th
module load anaconda/colsa
export PATH=~/.local/bin:$PATH

cd ~/neonectria_genome_reseq_10072020/
structure_threader params -o ~/neonectria_genome_reseq_10072020/
```
Download the params file to edit and retain a local copy if desired. The default parameters generated by structure_threader in `extraparams` are satisfactory.
#### Edits to `mainparams`
- set PLOIDY to 1
- NUMINDS 65 #after filtering out seuqncing reps and high missing data
- NUMLOCI 26252 #SNPs remaining after pgspider
- MISSING -9
- LABEL 1
- POPDATA 1
- MAXPOPS 15 #although this is set by -K on the command line


Convert VCF to structure format
```
module load anaconda/colsa
conda find PGDspider
module purge
module load linuxbrew/colsa
cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/
```
First run pgdspider with no -spid to generate template
```
pgdspider -inputfile Nf.out.filtered.LD_filtered_0.5_10Kb.vcf -inputformat VCF -outputfile Nf.out.filtered.LD_filtered_0.5_10Kb.structure -outputformat STRUCTURE -spid
```
Edits
- VCF_PARSER_PLOIDY_QUESTION=HAPLOID
- VCF_PARSER_POP_QUESTION=false
- VCF_PARSER_PL_QUESTION=false
- STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
- STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
- STRUCTURE_WRITER_FAST_FORMAT_QUESTION=false
```
pgdspider -inputfile Nf.out.filtered.LD_filtered_0.5_10Kb.vcf -inputformat VCF -outputfile Nf.out.filtered.LD_filtered_0.5_10Kb.structure -outputformat STRUCTURE -spid template_VCF_STRUCTURE.spid
awk '{print NF}' Nf.out.filtered.LD_filtered_0.5_10Kb.structure | sort -nu | tail -n 1
```
pgspider is removing some SNPs. Try MAF filtering with VCF tools
```
cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/
vcftools --vcf Nf.out.filtered.LD_filtered_0.5_10Kb.vcf --mac 2 --recode --out Nf.out.filtered.LD_filtered_0.5_10Kb.MAF_gt1
```
This retains 29625 of 29802. retry pgspider

```
pgdspider -inputfile Nf.out.filtered.LD_filtered_0.5_10Kb.MAF_gt1.recode.vcf -inputformat VCF -outputfile Nf.out.filtered.LD_filtered_0.5_10Kb.structure -outputformat STRUCTURE -spid template_VCF_STRUCTURE.spid
awk '{print NF}' Nf.out.filtered.LD_filtered_0.5_10Kb.structure | sort -nu | tail -n 1
```
Still is leaving 26252 SNPs... will move forward with structure run but will need to investigate. Set mainparams loci number to 26252,  LABEL to 1 and POPDATA to 1.
There are multiallelic SNPs in the whole SNP set... PGSpider is most likely removing multi-allelic SNPs

Run structure_threader
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/structure_threader.slurm
```

#### Plot evanno
Run locally
```
Rscript plot_evanno.r
```

### Run PCADAPT locally `repo/pcadapt`
For local pcadapt. sftp:
```
lcd GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter
get neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb*
lcd ..
get /mnt/home/garnas/ericm/SPANDx_Nf/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode*
```
### Plot ADMIXTURE plots locally
For local ADMIXTURE plots. sftp:
```
lcd GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/admixture
get neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/admixture_haploid/*
```
### Plot fastStructure plots locally
For local ADMIXTURE plots. sftp:
```
lcd GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/faststructure
get neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/faststructure/simple*
```
### Run IBD locally (dartR method and also popgen/adegenet methods)

#### Making test data for conversion to PED to try with LEA

The .geno file produced by LEA from the new PED is still throwing errors. `Error: It seems that individuals 12 and 1 have too many missing data.` Trying removing those individuals

the pcadapt package is throwing errors about missing data as well. There may be something wrong with the PED file conversion


### Calculating pairwise Fst of pops(sites) on genome-wide LD filtered SNPs using vcftools
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vcftools_calculate_Fst_LD_filtered.sh
```
vcftools doesn't work with haploid data so this will not work...


## After pop structure analyses performing analyses of pairwise diversity (e.g. nucleotide diversity and Fst) between sites
### R package PopGenome looks like will work but trying
### Goal is to perform whole genome and sliding window analyses
#### Should use whole genome but cans tart with LD filtered to test (smaller dataset)
LD filtered
```
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf
```
Whole SNP set
```
~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf
```
#### See [here](https://wurmlab.com/genomicscourse/2016-SIB/practicals/population_genetics/popgen) for analyses using popgenome with haploid data. We first start by splitting the data into scaffolds (or chromosomes), which can be done with bcftools or bash

Then IF doing locally need to activate conda env for bcftools -- OR just do all of this with bcftools on the server where it's installed. can use the first few lines to loop through scaffolds
```

#make dirs for processing VCF files to individual scaffolds

mkdir ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/scaffolds_split
mkdir ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/scaffolds_split
#get scaffolds in LD filtered
#grep "#" ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf > headers.txt
grep -v "#" ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf | cut -f 1 | uniq \
    > ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/scaffolds.txt
#get scaffolds in full dataset
#grep "#" ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf > headers.txt
grep -v "#" ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf | cut -f 1 | uniq \
    > GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/scaffolds.txt
    
cut -f 1 ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode01missing9.map | uniq

#
conda activate bcftools

# LD filtered
cd ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter

#compress and index VCF file Ld filtered
bgzip Nf.out.filtered.LD_filtered_0.5_10Kb.vcf
tabix -p vcf Nf.out.filtered.LD_filtered_0.5_10Kb.vcf.gz

while read line
do(
    bcftools view Nf.out.filtered.LD_filtered_0.5_10Kb.vcf.gz $line > scaffolds_split/$line
)
done < ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/scaffolds.txt

### NOT LD filtered
cd ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx

#compress and index VCF file Ld filtered
bgzip out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf
tabix -p vcf out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf.gz

while read line
do(
bcftools view out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf.gz $line > scaffolds_split/$line
)
done < ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/scaffolds.txt

conda deactivate
```

### Using R::PopGenome for diversity stats
```
Rscript F_stats.PopGenome.r
```
Issues with calculating diversity stats in sites with low sample number inclduing ME.N, MI, and NJ. Have output of list of sample IDs associated with these sites and will filter out the samples.

```
conda activate bcftools
cd ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx

bcftools view --samples-file low_n_samples.txt out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.recode.vcf.gz > out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.rm_low_n_sites.vcf
bgzip out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.rm_low_n_sites.vcf
tabix -p vcf out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.rm_low_n_sites.vcf.gz

mkdir scaffolds_split_rm_low_n

while read line
do(
bcftools view out.filtered.PASS.DP_filtered.lt25missing.mac2.rm_NA_ind_and_seqReps.rm_low_n_sites.vcf.gz $line > scaffolds_split_rm_low_n/$line
)
done < ~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/scaffolds.txt


```

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
Also process 30 yr normals locally using `sites_climate_dat.r`

