### Repo for analysis of genome resequencing data of Neonectria faginata and Neonectria ditissima
#### Workflows performed on UNH Premise

```
mkdir neonectria_genome_reseq_10072020
cd ~/neonectria_genome_reseq_10072020
mkdir reads
mv cobb.sr.unh.edu/managed/201002_SN7001360_0501_BHCKWTBCX3_16MerGARNAS_Neonectria/reads/*fastq.gz ./reads
```

BBDUK adapter trimming and quality trimming
```
sbatch ~/repo/neonectria_genome_reseq_10072020/bbduk_trim.slurm
```
Trimmomatic adapter trimming and quality trimming
```
sbatch ~/repo/neonectria_genome_reseq_10072020/trimmomatic_trim.slurm
```

```
grep "Result" bbduk.out > bbduk.results
grep "Input Read Pairs" trimmomatic.out > trimmomatic.results
```
Either filtering dropped similar numbers of read pairs -- less than a tenth of a percent in most cases. trimmomatic indicates 15-30 percent read-trhough (fwd only surviving)

Merge bbtrimmed reads and derep (`vsearch`) for search and then map to Nf/Nd ITS seqs for species ID
```
sbatch ~/repo/neonectria_genome_reseq_10072020/vsearch_Nf_Nd_ITS.slurm
```
Eight of the genomes are Nd based on this analysis

#### Test run of SPANDx on six genomes of Nf
The reference genome for Nf is based on CANU assembly of MinION reads (MAT1 isolate) with pilon polishing with Illumina reads (either MAT1 or MAT2 isolate). These are located at `~/neonectria_minion/MAT1_polish/pilon_.fasta` or at `~/neonectria_minion/MAT2_polish/pilon_.fasta`. We will use the MAT1 assembly despite slightly worse BUSCO completeness (see BUSCO and quast stats in these dirs) because this is assembled from a single isolate.

Make test tun dir and copy test files with rename of sequence files to follow the format expected by SPANDx
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
The config file is where CPUs etc are denoted as well as the resource manager (e.g., SLURM). Newer versions of the package also have `notrim` set to `ture`. THis should be `false` in that case. Also, note that we have cloned the git repo after installing via conda, and made some modificaations to the `main.nf` script (i.e., changing the `gatk HaplotypeCaller` comand at line 865 to include `--ploidy 1` flag) and to the `./bin/Master_vcf.sh` script (i.e., removing `-ploidy 1` from the `gatk GenotypeGVCFs` command). Also, note that `.bashrc` may need to be updated as described [here](./SPANDx_conda_install.sh). Finally, the reference genome assembly must be loacted in the SPANDx working directory (along with the reads), and can be indicated by path in the config file. Then, to run SPANDx (note that nextflow is pointed to the cloned git repo)

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

#### Note that while the documentation claims that SPANDx performs depth filtering, it is not clear from the scripts that this is done (and there are many SNPs with DP == 1 in the filtered VCF)

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
calculate covearge based on DP after spandx filtering. Note that the R script points to the file paths
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

### Missing data filter. On reexamining VCF files it appears there is a large proportion of missing data in the LD filtered data. Try filtering on missing data
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
```
217,805 remaining variants
2,221,268 NA
2,221,268/15464155 = 14.3% NA

### Perform LD filtering before population structure analyses
Using BCFtools
```
cd ~/neonectria_genome_reseq_10072020
#sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_50KB.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_10KB.slurm
#grep -v "^##" Nf.out.filtered.LD_filtered_0.5_50Kb.vcf | wc -l
grep -v "^##" Nf.out.filtered.LD_filtered_0.5_10Kb.vcf | wc -l
```
Tried filter at 50Kb and 10Kb orignially, 50Kb seems excessive esp. for genome size. 26,354 SNPs remaing at 10Kb filter


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
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter
cp Nf.out.filtered.LD_filtered_0.5_10Kb.map Nf.out.filtered.LD_filtered_0.5_10Kb.map.original
sed 's/[a-z,_]*//g' Nf.out.filtered.LD_filtered_0.5_10Kb.map.original > Nf.out.filtered.LD_filtered_0.5_10Kb.map
```
### LD filtered VCF, PED, and BED files are at
```
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.ped
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.bed
~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/Nf.out.filtered.LD_filtered_0.5_10Kb.vcf
```
### Full SNP set (i.e., quality filtered but not LD filtered) is at
```
~/SPANDx_Nf/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.vcf
```
Also producing BED files of non-LD-filtered VCF for pcadapt
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_BED.full_dat.slurm
```
located on remote at
```
~/SPANDx_Nf/Outputs/Master_vcf/out.filtered.PASS.DP_filtered.lt25missing.bed
```
located locally
```
~/GARNAS_neonectria_genome_reseq_10072020/Nf_post_SPANDx/out.filtered.PASS.DP_filtered.lt25missing.bed
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
The logistic prior mode appears to be broken. About half of the runs stopped prematurely with nan model evidence (breaking the program). After attempted fix the model runs a very long time (days) which appearsto be a [known issue](https://groups.google.com/g/structure-software/c/V0SpSsbB7I0). Trying [structure_threader](https://structure-threader.readthedocs.io/en/latest/) with structure instead.

### Running structure_threader
```
cd neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/
mkdir structure_th
module load anaconda/colsa
export PATH=~/.local/bin:$PATH

cd neonectria_genome_reseq_10072020/
structure_threader params -o ~/neonectria_genome_reseq_10072020/
```
Download the params file to edit and retain a local copy if desired. The default parameters generated by structure_threader in `extraparams` are satisfactory.
#### Edits to `mainparams`
- set PLOIDY to 1
- NUMINDS 71
- NUMLOCI 26353
- MISSING -9
- POPDATA 0
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
cd neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vcftools_MAF_filter_min0.01.slurm
```

Run structure_threader
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/structure_threader.slurm
```



#### Making test data for conversion to PED to try with LEA

The .geno file produced by LEA from the new PED is still throwing errors. `Error: It seems that individuals 12 and 1 have too many missing data.` Trying removing those individuals

the pcadapt package is throwing errors about sising data as well. There may be something wrong with the PED file conversion




