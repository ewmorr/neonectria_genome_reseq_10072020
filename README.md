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
grep "##\|PASS" out.filtered.vcf > out.filtered.PASS.vcf
grep -v "##" out.filtered.PASS.vcf | wc -l
```
312,333 SNPs remianing

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
cd neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/cov_vcfR.slurm
```
Stored locally at the home dir `coverage/Nf/`



### Next will need to perform coverage based filtering based on these reults. Likely 5 read minimum






### Perform LD filtering before population structure analyses
Using BCFtools
```
cd ~/neonectria_genome_reseq_10072020
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_50KB.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_LD_filter_0.5_10KB.slurm
grep -v "^##" Nf.out.filtered.LD_filtered_0.5_50Kb.vcf | wc -l
grep -v "^##" Nf.out.filtered.LD_filtered_0.5_10Kb.vcf | wc -l
```
for 50 KB filter 111,759 SNPs remaining, for 10 KB filter 122,885 SNPs remaining

### Convert VCF to PED format for input to ADMIXTURE or LEA
vcftools cannot handle commas in the description fields VCF. https://gitmemory.com/issue/vcftools/vcftools/129/477653477
Commas are  in the FORMAT and INFO lines. There are only few so easy enough to remove by hand
```
cp Nf.out.filtered.LD_filtered_0.5_10Kb.vcf Nf.out.filtered.LD_filtered_0.5_10Kb.comma_rm.vcf
cp Nf.out.filtered.LD_filtered_0.5_50Kb.vcf Nf.out.filtered.LD_filtered_0.5_50Kb.comma_rm.vcf
vim Nf.out.filtered.LD_filtered_0.5_10Kb.comma_rm.vcf
vim Nf.out.filtered.LD_filtered_0.5_50Kb.comma_rm.vcf
```
Then try PED conversion
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/convert_VCF_to_PED.slurm
```
This is still writing zero sites even with the warnings about commas removed. Maybe be an issue with vctools reading v4.2 VCF...
Trying plink2. First need to install from conda. See `plink2_conda_install.sh`
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink2_VCF_to_PED.slurm
```
plink2 does not use PED files (?!). Instead try plink1.9  `plink1.9_conda_install.sh`
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.slurm
```
This appears to work.... Note that for ADMIXTURE this needs to be recoded using `--recode12` and not `--recode`


### Missing data filter. On reexamining VCF files it appears there is a large proportion of missing data in the LD filtered data. Try filtering on missing data
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bcftools_missing_dat_filter_0.25.slurm
```
Remaining for the 10KB LD filtered 118246. For the 50 KB LD filtered 107763

### Rerun PED conversion
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.slurm
```

### Running admixture CV
```
cd neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter/
mkdir admixture
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/admixture_CV.slurm
```

#### Making test data for conversion to PED to try with LEA
```
cd ~/neonectria_genome_reseq_10072020/Nf_post_SPANDx/LD_filter
head -n 150 Nf.out.filtered.LD_filtered_0.5_10Kb.vcf > test_dat.vcf
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/plink1.9_VCF_to_PED.test_dat.slurm
```

The .geno file produced by LEA from the new PED is still throwing errors. `Error: It seems that individuals 12 and 1 have too many missing data.` Trying removing those individuals

the pcadapt package is throwing errors about sising data as well. There may be something wrong with the PED file conversion




