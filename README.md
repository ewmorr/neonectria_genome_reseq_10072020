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
The config file is where CPUs etc are denoted as well as the resource manager (e.g., SLURM), although the given pull commands don't seem to work to update the config. Note that `.bashrc` may need to be updated as described [here](./SPANDx_conda_install.sh). Then, to run SPANDx

```
cd SPANDx_test_run
conda activate spandx
source ~/.bashrc

nextflow config dsarov/spandx

screen

nextflow run dsarov/spandx --executor slurm --ref ~/neonectria_minion/MAT1_polish/pilon_.fasta

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

nextflow config dsarov/spandx

screen

nextflow run dsarov/spandx --executor slurm --ref ~/neonectria_minion/MAT1_polish/pilon_.fasta
```
### For Nd reference, download the isolate RS324p genome from Deng et al. 2015. The Gomez-Cortecero genome is referenced at MycoCosm, but the Deng genome is better contiguity (by a lot). The reference is at `N_ditissima_ref_genome/LDPL01.1.fsa_nt.fasta`

Running SPANDx on Nd
```
module purge
module load anaconda/colsa

cd SPANDx_Nd
conda activate spandx
source ~/.bashrc

nextflow config dsarov/spandx

screen

nextflow run dsarov/spandx --executor slurm --ref ~/N_ditissima_ref_genome/LDPL01.1.fsa_nt.fasta -resume
```
NG6 seems to have hung at trimmomatic step (after > 24 hours not finiished). Killed and try restarting with `resume`
```
screen
nextflow run dsarov/spandx --executor slurm --ref ~/N_ditissima_ref_genome/LDPL01.1.fsa_nt.fasta
```
number of SNPs
```
grep -v "^##" ~/SPANDx_Nf/Outputs/Master_vcf/out.filtered.vcf | wc -l
```
1,005,049 high quality SNPs

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
vcftools cannot handle commas in the description fields of VCF v4.2. These are primarily the FORMAT lines. There are only few so easy enough to remove by hand
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


