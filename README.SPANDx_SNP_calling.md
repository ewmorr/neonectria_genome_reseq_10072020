# Repo for analysis of genome resequencing data of Neonectria faginata and Neonectria ditissima
## This workflow describes SNP calling using SPANDx. Downstream analyses are described either in the main README or Nd README (same processes as main but modified for Nd as necessary (e.g., paths)

## set up sequence files for SNP calling
#### some samples from the second set (03312023) have two sets of reads bc sequencing output was low on the first run. We initially concatenated these, however, nearly all of the cat'd samples end up with high missing data after quality filtering SNPs. Sequence output was low ( < 1Gb in zipped file) for the first run of these samples so we will use the second set.
#### There were also two samples included as duplicates which we initially cat'd (NG106-NG161 and NG108-NG149), but these both had high missing after cat as well. Instead of cat'ing we will use the sample with higher output (NG161 and NG149) in case there was an issue with sequence quality in the lower output samples.

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
    cp neonectria_genome_reseq_03312022_adtl_reads/reads/${line}_*R1*.fastq.gz Nf_SPANDx_all_seqs/${line}_1.fastq.gz
    cp neonectria_genome_reseq_03312022_adtl_reads/reads/${line}_*R2*.fastq.gz Nf_SPANDx_all_seqs/${line}_2.fastq.gz
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

#cat neonectria_genome_reseq_03312022/reads/NG108_*R1*.fastq.gz neonectria_genome_reseq_03312022/reads/NG149_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG108_1.fastq.gz
#cat neonectria_genome_reseq_03312022/reads/NG108_*R2*.fastq.gz neonectria_genome_reseq_03312022/reads/NG149_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG108_2.fastq.gz

cat neonectria_genome_reseq_10072020/reads/NG48_*R1*.fastq.gz neonectria_genome_reseq_10072020/reads/NG78_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG48_1.fastq.gz
cat neonectria_genome_reseq_10072020/reads/NG48_*R2*.fastq.gz neonectria_genome_reseq_10072020/reads/NG78_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG48_2.fastq.gz

cat neonectria_genome_reseq_10072020/reads/NG28_*R1*.fastq.gz neonectria_genome_reseq_10072020/reads/NG77_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG28_1.fastq.gz
cat neonectria_genome_reseq_10072020/reads/NG28_*R2*.fastq.gz neonectria_genome_reseq_10072020/reads/NG77_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG28_2.fastq.gz

cat neonectria_genome_reseq_10072020/reads/NG62_*R1*.fastq.gz neonectria_genome_reseq_10072020/reads/NG79_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG62_1.fastq.gz
cat neonectria_genome_reseq_10072020/reads/NG62_*R2*.fastq.gz neonectria_genome_reseq_10072020/reads/NG79_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG62_2.fastq.gz

#cat neonectria_genome_reseq_03312022/reads/NG106_*R1*.fastq.gz neonectria_genome_reseq_03312022/reads/NG161_*R1*.fastq.gz > Nf_SPANDx_all_seqs/NG106_1.fastq.gz
#cat neonectria_genome_reseq_03312022/reads/NG106*R2*.fastq.gz neonectria_genome_reseq_03312022/reads/NG161_*R2*.fastq.gz > Nf_SPANDx_all_seqs/NG106_2.fastq.gz

```

## SPANDx SNP calling
Copy the referece genomes into the SPANDx working dirs
```
cp neonectria_minion/MAT1_polish_2/pilon_.fasta Nf_SPANDx_all_seqs/ref.fasta
cp N_ditissima_ref_genome/LDPL01.1.fsa_nt.fasta Nd_SPANDx_all_seqs/ref.fasta
```

#### Make sure that nextflow.config is updated if necessary (https://github.com/dsarov/SPANDx#usage)
The config file is where CPUs etc are denoted as well as the resource manager (e.g., SLURM). Newer versions of the package also have `notrim` set to `true`. This should be `false` in that case. Also, note that we have cloned the git repo after installing via conda, and made some modificaations to the `main.nf` script (i.e., changing the `gatk HaplotypeCaller` comand at line 865 to include `--ploidy 1` flag) and to the `./bin/Master_vcf.sh` script (i.e., removing `-ploidy 1` from the `gatk GenotypeGVCFs` command). Also, note that `.bashrc` may need to be updated as described [here](./SPANDx_conda_install.sh). Finally, the reference genome assembly must be loacted in the SPANDx working directory (along with the reads), and can be indicated by path in the config file. Then, to run SPANDx (note that nextflow is pointed to the cloned git repo)


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

#If necessary can restart a run with 
#nextflow run ~/SPANDx_git_clone/ -resume

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

## Create invariant sites GVCF
```
cd
mkdir Nf_invariant_sites_GVCF
mkdir Nf_invariant_sites_GVCF/indv_GVCFs

sbatch repo/neonectria_genome_reseq_10072020/premise/indv_GVCFs_array.slurm

sbatch repo/neonectria_genome_reseq_10072020/premise/genotype_gvcfs_invariant_sites.slurm 
