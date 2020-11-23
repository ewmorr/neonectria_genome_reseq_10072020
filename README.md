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

