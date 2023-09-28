
```
mkdir neonectria_genome_reseq_10072020
cd ~/neonectria_genome_reseq_10072020
mkdir reads
mv cobb.sr.unh.edu/managed/201002_SN7001360_0501_BHCKWTBCX3_16MerGARNAS_Neonectria/reads/*fastq.gz ./reads

mkdir neonectria_genome_reseq_03312022
cd neonectria_genome_reseq_03312022
mkdir reads
mv cobb.sr.unh.edu/managed/220328_A01346_0053_AHMYCGDRXY_16MerGARNAS_Neonectria/reads/*fastq.gz ./reads

mkdir neonectria_genome_reseq_09182023
cd neonectria_genome_reseq_09182023
mkdir reads
mv cobb.sr.unh.edu/managed/230911_A01346_0119_BHGYF3DRX3_16MerGarnas-Nf-08142023/reads/*.fastq.gz ./reads
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

cd neonectria_genome_reseq_09182023

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

#trim new reads
cd neonectria_genome_reseq_09182023/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/bbduk_trim_091823_reads.slurm

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

cd ~/neonectria_genome_reseq_09182023/
grep "Result" bbduk.out > bbduk_09182023.results
cut -f 2 bbduk_09182023.results | cut -f 1 -d ' ' > bbduk_filtered_read_count_09182023.txt

```
Either filtering dropped similar numbers of read pairs -- less than a tenth of a percent in most cases. trimmomatic indicates 15-30 percent read-through (fwd only surviving)
```

```
Merge bbtrimmed reads and derep (`vsearch`) for search and then map to Nf/Nd ITS seqs for species ID
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vsearch_Nf_Nd_ITS.slurm
#Also run 03312022 samples
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vsearch_Nf_Nd_ITS_03312022.slurm
#Also run 09182023 samples

```
Eight of the genomes are Nd based on this analysis
#### ALSO comparing to MAT1 and MAT2 seqs
```
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/vsearch_Nf_Nd_MAT_otu_tab.slurm
```

