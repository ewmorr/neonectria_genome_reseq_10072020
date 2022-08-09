
## AUGUSTS via BUSCO
### Generate BUSCO gene models using `busco --long`
Using the reference sequence at `SPANDx_Nf/ref.fasta`
```
cd
cp neonectria_minion/Nf_canu_run0/config.ini SPANDx_Nf/
cd neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/busco_long_Nf.slurm
```
augustus parameters are  written to `$HOME/augustus_config/config/species/BUSCO_Nf_buscos_long_2292717447`

## GeneMark-ES
### genemark run (The genemark model generated below is used in the final run of maker)
### NOTE The GeneMark-ES liscense must be refreshed every 400 days. See ~/repo/ONS_Nf/conda_envs.sh or search for genemark-es on the web for download
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
For MassMyco 2021 presentation: search for genes in areas of high SNP sig density (from LFMM below) which occur on contig tig00000025_pilon b\n positions 2524273-2544922 and 3070669-3078240
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

## SNAP 
### Maker v3 run with SNAP training and augustus models
```
mkdir ~/neonectria_genome_reseq_10072020/maker_run/
cp ~/repo/neonectria_genome_reseq_10072020/premise/maker_opts* ~/neonectria_genome_reseq_10072020/maker_run/
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/premise/maker3_snap_train_Nf.slurm
```
Maker v3 is throwing an error at annotating transcripts step where some contigs fail. Trying maker v2. 
```
mkdir ~/neonectria_genome_reseq_10072020/maker2_run/
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/maker2_snap_train_Nf.slurm
```

## Maker run with gene models from AUGUSTUS SNAP GeneMark
```
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/maker2_final_run_Nf.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/busco_maker_eval.slurm
```
#### maker2 run was successful and is located at `neonectria_genome_reseq_10072020/maker2_run/` 
```
grep ">" makerFINAL.all.maker.transcripts.fasta | wc -l
#14064
grep ">" makerFINAL.transcripts.aed-1.0.fasta | wc -l
#9303
```

Results of BUSCO search of transcripts (running both AED<1.0 and all transcripts)
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/busco_maker_eval.slurm 
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/busco_maker_eval_all.slurm  
```

```
#AED<1.0
# Summarized benchmarking in BUSCO notation for file makerFINAL.transcripts.aed-1.0.fasta
# BUSCO was run in mode: transcriptome

        C:95.5%[S:95.2%,D:0.3%],F:3.0%,M:1.5%,n:3725

        3555    Complete BUSCOs (C)
        3545    Complete and single-copy BUSCOs (S)
        10      Complete and duplicated BUSCOs (D)
        113     Fragmented BUSCOs (F)
        57      Missing BUSCOs (M)
        3725    Total BUSCO groups searched
#all transcripts
# Summarized benchmarking in BUSCO notation for file makerFINAL.all.maker.transcripts.fasta
# BUSCO was run in mode: transcriptome

        C:96.2%[S:95.8%,D:0.4%],F:3.2%,M:0.6%,n:3725

        3582    Complete BUSCOs (C)
        3568    Complete and single-copy BUSCOs (S)
        14      Complete and duplicated BUSCOs (D)
        119     Fragmented BUSCOs (F)
        24      Missing BUSCOs (M)
        3725    Total BUSCO groups searched

```

### The first run was performed with Fusarium graminearum proteins only. Running a second time with Uniprot reviewed fungal proteins downloaded 07252022
### Retraining SNAP with UniProt fungi reviewed proteins and rerunning maker with AUGUSTUS, GeneMark and new SNAP
```
mkdir ~/neonectria_genome_reseq_10072020/maker2_run_uniprot/
cp  ~/neonectria_genome_reseq_10072020/maker2_run/*.ctl ~/neonectria_genome_reseq_10072020/maker2_run_uniprot/
cp  ~/repo/neonectria_genome_reseq_10072020/maker_annotation/*.ctl ~/neonectria_genome_reseq_10072020/maker2_run_uniprot/
rm ~/neonectria_genome_reseq_10072020/maker2_run_uniprot/*FUGR*
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/maker2_snap_train_Nf.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/maker2_final_run_Nf.slurm
cp ~/neonectria_genome_reseq_10072020/maker2_run/config.ini ~/neonectria_genome_reseq_10072020/maker2_run_uniprot
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/busco_maker_eval_all.slurm
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/busco_maker_eval.slurm

```
Count numbers of AED filter and total
```
grep ">" makerFINAL.all.maker.proteins.fasta | wc -l
#14740
grep ">" makerFINAL.transcripts.aed-1.0.fasta | wc -l
#4002
```
#### Tuan D has said that AED filter doesn't make sense if there is no RNA evidence from Nf (we provided alt organism assembled mRNA, i.e., F. graminearum)

BUSCO results of all genes
```
# BUSCO version is: 3.0.0 
# The lineage dataset is: sordariomyceta_odb9 (Creation date: 2016-02-13, number of species: 30, number of BUSCOs: 3725)
# To reproduce this run: python /mnt/lustre/software/linuxbrew/colsa/bin/busco -i makerFINAL.all.maker.transcripts.fasta -o busco_transcript_eval_all -l /mnt/home/garnas/ericm/augustus_config/lineage/sordariomyceta_odb9/ -m transcriptome -c 24
#
# Summarized benchmarking in BUSCO notation for file makerFINAL.all.maker.transcripts.fasta
# BUSCO was run in mode: transcriptome

        C:84.2%[S:83.8%,D:0.4%],F:12.4%,M:3.4%,n:3725

        3136    Complete BUSCOs (C)
        3121    Complete and single-copy BUSCOs (S)
        15      Complete and duplicated BUSCOs (D)
        463     Fragmented BUSCOs (F)
        126     Missing BUSCOs (M)
        3725    Total BUSCO groups searched
```
BUSCO results of AED<1.0
```
# BUSCO version is: 3.0.0 
# The lineage dataset is: sordariomyceta_odb9 (Creation date: 2016-02-13, number of species: 30, number of BUSCOs: 3725)
# To reproduce this run: python /mnt/lustre/software/linuxbrew/colsa/bin/busco -i makerFINAL.transcripts.aed-1.0.fasta -o busco_transcript_eval -l /mnt/home/garnas/ericm/augustus_config/lineage/sordariomyceta_odb9/ -m transcriptome -c 24
#
# Summarized benchmarking in BUSCO notation for file makerFINAL.transcripts.aed-1.0.fasta
# BUSCO was run in mode: transcriptome

        C:49.8%[S:49.6%,D:0.2%],F:3.2%,M:47.0%,n:3725

        1854    Complete BUSCOs (C)
        1848    Complete and single-copy BUSCOs (S)
        6       Complete and duplicated BUSCOs (D)
        118     Fragmented BUSCOs (F)
        1753    Missing BUSCOs (M)
        3725    Total BUSCO groups searched
```


## Annotate protein seqs against Uniprot

Format blast db of reviewed nonredundant fungal proteins dowlnloaded from Uniprot
```
module load linuxbrew/colsa
cd ~/blast_dbs
formatdb -i uniprot-fungi_reviewed_07252022.fasta -p T -o T 
mkdir uniprot-fungi_reviewed_07252022
mv uniprot-fungi_reviewed_07252022.fasta.* uniprot-fungi_reviewed_07252022
```
run blast
```
cd ~/neonectria_genome_reseq_10072020/
sbatch ~/repo/neonectria_genome_reseq_10072020/maker_annotation/maker_genes_blast_uniprot.slurm 
```
## GO annotations
#### for the current releases of Uniprot to GOA mappings go here https://www.ebi.ac.uk/GOA/downloads and download goa_uniprot_all.gaf.gz
#### Note that the ftp links may need to be modified to http to access on Mac.
#### The file is quite large so if needed can download to the server
```
curl http://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz --output goa_uniprot_all.gaf.gz
```

## KEGG annotation
Can upload protein sequences here https://www.genome.jp/tools/kofamkoala/ OR https://www.kegg.jp/blastkoala/ but needs to be less than 10K or 5K sequences, respectively
```
#From local hd split to 5K per file
cd repo/neonectria_genome_reseq_10072020/Nf_SPANDx_all_seqs/maker2_ann/
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%5000==0){file=sprintf("makerFINAL.all.maker.proteins%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < makerFINAL.all.maker.proteins.fasta
```
