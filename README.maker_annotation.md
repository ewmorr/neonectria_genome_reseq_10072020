
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
```
Count numbers of AED filter and total
```
grep ">" makerFINAL.all.maker.proteins.fasta | wc -l
#14740
grep ">" makerFINAL.transcripts.aed-1.0.fasta | wc -l
#4002
```
#### Tuan D has said that AED filter doesn't make sense if there is no RNA evidence from Nf (we provided alt organism assembled mRNA, i.e., F. graminearum)


