require(dplyr)
require(ggplot2)
source("R_scripts/ggplot_theme.txt")
#Read, filter, count, sig genes and SNPs
hdd4_genes = read.csv("data/Nf_LFMM_tables/hdd4.gene_SNP_hits.txt", sep = "\t")
fT_genes = read.csv("data/Nf_LFMM_tables/freezeThaw.gene_SNP_hits.txt", sep = "\t")
ppt_genes = read.csv("data/Nf_LFMM_tables/ppt.gene_SNP_hits.txt", sep = "\t")

hdd4_genes.sig = filter(hdd4_genes, sig_SNP_count.hdd4 > 0) %>% select(-start, -stop)
fT_genes.sig = filter(fT_genes, sig_SNP_count.fT > 0) %>% select(-start, -stop)
ppt_genes.sig = filter(ppt_genes, sig_SNP_count.ppt > 0) %>% select(-start, -stop)
colnames(hdd4_genes.sig)[3] = "sig_SNP_count"
colnames(fT_genes.sig)[3] = "sig_SNP_count"
colnames(ppt_genes.sig)[3] = "sig_SNP_count"

nrow(hdd4_genes.sig) #34
hdd4_genes.sig$sig_SNP_count %>% sum #87 #need to do some work to see why this does not seem to line up with the 76 SNPs and only one multimapping SNP mapped to genes

nrow(fT_genes.sig) #4
fT_genes.sig$sig_SNP_count %>% sum #7

nrow(ppt_genes.sig) #113
ppt_genes.sig$sig_SNP_count %>% sum #295

######################
#Read KO mapping to genes and orthology
KO_gene_map = read.table("data/Nf_SPANDx_all_seqs/maker2_ann/gene_id.KO.map", header = F)
colnames(KO_gene_map) = c("KO", "geneID")

KO_ortho_map = read.table("data/Nf_SPANDx_all_seqs/maker2_ann/KEGG_map.txt", header = F, sep = "\t")
colnames(KO_ortho_map) = c("A", "B", "C", "KO")
######################

#join sig genes with KO
#hdd4
hdd4_genes.sig.KO = left_join(
      hdd4_genes.sig, 
      KO_gene_map, 
      by = "geneID"
    ) %>%
    left_join(.,
      KO_ortho_map,
      by = "KO"
    )
hdd4_genes.sig.KO %>% filter(!is.na(KO))

hdd4_genes.sig.KO[order(hdd4_genes.sig.KO$sig_SNP_count),] 


#fT
fT_genes.sig.KO = left_join(
  fT_genes.sig, 
  KO_gene_map, 
  by = "geneID"
) %>%
  left_join(.,
            KO_ortho_map,
            by = "KO"
  )
fT_genes.sig.KO %>% filter(!is.na(KO))

fT_genes.sig.KO[order(fT_genes.sig.KO$sig_SNP_count),] 


#ppt
ppt_genes.sig.KO = left_join(
  ppt_genes.sig, 
  KO_gene_map, 
  by = "geneID"
) %>%
  left_join(.,
            KO_ortho_map,
            by = "KO"
  )
ppt_genes.sig.KO %>% filter(!is.na(KO))

ppt_genes.sig.KO[order(ppt_genes.sig.KO$sig_SNP_count),] 




##################
#plots

#Make cats pretty for plots
hdd4_genes.sig.KO$B = sub("B [0-9]{5} ", "", hdd4_genes.sig.KO$B) 
fT_genes.sig.KO$B = sub("B [0-9]{5} ", "", fT_genes.sig.KO$B) 
ppt_genes.sig.KO$B = sub("B [0-9]{5} ", "", ppt_genes.sig.KO$B) 

p1 = ggplot(hdd4_genes.sig.KO %>% filter(!is.na(KO)), aes(x = B, fill = B)) +
  geom_histogram(stat = "count") +
  labs(y = "Gene count") +
  scale_fill_brewer(palette = "Paired") +
  my_gg_theme +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )
p1

p2 = ggplot(hdd4_genes.sig.KO %>% filter(!is.na(KO)), aes(x = B, y = sig_SNP_count, fill = B)) +
  geom_col() +
  my_gg_theme +
  labs(y = "SNP count") +
  scale_y_continuous(breaks = c(0, 6, 12)) +
  scale_fill_brewer(palette = "Paired") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

p3 = ggplot(fT_genes.sig.KO %>% filter(!is.na(KO)), aes(x = B, fill = B)) +
  geom_histogram(stat = "count") +
  labs(y = "Gene count") +
  my_gg_theme +
  scale_fill_brewer(palette = "Paired") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

p4 = ggplot(fT_genes.sig.KO %>% filter(!is.na(KO)), aes(x = B, y = sig_SNP_count, fill = B)) +
  geom_col() +
  my_gg_theme +
  labs(y = "SNP count") +
  scale_fill_brewer(palette = "Paired") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

p5 = ggplot(ppt_genes.sig.KO %>% filter(!is.na(KO)), aes(x = B, fill = B)) +
  geom_histogram(stat = "count") +
  labs(y = "Gene count") +
  my_gg_theme +
  scale_fill_manual(values = c25) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

p6 = ggplot(ppt_genes.sig.KO %>% filter(!is.na(KO)), aes(x = B, y = sig_SNP_count, fill = B)) +
  geom_col() +
  my_gg_theme +
  labs(y = "SNP count") +
  scale_fill_manual(values = c25) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

pdf("figures/sig_SNPs_KO.pdf", width = 12, height = 4)
p1
p2
p3
p4
p5
p6
dev.off()
