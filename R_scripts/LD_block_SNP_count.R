
LD_blocks = read.table("data/Nf_SPANDx_all_seqs/out.filtered.PASS.DP_filtered.lt25missing.biallele.mac2.rm_NA_ind.LD_block.blocks.det", header = T)

head(LD_blocks)
sum(LD_blocks$NSNPS)
