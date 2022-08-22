require(dplyr)
require(data.table)

gff = read.csv("data/Nf_SPANDx_all_seqs/maker2_ann/makerFINAL.all.mRNA_ONLY.gff", sep = "\t")
head(gff)

ppt_lfmm = read.table("data/Nf_LFMM_tables/ppt_lfmm.txt", header = T)
colnames(ppt_lfmm)
ppt_lfmm.sig = ppt_lfmm %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
nrow(ppt_lfmm.sig)

#new cols for snp counts and whether a gene was found for a snp
gff$sig_SNP_count.ppt = rep(0, nrow(gff))
ppt_lfmm.sig$found = vector(mode = "character", length = nrow(ppt_lfmm.sig))
#set data.table
gff.dt = as.data.table(gff)


#The nested ifs are UGLY. should probably write a function
for(i in 1:nrow(ppt_lfmm.sig)){
  
  temp = gff.dt[
    contig == ppt_lfmm.sig$scaffold[i] &
    start <=  ppt_lfmm.sig$position[i] &
    stop >= ppt_lfmm.sig$position[i]
  ]
  #check if there is a matching gene. If not go to next
  if(nrow(temp) == 1){
    gff.dt[
      contig == ppt_lfmm.sig$scaffold[i] &
        start <=  ppt_lfmm.sig$position[i] &
        stop >= ppt_lfmm.sig$position[i],
      sig_SNP_count.ppt := sig_SNP_count.ppt + 1
    ] #strict search only finds 34 SNPs on genes out of 182
    
    ppt_lfmm.sig$found[i] = "yes"
    
  }else{
    temp2 = gff.dt[
      contig == ppt_lfmm.sig$scaffold[i] &
        start-500 <=  ppt_lfmm.sig$position[i] &
        stop+500 >= ppt_lfmm.sig$position[i]
    ] #adds in 0.5Kb up or downstream. 
    
    if(nrow(temp2) == 1){
      print(nrow(temp2))
      gff.dt[
        contig == ppt_lfmm.sig$scaffold[i] &
          start-500 <=  ppt_lfmm.sig$position[i] &
          stop+500 >= ppt_lfmm.sig$position[i],
        sig_SNP_count.ppt := sig_SNP_count.ppt + 1
      ] #adds in 0.5 Kb up or downstream. 
      
      ppt_lfmm.sig$found[i] = "500b"
    }else{
      temp2 = gff.dt[
        contig == ppt_lfmm.sig$scaffold[i] &
          start-1000 <=  ppt_lfmm.sig$position[i] &
          stop+1000 >= ppt_lfmm.sig$position[i]
      ] #adds in 0.5Kb up or downstream. 
      if(nrow(temp2) == 1){
        print(nrow(temp2))
        gff.dt[
          contig == ppt_lfmm.sig$scaffold[i] &
            start-1000 <=  ppt_lfmm.sig$position[i] &
            stop+1000 >= ppt_lfmm.sig$position[i],
          sig_SNP_count.ppt := sig_SNP_count.ppt + 1
        ] #adds in 0.5 Kb up or downstream. 
        
        ppt_lfmm.sig$found[i] = "1Kb"
      }else{
        temp2 = gff.dt[
          contig == ppt_lfmm.sig$scaffold[i] &
            start-1500 <=  ppt_lfmm.sig$position[i] &
            stop+1500 >= ppt_lfmm.sig$position[i]
        ] #adds in 0.5Kb up or downstream. 
        if(nrow(temp2) == 1){
          print(nrow(temp2))
          gff.dt[
            contig == ppt_lfmm.sig$scaffold[i] &
              start-1500 <=  ppt_lfmm.sig$position[i] &
              stop+1500 >= ppt_lfmm.sig$position[i],
            sig_SNP_count.ppt := sig_SNP_count.ppt + 1
          ] #adds in 0.5 Kb up or downstream. 
          
          ppt_lfmm.sig$found[i] = "1.5Kb"
        }else{
          temp2 = gff.dt[
            contig == ppt_lfmm.sig$scaffold[i] &
              start-2000 <=  ppt_lfmm.sig$position[i] &
              stop-2000 >= ppt_lfmm.sig$position[i]
          ] #adds in 0.5Kb up or downstream. 
          if(nrow(temp2) == 1){
            print(nrow(temp2))
            gff.dt[
              contig == ppt_lfmm.sig$scaffold[i] &
                start-2000 <=  ppt_lfmm.sig$position[i] &
                stop+2000 >= ppt_lfmm.sig$position[i],
              sig_SNP_count.ppt := sig_SNP_count.ppt + 1
            ] #adds in 0.5 Kb up or downstream. 
            
            ppt_lfmm.sig$found[i] = "2Kb"
          }else{
            temp2 = gff.dt[
              contig == ppt_lfmm.sig$scaffold[i] &
                start-1000 <=  ppt_lfmm.sig$position[i] &
                stop+1000 >= ppt_lfmm.sig$position[i]
            ] #adds in 0.5Kb up or downstream. 
            if(nrow(temp2) > 0){
              gff.dt[
                contig == ppt_lfmm.sig$scaffold[i] &
                  start-1000 <=  ppt_lfmm.sig$position[i] &
                  stop+1000 >= ppt_lfmm.sig$position[i],
                sig_SNP_count.ppt := sig_SNP_count.ppt + 1
              ] #adds in 0.5 Kb up or downstream. 
              
              ppt_lfmm.sig$found[i] = "1KbMult"
            }else{
            ppt_lfmm.sig$found[i] = "no"
            }
          }
        }
      }  
    }
  }
}

sum(gff.dt$sig_SNP_count.ppt) #295
gff.dt[sig_SNP_count.ppt > 0,]

ppt_lfmm.sig[ppt_lfmm.sig$found == "yes",] |> nrow() #64
ppt_lfmm.sig[ppt_lfmm.sig$found == "500b",] |> nrow() #64
ppt_lfmm.sig[ppt_lfmm.sig$found == "1Kb",] |> nrow() #51
ppt_lfmm.sig[ppt_lfmm.sig$found == "1.5Kb",] |> nrow() #43
ppt_lfmm.sig[ppt_lfmm.sig$found == "2Kb",] |> nrow() #31
ppt_lfmm.sig[ppt_lfmm.sig$found == "1KbMult",] |> nrow() #9
ppt_lfmm.sig |> nrow() #548

#295 hits out of 548 total SNPs with 9 multi hits

write.table(gff.dt, "data/Nf_LFMM_tables/ppt.gene_SNP_hits.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(ppt_lfmm.sig, "data/Nf_LFMM_tables/ppt.SNPs.gene_found.txt", col.names = T, row.names = F, sep = "\t", quote = F)
