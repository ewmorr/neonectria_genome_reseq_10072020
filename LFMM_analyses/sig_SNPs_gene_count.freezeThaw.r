require(dplyr)
require(data.table)

gff = read.csv("data/Nf_SPANDx_all_seqs/maker2_ann/makerFINAL.all.mRNA_ONLY.gff", sep = "\t")
head(gff)

fT_lfmm = read.table("data/Nf_LFMM_tables/freezeThaw_lfmm.txt", header = T)
colnames(fT_lfmm)
fT_lfmm.sig = fT_lfmm %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
nrow(fT_lfmm.sig)

#new cols for snp counts and whether a gene was found for a snp
gff$sig_SNP_count.fT = rep(0, nrow(gff))
fT_lfmm.sig$found = vector(mode = "character", length = nrow(fT_lfmm.sig))
#set data.table
gff.dt = as.data.table(gff)


#The nested ifs are UGLY. should probably write a function
for(i in 1:nrow(fT_lfmm.sig)){
  
  temp = gff.dt[
    contig == fT_lfmm.sig$scaffold[i] &
    start <=  fT_lfmm.sig$position[i] &
    stop >= fT_lfmm.sig$position[i]
  ]
  #check if there is a matching gene. If not go to next
  if(nrow(temp) == 1){
    gff.dt[
      contig == fT_lfmm.sig$scaffold[i] &
        start <=  fT_lfmm.sig$position[i] &
        stop >= fT_lfmm.sig$position[i],
      sig_SNP_count.fT := sig_SNP_count.fT + 1
    ] #strict search only finds 34 SNPs on genes out of 182
    
    fT_lfmm.sig$found[i] = "yes"
    
  }else{
    temp2 = gff.dt[
      contig == fT_lfmm.sig$scaffold[i] &
        start-500 <=  fT_lfmm.sig$position[i] &
        stop+500 >= fT_lfmm.sig$position[i]
    ] #adds in 0.5Kb up or downstream. 
    
    if(nrow(temp2) == 1){
      print(nrow(temp2))
      gff.dt[
        contig == fT_lfmm.sig$scaffold[i] &
          start-500 <=  fT_lfmm.sig$position[i] &
          stop+500 >= fT_lfmm.sig$position[i],
        sig_SNP_count.fT := sig_SNP_count.fT + 1
      ] #adds in 0.5 Kb up or downstream. 
      
      fT_lfmm.sig$found[i] = "500b"
    }else{
      temp2 = gff.dt[
        contig == fT_lfmm.sig$scaffold[i] &
          start-1000 <=  fT_lfmm.sig$position[i] &
          stop+1000 >= fT_lfmm.sig$position[i]
      ] #adds in 0.5Kb up or downstream. 
      if(nrow(temp2) == 1){
        print(nrow(temp2))
        gff.dt[
          contig == fT_lfmm.sig$scaffold[i] &
            start-1000 <=  fT_lfmm.sig$position[i] &
            stop+1000 >= fT_lfmm.sig$position[i],
          sig_SNP_count.fT := sig_SNP_count.fT + 1
        ] #adds in 0.5 Kb up or downstream. 
        
        fT_lfmm.sig$found[i] = "1Kb"
      }else{
        temp2 = gff.dt[
          contig == fT_lfmm.sig$scaffold[i] &
            start-1500 <=  fT_lfmm.sig$position[i] &
            stop+1500 >= fT_lfmm.sig$position[i]
        ] #adds in 0.5Kb up or downstream. 
        if(nrow(temp2) == 1){
          print(nrow(temp2))
          gff.dt[
            contig == fT_lfmm.sig$scaffold[i] &
              start-1500 <=  fT_lfmm.sig$position[i] &
              stop+1500 >= fT_lfmm.sig$position[i],
            sig_SNP_count.fT := sig_SNP_count.fT + 1
          ] #adds in 0.5 Kb up or downstream. 
          
          fT_lfmm.sig$found[i] = "1.5Kb"
        }else{
          temp2 = gff.dt[
            contig == fT_lfmm.sig$scaffold[i] &
              start-2000 <=  fT_lfmm.sig$position[i] &
              stop-2000 >= fT_lfmm.sig$position[i]
          ] #adds in 0.5Kb up or downstream. 
          if(nrow(temp2) == 1){
            print(nrow(temp2))
            gff.dt[
              contig == fT_lfmm.sig$scaffold[i] &
                start-2000 <=  fT_lfmm.sig$position[i] &
                stop+2000 >= fT_lfmm.sig$position[i],
              sig_SNP_count.fT := sig_SNP_count.fT + 1
            ] #adds in 0.5 Kb up or downstream. 
            
            fT_lfmm.sig$found[i] = "2Kb"
          }else{
            temp2 = gff.dt[
              contig == fT_lfmm.sig$scaffold[i] &
                start-1000 <=  fT_lfmm.sig$position[i] &
                stop+1000 >= fT_lfmm.sig$position[i]
            ] #adds in 0.5Kb up or downstream. 
            if(nrow(temp2) > 0){
              gff.dt[
                contig == fT_lfmm.sig$scaffold[i] &
                  start-1000 <=  fT_lfmm.sig$position[i] &
                  stop+1000 >= fT_lfmm.sig$position[i],
                sig_SNP_count.fT := sig_SNP_count.fT + 1
              ] #adds in 0.5 Kb up or downstream. 
              
              fT_lfmm.sig$found[i] = "1KbMult"
            }else{
            fT_lfmm.sig$found[i] = "no"
            }
          }
        }
      }  
    }
  }
}

sum(gff.dt$sig_SNP_count.fT)
gff.dt[sig_SNP_count.fT > 0,]

fT_lfmm.sig[fT_lfmm.sig$found == "yes",] |> nrow()#1
fT_lfmm.sig[fT_lfmm.sig$found == "500b",] |> nrow()#0
fT_lfmm.sig[fT_lfmm.sig$found == "1Kb",] |> nrow()#1
fT_lfmm.sig[fT_lfmm.sig$found == "1.5Kb",] |> nrow()#3
fT_lfmm.sig[fT_lfmm.sig$found == "2Kb",] |> nrow()#2
fT_lfmm.sig[fT_lfmm.sig$found == "1KbMult",] |> nrow()#0
fT_lfmm.sig |> nrow()#11

#7 of 11 found

write.table(gff.dt, "data/Nf_LFMM_tables/freezeThaw.gene_SNP_hits.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(fT_lfmm.sig, "data/Nf_LFMM_tables/freezeThaw.SNPs.gene_found.txt", col.names = T, row.names = F, sep = "\t", quote = F)
