require(dplyr)
require(data.table)

gff = read.csv("data/Nf_SPANDx_all_seqs/maker2_ann/makerFINAL.all.mRNA_ONLY.gff", sep = "\t")
head(gff)

hdd4_lfmm = read.table("data/Nf_LFMM_tables/hdd4_lfmm.txt", header = T)
colnames(hdd4_lfmm)
hdd4_lfmm.sig = hdd4_lfmm %>% filter(FDR.sig == "sig") %>% select(scaffold, position)
nrow(hdd4_lfmm.sig)

#new cols for snp counts and whether a gene was found for a snp
gff$sig_SNP_count.hdd4 = rep(0, nrow(gff))
hdd4_lfmm.sig$found = vector(mode = "character", length = nrow(hdd4_lfmm.sig))
#set data.table
gff.dt = as.data.table(gff)


#The nested ifs are UGLY. should probably write a function
for(i in 1:nrow(hdd4_lfmm.sig)){
  
  temp = gff.dt[
    contig == hdd4_lfmm.sig$scaffold[i] &
    start <=  hdd4_lfmm.sig$position[i] &
    stop >= hdd4_lfmm.sig$position[i]
  ]
  #check if there is a matching gene. If not go to next
  if(nrow(temp) == 1){
    gff.dt[
      contig == hdd4_lfmm.sig$scaffold[i] &
        start <=  hdd4_lfmm.sig$position[i] &
        stop >= hdd4_lfmm.sig$position[i],
      sig_SNP_count.hdd4 := sig_SNP_count.hdd4 + 1
    ] #strict search only finds 34 SNPs on genes out of 182
    
    hdd4_lfmm.sig$found[i] = "yes"
    
  }else{
    temp2 = gff.dt[
      contig == hdd4_lfmm.sig$scaffold[i] &
        start-500 <=  hdd4_lfmm.sig$position[i] &
        stop+500 >= hdd4_lfmm.sig$position[i]
    ] #adds in 0.5Kb up or downstream. 
    
    if(nrow(temp2) == 1){
      print(nrow(temp2))
      gff.dt[
        contig == hdd4_lfmm.sig$scaffold[i] &
          start-500 <=  hdd4_lfmm.sig$position[i] &
          stop+500 >= hdd4_lfmm.sig$position[i],
        sig_SNP_count.hdd4 := sig_SNP_count.hdd4 + 1
      ] #adds in 0.5 Kb up or downstream. 
      
      hdd4_lfmm.sig$found[i] = "500b"
    }else{
      temp2 = gff.dt[
        contig == hdd4_lfmm.sig$scaffold[i] &
          start-1000 <=  hdd4_lfmm.sig$position[i] &
          stop+1000 >= hdd4_lfmm.sig$position[i]
      ] #adds in 0.5Kb up or downstream. 
      if(nrow(temp2) == 1){
        print(nrow(temp2))
        gff.dt[
          contig == hdd4_lfmm.sig$scaffold[i] &
            start-1000 <=  hdd4_lfmm.sig$position[i] &
            stop+1000 >= hdd4_lfmm.sig$position[i],
          sig_SNP_count.hdd4 := sig_SNP_count.hdd4 + 1
        ] #adds in 0.5 Kb up or downstream. 
        
        hdd4_lfmm.sig$found[i] = "1Kb"
      }else{
        temp2 = gff.dt[
          contig == hdd4_lfmm.sig$scaffold[i] &
            start-1500 <=  hdd4_lfmm.sig$position[i] &
            stop+1500 >= hdd4_lfmm.sig$position[i]
        ] #adds in 0.5Kb up or downstream. 
        if(nrow(temp2) == 1){
          print(nrow(temp2))
          gff.dt[
            contig == hdd4_lfmm.sig$scaffold[i] &
              start-1500 <=  hdd4_lfmm.sig$position[i] &
              stop+1500 >= hdd4_lfmm.sig$position[i],
            sig_SNP_count.hdd4 := sig_SNP_count.hdd4 + 1
          ] #adds in 0.5 Kb up or downstream. 
          
          hdd4_lfmm.sig$found[i] = "1.5Kb"
        }else{
          temp2 = gff.dt[
            contig == hdd4_lfmm.sig$scaffold[i] &
              start-2000 <=  hdd4_lfmm.sig$position[i] &
              stop-2000 >= hdd4_lfmm.sig$position[i]
          ] #adds in 0.5Kb up or downstream. 
          if(nrow(temp2) == 1){
            print(nrow(temp2))
            gff.dt[
              contig == hdd4_lfmm.sig$scaffold[i] &
                start-2000 <=  hdd4_lfmm.sig$position[i] &
                stop+2000 >= hdd4_lfmm.sig$position[i],
              sig_SNP_count.hdd4 := sig_SNP_count.hdd4 + 1
            ] #adds in 0.5 Kb up or downstream. 
            
            hdd4_lfmm.sig$found[i] = "2Kb"
          }else{
            temp2 = gff.dt[
              contig == hdd4_lfmm.sig$scaffold[i] &
                start-1000 <=  hdd4_lfmm.sig$position[i] &
                stop+1000 >= hdd4_lfmm.sig$position[i]
            ] #adds in 0.5Kb up or downstream. 
            if(nrow(temp2) > 0){
              gff.dt[
                contig == hdd4_lfmm.sig$scaffold[i] &
                  start-1000 <=  hdd4_lfmm.sig$position[i] &
                  stop+1000 >= hdd4_lfmm.sig$position[i],
                sig_SNP_count.hdd4 := sig_SNP_count.hdd4 + 1
              ] #adds in 0.5 Kb up or downstream. 
              
              hdd4_lfmm.sig$found[i] = "1KbMult"
            }else{
            hdd4_lfmm.sig$found[i] = "no"
            }
          }
        }
      }  
    }
  }
}

sum(gff.dt$sig_SNP_count.hdd4)
gff.dt[sig_SNP_count.hdd4 > 0,]

hdd4_lfmm.sig[hdd4_lfmm.sig$found == "yes",] |> nrow()#34
hdd4_lfmm.sig[hdd4_lfmm.sig$found == "500b",] |> nrow()#16
hdd4_lfmm.sig[hdd4_lfmm.sig$found == "1Kb",] |> nrow()#5
hdd4_lfmm.sig[hdd4_lfmm.sig$found == "1.5Kb",] |> nrow() #11
hdd4_lfmm.sig[hdd4_lfmm.sig$found == "2Kb",] |> nrow() #9
hdd4_lfmm.sig[hdd4_lfmm.sig$found == "1KbMult",] |> nrow() #1
hdd4_lfmm.sig |> nrow() #182

#setting the up/downstream to 1500 instead of 1kb picks up an additional 10 SNPs (60 over 50). Including the one Kb increases from 34 to 50. 

#Allowing SNPs to match multiple genes if within 1Kb up/down stream gives 75 hits over 62 SNPs (one SNP hits three genes at 1Kbp up/down)

#Allowing SNPs to match multiple genes within 1.5Kb gives 91 hits over 73

#We will want to go further with analyzing the SNPs that are not on genes but this is probably fine for now

#chaining from 500 to 1Kb to 1.5 Kb and then allowing multiple hits at 1Kb gives 87 hits, 34, 16, 5, 11, 9 at in gene, 500, 1K,1.5K,2K respetively, 1 mult
#182 total SNPs

write.table(gff.dt, "data/Nf_LFMM_tables/hdd4.gene_SNP_hits.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(hdd4_lfmm.sig, "data/Nf_LFMM_tables/hdd4.SNPs.gene_found.txt", col.names = T, row.names = F, sep = "\t", quote = F)
