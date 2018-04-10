#getting single and duplicated isoform from BUSCO result

trans <- read_delim("~/Desktop/full_table_BUSCO_output_plant.tsv", comment = "#", delim = "\t", col_names = c("Busco_id","Status","Sequence","Score","Length"))
dup.trans <- trans %>%
  filter(Status == "Duplicated"%>%
           group_by(Busco_id) %>%
           filter(Score == max(Score)) %>%
           ungroup()
         
         no_dup <- trans %>%
           filter(Status != "Duplicated")
         
         single_isoform_trans_complete <- rbind(no_dup, dup.trans) %>%
           arrange(Busco_id)
         
         bad.trans <- trans %>%
           filter(Status == "Duplicated") %>%
           group_by(Busco_id) %>%
           filter(Score != max(Score)) %>%
           ungroup()
         
         bad.trans2 <- trans %>%
           filter(Status == "Fragmented")
         
         bad.trans3 <- rbind(bad.trans,bad.trans2) %>%
           arrange(Busco_id)
         
         lowQ_trans_seqs <- as.data.frame(bad.trans3$Sequence)
         write.table(lowQ_trans_seqs, "Low_Quality_Transcript_ids.txt", row.names = F, quote = F)
         
         write.table(single_isoform_trans_complete, "single_isoform_trans_complete_Busco_ids.txt", row.names = F, quote = F)
         
         bad.trans3
         dim()
         lowQ_trans_seqs
         
##I want to see which contigs id of transrate result match with BUSCO, contigs.csv is the result of transrate 
         contigs <- read_csv("~/Desktop/contigs.csv")
         dim(contigs)
         dim(single_isoform_trans_complete)
         transrate.busco.contigs <- left_join(single_isoform_trans_complete,contigs, by = c("Sequence" = "contig_name"))
         
         write.table(transrate.busco.contigs, "single_isoform_trans_complete, row.names" = F, quote = F)

##I want to see the length distrubution of Trinity.fasta 
hist(tab$V3)
tab <- tab[tab$V3 <=500, ]
hist(tab$V3)
summary(tab)
tab<-read.delim("~/Desktop/Trinity_gene_lengths_Tessa", h=F, sep = "\t")
### Look at primary isoforms distribution
tab2 <-read.delim("~/Desktop/Trinity_isoform_header.txt", h=F, sep = "\t")
head(tab2)
sum(tab2[,2]<=1000)





  