---
  title: "KEGG_Analysis"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
#usefull website
#use this script to find the significant EC for pairwise comparision (KEGG enrichment)
#https://github.com/ctb/2017-ucdavis-igg201b/tree/master/lab9
#http://dibsi-rnaseq.readthedocs.io/en/latest/pathway_analysis.html
#https://github.com/dgrapov/TeachingDemos/blob/master/Demos/Pathway%20Analysis/KEGG%20Pathway%20Enrichment.md
#https://www.bioconductor.org/packages/3.7/bioc/vignettes/pathview/inst/doc/pathview.pdf
```{r}
#The result of the blastp in Blast_Arabidopsis folder 
#cut -f1,2,11 Drap_Oases_Plant6_blastp_outfmt6 | sed "s/::/\t/" | sed "s/::/\t/" | cut -f2,4,5 > Final_basltp_withEvalue_clean.txt
#sort -k3,3n -k1,1 -u Final_basltp_withEvalue_clean.txt > test.txt
#mv test.txt Final_basltp_bestEvalue.txt 

#cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.LvsR.csv | sed 's/,/	/g' | sed 's/"//g' | awk '{print $1}' > test1.txt
#grep -F -f test1.txt Final_basltp_bestEvalue.txt.sort  > test2.txt 
#less test2.txt
#Ruijuan command pick the best hit for gene name
#cat test2.txt | sort -k1,1 -k3,3 | sort -u -k1,1 --merge > test3.txt
#chack
#cat test2.txt | sort -k1,1 -k3,3 | sort -u -k1,1 --merge | awk '{print $1}' | sort | uniq | wc -l #1775
#cat test3.txt | awk '{print $2}' > LvsR.ath.name.txt
gene.data.LvsR <- read.csv(file="~/2017-ucdavis-igg201b/lab9/functional-analysis/Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.LvsR.csv", row.names=1) #dim 3779    5


sort -k3,3n -k1,1 -u Final_basltp_withEvalue_clean.txt > test.txt
mv test.txt Final_basltp_bestEvalue.txt have to find why I wrote this line


cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.LvsR.csv | sed 's/,/	/g' | sed 's/"//g' | awk '{print $1}' > test1.txt
grep -F -f test1.txt Final_basltp_bestEvalue.txt.sort  > test2.txt
Ruijuan command pick the best hit for gene name
cat test2.txt | sort -k1,1 -k3,3 | sort -u -k1,1 --merge > test3.txt
#chack
cat test2.txt | sort -k1,1 -k3,3 | sort -u -k1,1 --merge | awk '{print $1}' | sort | uniq | wc -l
1775
cat test3.txt | awk '{print $1,$2}' > test4.txt
join -1 1 -2 1 Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.LvsR.csv2 test4.txt | awk '{print $7,$2,$3,$4,$5,$6}' > gene.ath.LvsR
nano 
add "GeneName"	"logFC"	"logCPM"	"LR"	"PValue"	"FDR"
cat gene.ath.LvsR | sed 's/ / /'g > gene.ath.LvsR.v2
cat gene.ath.LvsR.v2 | sed 's/        / /g' > gene.ath.LvsR.v3


#Subset the data to identify the significantly DE genes as well as those that are increased in the leaf condition:
gene.ath.LvsR <- read.csv("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.LvsR.v3",row.names=NULL, )
#find the duplicate GeneName to use row.names in read.csv
#gene.ath.LvsR[!duplicated(gene.ath.LvsR$GeneName),] %>% dim() #remove duplicate but I need to pick the highest Fc so run following script in linux

remove the colnames line
nano gene.ath.LvsR

cat gene.ath.LvsR | sed "s/\./\t/" | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sed "s/\"//" | sed "s/\t/ /g" >  gene.ath.LvsR.v2
cat gene.ath.LvsR.v2 | awk '{if ($2<0) print $2*(-1); else print $2}' | sort -k1,1 -k2,2n > gene.ath.LvsR.Sorted

#Try to pick the highest FC by considering the up and down (+&-) Last version of edite
cat gene.ath.LvsR | sed "s/\./\t/" | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sed "s/\"//" | sed "s/\t/ /g" >  gene.ath.LvsR.v2
cat gene.ath.LvsR.v2 | sort -k1,1 > gene.ath.LvsR.v2.sortedName
g++ -o pickTop pickTop.cpp
./pickTop gene.ath.LvsR.v2.sortedName > gene.ath.LvsR.v3 #1520
nano gene.ath.LvsR.v3 to add the colnames again "GeneName"	"logFC"	"logCPM"	"LR"	"PValue"	"FDR"
gene.ath.LvsR.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.LvsR.v3", stringsAsFactors = F,row.names = 1, header = T) #now we can use GeneName as a row.names without duplicated GeneName
gene.ath.LvsR.v3%>% dim() #1519  5
#With any sort of RNAseq project, you will at some point (depending on your question), want to turn the lists of fold change and significantly differentially abundant genes into some sort of biologically meaningful information. In this lab we will be relying upon the KEGG database to help get a better idea of the biochemical consequences of the mutation. There are many different metabolic databases/frameworks you can choose from (e.g. GO), but I tend to gravitate towards KEGG for a few reasons: 1) the hierarchy of metabolic pathways is more logical than GO and 2) it can incorporate the co-analysis of genomic/transcriptomic/proteomic data and metabolics data (more on this later).
#R in terminal
#R
#Load the clusterProfiler library 
#library(clusterProfiler)
#library(DOSE)
#Read in the edgeR output
gene.ath.LvsR.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.LvsR.v3", stringsAsFactors = F,row.names = 1, header = T) #after edite the file
save(gene.ath.LvsR.v3, file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.LvsR.v3.RData")
#Subset the data to identify the significantly DE genes as well as those that are increased in the Leaf condition:
#de.genes.LvsR <- subset(gene.ath.LvsR.v3, FDR < 0.05) #do not need this it is already < 0.05
up.genes.LvsR <-subset(gene.ath.LvsR.v3, logFC > 1) #dim 807   5
up.genes.names.LvsR <-row.names(up.genes.LvsR)
#Run enrichKEGG (which calculates a p-value for enrichment using a hypergeometic distribution) to identify the significantly enriched KEGG pathways in the set of genes that are increased in abundance in the Leaf condition:
kegg.up.enrichKEGG.LvsR<-enrichKEGG(up.genes.names.LvsR, organism='ath')
save(kegg.up.enrichKEGG.LvsR, file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/kegg.up.enrichKEGG.LvsR.RData")
head(summary(kegg.up.enrichKEGG.LvsR))

x<- data.frame(kegg.up.enrichKEGG.LvsR) #now we can see the x$Description and other value or x$p.adjust
#x$Description
# [1] "Porphyrin and chlorophyll metabolism"       
#  [2] "Carbon fixation in photosynthetic organisms"
#  [3] "Photosynthesis"                             
#  [4] "Photosynthesis - antenna proteins"          
#  [5] "Carbon metabolism"                          
#  [6] "Glyoxylate and dicarboxylate metabolism"    
#  [7] "Carotenoid biosynthesis"                    
#  [8] "Glycolysis / Gluconeogenesis"               
#  [9] "Pentose phosphate pathway"                  
# [10] "alpha-Linolenic acid metabolism"            
# [11] "Arginine and proline metabolism" 
#ath00900     6 Terpenoid backbone biosynthesis but pvalue=0.18 I prefer to use this becuase the ath is not complete in terpenoid so the # count is not so high

#Visualize the enriched pathways:
#Load the pathview library:
#$library(pathview)
#pathview lets you take advantage of the genomes that are available in the KEGG database. Take a look at all of them here:
data(korg)
head(korg)
#Now, we are going to simulate some metabolite data just to show why we are using KEGG in particular and to show how pathview can handle both data types:
organism <- "Arabidopsis thaliana"
matches <- unlist(sapply(1:ncol(korg), function(i) {agrep(organism, korg[, i])}))
(kegg.code <- korg[matches, 1, drop = F])
cpd.data <- data.frame(FC = sim.mol.data(mol.type = "cpd", nmol = 4000))

#Now, let's plot a pathway that was identified as significantly enriched by the above analysis. We will start off looking at the pathway map for'ath01200' , which was signficantly enriched in the analysis with clusterProfiler:
#Now, grab the log fold change (logFC) data from the gene.data dataframe:
gene.ath.LvsR.v3.fc <-gene.ath.LvsR.v3[1]
save(gene.ath.LvsR.v3.fc,file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.LvsR.v3.fc.RData")

#https://www.rdocumentation.org/packages/pathview/versions/1.12.0/topics/pathview
map<-'ath01200'
pv.out <- pathview(gene.data = gene.ath.LvsR.v3.fc,  cpd.data = cpd.data, gene.idtype = "KEGG",
                   pathway.id = map, species = kegg.code, out.suffix = map,
                   keys.align = "y", kegg.native = T, match.data=T, key.pos = "topright")
plot.name<-paste(map,map,"png",sep=".")

#head(pv.out) to see compound and gene ID mapped

#Take a look at the pathways that are available for arabidopsis.

library(KEGGREST)
pathways<-keggList("pathway", kegg.code)
pathways
write.table(pathways,file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/pathways.ath")

# LvsR for down has no significant KEGG enrichment.
```


```{r}
#FvsR
cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsR.csv | sed 's/,/    /g' | sed 's/"//g' | awk '{print $1}' > GeneName_DEgenes.FvsR
grep -F -f GeneName_DEgenes.FvsR Final_basltp_bestEvalue.txt.sort > DEgenes_blastp.FvsR

Ruijuan command pick the best hit for gene name
cat DEgenes_blastp.FvsR | sort -k1,1 -k3,3 | sort -u -k1,1 --merge > DEgenes_blastp.FvsR.v2
#check
cat DEgenes_blastp.FvsR | sort -k1,1 -k3,3 | sort -u -k1,1 --merge | awk '{print $1}' | sort | uniq | wc -l
#427
cat DEgenes_blastp.FvsR.v2 | awk '{print $1,$2}' > DEgenes_blastp.FvsR.v3
less Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsR.csv2 | sort -k1,1 > Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsR.csv2.sort
join -1 1 -2 1 Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsR.csv2.sort DEgenes_blastp.FvsR.v3 | awk '{print $7,$2,$3,$4,$5,$6}' > gene.ath.FvsR

#remove after dot in GeneName
cat gene.ath.FvsR | sed "s/\./\t/" | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sed "s/\"//" | sed "s/\t/ /g" >  gene.ath.FvsR.v2

#Try to pick the highest FC by considering the up and down (+&-) Last version of edite to solve the " duplicate 'row.names' are not allowed" ERROR for read.csv
cat gene.ath.FvsR.v2 | sort -k1,1 > gene.ath.FvsR.v2.sortedName
g++ -o pickTop pickTop.cpp
./pickTop gene.ath.FvsR.v2.sortedName > gene.ath.FvsR.v3 #398
nano gene.ath.FvsR.v3 to add the colnames again "GeneName" "logFC" "logCPM" "LR" "PValue" "FDR"
gene.ath.FvsR.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsR.v3", stringsAsFactors = F,row.names = 1, header = T) #now we can use GeneName as a row.names without duplicated GeneName
gene.ath.FvsR.v3%>% dim() #397 5
#With any sort of RNAseq project, you will at some point (depending on your question), want to turn the lists of fold change and significantly differentially abundant genes into some sort of biologically meaningful information. In this lab we will be relying upon the KEGG database to help get a better idea of the biochemical consequences of the mutation. There are many different metabolic databases/frameworks you can choose from (e.g. GO), but I tend to gravitate towards KEGG for a few reasons: 1) the hierarchy of metabolic pathways is more logical than GO and 2) it can incorporate the co-analysis of genomic/transcriptomic/proteomic data and metabolics data (more on this later).
#R in terminal
#R
#Load the clusterProfiler library 
#library(clusterProfiler)
#library(DOSE)
#Read in the edgeR output
gene.ath.FvsR.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsR.v3", stringsAsFactors = F,row.names = 1, header = T) #after edite the file
save(gene.ath.FvsR.v3, file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsR.v3.RData")
#Subset the data to identify the significantly DE genes as well as those that are increased in the flower condition:
#de.genes.FvsR <- subset(gene.ath.FvsR.v3, FDR < 0.05) #do not need this it is already < 0.05
up.genes.FvsR <-subset(gene.ath.FvsR.v3, logFC > 1) #dim 337   5
up.genes.names.FvsR <-row.names(up.genes.FvsR)
#Run enrichKEGG (which calculates a p-value for enrichment using a hypergeometic distribution) to identify the significantly enriched KEGG pathways in the set of genes that are increased in abundance in the flower condition:
kegg.up.enrichKEGG.FvsR<-enrichKEGG(up.genes.names.FvsR, organism='ath')
save(kegg.up.enrichKEGG.FvsR, file =  "~/2017-ucdavis-igg201b/lab9/functional-analysis/kegg.up.enrichKEGG.FvsR.RData")
head(summary(kegg.up.enrichKEGG.FvsR))
# $ ID         : chr  "ath00196" "ath00195" "ath00941"
#  $ Description: chr  "Photosynthesis - antenna proteins" "Photosynthesis" "Flavonoid biosynthesis"
#  $ GeneRatio  : chr  "8/105" "9/105" "5/105"
#  $ BgRatio    : chr  "22/4916" "77/4916" "22/4916"
#  $ pvalue     : num  8.27e-09 3.16e-05 7.97e-05
#  $ p.adjust   : num  5.62e-07 1.07e-03 1.81e-03
#  $ qvalue     : num  5.22e-07 9.98e-04 1.68e-03
#  $ geneID     : chr  "AT1G29920/AT1G45474/AT1G61520/AT2G34430/AT3G54890/AT3G61470/AT4G10340/AT5G54270" "AT1G03130/AT1G08380/AT1G31330/AT1G55670/AT1G67740/AT1G79040/AT2G30570/AT2G30790/AT3G50820" "AT3G51240/AT4G22880/AT5G07990/AT5G13930/AT5G42800"
#  $ Count      : int  8 9 5

x<- data.frame(kegg.up.enrichKEGG.FvsR) #now we can see the x$Description and other value or x$p.adjust
#x$Description

#Visualize the enriched pathways:
#Load the pathview library:
#$library(pathview)
#pathview lets you take advantage of the genomes that are available in the KEGG database. Take a look at all of them here:
data(korg)
head(korg)
#Now, we are going to simulate some metabolite data just to show why we are using KEGG in particular and to show how pathview can handle both data types:
organism <- "Arabidopsis thaliana"
matches <- unlist(sapply(1:ncol(korg), function(i) {agrep(organism, korg[, i])}))
(kegg.code <- korg[matches, 1, drop = F])
#Note: "KEGG COMPOUND accession" has only 6599 unique IDs
cpd.data <- data.frame(FC = sim.mol.data(mol.type = "cpd", nmol = 4000))

#Now, let's plot a pathway that was identified as significantly enriched by the above analysis. We will start off looking at the pathway map for'ath01200' , which was signficantly enriched in the analysis with clusterProfiler:
#Now, grab the log fold change (logFC) data from the gene.data dataframe:
gene.ath.FvsR.v3.fc <-gene.ath.FvsR.v3[1]
save(gene.ath.FvsR.v3.fc,file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsR.v3.fc.RData")
map<-'ath00941'
pv.out <- pathview(gene.data = gene.ath.FvsR.v3.fc,  cpd.data = cpd.data, gene.idtype = "KEGG",
                   pathway.id = map, species = kegg.code, out.suffix = map,
                   keys.align = "y", kegg.native = T, match.data=T, key.pos = "topright")
plot.name<-paste(map,map,"png",sep=".")
```

```{r}
#Subset the data to identify the significantly DE genes as well as those that are decreased in the flower condition:
down.genes.FvsR <-subset(gene.ath.FvsR.v3, logFC < -1) #dim 60   5
down.genes.names.FvsR <-row.names(down.genes.FvsR)
#Run enrichKEGG (which calculates a p-value for enrichment using a hypergeometic distribution) to identify the significantly enriched KEGG pathways in the set of genes that are decreased in abundance in the flower condition:
kegg.down.enrichKEGG.FvsR<-enrichKEGG(down.genes.names.FvsR, organism='ath')
save(kegg.down.enrichKEGG.FvsR, file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/kegg.down.enrichKEGG.FvsR.RData")
head(summary(kegg.down.enrichKEGG.FvsR))
#                ID                             Description GeneRatio BgRatio
# ath01212 ath01212                   Fatty acid metabolism      3/16 67/4916
# ath01040 ath01040 Biosynthesis of unsaturated fatty acids      2/16 29/4916
# ath00061 ath00061                 Fatty acid biosynthesis      2/16 41/4916
# ath00592 ath00592         alpha-Linolenic acid metabolism      2/16 42/4916
#               pvalue   p.adjust     qvalue                        geneID Count
# ath01212 0.001193676 0.02387353 0.02010402 AT2G33150/AT2G43710/AT5G46290     3
# ath01040 0.003831197 0.03831197 0.03226271           AT2G33150/AT2G43710     2
# ath00061 0.007563853 0.03963507 0.03337690           AT2G43710/AT5G46290     2
# ath00592 0.007927013 0.03963507 0.03337690           AT1G17420/AT2G33150     2

#Visualize the enriched pathways:
#Load the pathview library:
#$library(pathview)
#pathview lets you take advantage of the genomes that are available in the KEGG database. Take a look at all of them here:
data(korg)
head(korg)
#Now, we are going to simulate some metabolite data just to show why we are using KEGG in particular and to show how pathview can handle both data types:
organism <- "Arabidopsis thaliana"
matches <- unlist(sapply(1:ncol(korg), function(i) {agrep(organism, korg[, i])}))
(kegg.code <- korg[matches, 1, drop = F])
cpd.data <- data.frame(FC = sim.mol.data(mol.type = "cpd", nmol = 4000))

#Now, let's plot a pathway that was identified as significantly enriched by the above analysis. We will start off looking at the pathway map , which was signficantly enriched in the analysis with clusterProfiler:
#Now, grab the log fold change (logFC) data from the gene.data dataframe:
gene.ath.FvsR.v3.fc <-gene.ath.FvsR.v3[1]
save(gene.ath.FvsR.v3.fc,file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsR.v3.fc.RData")
map<-'ath01212'
pv.out <- pathview(gene.data = gene.ath.FvsR.v3.fc,  cpd.data = cpd.data, gene.idtype = "KEGG",
                   pathway.id = map, species = kegg.code, out.suffix = map,
                   keys.align = "y", kegg.native = T, match.data=T, key.pos = "topright")
plot.name<-paste(map,map,"png",sep=".")
```


```{r}
#FvsL
cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsL.csv | sed 's/,/    /g' | sed 's/"//g' | awk '{print $1}' > GeneName_DEgenes.FvsL
grep -F -f GeneName_DEgenes.FvsL Final_basltp_bestEvalue.txt.sort > DEgenes_blastp.FvsL

Ruijuan command pick the best hit for gene name
cat DEgenes_blastp.FvsL | sort -k1,1 -k3,3 | sort -u -k1,1 --merge > DEgenes_blastp.FvsL.v2
#check
cat DEgenes_blastp.FvsL | sort -k1,1 -k3,3 | sort -u -k1,1 --merge | awk '{print $1}' | sort | uniq | wc -l
#556
cat DEgenes_blastp.FvsL.v2 | awk '{print $1,$2}' > DEgenes_blastp.FvsL.v3
cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsL.csv | sed 's/,/    /g' | sed 's/"//g' >  Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsL.csv2
cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsL.csv2 | sort -k1,1 > Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsL.csv2.sort
join -1 1 -2 1 Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsL.csv2.sort DEgenes_blastp.FvsL.v3 | awk '{print $7,$2,$3,$4,$5,$6}' > gene.ath.FvsL

#remove after dot in GeneName
cat gene.ath.FvsL | sed "s/\./\t/" | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sed "s/\"//" | sed "s/\t/ /g" >  gene.ath.FvsL.v2

#Try to pick the highest FC by considering the up and down (+&-) Last version of edite to solve the " duplicate 'row.names' are not allowed" ERROR for read.csv
cat gene.ath.FvsL.v2 | sort -k1,1 > gene.ath.FvsL.v2.sortedName
g++ -o pickTop pickTop.cpp
./pickTop gene.ath.FvsL.v2.sortedName > gene.ath.FvsL.v3 #517
nano gene.ath.FvsR.v3 to add the colnames again "GeneName" "logFC" "logCPM" "LR" "PValue" "FDR"
gene.ath.FvsL.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsL.v3", stringsAsFactors = F,row.names = 1, header = T) #now we can use GeneName as a row.names without duplicated GeneName
gene.ath.FvsL.v3%>% dim() #517 5
#With any sort of RNAseq project, you will at some point (depending on your question), want to turn the lists of fold change and significantly differentially abundant genes into some sort of biologically meaningful information. In this lab we will be relying upon the KEGG database to help get a better idea of the biochemical consequences of the mutation. There are many different metabolic databases/frameworks you can choose from (e.g. GO), but I tend to gravitate towards KEGG for a few reasons: 1) the hierarchy of metabolic pathways is more logical than GO and 2) it can incorporate the co-analysis of genomic/transcriptomic/proteomic data and metabolics data (more on this later).
#R in terminal
#R
#Load the clusterProfiler library 
#library(clusterProfiler)
#library(DOSE)
#Read in the edgeR output
gene.ath.FvsL.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsL.v3", stringsAsFactors = F,row.names = 1, header = T) #after edite the file
save(gene.ath.FvsL.v3, file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsL.v3.RData")
#Subset the data to identify the significantly DE genes as well as those that are increased in the flower condition:
#de.genes.FvsR <- subset(gene.ath.FvsR.v3, FDR < 0.05) #do not need this it is already < 0.05
up.genes.FvsL <-subset(gene.ath.FvsL.v3, logFC > 1) #dim 468  5
up.genes.names.FvsL <-row.names(up.genes.FvsL)
#Run enrichKEGG (which calculates a p-value for enrichment using a hypergeometic distribution) to identify the significantly enriched KEGG pathways in the set of genes that are increased in abundance in the flower condition:
kegg.up.enrichKEGG.FvsL<-enrichKEGG(up.genes.names.FvsL, organism='ath') #NO enriched
```


```{r}
#FvsL decreased in Flower
down.genes.FvsL <-subset(gene.ath.FvsL.v3, logFC < -1)
down.genes.names.FvsL <-row.names(down.genes.FvsL)
kegg.down.enrichKEGG.FvsL<-enrichKEGG(down.genes.names.FvsL, organism='ath') 
#                ID                           Description GeneRatio  BgRatio
# ath04626 ath04626            Plant-pathogen interaction       3/8 170/4916
# ath00562 ath00562         Inositol phosphate metabolism       2/8  73/4916
# ath04070 ath04070 Phosphatidylinositol signaling system       2/8  76/4916
# ath00740 ath00740                 Riboflavin metabolism       1/8  15/4916
#               pvalue   p.adjust      qvalue                        geneID Count
# ath04626 0.002001552 0.01401087 0.006320692 AT4G04720/AT4G38230/AT5G21274     3
# ath00562 0.005748101 0.01450962 0.006545691           AT2G22240/AT4G08170     2
# ath04070 0.006218407 0.01450962 0.006545691           AT4G08170/AT5G21274     2
# ath00740 0.024168017 0.04229403 0.019080014                     AT2G01890     1
map<-'ath04626'
pv.out <- pathview(gene.data = gene.ath.FvsL.v3.fc,  cpd.data = cpd.data, gene.idtype = "KEGG",
                   pathway.id = map, species = kegg.code, out.suffix = map,
                   keys.align = "y", kegg.native = T, match.data=T, key.pos = "topright")
plot.name<-paste(map,map,"png",sep=".") #NOT too much gene mapped
```

```{r}
#FvsS
cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsS.csv | sed 's/,/    /g' | sed 's/"//g' | awk '{print $1}' > GeneName_DEgenes.FvsS
grep -F -f GeneName_DEgenes.FvsS Final_basltp_bestEvalue.txt.sort > DEgenes_blastp.FvsS

#Ruijuan command pick the best hit for gene name
cat DEgenes_blastp.FvsS | sort -k1,1 -k3,3 | sort -u -k1,1 --merge > DEgenes_blastp.FvsS.v2
#check
cat DEgenes_blastp.FvsS | sort -k1,1 -k3,3 | sort -u -k1,1 --merge | awk '{print $1}' | sort | uniq | wc -l
#385
cat DEgenes_blastp.FvsS.v2 | awk '{print $1,$2}' > DEgenes_blastp.FvsS.v3
cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsS.csv | sed 's/,/        /g' | sed 's/"//g' > Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsS.csv2
cat Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsS.csv2 | sort -k1,1 > Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsS.csv2.sort
join -1 1 -2 1 Ferula_RNAseq_drap_oases_plant6_Davis_DEgenes.FvsS.csv2.sort DEgenes_blastp.FvsS.v3 | awk '{print $7,$2,$3,$4,$5,$6}' > gene.ath.FvsS

#remove after dot in GeneName
cat gene.ath.FvsS | sed "s/\./\t/" | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | sed "s/\"//" | sed "s/\t/ /g" >  gene.ath.FvsS.v2

#Try to pick the highest FC by considering the up and down (+&-) Last version of edite to solve the " duplicate 'row.names' are not allowed" ERROR for read.csv
cat gene.ath.FvsS.v2 | sort -k1,1 > gene.ath.FvsS.v2.sortedName
g++ -o pickTop pickTop.cpp
./pickTop gene.ath.FvsS.v2.sortedName > gene.ath.FvsS.v3 #357
nano gene.ath.FvsS.v3 to add the colnames again "GeneName" "logFC" "logCPM" "LR" "PValue" "FDR"
gene.ath.FvsS.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsS.v3", stringsAsFactors = F,row.names = 1, header = T) #now we can use GeneName as a row.names without duplicated GeneName
gene.ath.FvsS.v3%>% dim() #356 5
#library(clusterProfiler)
#library(DOSE)
#Read in the edgeR output
gene.ath.FvsS.v3 <- read.table("~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsS.v3", stringsAsFactors = F,row.names = 1, header = T) #after edite the file
save(gene.ath.FvsS.v3, file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsS.v3.RData")
#Subset the data to identify the significantly DE genes as well as those that are increased in the flower condition:
up.genes.FvsS <-subset(gene.ath.FvsS.v3, logFC > 1) #dim 299   5
up.genes.names.FvsS <-row.names(up.genes.FvsS)
#Run enrichKEGG (which calculates a p-value for enrichment using a hypergeometic distribution) to identify the significantly enriched KEGG pathways in the set of genes that are increased in abundance in the flower condition:
kegg.up.enrichKEGG.FvsS<-enrichKEGG(up.genes.names.FvsS, organism='ath')
save(kegg.up.enrichKEGG.FvsS, file =  "~/2017-ucdavis-igg201b/lab9/functional-analysis/kegg.up.enrichKEGG.FvsS.RData")
summary(kegg.up.enrichKEGG.FvsS)
#               ID                                   Description GeneRatio
# ath00941 ath00941                        Flavonoid biosynthesis      7/97
# ath00040 ath00040      Pentose and glucuronate interconversions      7/97
# ath00909 ath00909 Sesquiterpenoid and triterpenoid biosynthesis      4/97
# ath00073 ath00073          Cutin, suberine and wax biosynthesis      4/97
#          BgRatio       pvalue     p.adjust       qvalue
# ath00941 22/4916 1.254081e-07 8.653160e-06 8.052521e-06
# ath00040 76/4916 6.737662e-04 2.175253e-02 2.024263e-02
# ath00909 23/4916 9.457621e-04 2.175253e-02 2.024263e-02
# ath00073 28/4916 2.027798e-03 3.497952e-02 3.255150e-02
#                                                                         geneID
# ath00941 AT3G51240/AT3G55120/AT4G22880/AT5G07990/AT5G08640/AT5G13930/AT5G42800
# ath00040 AT1G11580/AT3G01270/AT3G05610/AT3G14310/AT3G43270/AT4G13710/AT5G51500
# ath00909                               AT1G78960/AT4G34640/AT4G34650/AT5G23960
# ath00073                               AT1G02205/AT1G57750/AT3G44540/AT4G00360
#          Count
# ath00941     7
# ath00040     7
# ath00909     4
# ath00073     4

organism <- "Arabidopsis thaliana"
matches <- unlist(sapply(1:ncol(korg), function(i) {agrep(organism, korg[, i])}))
(kegg.code <- korg[matches, 1, drop = F])
#Note: "KEGG COMPOUND accession" has only 6599 unique IDs
cpd.data <- data.frame(FC = sim.mol.data(mol.type = "cpd", nmol = 4000))

#Now, grab the log fold change (logFC) data from the gene.data dataframe:
gene.ath.FvsS.v3.fc <-gene.ath.FvsS.v3[1]
save(gene.ath.FvsS.v3.fc,file = "~/2017-ucdavis-igg201b/lab9/functional-analysis/gene.ath.FvsS.v3.fc.RData")
map <- "ath00909"
map <-"ath00073"
pv.out <- pathview(gene.data = gene.ath.FvsS.v3.fc,  cpd.data = cpd.data, gene.idtype = "KEGG",
                   pathway.id = map, species = kegg.code, out.suffix = map,
                   keys.align = "y", kegg.native = T, match.data=T, key.pos = "topright")
plot.name<-paste(map,map,"png",sep=".")
```

```{r}
down.genes.FvsS <-subset(gene.ath.FvsS.v3, logFC < 1) #dim 299   5
down.genes.names.FvsS <-row.names(down.genes.FvsS)
kegg.down.enrichKEGG.FvsS<-enrichKEGG(down.genes.names.FvsS, organism='ath') #NOT enriched
```

