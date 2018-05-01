#drow plot for amount of three matabolite based on Mississippi result
metabolites <-read.csv("~/Desktop/Mississippi_Metabolites.tsv",row.names = 1)
colnames(metabolites)

metabolites <- data.frame(file=colnames(metabolites),
                          batch=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\1",colnames(metabolites))),
                          trt=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\2",colnames(metabolites)))
)
write.csv(metabolites, file="~/Desktop/Missi_Metabolites.csv")


#
metabolites <-read.csv("~/Miss_Result.csv", header=T, #check.names = F, row.names = 1)

rownames(metabolites) <- metabolites$Type
rownames(metabolites) 

metabolites.samples <- data.frame(file=colnames(metabolites),
                                  batch=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\1",colnames(metabolites))),
                                  trt=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\2",colnames(metabolites)))
)

#####
metabolites <-read.csv("~/Desktop/Miss_Result.csv",row.names = 1)

head(metabolites)

metabolites.samples <- data.frame(file=colnames(metabolites),
                                  batch=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\1",colnames(metabolites))),
                                  trt=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\2",colnames(metabolites)))
)


######
t.metabolites <- as.data.frame(t(metabolites))
t.metabolites$gt <- gsub("(A|B|C|D|E|F)(F|R|S|L|O)","\\1",rownames(t.metabolites))
t.metabolites$trt <- gsub("(A|B|C|D|E|F)(F|R|S|L|O)","\\2",rownames(t.metabolites))

head(t.metabolites)

#visulize 
library(ggplot2)
mt.t.metabolites <- melt(t.metabolites)
mt.t.metabolites
x <-ggplot(mt.t.metabolites, aes(gt, value)) +
  geom_boxplot(aes(color = trt)) 
  facet_grid(gt ~ variable))
ggsave("~/Miss.metabolite.png", width = 8, height = 8)


####NEW version for Figure
metabolites <-read.csv("~/Desktop/Missi_Metabolites.csv")
rownames(metabolites)
head(metabolites)
metabolites$sample <- gsub("(Flower|Root|Stem|Leaf)(6|2|3)","\\2",metabolites$Sample)
metabolites$tissue <- gsub("(Flower|Root|Stem|Leaf)(6|2|3)","\\1",metabolites$Sample)
metabolites<-metabolites[,-2]
t.metabolites <- as.data.frame(t(metabolites))
mt.t.metabolites <- melt(metabolites)
mt.t.metabolites
umbelliferone <- mt.t.metabolites[c(1:12),]
g <- ggplot(umbelliferone, aes(tissue, value, color=tissue))+ geom_boxplot()

Luteolin <- mt.t.metabolites[c(27:38),]
q <- ggplot(Luteolin, aes(tissue, value, color=tissue)) + geom_boxplot()

Oleo <- mt.t.metabolites[c(13,26,39),]
p <- ggplot(Oleo, mapping=aes(variable, value, fill=variable))+geom_bar(stat="identity") 


library(cowplot)
R<- plot_grid(g,q,p, labels = c("A", "B", "C"), nrow = 2)
ggsave(R,filename = "~/Desktop/metabolite.pdf",width = 12,height = 8)
