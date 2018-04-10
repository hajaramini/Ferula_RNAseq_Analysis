#drow plot for amount of three matabolite based on Mississippi result
metabolites <-read.csv("~/Desktop/Mississippi_Metabolites.tsv",row.names = 1)
colnames(metabolites)

metabolites <- data.frame(file=colnames(metabolites),
                          batch=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\1",colnames(metabolites))),
                          trt=factor(sub("(A|B|C|D|E|F)(F|R|S|L|O)","\\2",colnames(metabolites)))
)
write.csv(metabolites, file="~/Desktop/Missi_Metabolites.csv")


#
metabolites <-read.csv("~/Desktop/Miss_Result.csv", header=T, check.names = F, row.names = 1)

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
ggplot(mt.t.metabolites, aes(gt, value)) +
  geom_boxplot(aes(color = trt)) +
  facet_grid(gt ~ variable) 
