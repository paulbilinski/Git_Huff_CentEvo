setwd("~/Documents/Projects/Huff_CentromereEvo/Git_Huff_CentEvo/playdir/")

Orig <- read.csv("anepalRFM_assemblynames_forsubset.txt", header=FALSE)
remain <- read.csv("assemblynamesinDB_NEWknobs_vs_anepaltrf.txt", header=FALSE)
notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
write.csv(notknownrepeats,file = "nonNEWknob_anepal.csv")

###Start analyses for kinetichoreproject

setwd("~/Documents/Projects/Huff_CentromereEvo/Git_Huff_CentEvo/")

data<- read.csv("DataSheet_foranalysis.csv")
plot(data$Bignet_mb ~ data$genome_size, xlab="Genome Size (GB)", ylab="Cent Repeat Content (MB)")
abline(lm(data$Bignet_mb ~ data$genome_size), col="red")
cor.test(data$Bignet_mb, data$genome_size)

plot(data$genome_size, data$Length)

install.packages("Hmisc")
library(Hmisc)
rcorr(as.matrix(data))
hist(data$genome_size,breaks=10, col="lightblue", main="Histogram of Genome Size", xlab="Genome Size (GB)", ylab="")
hist(data$Bignet_mb)





