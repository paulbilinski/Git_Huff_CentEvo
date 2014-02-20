setwd("~/Documents/Projects/Huff_CentromereEvo/")

data <- read.csv("Master_justdata.csv")
data2 <- subset(data, data$FWD_read.>5000)

par(mfrow=c(1,2))
plot(data2$genome_size ~ data2$Maize_CRMUTE_.,xlab="% CRM",ylab="Genome Size")
abline(lm(data2$genome_size ~ data2$Maize_CRMUTE_.),col="red")
plot(data2$genome_size ~ data2$Maize_CRMUTE_mb,xlab="MB CRM",ylab="Genome Size")
abline(lm(data2$genome_size ~ data2$Maize_CRMUTE_mb),col="red")
