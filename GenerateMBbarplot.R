setwd("~/Documents/Projects/Huff_CentromereEvo/Git_Huff_CentEvo")

data <- read.csv("DataSheet_foranalysis.csv")
barplot(data$Bignet_mb, names.arg=data$Species)
?barplot

tempxlab <- barplot(data$Bignet_mb, space=0.4, xaxt='n', ann=FALSE, col="peru", ylab="Mb of Centromere Repeat")
axis(1, cex.axis=0.45, las=2, at=tempxlab, space=0.4, labels=data$Species)
