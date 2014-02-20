setwd("~/Documents/Projects/Huff_CentromereEvo/Git_Huff_CentEvo/playdir/")

Orig <- read.csv("Original_assemblynames_forsubset_tdact.txt", header=FALSE)
remain <- read.csv("assemblynamesinDB_knobs_vs_tdacttrf.txt", header=FALSE)
notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
write.csv(notknownrepeats,file = "nonrepeatsknob_tdact.csv")

Orig <- read.csv("Tander_assemblynames_forsubset.txt", header=FALSE)
remain <- read.csv("assemblynamesinDB_knobs_vs_tandertrf.txt", header=FALSE)
notknownrepeats <- as.data.frame(setdiff(Orig$V1, remain$V1))
write.csv(notknownrepeats,file = "Tandernonrepeat.csv")

