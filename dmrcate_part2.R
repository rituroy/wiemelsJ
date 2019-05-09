## Wiemel project

cat("\n\n============================ dmrcate_part2.R ===========================\n\n",sep="")

## ---------------------------------

computerFlag="cluster"
computerFlag=""

##############################################

if (computerFlag=="cluster") {
	setwd("/home/royr/project/JoeWiemels")
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
	setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

##############################################

library(DMRcate)

varThis="_caco"; subsetFlag="_dsal"

if (computerFlag=="cluster") {
    datadir=""
} else {
    datadir="epic/results/dmrcate/"
}
load(paste(datadir,"tmp3_dmrcate_methResp",varThis,subsetFlag,"Subset_covSexPlateCelltype_periDsal_bmiq_mVal.RData",sep=""))

varList1=sub("_","",varThis)
groups <- c(Case="magenta", Ctrl="forestgreen")
cols <- groups[as.character(phen[,varList1])]
samps <- c(1:6, 38+(1:6))
samps=1:nrow(phen)
png(paste("plot_dmrcate_",varThis,subsetFlag,sep=""))
DMR.plot(ranges=results.ranges, dmr=1, CpGs=meth, what="M", arraytype = "EPIC",phen.col=cols, genome="hg19", samps=samps)
dev.off()

