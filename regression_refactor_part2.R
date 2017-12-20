## Goodman test not used. Gives negative variance sometimes

computerFlag="cluster"
computerFlag=""

## ----------------------------------------------
if (computerFlag=="cluster") {
	setwd("/home/royr/project/JoeWiemels")
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
	setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

## ----------------------------------------------

covESFlag=""
covESFlag="_covEpStr"

## ----------------------------------------------
if (computerFlag=="cluster") {
	datadir1=datadir2=""
} else {
	datadir1="results/set1/guthrie/comparison/"
	datadir2="results/set2/guthrie/comparison/"
}

if (F) {
	snp1 <- read.delim(paste("list_to_exclude.txt",sep=""),header=T, sep=" ",quote="",comment.char="",as.is=T,fill=T)
	snp2 <- read.delim(paste("CpGs to exclude_FINAL.txt",sep=""),header=F, sep="\t",quote="",comment.char="",as.is=T,fill=T)
	snp3 <- read.delim(paste("list_to_exclude_Sept_24.txt",sep=""),header=T, sep=" ",quote="",comment.char="",as.is=T,fill=T)
	snp1[,1]=gsub("\"","",snp1[,1])
	write.table(c("cpgId",snp1[,1]),file="list_to_exclude_Sept_24.txt", sep="\t", col.names=F, row.names=F, quote=F, append=F)
}

## ----------------------------------------------
if (computerFlag=="cluster") {
	ann <- read.delim(paste("data/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
	snpVec <- read.table(paste("data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
} else {
	ann <- read.delim(paste("docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
	snpVec <- read.table(paste("docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	dmr <- read.table(paste("docs/SeungTae/leukemia.DMRs/leukemia.DMRs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
}
ann[which(ann[,"CHR"]=="X"),"CHR"]="23"
ann[which(ann[,"CHR"]=="Y"),"CHR"]="24"
ann[,"CHR"]=as.integer(ann[,"CHR"])
ann <- ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
for (k in 1:ncol(ann)) if (class(ann[,k])=="factor") ann[,k]=as.character(ann[,k])

snpVec=snpVec[,1]
ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
ann$keep=0
ann$keep[which(ann$snp==0 & ann$CHR%in%1:22)]=1

if (covESFlag=="_covEpStr") {
    dirEpistructure="docs/SemiraGonsethNussle/epistructure/"
    prId=read.table(paste(dirEpistructure,"epistructure_reference_sites.txt",sep=""), sep="\t", h=F, quote="", comment.char="",as.is=T,fill=T)
    prId=prId[,1]
    ann$keep[which(ann$keep==1 & ann$IlmnID%in%prId)]=0
}


ann$geneSym=sapply(toupper(ann$UCSC_RefGene_Name),function(x) {
	strsplit(x,";")[[1]][1]
},USE.NAMES=F)
ann$geneSym[is.na(ann$geneSym)]=""

save.image("tmp2_1.RData")

## -----------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------

## ----------------------------------------------
## Subsets

library(qvalue)

transformFlag=""
transformFlag="_mVal"

datadir="results/comparison/pcb/"

colIdPV="qv"; pvName="qv"
pThres=0.1
pThres=0.05
pThres=0.2

colIdPV="pv"; pvName="pv"
pThres=0.1
pThres=0.05
pThres=0.001
pThres=0.01
pThres=0.2

fileList=dir(datadir)
f1=read.table(paste(datadir,fileList[1],sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl=c()
tbl2=NULL
for (fId in 1:length(fileList)) {
    cat("\n\n================= ",fileList[fId],"\n",sep="")
    f2=read.table(paste(datadir,fileList[fId],sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    nm=sub("coef_","",names(f2)[grep("coef",names(f2))])
    names(f2)=c("cpgId","coef","se","pv")
    if (any(f1$cpgId!=f2$cpgId,na.rm=T) | any(is.na(f1$cpgId)!=is.na(f2$cpgId))) print(fId)
    i=which(!is.na(f2$pv))
    f2$qv[i]=qvalue(f1$pv[i])$qvalues
    i=which(f2[,colIdPV]<pThres)
    print(length(i))
    if (length(i)!=0) {
        tbl2=rbind(tbl2,cbind(comparison=rep(nm,length(i)),f2[i,]))
        tbl=c(tbl,f2$cpgId[i])
    }
}
i=match(tbl,ann$IlmnID)
x=as.character(ann$CHR[i])
x[which(x=="23")]="X"
x[which(x=="24")]="Y"
x=paste("chr",x,":",ann$MAPINFO[i],"-",ann$MAPINFO[i],sep="")
sum(!duplicated(x))
write.table(x[!duplicated(x)],file=paste("ann_pcb_",colIdPV,pThres,"_forUCSCLiftover_hg19.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)
write.table(tbl2,file=paste("stat_pcb_",colIdPV,pThres,".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)

i=which(tbl2$comparison=="logged_PCB_105_SRS")
i=which(tbl2$comparison=="logged_PCB_aroclor1260")
for (compFlag in c("logged_PCB_105_SRS","logged_PCB_138_SRS","logged_PCB_aroclor1260")) {
    i=which(tbl2$comparison==compFlag)
    write.table(x[i],file=paste("ann_pcb_",compFlag,"_pv",pThres,"_forUCSCLiftover_hg19.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)
    write.table(tbl2[i,],file=paste("stat_pcb_",compFlag,"_pv",pThres,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

if (F) {
    varFlag="meth"
    stat_cc_s1=read.table(paste(datadir,"stat_",varFlag,"Resp_logged_PCB_105_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_cc_s1=read.table(paste(datadir,"stat_",varFlag,"Resp_logged_PCB_105_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
}

## --------------

transformFlag=""
transformFlag="_mVal"

datadir="results/comparison/pcb/"

varFlag=""
stat_pcb105=read.table(paste(datadir,"stat_methResp_logged_PCB_105_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb118=read.table(paste(datadir,"stat_methResp_logged_PCB_118_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb138=read.table(paste(datadir,"stat_methResp_logged_PCB_138_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb153=read.table(paste(datadir,"stat_methResp_logged_PCB_153_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb170=read.table(paste(datadir,"stat_methResp_logged_PCB_170_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb180=read.table(paste(datadir,"stat_methResp_logged_PCB_180_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb1260=read.table(paste(datadir,"stat_methResp_logged_PCB_aroclor1260_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

varFlag=""
stat_pcb105=read.table(paste(datadir,"stat_methResp_logged_PCB_105_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb118=read.table(paste(datadir,"stat_methResp_logged_PCB_118_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb138=read.table(paste(datadir,"stat_methResp_logged_PCB_138_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb153=read.table(paste(datadir,"stat_methResp_logged_PCB_153_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb170=read.table(paste(datadir,"stat_methResp_logged_PCB_170_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb180=read.table(paste(datadir,"stat_methResp_logged_PCB_180_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb1260=read.table(paste(datadir,"stat_methResp_logged_PCB_aroclor1260_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

varFlag=""
stat_pcb105=read.table(paste(datadir,"stat_methResp_logged_PCB_105_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb118=read.table(paste(datadir,"stat_methResp_logged_PCB_118_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb138=read.table(paste(datadir,"stat_methResp_logged_PCB_138_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb153=read.table(paste(datadir,"stat_methResp_logged_PCB_153_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb170=read.table(paste(datadir,"stat_methResp_logged_PCB_170_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb180=read.table(paste(datadir,"stat_methResp_logged_PCB_180_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_pcb1260=read.table(paste(datadir,"stat_methResp_logged_PCB_aroclor1260_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## ----------------------------------------------
## Subsets

transformFlag=""
transformFlag="_mVal"

datadir="results/comparison/caco/"

varFlag="caco"
stat_cc_s1=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsem_hispSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s2=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSem_hispSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s3=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsemNoSexChr_hispSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s4=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSemNoSexChr_hispSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s5=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsemNoSnpNoSexChr_hispSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s6=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSemNoSnpNoSexChr_hispSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

varFlag="caco"
stat_cc_s1=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsem_noHispWtSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s2=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSem_noHispWtSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s3=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsemNoSexChr_noHispWtSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s4=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSemNoSexChr_noHispWtSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s5=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsemNoSnpNoSexChr_noHispWtSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s6=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSemNoSnpNoSexChr_noHispWtSubset_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

varFlag="caco"
stat_cc_s1=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsem_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s2=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSem_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s3=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsemNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s4=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSemNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s5=read.table(paste(datadir,"stat_",varFlag,"Resp_methXsemNoSnpNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_s6=read.table(paste(datadir,"stat_",varFlag,"Resp_methXlogSemNoSnpNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## ----------------------------------------------
## Subsets

transformFlag=""
transformFlag="_mVal"

datadir="results/comparison/caco/"

varFlag="caco"
stat_cc_2=read.table(paste(datadir,"stat_methResp_",varFlag,"_covSexGestage_covEpStr_allGuthSet2_bmiq_mVal.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_fs_2=read.table(paste(datadir,"stat_methResp_",varFlag,"_femaleSubset_covGestage_covEpStr_allGuthSet2_bmiq_mVal.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_ms_2=read.table(paste(datadir,"stat_methResp_",varFlag,"_maleSubset_covGestage_covEpStr_allGuthSet2_bmiq_mVal.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_he_2=read.table(paste(datadir,"stat_methResp_",varFlag,"_hispSubset_covSexGestage_covPrinComp1234_allGuthSet2_bmiq_mVal.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_cc_we_2=read.table(paste(datadir,"stat_methResp_",varFlag,"_noHispWtSubset_covSexGestage_covPrinComp1234_allGuthSet2_bmiq_mVal.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

stat_cc_3=read.table(paste(datadir,"stat_methResp_",varFlag,"_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq_mVal.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

png("tmp.png")
par(mfrow=c(2,2))
lim=range(c(stat_cc_2$coef_caco,stat_cc_3$coef_caco),na.rm=T)
plot(stat_cc_2$coef_caco, stat_cc_3$coef_caco,xlim=lim,ylim=lim); abline(c(0,1))
lim=range(c(stat_cc_2$pv_caco,stat_cc_3$pv_caco),na.rm=T)
plot(stat_cc_2$pv_caco, stat_cc_3$pv_caco,xlim=lim,ylim=lim); abline(c(0,1))
x1=order(stat_cc_2$pv_caco); x2=order(stat_cc_3$pv_caco)
lim=range(c(x1,x2),na.rm=T)
plot(x1, x2,xlim=lim,ylim=lim,xlab="stat_cc_2: order(pv)",ylab="stat_cc_3: order(pv)"); abline(c(0,1))
dev.off()

pThres=0.05
table(stat_cc_2$pv_caco<pThres,stat_cc_3$pv_caco<pThres)
table(stat_cc_2$qv_caco<pThres,stat_cc_3$qv_caco<pThres)


tbl=read.table(paste("/Users/royr/Downloads/ann_methResp_caco_covSexGestage_covEpStr_allGuthSet2_bmiq_mVal.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)


## ----------------------------------------------
## Subsets


transformFlag=""
modelFlag=""

datadir="results/comparison/aml/"
stat_aml_1_0=read.table(paste(datadir,"stat_methResp_caco_covSexEthnGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_aml_1_00=read.table(paste(datadir,"stat_methResp_caco_covSexGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_aml_1_1=read.table(paste(datadir,"stat_methResp_caco_hispSubset_covSexGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_aml_1_2=read.table(paste(datadir,"stat_methResp_caco_noHispWtSubset_covSexGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_aml_2_0=read.table(paste(datadir,"stat_cacoResp_meth_covSexEthnGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_aml_2_00=read.table(paste(datadir,"stat_cacoResp_meth_covSexGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_aml_2_1=read.table(paste(datadir,"stat_cacoResp_meth_hispSubset_covSexGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_aml_2_2=read.table(paste(datadir,"stat_cacoResp_meth_noHispWtSubset_covSexGestage_covEpStr_aml_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
## ----------------------------------------------
## Subsets

datadir="results/comparison/birthWt/refactor/"


modelFlag="pobw->meth->caco"
stat_bw_1=read.table(paste(datadir,"stat_mediation_cacoResp_pobw_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_bw_2=read.table(paste(datadir,"stat_mediation_methResp_pobw_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_bw_3=read.table(paste(datadir,"stat_mediation_cacoResp_meth_pobw_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_bw_c=read.table(paste(datadir,"stat_mediation_covariance_cacoResp_meth_pobw_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

modelFlag="meth->pobw->caco"
stat_bw_1=read.table(paste(datadir,"stat_mediation_cacoResp_meth_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_bw_2=read.table(paste(datadir,"stat_mediation_pobwResp_meth_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_bw_3=read.table(paste(datadir,"stat_mediation_cacoResp_meth_pobw_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat_bw_c=read.table(paste(datadir,"stat_mediation_covariance_cacoResp_meth_pobw_covSexHisp_covPrinComp1234_allGuthSet2_bmiq.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)


load("tmp2.RData")

transformFlag=""
transformFlag="_mVal"

if (F) {
    stat_bw_i_2=read.table(paste(datadir,"stat_refactor_cacoXpobw_allGuthSet2_bmiq_winsorTailP0.05",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_1_0=read.table(paste(datadir,"stat_refactor_pobw_ctrlSubset_allGuthSet1_bmiq_winsorTailP0.05",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_2_0=read.table(paste(datadir,"stat_refactor_pobw_ctrlSubset_allGuthSet2_bmiq_winsorTailP0.05",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
}

#if (F) {
    varFlag="birthWt"
    varFlag="pobw"
    stat_bw_1_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_1=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_1_2=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"Bi_covPrinComp1234_allGuthSet1_noNonRndChip_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_2_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_2=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    stat_bw_3_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_3=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_4_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_4=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    if (transformFlag=="_mVal") {
        stat_bw_i_2_m_co=stat_bw_i_2
        stat_bw_1_0_m_co=stat_bw_1_0
        stat_bw_2_0_m_co=stat_bw_2_0
        
        stat_bw_i_4_m_co=stat_bw_i_4
        stat_bw_3_0_m_co=stat_bw_3_0
        stat_bw_4_0_m_co=stat_bw_4_0
    } else {
        stat_bw_i_2_b_co=stat_bw_i_2
        stat_bw_1_0_b_co=stat_bw_1_0
        stat_bw_2_0_b_co=stat_bw_2_0

        stat_bw_i_4_b_co=stat_bw_i_4
        stat_bw_3_0_b_co=stat_bw_3_0
        stat_bw_4_0_b_co=stat_bw_4_0
    }
#}

#load("tmp_all.RData")
#load("tmp2.RData")
#if (F) {
    varFlag="birthWtBi"; cutoff=3500
    varFlag="pobwBi"; cutoff=1
    stat_bw_1_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_1=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_1_2=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet1_noNonRndChip_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_2_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_2=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    stat_bw_3_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_3=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covEpStr_allGuthSet1_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_4_0=read.table(paste(datadir,"stat_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    stat_bw_i_4=read.table(paste(datadir,"stat_cacoResp_methX",varFlag,"_covEpStr_allGuthSet2_bmiq",transformFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    if (transformFlag=="_mVal") {
        stat_bw_i_2_m_ca=stat_bw_i_2
        stat_bw_1_0_m_ca=stat_bw_1_0
        stat_bw_2_0_m_ca=stat_bw_2_0
        
        stat_bw_i_4_m_ca=stat_bw_i_4
        stat_bw_3_0_m_ca=stat_bw_3_0
        stat_bw_4_0_m_ca=stat_bw_4_0
    } else {
        stat_bw_i_2_b_ca=stat_bw_i_2
        stat_bw_1_0_b_ca=stat_bw_1_0
        stat_bw_2_0_b_ca=stat_bw_2_0
        
        stat_bw_i_4_b_ca=stat_bw_i_4
        stat_bw_3_0_b_ca=stat_bw_3_0
        stat_bw_4_0_b_ca=stat_bw_4_0
    }
#}


## ------------
stat_bw_i_2=stat_bw_i_2_b_ca
stat_bw_1_0=stat_bw_1_0_b_ca
stat_bw_2_0=stat_bw_2_0_b_ca

## ------------
stat_bw_i_2=stat_bw_i_2_m_co
stat_bw_1_0=stat_bw_1_0_m_co
stat_bw_2_0=stat_bw_2_0_m_co

stat_bw_i_2=stat_bw_i_2_m_ca
stat_bw_1_0=stat_bw_1_0_m_ca
stat_bw_2_0=stat_bw_2_0_m_ca

## ------------
stat_bw_i_2=stat_bw_i_2_b_co
stat_bw_1_0=stat_bw_1_0_b_co
stat_bw_2_0=stat_bw_2_0_b_co

stat_bw_i_4=stat_bw_i_4_b_co
stat_bw_3_0=stat_bw_3_0_b_co
stat_bw_4_0=stat_bw_4_0_b_co

## ------------
stat_bw_i_2=stat_bw_i_2_b_ca
stat_bw_1_0=stat_bw_1_0_b_ca
stat_bw_2_0=stat_bw_2_0_b_ca

stat_bw_i_4=stat_bw_i_4_b_ca
stat_bw_3_0=stat_bw_3_0_b_ca
stat_bw_4_0=stat_bw_4_0_b_ca
## ------------


#save.image("tmp2_1.RData")

stat=stat_bw_c
for (k in 1:ncol(stat)) {
    if (any(is.na(stat[,k]))) print(colnames(stat)[k])
}

compInfo=data.frame(comp1=c("m_co","m_co","m_ca","b_co"),comp2=c("m_ca","b_co","b_ca","b_ca"),stringsAsFactors=F)
png("plotCoefPair_%1d.png")
par(mfcol=c(2,2))
#par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(5, 5, 4, 2) + 0.1)
for (k in 1:nrow(compInfo)) {
    colId=match(c("comp1","comp2"),names(compInfo))
    lab=rep("",length(colId))
    for (p in 1:length(colId)) {
        switch(compInfo[k,colId[p]],
        "m_co"={
            stat=stat_bw_i_2_m_co
            ttl="Methylation-pobw interaction coefficient\nContinuous pobw, m-value based"
        },
        "m_ca"={
            stat=stat_bw_i_2_m_ca
            ttl="Methylation-pobw interaction coefficient\nCategorical pobw, m-value based"
        },
        "b_co"={
            stat=stat_bw_i_2_b_co
            ttl="Methylation-pobw interaction coefficient\nContinuous pobw, beta-value based"
        },
        "b_ca"={
            stat=stat_bw_i_2_b_ca
            ttl="Methylation-pobw interaction coefficient\nCategorical pobw, beta-value based"
        }
        )
        names(stat)=sub("_meth.pobwBi","",names(stat),fixed=T)
        names(stat)=sub("_meth.pobw","",names(stat),fixed=T)
        lab[p]=ttl
        if (p==1) {
            stat1=stat
        } else {
            stat2=stat
        }
        
    }
    i=1:10
    i=1:nrow(stat1)
    plot(stat1$coef[i],stat2$coef[i],xlab=lab[1],ylab=lab[2])
}
dev.off()

png("tmp_pv.png")
plot(stat_bw_i_2_co$pv_meth.pobw,stat_bw_i_2_ca$pv_meth.pobwBi)
dev.off()


stat2=stat_bw_i_2
iA2=match(stat2$cpgId,ann$IlmnID)
i=grep("ARID5B",ann$geneSym[iA2])
summary(stat_bw_i_2$pv_meth.pobwBi[i])
summary(stat_bw_2_0$pv_pobwBi[i])
stat2=stat_bw_1_0
iA2=match(stat2$cpgId,ann$IlmnID)
i=grep("ARID5B",ann$geneSym[iA2])
summary(stat_bw_1_0$pv_pobwBi[i])
i1=match(stat_bw_1_0$cpgId[i],stat_bw_2_0$cpgId)
i2=
plot(stat_bw_1_0$pv_pobwBi[i],

## ----------------------------------------------
## Check if there are any significant interaction effects for a short list of genes based on other results
## Implemented for
##       Set2: caco ~ methylation * pobw, where pobw is a categorical variable, pobw <= 1 vs. pobw > 1. Median pobw is ~ 1

library(qvalue)

id=c()

if (transformFlag=="_mVal") {
    pThres=0.01
} else {
    pThres=0.005
}
fName=paste("_cacoResp_methXpobwBi_qv0.05_for_methResp_pobwBi_ctrlSubset_pv",pThres,transformFlag,sep="")
i1=which(stat_bw_2_0$pv_pobw<pThres)
stat2=stat_bw_i_2[i1,]
nm=names(stat2)
names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
stat2$qv=NA
iA2=match(stat2$cpgId,ann$IlmnID)
i=which(ann$snp[iA2]==0 & ann$CHR[iA2]%in%1:22)
stat2$qv[i]=getAdjustedPvalue(stat2$pv[i],method="qv",strict=T)
i=which(stat2$qv<0.05)
id=stat2$cpgId[i]
x=sapply(ann$UCSC_RefGene_Name[iA2][i][order(stat2$pv[i])],function(x) {strsplit(x,";")[[1]][1]},USE.NAMES=F)
print(unique(x))
names(stat2)=nm
write.table(paste("Set2: caco ~ methylation * pobw, where pobw is a categorical variable, pobw <= 1 vs. pobw > 1. Median pobw is ~ 1\nLoci with q-value < 0.05 for methylation-pobw interaction term\nConsidered only ",length(i1)," loci where p-value < ",pThres," for methylation ~ pobw (categorical) in set2 controls\n",ifelse(transformFlag=="_mVal","M-value","Beta-value")," based",sep=""), file=paste("stat",fName,".txt",sep=""), append=F,col.names=F,row.names=F, sep="\t",quote=F)
write.table(cbind(ann[iA2,][i,],stat2[i,]), file=paste("stat",fName,".txt",sep=""), append=T,col.names=T,row.names=F, sep="\t",quote=F)

i1=match(stat2$cpgId[i],stat_bw_i_1$cpgId)
table(is.na(i1))
#plot(stat_bw_i_1$coef_meth.pobwBi[i1],stat2$coef[i])
cor(stat_bw_i_1$coef_meth.pobwBi[i1],stat2$coef[i],use="complete.obs",method="pearson")
cor(stat_bw_i_1$coef_meth[i1],stat2$coef_meth[i],use="complete.obs",method="pearson")
cor(stat_bw_i_1$coef_pobwBi[i1],stat2$coef_pobwBi[i],use="complete.obs",method="pearson")
cor(stat_bw_i_1_2$coef_meth.pobwBi[i1],stat2$coef[i],use="complete.obs",method="pearson")
cor(stat_bw_i_1_2$coef_meth[i1],stat2$coef_meth[i],use="complete.obs",method="pearson")
cor(stat_bw_i_1_2$coef_pobwBi[i1],stat2$coef_pobwBi[i],use="complete.obs",method="pearson")
cor(stat_bw_i_1$coef_meth.pobwBi[i1],stat_bw_i_1_2$coef_meth.pobwBi[i1],use="complete.obs",method="pearson")
cor(stat_bw_i_1$coef_meth[i1],stat_bw_i_1_2$coef_meth[i1],use="complete.obs",method="pearson")
cor(stat_bw_i_1$coef_pobwBi[i1],stat_bw_i_1_2$coef_pobwBi[i1],use="complete.obs",method="pearson")

i1=match(stat2$cpgId[i],stat_bw_i_1$cpgId)
table(is.na(i1))

cutoff=1
phen$pobwBi=as.integer(phen$pobw>cutoff)
for (genesetFlag in c("IRS2","FANC")) {
    switch(genesetFlag,
        "IRS2"={
            geneDesc="Insulin receptor substrate 2"
            id2=id[grep("IGF1R|PLCG1|UBTF",ann$UCSC_RefGene_Name[iA2][i])]
        },
        "FANC"={
         geneDesc="Fanconi anemia"
         id2=id[grep("FANC",ann$UCSC_RefGene_Name[iA2][i])]
        }
    )
    cbind(geneSym=ann$UCSC_RefGene_Name[match(id2,ann$IlmnID)],stat_bw_i_2[match(id2,stat_bw_i_2$cpgId),])
    cbind(geneSym=ann$UCSC_RefGene_Name[match(id2,ann$IlmnID)],stat_bw_1_0[match(id2,stat_bw_1_0$cpgId),])
    cbind(geneSym=ann$UCSC_RefGene_Name[match(id2,ann$IlmnID)],stat_bw_2_0[match(id2,stat_bw_2_0$cpgId),])
    i2=match(id2,rownames(meth))
    png(paste("interactionPlot",fName,"_protein",genesetFlag,".png",sep=""))
    par(mfcol=c(2,2))
    for (ii in i2) {
        interaction.plot(x.factor=phen$pobwBi,trace.factor=phen$caco,response=meth[ii,])
    }
    par(mfcol=c(1,1))
    title(sub=geneDesc,cex=2)
    dev.off()
    grp=as.character(phen$pobwBi)
    grp[which(phen$pobwBi==0)]="low pobw"
    grp[which(phen$pobwBi==1)]="high pobw"
    j=which(phen$caco==0)
    grp[j]=paste(grp[j],"ctrl",sep="\n")
    j=which(phen$caco==1)
    grp[j]=paste(grp[j],"case",sep="\n")
    grp=paste(phen$pobwBi,phen$caco,sep="/")
    grpUniq=sort(unique(grp))
    ttl=paste(rep(c("low pobw","high pobw"),each=2),c("ctrl","case"),sep="\n")
    png(paste("boxplot",fName,"_protein",genesetFlag,".png",sep=""))
    par(mfcol=c(2,2))
    lim=quantile(c(meth[i2,]),probs=c(.1,.9),na.rm=T)
    for (ii in 1:length(id2)) {
        lim=quantile(c(meth[i2[ii],]),probs=c(.1,.9),na.rm=T)
        boxplot(meth[i2[ii],]~grp,names=ttl,ylim=lim,main=paste(id2[ii],": ",ann$geneSym[match(id2[ii],ann$IlmnID)],sep=""),ylab="Methylation",las=2)
        for (gId in 1:length(grpUniq)) {
            x=mean(meth[i2[ii],which(grp==grpUniq[gId])],na.rm=T)
            lines(gId+c(-.5,+.5),rep(x,2),col="red")
        }
    }
    plot(0:2,0:2,type="n",axes=F,xlab="",ylab="")
    title(sub=geneDesc,cex=2,lty="solid")
    legend(1,1,"mean",lty="solid",col="red")
    par(mfcol=c(1,1))
    dev.off()
}

for (cohortFlag in c("set1","set2","set1set2")) {
    #for (cohortFlag in c("set2")) {
    switch(cohortFlag,
    "set1"={
        stat2=stat_bw_1_0
    },
    "set2"={
        stat2=stat_bw_2_0
    },
    "set1set2"={
        stat1=stat_bw_1_0
        colnames(stat1)=sub("_","Set1_",colnames(stat1))
        stat2=stat_bw_2_0
        i=match(stat1$cpgId,stat2$cpgId); i1=which(!is.na(i)); i2=i[i1]
        stat2=cbind(stat1[i1,],stat2[i2,which(!colnames(stat2)%in%colnames(stat1))])
    }
    )
    names(stat2)=sapply(names(stat2),function(x) strsplit(x,"_")[[1]][1],USE.NAMES=F)
    colIdPV="pv"; pThres=0.05
    stat2$qv=NA
    i=which(!is.na(stat2[i,colIdPV]))
    iA2=match(stat2$cpgId,ann$IlmnID)
    i=which(ann$snp[iA2]==0 & ann$CHR[iA2]%in%1:22)
    stat2$qv[i]=getAdjustedPvalue(stat2$pv[i],method="qv",strict=T)
    if (cohortFlag=="set1set2") {
        stat2$qvSet1=NA
        stat2$qvSet1[i]=getAdjustedPvalue(stat2$pvSet1[i],method="qv",strict=T)
    }
    stat1=stat2
    i=match(stat_bw_i_2$cpgId,stat1$cpgId); i1=which(!is.na(i)); i2=i[i1]
    colIdPV="qv"; pThres=0.2
    colIdPV="qv"; pThres=0.1
    colIdPV="qv"; pThres=0.05
    colIdPV="pv"; pThres=0.2
    colIdPV="pv"; pThres=0.1
    colIdPV="pv"; pThres=0.05
    colIdPV="pv"; pThres=0.001
    colIdPV="pv"; pThres=0.0005
    colIdPV="pv"; pThres=0.0001
    thInfo=data.frame(colIdPV=c(rep("qv",3),rep("pv",8)),pThres=c(0.2,0.1,0.05,0.2,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001),stringsAsFactors=F)
    for (tId in 1:nrow(thInfo)) {
        colIdPV=thInfo$colIdPV[tId]
        pThres=thInfo$pThres[tId]
        cThres=0.1
        cThres=0
        cat("\n\n================== ",cohortFlag," ",colIdPV,"<",pThres,", abs(coef)>=",cThres," =============\n",sep="")
        if (cohortFlag=="set1set2") {
            i=which(stat1[i2,colIdPV]<pThres & abs(stat1$coef[i2])>=cThres & stat1[i2,paste(colIdPV,"Set1",sep="")]<pThres & abs(stat1$coefSet1[i2])>=cThres)
        } else {
            i=which(stat1[i2,colIdPV]<pThres & abs(stat1$coef[i2])>=cThres)
        }
        if (length(i)!=0) {
            stat2=stat_bw_i_2[i1,][i,]
            #names(stat2)=sub("_interaction","",names(stat2))
            names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
            print("summary(stat2$pv)")
            print(summary(stat2$pv))
            stat2$qv=NA
            iA2=match(stat2$cpgId,ann$IlmnID)
            i=which(ann$snp[iA2]==0 & ann$CHR[iA2]%in%1:22)
            stat2$qv[i]=getAdjustedPvalue(stat2$pv[i],method="qv",strict=T)
            iA2=match(stat2$cpgId,ann$IlmnID)
            print("summary(stat2$qv)")
            print(summary(stat2$qv))
            print("table(stat2$pv<0.05)")
            print(table(stat2$pv<0.05))
            print("table(stat2$qv<0.05)")
            print(table(stat2$qv<0.05))
            i=which(stat2$qv<0.05)
            if (length(i)!=0) {
                x=sapply(ann$UCSC_RefGene_Name[iA2][i][order(stat2$pv[i])],function(x) {strsplit(x,";")[[1]][1]},USE.NAMES=F)
                print(unique(x))
                print("table(stat2$cpgId[i]%in%id)")
                print(table(stat2$cpgId[i]%in%id))
            }
            i=which(stat2$qv<0.2)
            if (length(i)!=0) {
                #print("qv<0.2")
                #print(cbind(ann[iA2,c("IlmnID","UCSC_RefGene_Name")][i,],coef=round(stat2$coef[i],2),signif(stat2[i,c("pv","qv")],2)))
            }
            x=p.adjust(stat2$pv,method="BH")
            cat("\ntable(p.adjust(stat2$pv,method='BH')<0.05)\n")
            print(table(x<0.05))
            i=which(x<0.2)
            if (length(i)!=0) {
                #cat("\nfdr<0.2\n")
                #print(cbind(ann[iA2,c("IlmnID","UCSC_RefGene_Name")][i,],stat2[i,c("coef","pv","qv")]))
            }
        }
    }
}

lim=c(0,1)
i=i2[which(stat1[i2,colIdPV]<pThres & abs(stat1$coef[i2])>=cThres)]
plot(stat1$pv[i],stat2$pv,xlim=lim,ylim=lim)

#library(multtest)

## -------------------
library(bumphunter)
library(methyAnalysis)

datadir="tmp/"

if (T) {
fName="bumphunter_methResp_pobw_ctrlSubset_covPrinComp1234_allGuthSet2_bmiq_mVal_250maxGap_1000perms_2cores.RData"
fName="bumphunter_methResp_pobw_ctrlSubset_covPrinComp1234_allGuthSet1_bmiq_mVal_250maxGap_1000perms_2cores.RData"

fName="bumphunter_methResp_pobwBictrl_covPrinComp1234_allGuthSet1_bmiq_mVal_300maxGap_1000perms_0cores.RData"
fName="bumphunter_methResp_cacoXpobwBi_covPrinComp1234_allGuthSet2_bmiq_mVal_400maxGap_1000perms_0cores.RData"
fName="bumphunter_methResp_pobwBictrl_covPrinComp1234_allGuthSet2_bmiq_mVal_400maxGap_1000perms_0cores.RData"

cat("\n\n============================ ",fName," ===========================\n\n",sep="")

load(file=paste(datadir,fName,sep=""))
tab=dmrs$tab
GR_tab=as(tab[1:1000,],"GRanges")                            # If you get an error here reduce the number of to  tab[1:100]
tab.ann = annotateDMRInfo(GR_tab,'TxDb.Hsapiens.UCSC.hg19.knownGene')
dmrInfo=as.data.frame(tab.ann$sigDMRInfo)
dmrInfo$chr=as.integer(sub("chr","",dmrInfo$seqnames))
rownames(dmrInfo)=NULL

if (length(grep("allGuthSet1",fName))==1) {
    fName1="allGuthSet1"
} else {
    fName1="allGuthSet2"
}

if (length(grep("X",fName))==1) {
    stat2=stat_bw_i_2
    i2=1
    i2=which(dmrInfo$GeneSymbol=="ZNF148")
} else {
    if (length(grep("allGuthSet1",fName))==1) {
        stat2=stat_bw_1_0
        i2=1:3
        i2=which(dmrInfo$GeneSymbol=="ZNF148")
    } else {
        stat2=stat_bw_2_0
        i2=1
        i2=which(dmrInfo$GeneSymbol=="ZNF148")
    }
}
i21=i2[1]
iA2=match(stat2$cpgId,ann$IlmnID)
#i1=dmrInfo$indexStart[i2]:dmrInfo$indexEnd[i2]
#i1=which(ann$CHR[iA2]==dmrInfo$chr[i2] & ann$MAPINFO[iA2]>=dmrInfo$start[i2] &  ann$MAPINFO[iA2]<=dmrInfo$end[i2])
i11=which(ann$CHR[iA2]==dmrInfo$chr[i21] & ann$MAPINFO[iA2]>=dmrInfo$start[i21] &  ann$MAPINFO[iA2]<=dmrInfo$end[i21])
i1=c()
for (i in 1:length(i2)) {
    i1=c(i1,which(ann$CHR[iA2]==dmrInfo$chr[i2[i]] & ann$MAPINFO[iA2]>=dmrInfo$start[i2[i]] &  ann$MAPINFO[iA2]<=dmrInfo$end[i2[i]]))
}
i1=unique(i1)
x=strsplit(fName,"_")[[1]]
cat(paste(fName1,"\n\n"))
cat(paste("From Bumphunter: maxGap ",sub("maxGap","",x[grep("maxGap",x)]),", nResamples ",sub("perms","",x[grep("perms",x)]),"\n",sep=""))
x=rep("",length(i2)); x[which(i2==i21)]="***"
tbl=cbind(dmrInfo[i2,c("GeneSymbol","seqnames","start","end","value","p.value","fwer")],x)
names(tbl)[match(c("x"),names(tbl))]=c(" ")
print(tbl)
cat("\nFrom locus-level analysis:\n")
tbl=stat2[i1,-(grep("cpgId|se_",names(stat2)))]
k=grep("coef_",names(tbl)); tbl[,k]=round(tbl[,k],2)
k=grep("pv_",names(tbl)); tbl[,k]=signif(tbl[,k],2)
k=grep("qv_",names(tbl)); if (length(k)!=0) tbl[,k]=signif(tbl[,k],2)
x=rep("",length(i1)); x[which(i1%in%i11)]="***"
tbl=cbind(ann[iA2,c("IlmnID","geneSym","CHR","MAPINFO","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")][i1,],tbl,x)
names(tbl)[match(c("UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island","x"),names(tbl))]=c("geneGroup","RelToCpG"," ")
rownames(tbl)=NULL
print(tbl)
#ann[iA2,][i1,]
#stat2[i1,]
}

i1=which(ann$geneSym[iA2]=="ZNF148")
tbl=stat2[i1,-(grep("cpgId|se_",names(stat2)))]
tbl=cbind(ann[iA2,c("IlmnID","geneSym","CHR","MAPINFO","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")][i1,],tbl)
names(tbl)[match(c("UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island","x"),names(tbl))]=c("geneGroup","RelToCpG"," ")
rownames(tbl)=NULL
print(tbl)

"

============================ bumphunter_methResp_pobw_ctrlSubset_covPrinComp1234_allGuthSet2_bmiq_mVal_250maxGap_1000perms_2cores.RData ===========================

allGuthSet2

From Bumphunter: maxGap 250, nResamples 1000
GeneSymbol seqnames     start       end     value      p.value  fwer
1       ZNF148     chr3 125093764 125094316 1.8734191 1.408238e-05 0.043 ***
13      ZNF148     chr3 125076233 125076480 1.2062644 1.004761e-03 0.958
70      ZNF148     chr3 125090636 125090636 2.0490720 9.968688e-03 1.000
125     ZNF148     chr3 125077434 125077434 1.7728251 1.855698e-02 1.000
144     ZNF148     chr3 125093330 125093330 1.6961477 2.228553e-02 1.000
833     ZNF148     chr3 125085322 125085322 0.9680929 2.204738e-01 1.000

From locus-level analysis:
IlmnID geneSym CHR   MAPINFO     geneGroup RelToCpG coef_pobw pv_pobw
1 cg17563338  ZNF148   3 125093799         5'UTR   Island     -0.12 0.45000 ***
2 cg15705999  ZNF148   3 125093863         5'UTR   Island     -0.14 0.29000 ***
3 cg24591090  ZNF148   3 125094085 1stExon;5'UTR   Island     -0.26 0.07900 ***
4 cg26994334  ZNF148   3 125094183 1stExon;5'UTR   Island     -0.03 0.85000 ***
5 cg24928100  ZNF148   3 125094316        TSS200   Island      0.57 0.00300 ***
6 cg12508343  ZNF148   3 125076233         5'UTR   Island      0.08 0.74000
7 cg01289769  ZNF148   3 125077434         5'UTR  S_Shore     -0.03 0.89000
8 cg15688767  ZNF148   3 125093330         5'UTR  N_Shore     -0.22 0.22000
9 cg10687936  ZNF148   3 125085322         5'UTR              -0.80 0.00072


============================ bumphunter_methResp_pobw_ctrlSubset_covPrinComp1234_allGuthSet1_bmiq_mVal_250maxGap_1000perms_2cores.RData ===========================

allGuthSet1

From Bumphunter: maxGap 250, nResamples 1000
GeneSymbol seqnames     start       end      value      p.value  fwer
3       ZNF148     chr3 125093764 125094322 -1.1145081 5.432449e-05 0.167 ***
41      ZNF148     chr3 125076233 125076480 -0.8513187 2.700936e-03 1.000
240     ZNF148     chr3 125085322 125085322 -1.2987742 4.403472e-02 1.000
254     ZNF148     chr3 125093330 125093330 -1.2619692 4.963372e-02 1.000

From locus-level analysis:
IlmnID geneSym CHR   MAPINFO     geneGroup RelToCpG coef_pobw pv_pobw
1 cg17563338  ZNF148   3 125093799         5'UTR   Island      0.05 7.3e-01 ***
2 cg15705999  ZNF148   3 125093863         5'UTR   Island     -0.04 7.7e-01 ***
3 cg24591090  ZNF148   3 125094085 1stExon;5'UTR   Island      0.11 4.2e-01 ***
4 cg26994334  ZNF148   3 125094183 1stExon;5'UTR   Island      0.04 7.4e-01 ***
5 cg24928100  ZNF148   3 125094316        TSS200   Island     -0.05 7.9e-01 ***
6 cg23506350  ZNF148   3 125094322        TSS200   Island     -0.26 3.7e-01 ***
7 cg12508343  ZNF148   3 125076233         5'UTR   Island      0.04 8.2e-01
8 cg10687936  ZNF148   3 125085322         5'UTR              -0.67 9.5e-05
9 cg15688767  ZNF148   3 125093330         5'UTR  N_Shore     -0.27 2.5e-02





===================================================================================

============================ bumphunter_methResp_pobwBictrl_covPrinComp1234_allGuthSet2_bmiq_mVal_400maxGap_1000perms_0cores.RData ===========================

allGuthSet2

From Bumphunter: maxGap 400, nResamples 1000
GeneSymbol seqnames     start       end     value      p.value  fwer
1       ZNF148     chr3 125093764 125094316 0.6684035 3.619054e-06 0.011 ***
6       ZNF148     chr3 125076067 125076480 0.4100537 3.928319e-04 0.722
30      ZNF148     chr3 125090636 125090636 0.7146316 3.063365e-03 0.998
52      ZNF148     chr3 125077434 125077434 0.6409514 5.302902e-03 1.000
60      ZNF148     chr3 125093330 125093330 0.6176307 6.318211e-03 1.000
288     ZNF148     chr3 125085322 125085322 0.3748261 5.401768e-02 1.000

From locus-level analysis:
IlmnID geneSym CHR   MAPINFO     geneGroup RelToCpG coef_pobwBi pv_pobwBi qv_pobwBi
1  cg17563338  ZNF148   3 125093799         5'UTR   Island       -0.03      0.49         1 ***
2  cg15705999  ZNF148   3 125093863         5'UTR   Island       -0.05      0.20         1 ***
3  cg24591090  ZNF148   3 125094085 1stExon;5'UTR   Island       -0.05      0.18         1 ***
4  cg26994334  ZNF148   3 125094183 1stExon;5'UTR   Island        0.04      0.42         1 ***
5  cg24928100  ZNF148   3 125094316        TSS200   Island        0.08      0.16         1 ***
6  cg16704739  ZNF148   3 125076067         5'UTR   Island       -0.03      0.32         1
7  cg12508343  ZNF148   3 125076233         5'UTR   Island        0.08      0.23         1
8  cg01289769  ZNF148   3 125077434         5'UTR  S_Shore        0.02      0.72         1
9  cg15688767  ZNF148   3 125093330         5'UTR  N_Shore       -0.03      0.54         1
10 cg10687936  ZNF148   3 125085322         5'UTR                -0.14      0.03         1


============================ bumphunter_methResp_pobwBictrl_covPrinComp1234_allGuthSet1_bmiq_mVal_300maxGap_1000perms_0cores.RData ===========================

allGuthSet1

From Bumphunter: maxGap 300, nResamples 1000
GeneSymbol seqnames     start       end      value    p.value fwer
120     ZNF148     chr3 125093863 125094085 -0.2425693 0.01576521    1 ***
962     ZNF148     chr3 125085322 125085322 -0.2622196 0.31379861    1

From locus-level analysis:
IlmnID geneSym CHR   MAPINFO     geneGroup RelToCpG coef_pobwBi pv_pobwBi qv_pobwBi
1 cg15705999  ZNF148   3 125093863         5'UTR   Island       -0.03    0.4100      0.98 ***
2 cg24591090  ZNF148   3 125094085 1stExon;5'UTR   Island        0.06    0.1500      0.98 ***
3 cg10687936  ZNF148   3 125085322         5'UTR                -0.18    0.0012      0.91


============================ bumphunter_methResp_cacoXpobwBi_covPrinComp1234_allGuthSet2_bmiq_mVal_400maxGap_1000perms_0cores.RData ===========================

allGuthSet2

From Bumphunter: maxGap 400, nResamples 1000
GeneSymbol seqnames     start       end      value      p.value  fwer
1       ZNF148     chr3 125093764 125094316 -0.6215160 3.723558e-05 0.108 ***
16      ZNF148     chr3 125076233 125076480 -0.4417863 1.114309e-03 0.961
131     ZNF148     chr3 125090636 125090636 -0.6014524 1.330344e-02 1.000
217     ZNF148     chr3 125077434 125077434 -0.5136083 2.647036e-02 1.000
234     ZNF148     chr3 125093330 125093330 -0.5052814 2.838868e-02 1.000

From locus-level analysis:
IlmnID geneSym CHR   MAPINFO     geneGroup RelToCpG coef_meth pv_meth coef_pobwBi pv_pobwBi
1 cg17563338  ZNF148   3 125093799         5'UTR   Island     -0.70   0.100        1.12     0.530
2 cg15705999  ZNF148   3 125093863         5'UTR   Island      0.19   0.710       -0.06     0.910
3 cg24591090  ZNF148   3 125094085 1stExon;5'UTR   Island     -0.98   0.049        3.10     0.270
4 cg26994334  ZNF148   3 125094183 1stExon;5'UTR   Island     -0.10   0.820       -2.75     0.260
5 cg24928100  ZNF148   3 125094316        TSS200   Island      0.41   0.280       -4.79     0.100
6 cg12508343  ZNF148   3 125076233         5'UTR   Island      0.75   0.018       -3.11     0.094
7 cg01289769  ZNF148   3 125077434         5'UTR  S_Shore      0.45   0.260        0.68     0.640
8 cg15688767  ZNF148   3 125093330         5'UTR  N_Shore      0.24   0.550        0.28     0.920
coef_meth.pobwBi pv_meth.pobwBi qv_meth.pobwBi
1             0.44           0.46           0.96 ***
2             0.20           0.79           0.97 ***
3             0.77           0.24           0.95 ***
4            -0.61           0.29           0.95 ***
5            -0.89           0.11           0.95 ***
6            -0.67           0.11           0.95
7            -0.31           0.54           0.96
8             0.09           0.87           0.97


============================ methResp_cacoXpobwBi_covPrinComp1234
allGuthSet2

IlmnID geneSym CHR   MAPINFO UCSC_RefGene_Group Relation_to_UCSC_CpG_Island   coef_meth    pv_meth
1  cg13483431  ZNF148   3 124949647              3'UTR                              0.06583663 0.66936078
2  cg00710736  ZNF148   3 125000990               Body                             -0.26623876 0.15643027
3  cg21401200  ZNF148   3 125022192               Body                             -0.75114110 0.04569460
4  cg26313511  ZNF148   3 125053815              5'UTR                             -0.12668696 0.80089271
5  cg16704739  ZNF148   3 125076067              5'UTR                      Island  0.53901158 0.32871577
6  cg12508343  ZNF148   3 125076233              5'UTR                      Island  0.75391135 0.01786379
7  cg01289769  ZNF148   3 125077434              5'UTR                     S_Shore  0.45387651 0.25552188
8  cg10687936  ZNF148   3 125085322              5'UTR                             -0.18797876 0.47585931
9  cg15688767  ZNF148   3 125093330              5'UTR                     N_Shore  0.23568175 0.55235160
10 cg17563338  ZNF148   3 125093799              5'UTR                      Island -0.69971080 0.10394398
11 cg15705999  ZNF148   3 125093863              5'UTR                      Island  0.19447424 0.70941897
12 cg24591090  ZNF148   3 125094085      1stExon;5'UTR                      Island -0.97620154 0.04867693
13 cg26994334  ZNF148   3 125094183      1stExon;5'UTR                      Island -0.09888253 0.82140309
14 cg24928100  ZNF148   3 125094316             TSS200                      Island  0.41359973 0.28106553
15 cg23506350  ZNF148   3 125094322             TSS200                      Island -0.16265362 0.63127784
16 cg15094636  ZNF148   3 125094397             TSS200                     S_Shore  0.45932616 0.33801285
17 cg21664828  ZNF148   3 125094611            TSS1500                     S_Shore -0.84281303 0.07982023
18 cg18107105  ZNF148   3 125094680            TSS1500                     S_Shore  0.02102521 0.93946385
coef_pobwBi  pv_pobwBi coef_meth.pobwBi pv_meth.pobwBi qv_meth.pobwBi
1  -0.73342518 0.28829182      0.150481435     0.41557932      0.9574668
2  -0.17917032 0.72462497     -0.003023948     0.98779789      0.9759100
3  -3.52445308 0.02480313      1.061348408     0.03220221      0.9517068
4  -0.47200285 0.13010056     -0.472684741     0.25868343      0.9524306
5   1.37873320 0.60553845      0.431231561     0.55651048      0.9631672
6  -3.11118167 0.09387973     -0.668043990     0.11425598      0.9517068
7   0.67723545 0.63968132     -0.307781130     0.54249750      0.9622882
8   0.51355559 0.59430836     -0.289775610     0.41935501      0.9576622
9   0.27775814 0.92330426      0.091849171     0.87077560      0.9721882
10  1.12452979 0.52928746      0.439783532     0.45936962      0.9592739
11 -0.05503444 0.91423684      0.204011429     0.78827558      0.9703705
12  3.10113121 0.26775696      0.771661157     0.23609147      0.9519707
13 -2.75126948 0.26146557     -0.612885284     0.29482156      0.9534980
14 -4.78875080 0.10054105     -0.888969942     0.11393409      0.9517068
15  0.97957115 0.69156173      0.197627160     0.63244857      0.9657728
16 -1.20637301 0.72029683     -0.202895562     0.76352119      0.9697980
17  2.42533435 0.39446717      0.557259438     0.35399585      0.9550252
18 -0.32843530 0.87368967     -0.024817147     0.94864702      0.9751662

"

load("tmp2_1.RData")
## Run following section to get clin1 & clin2
## misc.R: Compare set1 & set2 clinical information & beta values"

datadir="docs/snp_arid5b/"
snpH=read.table(paste(datadir,"summary_birthweight_hisp_int.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
snpW=read.table(paste(datadir,"summary_birthweight_white_int.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(snpH)=names(snpW)=c("snpId","chr","pos","pv","logpv")

gene="ZNF148"
gene="ARID5B"

addFlag="snp"

typeFlag="_allcpg"
typeFlag=""

tbl=read.table("docs/probeId.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
iA2=match(tbl$cpgId,ann$IlmnID)
#stat2=stat_bw_i_2
#iA2=match(stat2$cpgId,ann$IlmnID)
#x=ann$CHR+10^9*ann$MAPINFO
#summary(diff(x))
#iA2=order(x)
#i=grep(gene,ann$geneSym[iA2]); i1=range(i)
i=grep(gene,ann$geneSym[iA2])
#ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]<=(min(ann$MAPINFO[iA2][i])-500000)); ii=which.max(ann$MAPINFO[iA2][ii]); i=c(i,ii)
#ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]>=(max(ann$MAPINFO[iA2][i])+500000)); ii=which.min(ann$MAPINFO[iA2][ii]); i=c(i,ii)
ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]<=(min(ann$MAPINFO[iA2][i])-500000)); ii=iA2[ii][which.max(ann$MAPINFO[iA2][ii])]; i=c(i,ii)
ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]>=(max(ann$MAPINFO[iA2][i])+500000)); ii=iA2[ii][which.min(ann$MAPINFO[iA2][ii])]; i=c(i,ii)
i=sort(unique(i))
i1=range(i)

if (F) {
    datadir="docs/all/set1/"; fName="beta_bmiq_allGuthSet1"
    tmp=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
    beta=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,skip=i1[1],nrow=diff(i1)+1)
    colnames(beta)=colnames(tmp)
    table(beta[,1]%in%ann[iA2,1][i])
    table(ann[iA2,1][i]%in%beta[,1])
    rownames(beta)=beta$probeId
    beta=as.matrix(beta[-1])
    iA2=match(rownames(beta),ann$IlmnID)
    beta=beta[order(ann$CHR[iA2]+10^9*ann$MAPINFO[iA2]),]
    beta1=beta
    clin1=clin1[match(colnames(beta1),clin1$id),]
}

datadir="docs/all/set2/"; fName="beta_bmiq_allGuthSet2"
tmp=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
beta=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,skip=i1[1],nrow=diff(i1)+1)
colnames(beta)=colnames(tmp)
table(beta[,1]%in%ann[iA2,1][i])
table(ann[iA2,1][i]%in%beta[,1])
rownames(beta)=beta$probeId
beta=as.matrix(beta[-1])
iA2=match(rownames(beta),ann$IlmnID)
beta=beta[order(ann$CHR[iA2]+10^9*ann$MAPINFO[iA2]),]
stat2=stat_bw_1_0
i=match(rownames(beta),stat2$cpgId); i1=which(!is.na(i)); i2=i[i1]; i12=which(is.na(i))
stat22=stat2[1:length(i12),]
for (k in 1:ncol(stat22)) stat22[,k]=NA
stat22$cpgId=rownames(beta)[i12]
stat2=rbind(stat2,stat22)
stat2=stat2[match(rownames(beta),stat2$cpgId),]
stat10=stat2
stat2=stat_bw_2_0
i=match(rownames(beta),stat2$cpgId); i1=which(!is.na(i)); i2=i[i1]; i12=which(is.na(i))
stat22=stat2[1:length(i12),]
for (k in 1:ncol(stat22)) stat22[,k]=NA
stat22$cpgId=rownames(beta)[i12]
stat2=rbind(stat2,stat22)
stat2=stat2[match(rownames(beta),stat2$cpgId),]
stat20=stat2
stat2=stat_bw_i_2
i=match(rownames(beta),stat2$cpgId); i1=which(!is.na(i)); i2=i[i1]; i12=which(is.na(i))
stat22=stat2[1:length(i12),]
for (k in 1:ncol(stat22)) stat22[,k]=NA
stat22$cpgId=rownames(beta)[i12]
stat2=rbind(stat2,stat22)
stat2=stat2[match(rownames(beta),stat2$cpgId),]
stat2i=stat2
beta2=beta
clin2=clin2[match(colnames(beta2),clin2$id),]

## Scatter plot
levelInfo=data.frame(caco=c(0,0,1,1),pobwBi=c(0,1,0,1),stringsAsFactors=F)
colList=c("skyblue","skyblue4","orange","red")
colList=rep("grey",4)
ltyList=c("solid","dashed","dotted","dotdash")
cutoff=1
pThres=0.05
for (compFlag in c("","_zoom")) {
    fName=paste(compFlag,"_",gene,sep="")
    iA2=match(rownames(beta2),ann$IlmnID)
    if (compFlag=="_zoom") {
        header=paste("Zoomed plot: ",gene,sep="")
        ii=which(ann$geneSym[iA2]==gene & stat2i$pv_meth.pobwBi<pThres)
        i=(min(ii)-1):(max(ii)+1)
        if (F) {
            i=which(ann$geneSym[iA2]==gene & stat2i$pv_meth.pobwBi<pThres)
            ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]<=(min(ann$MAPINFO[iA2][i])-500000)); ii=which.max(ann$MAPINFO[iA2][ii]); i=c(i,ii)
            ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]>=(max(ann$MAPINFO[iA2][i])+500000)); ii=which.min(ann$MAPINFO[iA2][ii]); i=c(i,ii)
            i=sort(unique(i))
            i1=range(i)
            i=i1[1]:i1[2]
        }
    } else {
        header=gene
        i=which(ann$geneSym[iA2]==gene)
    }
    dat=beta2[i,]; iA2=iA2[i]; statThis10=stat10[i,]; statThis20=stat20[i,]; statThis2i=stat2i[i,]
    phen=clin2
    phen$pobwBi=as.integer(phen$pobw>cutoff)
    corMat=cor(t(dat),use="complete.obs",method="spearman")
    xThis=ann$MAPINFO[iA2]/10^3
    xlim=range(xThis,na.rm=T)
    if (compFlag=="_zoom") {
        ylim=c(0,0.3)
    } else {
        ylim=range(c(dat),na.rm=T)
    }
    #png(paste("plot",fName,".png",sep=""))
    png(paste("plot",fName,".png",sep=""),width=8*240,height=4*240)
    #par(mfrow=c(2,1))
    plot(xThis,dat[,1],xlim=xlim,ylim=ylim,main=header,xlab="KB",ylab="Set2: Beta-value",las=2)
    for (grpId in 1:nrow(levelInfo)) {
        j=which(phen$caco==levelInfo$caco[grpId] & phen$pobwBi==levelInfo$pobwBi[grpId])
        x=rep(NA,nrow(dat))
        for (i in 1:nrow(dat)) {
            points(rep(xThis[i],length(j)),dat[i,j],col=colList[grpId])
            #x[i]=median(dat[i,j],na.rm=T)
            x[i]=mean(dat[i,j],na.rm=T)
            if (i!=1) {
                lines(xThis[(i-1):i],x[(i-1):i],lty=ltyList[grpId])
            }
        }
    }
    legend(xThis[length(xThis)*3/4],0.29,legend=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw"),lty=ltyList)
    if (compFlag=="_zoom") {
        legend(xThis[length(xThis)*3/4],0.2,legend=paste(c("set2: caco~meth:pobw"),", pv<",pThres,sep=""),fill=c("red"))
    } else {
        legend(xThis[length(xThis)*3/4],0.1,legend=paste(c("set1 ctrl: meth~pobw","set2 ctrl: meth~pobw","set2: caco~meth:pobw"),", pv<",pThres,sep=""),fill=c("skyblue","skyblue4","red"))
    }
    axis(side=3,at=xThis[which(statThis10$pv_pobwBi<pThres)],col.ticks="skyblue",lwd.ticks=3,las=2)
    axis(side=1,at=xThis[which(statThis20$pv_pobwBi<pThres)],col.ticks="skyblue4",lwd.ticks=3,las=2)
    axis(side=1,at=xThis[which(statThis2i$pv_meth.pobwBi<pThres)],col.ticks="red",las=2)
    if ("snp"%in%addFlag) {
        axis(side=3,at=xThis[iList$snpH],col.ticks="palegreen",lwd.ticks=3,las=2)
        axis(side=3,at=xThis[iList$snpW],col.ticks="seagreen",lwd.ticks=3,las=2)
    }
    #plot(range(ann$MAPINFO[i1[i]]),range(c(corMat)),type="n")
    #plot(1:length(i),1:length(i),type="n")
    #rect(1,corMat[1,2]-.5,2,corMat[1,2]+.5,col="red")
    dev.off()
}

n=nrow(dat)*(nrow(dat)-1)/2
plot(1:n,rep(1,n),ylim=c(-1,1),type="n")
x=c()
for (grpId in 1:nrow(levelInfo)) {
    j=which(phen$caco==levelInfo$caco[grpId] & phen$pobwBi==levelInfo$pobwBi[grpId])
    corMat=cor(t(dat[,j]),use="complete.obs",method="spearman")
    print(summary(abs(corMat[lower.tri(corMat)])))
    lines(1:n,corMat[lower.tri(corMat)],lty=ltyList[grpId])
    x=cbind(x,corMat[lower.tri(corMat)])
}
j=apply(x,1,function(x) {any(abs(x)>=.4)})
sdThis=apply(x,1,sd,na.rm=T)
table(j)
summary(sdThis)
summary(sdThis[j])


## Boxplot
levelInfo=data.frame(caco=c(0,0,1,1),pobwBi=c(0,1,0,1),stringsAsFactors=F)
colList=rep("grey",4)
colList=c("skyblue","skyblue4","orange","red")
ltyList=c("solid","dashed","dotted","dotdash")
cutoff=1
pThres=0.05
pThres2=10^-5
pThres2=0.05
for (compFlag in c("","_zoomSignifInteraction","_zoomSignifPobw")) {
    if (gene=="ZNF148") pos=c(125093764,125094316) else pos=NULL
    fName=paste(compFlag,"_",gene,sep="")
    iA2=match(rownames(beta2),ann$IlmnID)
    if (compFlag=="_zoomSignifInteraction") {
        header=paste("Region of significant interaction: ",gene,sep="")
        if (typeFlag=="_allcpg") {
            ii=which(ann$geneSym[iA2]==gene & stat2i$pv_meth.pobwBi<pThres)
            i=(min(ii)-1):(max(ii)+1)
        } else {
            i=which(ann$geneSym[iA2]==gene & stat2i$pv_meth.pobwBi<pThres)
            if (length(i)==1) i=(min(i)-1):(max(i)+1)
        }
        if (F) {
            i=which(ann$geneSym[iA2]==gene & stat2i$pv_meth.pobwBi<pThres)
            ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]<=(min(ann$MAPINFO[iA2][i])-500000)); ii=which.max(ann$MAPINFO[iA2][ii]); i=c(i,ii)
            ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]>=(max(ann$MAPINFO[iA2][i])+500000)); ii=which.min(ann$MAPINFO[iA2][ii]); i=c(i,ii)
            i=sort(unique(i))
            i1=range(i)
            i=i1[1]:i1[2]
        }
    } else if (compFlag=="_zoomSignifPobw") {
        header=paste("Region of significant meth~pobw: ",gene,sep="")
        if (typeFlag=="_allcpg") {
            ii=which(ann$geneSym[iA2]==gene & (stat10$pv_pobwBi<pThres | stat20$pv_pobwBi<pThres))
            i=(min(ii)-1):(max(ii)+1)
        } else {
            i=which(ann$geneSym[iA2]==gene & (stat10$pv_pobwBi<pThres | stat20$pv_pobwBi<pThres))
            if (length(i)==1) i=(min(i)-1):(max(i)+1)
        }
        if (F) {
            i=which(ann$geneSym[iA2]==gene & stat2i$pv_meth.pobwBi<pThres)
            ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]<=(min(ann$MAPINFO[iA2][i])-500000)); ii=which.max(ann$MAPINFO[iA2][ii]); i=c(i,ii)
            ii=which(ann$CHR[iA2]==ann$CHR[iA2[i[1]]] & ann$MAPINFO[iA2]>=(max(ann$MAPINFO[iA2][i])+500000)); ii=which.min(ann$MAPINFO[iA2][ii]); i=c(i,ii)
            i=sort(unique(i))
            i1=range(i)
            i=i1[1]:i1[2]
        }
    } else {
        header=gene
        if (typeFlag=="_allcpg") {
            ii=which(ann$geneSym[iA2]==gene)
            i=(min(ii)-1):(max(ii)+1)
        } else {
            i=which(ann$geneSym[iA2]==gene)
            if (length(i)==1) i=(min(i)-1):(max(i)+1)
        }
    }
    dat=beta2[i,]; iA2=iA2[i]; statThis10=stat10[i,]; statThis20=stat20[i,]; statThis2i=stat2i[i,]
    iList=list(snpH=NULL,snpW=NULL)
    posTh=100
    snp=snpH
    ii=which(snp$chr==ann$CHR[iA2][1] & snp$pos>=min(ann$MAPINFO[iA2]) & snp$pos<=max(ann$MAPINFO[iA2])  & snp$pv<pThres2)
    i=c()
    for (k in 1:length(ii)) {
        i=c(i,which(ann$MAPINFO[iA2]>=(snp$pos[ii[k]]-posTh) & ann$MAPINFO[iA2]<=(snp$pos[ii[k]]+posTh)))
    }
    iList$snpH=i
    snp=snpW
    ii=which(snp$chr==ann$CHR[iA2][1] & snp$pos>=min(ann$MAPINFO[iA2]) & snp$pos<=max(ann$MAPINFO[iA2])  & snp$pv<pThres2)
    i=c()
    for (k in 1:length(ii)) {
        i=c(i,which(ann$MAPINFO[iA2]>=(snp$pos[ii[k]]-posTh) & ann$MAPINFO[iA2]<=(snp$pos[ii[k]]+posTh)))
    }
    iList$snpW=i
    phen=clin2
    phen$pobwBi=as.integer(phen$pobw>cutoff)
    dat2=xThis=c()
    k=1
    for (i in 1:nrow(dat)) {
        for (grpId in 1:nrow(levelInfo)) {
            j=which(phen$caco==levelInfo$caco[grpId] & phen$pobwBi==levelInfo$pobwBi[grpId])
            x=rep("",nrow(dat))
            dat2=c(dat2,dat[i,j])
            xThis=c(xThis,rep(k,length(j)))
            k=k+1
        }
    }
    ttl=rep("",sum(!duplicated(xThis)))
    ylim=NULL
    if (compFlag=="_zoomSignifInteraction") {
        if (gene=="ARID5B") {
            ylim=c(0,0.17)
        } else if (gene=="ZNF148") {
            ylim=c(0.7,1)
        }
    } else {
        ylim=range(c(dat2),na.rm=T)
    }
    png(paste("boxplot",fName,".png",sep=""),width=4*240,height=4*240)
    boxplot(dat2~xThis,names=ttl,ylim=ylim,main=header,xlab="MB",ylab="Set2: Beta-value",col=colList)
    if (!is.null(pos)) {
        i=which(ann$MAPINFO[iA2]>=pos[1] & ann$MAPINFO[iA2]<=pos[2])
        if (length(i)==0) {
            pos=NULL
        } else {
            i=range(i)
            #lines(i*4-1.5,rep(ylim[2],2),col="green",lwd=3)
        }
    }
    for (grpId in 1:nrow(levelInfo)) {
        j=which(phen$caco==levelInfo$caco[grpId] & phen$pobwBi==levelInfo$pobwBi[grpId])
        x=rep(NA,nrow(dat))
        for (i in 1:nrow(dat)) {
            x[i]=mean(dat[i,j],na.rm=T)
            if (i!=1) {
                #lines(nrow(levelInfo)*((i-1):i)-1.5,x[(i-1):i],lty=ltyList[grpId])
            }
        }
    }
    if (compFlag=="_zoomSignifInteraction") {
        leg=paste(c("set2: caco~meth:pobw"),", pv<",pThres,sep="")
        colThis=c("red")
        if ("snp"%in%addFlag) {
            leg=c(leg,paste(c("hisp","white")," int snp within ",posTh,"bp, pv<",pThres2,sep=""))
            colThis=c(colThis,"palegreen","seagreen")
        }
        if (gene=="ARID5B") {
            legend(nrow(dat)*3,0.169,legend=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw"),fill=colList)
            legend(nrow(dat)*3,0.15,legend=leg,lty="solid",lwd=3,col=colThis)
        } else {
            legend(nrow(dat)*3,0.8,legend=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw"),fill=colList)
            legend(nrow(dat)*3,0.75,legend=leg,lty="solid",lwd=3,col=colThis)
        }
    } else if (compFlag=="_zoomSignifPobw") {
        leg=paste(c("set1 ctrl: meth~pobw","set2 ctrl: meth~pobw"),", pv<",pThres,sep="")
        colThis=c("skyblue","skyblue4")
        if (!is.null(pos)) {
            leg=c(leg,paste("set2 ctrl: meth~pobw BH region fwer<",pThres,sep=""))
            colThis=c(colThis,"green")
        }
        if ("snp"%in%addFlag) {
            leg=c(leg,paste(c("hisp","white")," int snp within ",posTh,"bp, pv<",pThres2,sep=""))
            colThis=c(colThis,"palegreen","seagreen")
        }
        if (gene=="ARID5B") {
            #legend(nrow(dat)*3+1,0.19,legend=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw"),fill=colList)
            #legend(nrow(dat)*3+1,0.1,legend=leg,lty="solid",lwd=3,col=colThis)
            legend(nrow(dat)*3-4,0.29,legend=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw"),fill=colList)
            legend(nrow(dat)*3-4,0.2,legend=leg,lty="solid",lwd=3,col=colThis)
        } else {
            legend(nrow(dat)*3-4,0.29,legend=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw"),fill=colList)
            legend(nrow(dat)*3-4,0.2,legend=leg,lty="solid",lwd=3,col=colThis)
        }
    } else {
        leg=paste(c("set1 ctrl: meth~pobw","set2 ctrl: meth~pobw","set2: caco~meth:pobw"),", pv<",pThres,sep="")
        colThis=c("skyblue","skyblue4","red")
        if (!is.null(pos)) {
            leg=c(leg,paste("set2 ctrl: meth~pobw BH region fwer<",pThres,sep=""))
            colThis=c(colThis,"green")
        }
        if ("snp"%in%addFlag) {
            leg=c(leg,paste(c("hisp","white")," int snp within ",posTh,"bp, pv<",pThres2,sep=""))
            colThis=c(colThis,"palegreen","seagreen")
        }
        if (gene=="ARID5B") {
            leg2=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw")
            col2This=colList
            legend(nrow(dat)*3+1,0.29,legend=leg2,fill=col2This)
            legend(nrow(dat)*3+1,0.2,legend=leg,lty="solid",lwd=3,col=colThis)
        } else {
            legend(nrow(dat)*3-4,0.29,legend=c("ctrl/low pobw","ctrl/high pobw","case/low pobw","case/high pobw"),fill=colList)
            legend(nrow(dat)*3-4,0.2,legend=leg,lty="solid",lwd=3,col=colThis)
        }
    }
    ii=which(statThis10$pv_pobwBi<pThres)
    axis(side=3,at=nrow(levelInfo)*ii-1.5,labels=round(ann$MAPINFO[iA2][ii]/10^3),col.ticks="skyblue",lwd.ticks=4,las=2)
    if (!is.null(pos)) {
        ii=which(ann$MAPINFO[iA2]>=pos[1] & ann$MAPINFO[iA2]<=pos[2]); ii=range(ii)
        axis(side=3,at=nrow(levelInfo)*ii-1.5,labels=round(ann$MAPINFO[iA2][ii]/10^3),col.ticks="green",lwd.ticks=4,las=2)
    }
    ii=which(statThis20$pv_pobwBi<pThres)
    axis(side=1,at=nrow(levelInfo)*ii-1.5,labels=round(ann$MAPINFO[iA2][ii]/10^3),col.ticks="skyblue4",lwd.ticks=4,las=2)
    ii=which(statThis2i$pv_meth.pobwBi<pThres)
    axis(side=1,at=nrow(levelInfo)*ii-1.5,labels=round(ann$MAPINFO[iA2][ii]/10^3),col.ticks="red",lwd.ticks=4,las=2)
    if ("snp"%in%addFlag) {
        ii=iList$snpH
        axis(side=3,at=nrow(levelInfo)*ii-1.5,labels=round(ann$MAPINFO[iA2][ii]/10^3),col.ticks="palegreen",lwd.ticks=4,las=2)
        ii=iList$snpW
        axis(side=3,at=nrow(levelInfo)*ii-1.5,labels=round(ann$MAPINFO[iA2][ii]/10^3),col.ticks="seagreen",lwd.ticks=4,las=2)
    }
    dev.off()
}

#save.image(paste("tmp_",gene,".RData",sep=""))

## -------------------
"
Model 1: Y’ = cX + E1
Model 2: M’ = aX + E2
Model 3: Y’’ = bM + c’X +E3

c'^2 * Var(X) + b^2 * Var(M) + 2*b*c'*Cov(X,M) + p^2/3
sqrt((comp a)^2* SE(comp b)^2 + (comp b)^2* SE(comp a)^2)
"

fName1="_mediation_pobw"
heading1="Mediation model: pobw"

testFlag=c("Sobel","Aroian","Goodman")
testFlag=c("Sobel","Aroian")

if (modelFlag=="pobw->meth->caco") {
    sdX=sqrt(stat_bw_c$var_pobw)
    sdM=sqrt(stat_bw_c$var_meth)
    sdY1=sqrt(stat_bw_1$coef_pobw^2 * stat_bw_c$var_pobw + pi^2/3)
    sdY11=sqrt(stat_bw_3$coef_pobw^2 * stat_bw_c$var_pobw + stat_bw_3$coef_meth^2 * stat_bw_c$var_meth + 2 * stat_bw_3$coef_pobw * stat_bw_3$coef_meth * abs(stat_bw_c$cov) + pi^2/3)

    cA=stat_bw_2$coef_pobw
    cB=stat_bw_3$coef_meth*sdM/sdY11
    cC=stat_bw_1$coef_pobw*sdX/sdY1
    cC1=stat_bw_3$coef_pobw*sdX/sdY11

    seA=stat_bw_2$se_pobw
    seB=stat_bw_3$se_meth*sdM/sdY11
    seC=stat_bw_1$se_pobw*sdX/sdY1
    seC1=stat_bw_3$se_pobw*sdX/sdY11
} else {
    sdX=sqrt(stat_bw_c$var_meth)
    sdM=sqrt(stat_bw_c$var_pobw)
    sdY1=sqrt(stat_bw_1$coef_meth^2 * stat_bw_c$var_meth + pi^2/3)
    sdY11=sqrt(stat_bw_3$coef_meth^2 * stat_bw_c$var_meth + stat_bw_3$coef_pobw^2 * stat_bw_c$var_pobw + 2 * stat_bw_3$coef_meth * stat_bw_3$coef_pobw * abs(stat_bw_c$cov) + pi^2/3)
    
    cA=stat_bw_2$coef_meth
    cB=stat_bw_3$coef_pobw*sdM/sdY11
    cC=stat_bw_1$coef_meth*sdX/sdY1
    cC1=stat_bw_3$coef_meth*sdX/sdY11
    
    seA=stat_bw_2$se_meth
    seB=stat_bw_3$se_pobw*sdM/sdY11
    seC=stat_bw_1$se_meth*sdX/sdY1
    seC1=stat_bw_3$se_meth*sdX/sdY11
}

fName=fName1
heading=heading1

png(paste("plot_abVcc1",fName,".png",sep=""))
x1=cA*cB; x2=cC-cC1
lim=range(c(x1,x2),na.rm=T)
plot(x1,x2,xlim=lim,ylim=lim,main=heading,xlab="a*b",ylab="c-c1")
abline(c(0,1),col="red")
dev.off()

z_S=cA*cB/sqrt(cA^2*seB^2 + cB^2*seA^2)
z_A=cA*cB/sqrt(cA^2*seB^2 + cB^2*seA^2 + seA^2*seB^2)
z_G=cA*cB/sqrt(cA^2*seB^2 + cB^2*seA^2 - seA^2*seB^2)

png(paste("histogram_zValue",fName,".png",sep=""))
if ("Goodman"%in%testFlag) par(mfrow=c(3,1)) else par(mfrow=c(2,1))
hist(z_S,main="Sobel test",xlab="z-value")
hist(z_A,main="Aroian test",xlab="z-value")
if ("Goodman"%in%testFlag) {
    hist(z_G,main="Goodman test",xlab="z-value")
}
dev.off()

pv_S=2*(1-pnorm(abs(z_S),mean=0,sd=1))
pv_A=2*(1-pnorm(abs(z_A),mean=0,sd=1))
pv_G=2*(1-pnorm(abs(z_G),mean=0,sd=1))

png(paste("plot_pValue",fName,".png",sep=""))
x1=pv_S; x2=pv_A
x3=pv_G
lim=range(c(x1,x2),na.rm=T)
lim=range(c(x1,x2,x3),na.rm=T)
if ("Goodman"%in%testFlag) par(mfrow=c(2,2))
plot(x1,x2,xlim=lim,ylim=lim,main=heading,xlab="Sobel test: P-value",ylab="Aroian test: P-value"); abline(c(0,1),col="red")
if ("Goodman"%in%testFlag) {
    plot(x1,x3,xlim=lim,ylim=lim,main=heading,xlab="Sobel test: P-value",ylab="Goodman test: P-value"); abline(c(0,1),col="red")
    plot(x2,x3,xlim=lim,ylim=lim,main=heading,xlab="Aroian test: P-value",ylab="Goodman test: P-value"); abline(c(0,1),col="red")
}
dev.off()

png("tmp.png"); plot(z_S,pv_S); dev.off()

summary(cC)
summary(seC)
summary(stat_bw_1[,paste("pv_",strsplit(modelFlag,"->")[[1]][1],sep="")])

fName=paste("_",strsplit(modelFlag,"->")[[1]][2],"Mediated_cacoResp_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2_bmiq",sep="")
heading=paste("Set2: caco ~ ",strsplit(modelFlag,"->")[[1]][1],", ",strsplit(modelFlag,"->")[[1]][2]," mediated. Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
par(mfcol=c(1,3))
hist(cC,main="",xlab="Coefficient c")
hist(seC,main=heading,xlab="SE(c)")
hist(stat_bw_1[,paste("pv_",strsplit(modelFlag,"->")[[1]][1],sep="")],main="",xlab="P-value")
dev.off()

if ("Goodman"%in%testFlag) {
    stat_bw_z=data.frame(cpgId=stat_bw_1$cpgId,coef_a=cA,coef_b=cB,coef_c=cC,coef_c1=cC1,coef_cMinusC1=cC-cC1,pv_S,pv_A,pv_G)
} else {
    stat_bw_z=data.frame(cpgId=stat_bw_1$cpgId,coef_a=cA,coef_b=cB,coef_c=cC,coef_c1=cC1,coef_cMinusC1=cC-cC1,pv_S,pv_A)
}


## -------------------
library(qvalue)
#source(paste(dirSrc,"functions/TTest.9.1.6.R",sep=""))
source(paste(dirSrc,"functions/TTest.9.1.7.R",sep=""))

#for (compId in paste("med_",c("zS","zA","zG"),sep="")) {
#for (compId in c("p_co_1","p_co_2","mp_2")) {
#for (compId in c("p2_co_1","p2_co_2","mp2_2")) {
#for (compId in c("p_co_m_1","p_co_m_2","mp_m_2")) {
#for (compId in c("p2_co_m_1","p2_co_m_2","mp2_m_2")) {
#for (compId in c("p3_co_3","p3_co_4","mp3_4")) {
#for (compId in c("p4_co_3","p4_co_4","mp4_4")) {
#for (compId in c("aml_cc_1_00","aml_cc_1_0","aml_cc_1_1","aml_cc_1_2","aml_cc_2_00","aml_cc_2_0","aml_cc_2_1","aml_cc_2_2")) {
#for (compId in c("cc_2","cc_fs_2","cc_ms_2","cc_he_2","cc_we_2","cc_3")) {
#for (compId in c("cc_s1","cc_s2","cc_s3","cc_s4","cc_s5","cc_s6")) {
for (compId in c("pcb105","pcb118","pcb138","pcb153","pcb170","pcb180","pcb1260")) {
    colIdPV="pv"
	cat("\n\n==================",compId,"==================\n")
	switch(compId,
            "pcb105"={
                stat2=stat_pcb105
                colIdPV=names(stat2)[grep("pv_",names(stat2))]
            },
            "pcb118"={
                stat2=stat_pcb118
                colIdPV=names(stat2)[grep("pv_",names(stat2))]
            },
            "pcb138"={
                stat2=stat_pcb138
                colIdPV=names(stat2)[grep("pv_",names(stat2))]
            },
            "pcb153"={
                stat2=stat_pcb153
                colIdPV=names(stat2)[grep("pv_",names(stat2))]
            },
            "pcb170"={
                stat2=stat_pcb170
                colIdPV=names(stat2)[grep("pv_",names(stat2))]
            },
            "pcb180"={
                stat2=stat_pcb180
                colIdPV=names(stat2)[grep("pv_",names(stat2))]
            },
            "pcb1260"={
                stat2=stat_pcb1260
                colIdPV=names(stat2)[grep("pv_",names(stat2))]
            },
            "cc_s1"={
                stat2=stat_cc_s1
                colIdPV=names(stat2)[grep("pv_meth.",names(stat2))]
            },
            "cc_s2"={
                stat2=stat_cc_s2
                colIdPV=names(stat2)[grep("pv_meth.",names(stat2))]
            },
            "cc_s3"={
                stat2=stat_cc_s3
                colIdPV=names(stat2)[grep("pv_meth.",names(stat2))]
            },
            "cc_s4"={
                stat2=stat_cc_s4
                colIdPV=names(stat2)[grep("pv_meth.",names(stat2))]
            },
            "cc_s5"={
                stat2=stat_cc_s5
                colIdPV=names(stat2)[grep("pv_meth.",names(stat2))]
            },
            "cc_s6"={
                stat2=stat_cc_s6
                colIdPV=names(stat2)[grep("pv_meth.",names(stat2))]
            },
            "cc_2"={
                stat2=stat_cc_2
                colIdPV=paste("pv_",varFlag,sep="")
            },
            "cc_fs_2"={
                stat2=stat_cc_fs_2
                colIdPV=paste("pv_",varFlag,sep="")
            },
            "cc_ms_2"={
                stat2=stat_cc_ms_2
                colIdPV=paste("pv_",varFlag,sep="")
            },
            "cc_he_2"={
                stat2=stat_cc_he_2
                colIdPV=paste("pv_",varFlag,sep="")
            },
            "cc_we_2"={
                stat2=stat_cc_we_2
                colIdPV=paste("pv_",varFlag,sep="")
            },
            "cc_3"={
                stat2=stat_cc_3
                colIdPV=paste("pv_",varFlag,sep="")
            },
            "med_zS"={
               stat2=stat_bw_z
               colIdPV="pv_S"
		   },
           "med_zA"={
               stat2=stat_bw_z
               colIdPV="pv_A"
           },
           "med_zG"={
               stat2=stat_bw_z
               colIdPV="pv_G"
           },
           "p_co_1"={
               stat2=stat_bw_1_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "p_co_2"={
               stat2=stat_bw_2_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "mp_2"={
               stat2=stat_bw_i_2
               colIdPV=paste("pv_meth.",varFlag,sep="")
           },
           "p2_co_1"={
               stat2=stat_bw_1_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "p2_co_2"={
               stat2=stat_bw_2_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "mp2_2"={
               stat2=stat_bw_i_2
               colIdPV=paste("pv_meth.",varFlag,sep="")
           },
           "p_co_m_1"={
               stat2=stat_bw_1_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "p_co_m_2"={
               stat2=stat_bw_2_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "mp_m_2"={
               stat2=stat_bw_i_2
               colIdPV=paste("pv_meth.",varFlag,sep="")
           },
           "p2_co_m_1"={
               stat2=stat_bw_1_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "p2_co_m_2"={
               stat2=stat_bw_2_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "mp2_m_2"={
               stat2=stat_bw_i_2
               colIdPV=paste("pv_meth.",varFlag,sep="")
           },
           "p3_co_3"={
               stat2=stat_bw_3_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "p3_co_4"={
               stat2=stat_bw_4_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "mp3_4"={
               stat2=stat_bw_i_4
               colIdPV=paste("pv_meth.",varFlag,sep="")
           },
           "p4_co_3"={
               stat2=stat_bw_3_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "p4_co_4"={
               stat2=stat_bw_4_0
               colIdPV=paste("pv_",varFlag,sep="")
           },
           "mp4_4"={
               stat2=stat_bw_i_4
               colIdPV=paste("pv_meth.",varFlag,sep="")
           },
           "aml_cc_1_00"={
               stat2=stat_aml_1_00
               colIdPV="pv_caco"
           },
           "aml_cc_1_0"={
               stat2=stat_aml_1_0
               colIdPV="pv_caco"
           },
           "aml_cc_1_1"={
               stat2=stat_aml_1_1
               colIdPV="pv_caco"
           },
           "aml_cc_1_2"={
               stat2=stat_aml_1_2
               colIdPV="pv_caco"
           },
           "aml_cc_2_00"={
               stat2=stat_aml_2_00
               colIdPV="pv_meth"
           },
           "aml_cc_2_0"={
               stat2=stat_aml_2_0
               colIdPV="pv_meth"
           },
           "aml_cc_2_1"={
               stat2=stat_aml_2_1
               colIdPV="pv_meth"
           },
           "aml_cc_2_2"={
               stat2=stat_aml_2_2
               colIdPV="pv_meth"
           }
		   )
	iA2=match(stat2$cpgId,ann$IlmnID)
	stat2$qv=NA
    #i=which(ann$snp[iA2]==0 & ann$CHR[iA2]%in%1:22)
    i=which(ann$keep[iA2]==1)
    stat2$qv[i]=getAdjustedPvalue(stat2[i,colIdPV],method="qv",strict=T)
	if (all(is.na(stat2$qv))) {
		cat("Using BH method !!!\n")
        stat2$qv[i]=getAdjustedPvalue(stat2[i,colIdPV],method="qv",strict=F)
        #stat2$qv[i]=p.adjust(stat2[i,colIdPV],method="BH")
    }
	i=which(stat2$qv<.05)
	if (any(stat2$qv[i]<stat2[i,colIdPV])) {cat("Q-value < p-value !!!\n")}
    names(stat2)[which(names(stat2)=="qv")]=sub("pv","qv",colIdPV)
    i=which(ann$keep[iA2]==1)
    stat2$lambda=stat2$qvGI=stat2$pvGI=NA
    stat2$lambda[i]=qchisq(median(stat2[i,colIdPV],na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
    stat2$pvGI[i]=stat2[i,colIdPV]*stat2$lambda[i]
    stat2$qvGI[i]=getAdjustedPvalue(stat2[i,"pvGI"],method="qv",strict=T)
    if (all(is.na(stat2$qvGI))) {
        cat("Using BH method !!!\n")
        stat2$qvGI[i]=getAdjustedPvalue(stat2[i,"pvGI"],method="qv",strict=F)
        #stat2$qvGI[i]=p.adjust(stat2[i,"pvGI"],method="BH")
    }
    i=which(stat2$qvGI<.05)
    if (any(stat2$qvGI[i]<stat2[i,"pvGI"])) {cat("Q-value < p-value !!!\n")}
	switch(compId,
            "pcb105"={
                stat_pcb105=stat2
            },
            "pcb118"={
                stat_pcb118=stat2
            },
            "pcb138"={
                stat_pcb138=stat2
            },
            "pcb153"={
                stat_pcb153=stat2
            },
            "pcb170"={
                stat_pcb170=stat2
            },
            "pcb180"={
                stat_pcb180=stat2
            },
            "pcb1260"={
                stat_pcb1260=stat2
            },
            "cc_s1"={
                stat_cc_s1=stat2
            },
            "cc_s2"={
                stat_cc_s2=stat2
            },
            "cc_s3"={
                stat_cc_s3=stat2
            },
            "cc_s4"={
                stat_cc_s4=stat2
            },
            "cc_s5"={
                stat_cc_s5=stat2
            },
            "cc_s6"={
                stat_cc_s6=stat2
            },
            "cc_2"={
                stat_cc_2=stat2
            },
            "cc_fs_2"={
                stat_cc_fs_2=stat2
            },
            "cc_ms_2"={
                stat_cc_ms_2=stat2
            },
            "cc_he_2"={
                stat_cc_he_2=stat2
            },
            "cc_we_2"={
                stat_cc_we_2=stat2
            },
            "cc_3"={
                stat_cc_3=stat2
            },
            "med_zS"={
               stat_bw_z=stat2
		   },
           "med_zA"={
               stat_bw_z=stat2
           },
           "med_zG"={
               stat_bw_z=stat2
           },
           "p_co_1"={
               stat_bw_1_0=stat2
           },
           "p_co_2"={
               stat_bw_2_0=stat2
           },
           "mp_2"={
               stat_bw_i_2=stat2
           },
           "p2_co_1"={
               stat_bw_1_0=stat2
           },
           "p2_co_2"={
               stat_bw_2_0=stat2
           },
           "mp2_2"={
               stat_bw_i_2=stat2
           },
           "p_co_m_1"={
               stat_bw_1_0=stat2
           },
           "p_co_m_2"={
               stat_bw_2_0=stat2
           },
           "mp_m_2"={
               stat_bw_i_2=stat2
           },
           "p2_co_m_1"={
               stat_bw_1_0=stat2
           },
           "p2_co_m_2"={
               stat_bw_2_0=stat2
           },
           "mp2_m_2"={
               stat_bw_i_2=stat2
           },
           "p3_co_3"={
               stat_bw_3_0=stat2
           },
           "p3_co_4"={
               stat_bw_4_0=stat2
           },
           "mp3_4"={
               stat_bw_i_4=stat2
           },
           "p4_co_3"={
               stat_bw_3_0=stat2
           },
           "p4_co_4"={
               stat_bw_4_0=stat2
           },
           "mp4_4"={
               stat_bw_i_4=stat2
           },
           "aml_cc_1_00"={
               stat_aml_1_00=stat2
           },
           "aml_cc_1_0"={
               stat_aml_1_0=stat2
           },
           "aml_cc_1_1"={
               stat_aml_1_1=stat2
           },
           "aml_cc_1_2"={
               stat_aml_1_2=stat2
           },
           "aml_cc_2_00"={
               stat_aml_2_00=stat2
           },
           "aml_cc_2_0"={
               stat_aml_2_0=stat2
           },
           "aml_cc_2_1"={
               stat_aml_2_1=stat2
           },
           "aml_cc_2_2"={
               stat_aml_2_2=stat2
           }
		   )
}

save.image("tmp2.RData")


png("tmp_%01d.png")
for (k in grep("coef_",names(stat_bw_z))) {
    lim=range(c(stat_bw_i_2$coef_interaction,stat_bw_z[,k]),na.rm=T)
    plot(stat_bw_i_2$coef_interaction,stat_bw_z[,k],xlim=lim,ylim=lim,ylab=names(stat_bw_z)[k])
    abline(c(0,1))
}
dev.off()


## --------------------------------------------------------

pThres=0.05
colIdPV="pv"

stat2=stat_bw_z

iA2=match(stat2$cpgId,ann$IlmnID)
ii=which(stat_bw_z$pv_S<pThres | stat_bw_z$pv_A<pThres)
if (F) {
    tbl=stat_bw_1[ii,which(!names(stat_bw_1)%in%c("cpgId"))]
    names(tbl)=paste(names(tbl),"_c",sep="")
    tbl2=stat_bw_2[ii,which(!names(stat_bw_2)%in%c("cpgId"))]
    names(tbl2)=paste(names(tbl2),"_a",sep="")
    tbl=cbind(tbl,tbl2)
    tbl2=stat_bw_3[ii,which(!names(stat_bw_3)%in%c("cpgId"))]
    names(tbl2)=paste(names(tbl2),rep(c("_b","_c1"),each=3),sep="")
    tbl=cbind(tbl,tbl2)
    tbl2=stat_bw_z[ii,which(!names(stat_bw_z)%in%c("cpgId"))]
    names(tbl2)=c("coef_cMinusC1","pv_Sobel","pv_Aroian","qv_Sobel","qv_Aroian")
    tbl=cbind(tbl,tbl2)
}
if (length(ii)!=0) {
    tbl=stat_bw_z[ii,which(!names(stat_bw_z)%in%c("cpgId"))]
    names(tbl)=c("coef_a","coef_b","coef_c","coef_c1","coef_cMinusC1","pv_Sobel","pv_Aroian","qv_Sobel","qv_Aroian")
    tbl=cbind(ann[iA2,][ii,],tbl)
    write.table(paste("Loci with Sobel or Aroian test p-value < ",pThres,sep=""), file=paste("stat",fName,"_",colIdPV,pThres,".txt",sep=""), append=F,col.names=F,row.names=F, sep="\t",quote=F)
    write.table(tbl, file=paste("stat",fName,"_",colIdPV,pThres,".txt",sep=""), append=T,col.names=T,row.names=F, sep="\t",quote=F)
}

## -----------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------


####################################################################
####################################################################
####################################################################
####################################################################
## ----------------------------------------------

plotFlag=""
plotFlag="_onePlot"

outlierFlag=F
outlierFlag=T

## -------------------
#for (compId in paste("med_",c("1","2","3X","3M","zS","zA"),sep="")) {
#for (compId in c("p_co_1","p_co_2","mp_2")) {
#for (compId in c("p2_co_1","p2_co_2","mp2_2")) {
#for (compId in c("p_co_m_1","p_co_m_2","mp_m_2")) {
#for (compId in c("p2_co_m_1","p2_co_m_2","mp2_m_2")) {
#for (compId in c("p3_co_3","p3_co_4","mp3_4")) {
#for (compId in c("p4_co_3","p4_co_4","mp4_4")) {
#for (compId in c("aml_cc_1_0","aml_cc_1_1","aml_cc_1_2","aml_cc_2_00","aml_cc_2_0","aml_cc_2_1","aml_cc_2_2")) {
#for (compId in c("cc_2","cc_fs_2","cc_ms_2","cc_he_2","cc_we_2","cc_3")) {
#for (compId in c("cc_s1","cc_s2","cc_s3","cc_s4","cc_s5","cc_s6")) {
for (compId in c("pcb105","pcb118","pcb138","pcb153","pcb170","pcb180","pcb1260")) {
    colIdEst="coef"; colIdPV=c("pv","pv"); 	pThres=0.001
    colIdEst="coef"; colIdPV=c("pv","qv"); 	pThres=0.1
    colIdEst="coef"; colIdPV=c("pv","qv"); 	pThres=0.01
    colIdEst="coef"; colIdPV=c("pv","qv"); 	pThres=10^-28
    colIdEst="coef"; colIdPV=c("pv","pv"); 	pThres=10^-28
    colIdEst="coef"; colIdPV=c("pv","pv"); 	pThres=0.05
    colIdEst="coef"; colIdPV=c("pvGI","qvGI"); 	pThres=0.05
    colIdEst="coef"; colIdPV=c("pv","qv"); 	pThres=0.05
    
	switch(compId,
            "pcb105"={
                stat2=stat_pcb105; fName1=paste("_methResp_logged_PCB_105_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,sep=""); compName=paste("M-value based\nSet1+Set2 ctrl: meth ~ pcb105\nCov: set, ReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "pcb118"={
                stat2=stat_pcb118; fName1=paste("_methResp_logged_PCB_118_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,sep=""); compName=paste("M-value based\nSet1+Set2 ctrl: meth ~ pcb118\nCov: set, ReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "pcb138"={
                stat2=stat_pcb138; fName1=paste("_methResp_logged_PCB_138_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,sep=""); compName=paste("M-value based\nSet1+Set2 ctrl: meth ~ pcb138\nCov: set, ReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "pcb153"={
                stat2=stat_pcb153; fName1=paste("_methResp_logged_PCB_153_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,sep=""); compName=paste("M-value based\nSet1+Set2 ctrl: meth ~ pcb153\nCov: set, ReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "pcb170"={
                stat2=stat_pcb170; fName1=paste("_methResp_logged_PCB_170_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,sep=""); compName=paste("M-value based\nSet1+Set2 ctrl: meth ~ pcb170\nCov: set, ReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "pcb180"={
                stat2=stat_pcb180; fName1=paste("_methResp_logged_PCB_180_SRS_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,sep=""); compName=paste("M-value based\nSet1+Set2 ctrl: meth ~ pcb180\nCov: set, ReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "pcb1260"={
                stat2=stat_pcb1260; fName1=paste("_methResp_logged_PCB_aroclor1260_ctrlSubset_covSet_covPrinComp1234_covEpStr_allGuthSet1Set2_bmiq",transformFlag,sep=""); compName=paste("M-value based\nSet1+Set2 ctrl: meth ~ aroclor1260\nCov: set, ReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "cc_s1"={
                stat2=stat_cc_s1; fName1=paste("_",varFlag,"Resp_methXsem_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,sep=""); compName=paste("M-value based, coefficient - interaction\nSet2: ",varFlag," ~ meth * sem\nCov: Sex, gestage, ReFACTor comp 1,2,3,4, epistructure",sep="")
                k=grep("pv_meth.",names(stat2))
                x=sub("pv","",names(stat2)[k])
                names(stat2)=sub(x,"",names(stat2),fixed=T)
            },
            "cc_s2"={
                stat2=stat_cc_s2; fName1=paste("_",varFlag,"Resp_methXlogSem_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,sep=""); compName=paste("M-value based, coefficient - interaction\nSet2: ",varFlag," ~ meth * log(sem)\nCov: Sex, gestage, ReFACTor comp 1,2,3,4, epistructure",sep="")
                k=grep("pv_meth.",names(stat2))
                x=sub("pv","",names(stat2)[k])
                names(stat2)=sub(x,"",names(stat2),fixed=T)
            },
            "cc_s3"={
                stat2=stat_cc_s3; fName1=paste("_",varFlag,"Resp_methXsemNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,sep=""); compName=paste("M-value based, coefficient - interaction\nSet2: ",varFlag," ~ meth * sem\nCov: Sex, gestage, ReFACTor comp 1,2,3,4, epistructure. No sex chroms for SEM",sep="")
                k=grep("pv_meth.",names(stat2))
                x=sub("pv","",names(stat2)[k])
                names(stat2)=sub(x,"",names(stat2),fixed=T)
            },
            "cc_s4"={
                stat2=stat_cc_s4; fName1=paste("_",varFlag,"Resp_methXlogSemNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,sep=""); compName=paste("M-value based, coefficient - interaction\nSet2: ",varFlag," ~ meth * log(sem)\nCov: Sex, gestage, ReFACTor comp 1,2,3,4, epistructure. No sex chroms for SEM",sep="")
                k=grep("pv_meth.",names(stat2))
                x=sub("pv","",names(stat2)[k])
                names(stat2)=sub(x,"",names(stat2),fixed=T)
            },
            "cc_s5"={
                stat2=stat_cc_s5; fName1=paste("_",varFlag,"Resp_methXsemNoSnpNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,sep=""); compName=paste("M-value based, coefficient - interaction\nSet2: ",varFlag," ~ meth * sem\nCov: Sex, gestage, ReFACTor comp 1,2,3,4, epistructure. No SNP, no sex chroms for SEM",sep="")
                k=grep("pv_meth.",names(stat2))
                x=sub("pv","",names(stat2)[k])
                names(stat2)=sub(x,"",names(stat2),fixed=T)
            },
            "cc_s6"={
                stat2=stat_cc_s6; fName1=paste("_",varFlag,"Resp_methXlogSemNoSnpNoSexChr_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2_bmiq",transformFlag,sep=""); compName=paste("M-value based, coefficient - interaction\nSet2: ",varFlag," ~ meth * log(sem)\nCov: Sex, gestage, ReFACTor comp 1,2,3,4, epistructure. No SNP, no sex chroms for SEM",sep="")
                k=grep("pv_meth.",names(stat2))
                x=sub("pv","",names(stat2)[k])
                names(stat2)=sub(x,"",names(stat2),fixed=T)
            },
            "cc_2"={
                stat2=stat_cc_2; fName1=paste("_methResp_",varFlag,"_covSexGestage_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2: meth ~ ",varFlag,"\nCov: Sex, gestage, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "cc_fs_2"={
                stat2=stat_cc_fs_2; fName1=paste("_methResp_",varFlag,"_femaleSubset_covGestage_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2 female: meth ~ ",varFlag,"\nCov: Gestage, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "cc_ms_2"={
                stat2=stat_cc_ms_2; fName1=paste("_methResp_",varFlag,"_maleSubset_covGestage_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2 male: meth ~ ",varFlag,"\nCov: Gestage, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "cc_he_2"={
                stat2=stat_cc_he_2; fName1=paste("_methResp_",varFlag,"_hispSubset_covSexGestage_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2 hispanic: meth ~ ",varFlag,"\nCov: Sex, gestage, ReFACTor comp 1,2,3,4",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "cc_we_2"={
                stat2=stat_cc_we_2; fName1=paste("_methResp_",varFlag,"_noHistWhiteSubset_covSexGestage_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2 non-hispanic white: meth ~ ",varFlag,"\nCov: Sex, gestage, ReFACTor comp 1,2,3,4",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "cc_3"={
                stat2=stat_cc_3; fName1=paste("_methResp_",varFlag,"_covSexGestage_covPrinComp1234_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2: meth ~ ",varFlag,", cov: sex, gestage,\nReFACTor comp 1,2,3,4, epistructure",sep="")
                names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
            },
            "med_1"={
               stat2=stat_bw_1; fName1=paste("_mediation_cacoResp_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2",sep=""); compName=paste("Coefficient: c (",strsplit(modelFlag,"->")[[1]][1],")\nSet2: caco ~ ",strsplit(modelFlag,"->")[[1]][1],". Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_",strsplit(modelFlag,"->")[[1]][1],sep=""),"",names(stat2))
		   },
           "med_2"={
               stat2=stat_bw_2; fName1=paste("_mediation_",strsplit(modelFlag,"->")[[1]][2],"Resp_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2",sep=""); compName=paste("Coefficient: a (",strsplit(modelFlag,"->")[[1]][1],")\nSet2: ",strsplit(modelFlag,"->")[[1]][2]," ~ ",strsplit(modelFlag,"->")[[1]][1],". Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_",strsplit(modelFlag,"->")[[1]][1],sep=""),"",names(stat2))
           },
           "med_3X"={
               stat2=stat_bw_3; fName1=paste("_mediation_cacoResp_",strsplit(modelFlag,"->")[[1]][2],"_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2",sep=""); compName=paste("Coefficient: c' (",strsplit(modelFlag,"->")[[1]][1],")\nSet2 : caco ~ ",strsplit(modelFlag,"->")[[1]][2]," + ",strsplit(modelFlag,"->")[[1]][1],". Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_",strsplit(modelFlag,"->")[[1]][1],sep=""),"",names(stat2))
           },
           "med_3M"={
               stat2=stat_bw_3; fName1=paste("_mediation_cacoResp_",strsplit(modelFlag,"->")[[1]][2],"_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2",sep=""); compName=paste("Coefficient: b (",strsplit(modelFlag,"->")[[1]][2],")\nSet2 : caco ~ ",strsplit(modelFlag,"->")[[1]][2]," + ",strsplit(modelFlag,"->")[[1]][1],". Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_",strsplit(modelFlag,"->")[[1]][2],sep=""),"",names(stat2))
           },
           "med_zS"={
               stat2=stat_bw_z; fName1=paste("_mediation_ztestSobel_cacoResp_",strsplit(modelFlag,"->")[[1]][2],"_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2",sep=""); compName=paste("Coefficient: c-c', Sobel test\nSet2 : caco ~ ",strsplit(modelFlag,"->")[[1]][2]," + ",strsplit(modelFlag,"->")[[1]][1],". Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
               names(stat2)=gsub("_cMinusC1|_S","",names(stat2))
           },
           "med_zA"={
               stat2=stat_bw_z; fName1=paste("_mediation_ztestAroian_cacoResp_",strsplit(modelFlag,"->")[[1]][2],"_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2",sep=""); compName=paste("Coefficient: c-c', Aroian test\nSet2 : caco ~ ",strsplit(modelFlag,"->")[[1]][2]," + ",strsplit(modelFlag,"->")[[1]][1],". Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
               names(stat2)=gsub("_cMinusC1|_A","",names(stat2))
           },
           "med_zG"={
               stat2=stat_bw_z; fName1=paste("_mediation_ztestGoodman_cacoResp_",strsplit(modelFlag,"->")[[1]][2],"_",strsplit(modelFlag,"->")[[1]][1],"_covSexHisp_covPrinComp1234_allGuthSet2",sep=""); compName=paste("Coefficient: c-c', Goodman test\nSet2 : caco ~ ",strsplit(modelFlag,"->")[[1]][2]," + ",strsplit(modelFlag,"->")[[1]][1],". Cov: Sex, hisp,\nReFACTor comp 1,2,3,4",sep="")
               names(stat2)=gsub("_cMinusC1|_G","",names(stat2))
           },
           "p_co_1"={
               stat2=stat_bw_1_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1",transformFlag,sep=""); compName=paste("Beta-value based\nSet1: meth ~ ",varFlag,"\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p_co_2"={
               stat2=stat_bw_2_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nSet2: meth ~ ",varFlag,"\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp_2"={
               stat2=stat_bw_i_2; fName1=paste("_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag,"\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "p2_co_1"={
               stat2=stat_bw_1_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1",transformFlag,sep=""); compName=paste("Beta-value based\nSet1 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p2_co_2"={
               stat2=stat_bw_2_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nSet2 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp2_2"={
               stat2=stat_bw_i_2; fName1=paste("_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "p3_co_3"={
               stat2=stat_bw_3_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet1",transformFlag,sep=""); compName=paste("Beta-value based\nSet1: meth ~ ",varFlag,"\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p3_co_4"={
               stat2=stat_bw_4_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nSet2: meth ~ ",varFlag,"\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp3_4"={
               stat2=stat_bw_i_4; fName1=paste("_cacoResp_methX",varFlag,"_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag,"\nCov: Epistructure",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "p4_co_3"={
               stat2=stat_bw_3_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet1",transformFlag,sep=""); compName=paste("Beta-value based\nSet1 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p4_co_4"={
               stat2=stat_bw_4_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nSet2 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp4_4"={
               stat2=stat_bw_i_4; fName1=paste("_cacoResp_methX",varFlag,"_covEpStr_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: Epistructure",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "p_co_m_1"={
               stat2=stat_bw_1_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1",transformFlag,sep=""); compName=paste("M-value based\nSet1: meth ~ ",varFlag,"\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p_co_m_2"={
               stat2=stat_bw_2_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2: meth ~ ",varFlag,"\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp_m_2"={
               stat2=stat_bw_i_2; fName1=paste("_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag,"\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "p2_co_m_1"={
               stat2=stat_bw_1_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1",transformFlag,sep=""); compName=paste("M-value based\nSet1 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p2_co_m_2"={
               stat2=stat_bw_2_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nSet2 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp2_m_2"={
               stat2=stat_bw_i_2; fName1=paste("_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("M-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: ReFACTor comp 1,2,3,4",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "p3_co_1"={
               stat2=stat_bw_3_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1",transformFlag,sep=""); compName=paste("Beta-value based\nSet1: meth ~ ",varFlag,"\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p3_co_2"={
               stat2=stat_bw_4_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nSet2: meth ~ ",varFlag,"\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp3_2"={
               stat2=stat_bw_i_4; fName1=paste("_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag,"\nCov: Epistructure",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "p4_co_1"={
               stat2=stat_bw_3_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet1",transformFlag,sep=""); compName=paste("Beta-value based\nSet1 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "p4_co_2"={
               stat2=stat_bw_4_0; fName1=paste("_methResp_",varFlag,"_ctrlSubset_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nSet2 ctrl: meth ~ ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: Epistructure",sep="")
               names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
           },
           "mp4_2"={
               stat2=stat_bw_i_4; fName1=paste("_cacoResp_methX",varFlag,"_covPrinComp1234_allGuthSet2",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: interaction\nSet2: caco ~ meth * ",varFlag," (",varFlag,"<=",cutoff," vs. ",varFlag,">",cutoff,")\nCov: Epistructure",sep="")
               names(stat2)=sub(paste("_meth.",varFlag,sep=""),"",names(stat2),fixed=T)
           },
           "aml_cc_1_0"={
               stat2=stat_aml_1_0; fName1=paste("_methResp_caco_covSexEthnGestage_covEpStr_aml_bmiq",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: caco\nAML: meth ~ caco\nCov: sex, ethn, gestage, epistr",sep="")
               names(stat2)=sub("_caco","",names(stat2),fixed=T)
           },
           "aml_cc_1_1"={
               stat2=stat_aml_1_1; fName1=paste("_methResp_caco_hispSubset_covSexGestage_covEpStr_aml_bmiq",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: caco\nAML hisp: meth ~ caco\nCov: sex, gestage, epistr",sep="")
               names(stat2)=sub("_caco","",names(stat2),fixed=T)
           },
           "aml_cc_1_2"={
               stat2=stat_aml_1_2; fName1=paste("_methResp_caco_noHispWtSubset_covSexGestage_covEpStr_aml_bmiq",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: caco\nAML noHispWt: meth ~ caco\nCov: sex, gestage, epistr",sep="")
               names(stat2)=sub("_caco","",names(stat2),fixed=T)
           },
           "aml_cc_2_00"={
               stat2=stat_aml_2_00; fName1=paste("_cacoResp_meth_covSexGestage_covEpStr_aml_bmiq",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: meth\nAML: caco ~ meth\nCov: sex, gestage, epistr",sep="")
               names(stat2)=sub("_meth","",names(stat2),fixed=T)
           },
           "aml_cc_2_0"={
               stat2=stat_aml_2_0; fName1=paste("_cacoResp_meth_covSexEthnGestage_covEpStr_aml_bmiq",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: meth\nAML: caco ~ meth\nCov: sex, ethn, gestage, epistr",sep="")
               names(stat2)=sub("_meth","",names(stat2),fixed=T)
           },
           "aml_cc_2_1"={
               stat2=stat_aml_2_1; fName1=paste("_cacoResp_meth_hispSubset_covSexGestage_covEpStr_aml_bmiq",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: meth\nAML hisp: caco ~ meth\nCov: sex, gestage, epistr",sep="")
               names(stat2)=sub("_meth","",names(stat2),fixed=T)
           },
           "aml_cc_2_2"={
               stat2=stat_aml_2_2; fName1=paste("_cacoResp_meth_noHispWtSubset_covSexGestage_covEpStr_aml_bmiq",transformFlag,sep=""); compName=paste("Beta-value based\nCoefficient: meth\nAML noHispWt: caco ~ meth\nCov: sex, gestage, epistr",sep="")
               names(stat2)=sub("_meth","",names(stat2),fixed=T)
           }
	)
	i=match(stat2$cpgId,ann$IlmnID)
    #stat2=stat2[which(ann$snp[i]==0),]
    stat2=stat2[which(ann$keep[i]==1),]
    iA2=match(stat2$cpgId,ann$IlmnID)
    #iL=match(stat2$cpgId,dmr$CpG)
    #table(is.na(iL))
	
	ann2=ann[iA2,]
	
####################################################################
	
	cat("\n\n",compName,"\n",sep="")
	i=which(ann$keep[iA2]==1)
	stat=stat2
	fName=fName1
	
	x1=matrix(0,nrow=2,ncol=2,dimnames=list(c("dn","up"),paste(colIdPV[2],c(">=","<"),pThres,sep="")))
	ii=i[which(stat[i,colIdEst]!=0)]
	x2=table(stat[ii,colIdEst]>0,stat[ii,colIdPV[2]]<pThres)
	x1[match(rownames(x2),c("FALSE","TRUE")),match(colnames(x2),c("FALSE","TRUE"))]=x2
	print(x1)
	
	if (plotFlag=="_onePlot") {
		png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
		par(mfcol=c(1,3))
	}
	
	if (plotFlag=="") {
		png(paste("qqplot",fName,".png",sep=""))
		header=compName
	} else {
		header=""
	}
	pvs <- sort(na.exclude(stat[i,colIdPV[1]]))
	qqplot(-log10(runif(length(pvs),0,1)),-log10(pvs),xlab="Expected -log10(p-values) by random",ylab="Observed -log10(p-values)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
	abline(0,1)
	if (plotFlag=="") {
		dev.off()
	}
	
	if (plotFlag=="") {
		png(paste("histogram",fName,".png",sep=""))
		header=compName
	} else {
		header=""
	}
	hist(stat[i,colIdPV[1]],xlab="P-value",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
	if (plotFlag=="") {
		dev.off()
	}
#	title(main=sub(" ReFACTor","\nReFACTor",compName))
	title(main=compName)

	if (plotFlag=="") {
		png(paste("volcanoPlot",fName,"_",colIdPV[1],pThres,".png",sep=""))
		header=compName
	} else {
		header=""
	}
    iThis=i
    if (!outlierFlag) {
        x=quantile(abs(stat[iThis,colIdEst]),probs=.95,na.rm=T)
        iThis=iThis[which(abs(stat[iThis,colIdEst])<=x)]
    }
	plot(stat[iThis,colIdEst],-log10(stat[iThis,colIdPV[1]]),xlab="Estimate",ylab="-log10(p-value)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
	ii=iThis[which(stat[iThis,colIdPV[2]]<pThres)]
	points(stat[ii,colIdEst],-log10(stat[ii,colIdPV[1]]),col="red")
    if (plotFlag=="") {
		dev.off()
	}
    
    #plot(stat2$coef,-log10(stat2$pv)); iii=which(stat2$qv<.05); points(stat2$coef[iii],-log10(stat2$pv)[iii],col="red")
    
    if (plotFlag=="_onePlot") {
#		par(mfrow=c(1,1))
#		title(main=compName)
		dev.off()
	}
	
	if (F) {
	
	png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
	par(mfcol=c(1,3))
	plot(1:6)
	plot(1:6)
	title(main=sub(" ReFACTor","\nReFACTor",compName))
	plot(1:6)
	dev.off()
	}

ii=i[order(stat[i,colIdPV[2]])]
ii=ii[stat[ii,colIdPV[2]]<pThres]
tbl=cbind(ann[iA2,][ii,],stat[ii,c("coef","pv")])
write.table(tbl, file=paste("stat",fName,"_",colIdPV[2],pThres,".txt",sep=""), append=F,col.names=T,row.names=F, sep="\t",quote=F)

####################################################################
## Gene level summarization of p-values
####################################################################
	
	pvFlag=""
	
	k=which(colnames(stat2)==colIdPV[2])
	kk=grep("Mean",colnames(stat2)[k])
	if (length(kk)!=0) k=k[-kk]
	
	fName2=fName
    # pThres=.001; prId=which(ann2$CHR%in%1:22 & !ann2$IlmnID%in%snpVec & !is.na(stat2[,k]) & stat2[,k]<pThres)
    # prId=which(ann2$CHR%in%1:22 & !ann2$IlmnID%in%snpVec & !is.na(stat2[,k]) & stat2[,k]<pThres)
    prId=which(ann2$keep==1 & !is.na(stat2[,k]) & stat2[,k]<pThres)
    
	if (length(prId)!=0) {
		annotSel=ann2[prId,]
		pval=stat2[,k][prId]
		coef=stat2[,colIdEst[1]][prId]
		signGiven=stat2$sLEU[prId]
		
		
		if (F) {
            # Assign polycomb status to each CpG
			tmp1 = strsplit(annotSel$UCSC_RefGene_Accession,";")
			names(tmp1)=paste(1:length(tmp1),":",sep="")
			tmp2 = unlist(tmp1)
			tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
			tmp4 = data.frame(UCSC_REFGENE_ACCESSION=tmp2,
							  rowid=tmp3, stringsAsFactors=FALSE)
			load("PolycombComplete-120109.RData")
			tmp5 = merge(polycombTab,tmp4)
			tmp6 = unique(tmp5$rowid)
			annotSel$PcG = rep(0, dim(annotSel)[1])
			annotSel$PcG[tmp6] = 1 
		}
		
        # Get an index of CpGs by gene region
		tmp1 = strsplit(annotSel$UCSC_RefGene_Name,";")
		names(tmp1)=paste(1:length(tmp1),":",sep="")
		tmp2 = unlist(tmp1)
		tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
		tmp4 = strsplit(annotSel$UCSC_RefGene_Group,";")
		names(tmp4)=paste(1:length(tmp4),":",sep="")
		tmp5 = unlist(tmp4)
		tmp6 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp5)))
		all(tmp3==tmp6)
		
		if (F) {
			library(qvalue)
			qval = qvalue(pval)
			qval$pi0
			qThresh = max(pval[qval$qv<=0.05])
			qThresh
		}
		
		GeneAnnotation = unique(data.frame(UCSC_REFGENE_NAME=tmp2,UCSC_REFGENE_GROUP=tmp5,rowid=tmp3, stringsAsFactors=FALSE))
		GeneIndex = split(GeneAnnotation$rowid,with(GeneAnnotation,paste(UCSC_REFGENE_NAME,UCSC_REFGENE_GROUP,sep=":")))
		GeneIndexN = sapply(GeneIndex, length)
		
		if (length(GeneIndexN)!=0) {
			medPval = sapply(GeneIndex, function(u) median(pval[u],na.rm=T))
			propHit = sapply(GeneIndex, function(u) mean(pval[u]<pThres,na.rm=T))
			
			#isPcG = sapply(GeneIndex, function(u) min(annotSel$PcG[u]))
			#isNearPcG = sapply(GeneIndex, function(u) max(annotSel$PcG[u]))
			
			tmp1 = sapply(GeneIndex, function(u) u[which.min(annotSel$MAPINFO[u])])
			tmp2 = sapply(GeneIndex, function(u) u[which.max(annotSel$MAPINFO[u])])
			GeneSym1 = annotSel$UCSC_RefGene_Name[tmp1]
			GeneSym2 = annotSel$UCSC_RefGene_Name[tmp2]
			
			annotSelMap <- as.numeric(annotSel$MAPINFO)
			tmpF <- function(u) {
				sgn = sign(coef[u])
				pv = pval[u]
				sgnChar = ifelse(pv>pThres, ".", ifelse(sgn<0,"-","+"))
				paste(sgnChar[order(annotSel$CHR[u], annotSelMap[u])],collapse="")
			}
			#Check
            #if (length(GeneIndex)>2) print(tmpF(GeneIndex[[3]]))
			
			hypohyper =sapply(GeneIndex, tmpF)
			
			combineSigns <- function(u) {
				sgn = signGiven[u]
				sgnChar = ifelse(sgn==0, ".", ifelse(sgn<0,"-","+"))
				paste(sgnChar[order(annotSel$CHR[u], annotSelMap[u])],collapse="")
			}
			hypohyper2 =sapply(GeneIndex, combineSigns)
			
			tmp = strsplit(names(GeneIndexN),":")
            #GeneResults = data.frame(Gene=sapply(tmp,function(u)u[1]),Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN, Sign=hypohyper, SignLeuk=hypohyper2, medPval, propHit, stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)
            GeneResults = data.frame(Gene=sapply(tmp,function(u)u[1]),Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN, Sign=hypohyper, medPval, propHit, stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)
			
			rownames(GeneResults)=NULL
			ord = order(medPval)
			
			if (pThres>1) {
				fName=paste("geneSummary_top500_",ifelse(pvFlag=="",colIdPV[2],"pvPerm"),fName2,".txt",sep="")
				write.table(GeneResults[ord[1:500],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
				
				fName=paste("geneSummary_med",ifelse(pvFlag=="",colIdPV[2],"PvPerm"),".001",fName2,".txt",sep="")
				write.table(GeneResults[ord[medPval[ord]<.001],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
				
				fName=paste("geneSummary_med",ifelse(pvFlag=="",colIdPV[2],"PvPerm"),".05",fName2,".txt",sep="")
				write.table(GeneResults[ord[medPval[ord]<.05],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
			} else {
				fName=paste("geneSummary_",ifelse(pvFlag=="",colIdPV[2],"pvPerm"),pThres,fName2,".txt",sep="")
				write.table(GeneResults[ord,which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
			}
		}
		
		tbl=GeneResults[ord,which(!names(GeneResults)%in%c("isPcG","isNearPcG"))]
		
		i=prId[grep("FAM5C",ann2$UCSC_RefGene_Name[prId])]
		i=i[grep("5'UTR",ann2$UCSC_RefGene_Group[i])]
		tbl[which(tbl$Gene=="FAM5C"),]
		stat2[i,]
	}
	
####################################################################
	
}

####################################################################
####################################################################
####################################################################
####################################################################

pThres=c(0.001,0.05)
iA2=match(stat2$cpgId,ann$IlmnID)

i=c()
for (compId in c("pcb105","pcb118","pcb138","pcb153","pcb170","pcb180","pcb1260")) {
    colIdPV="pv"
    cat("\n\n==================",compId,"==================\n")
    switch(compId,
    "pcb105"={
        stat2=stat_pcb105
        colIdPV=names(stat2)[grep("pv_",names(stat2))]
    },
    "pcb118"={
        stat2=stat_pcb118
        colIdPV=names(stat2)[grep("pv_",names(stat2))]
    },
    "pcb138"={
        stat2=stat_pcb138
        colIdPV=names(stat2)[grep("pv_",names(stat2))]
    },
    "pcb153"={
        stat2=stat_pcb153
        colIdPV=names(stat2)[grep("pv_",names(stat2))]
    },
    "pcb170"={
        stat2=stat_pcb170
        colIdPV=names(stat2)[grep("pv_",names(stat2))]
    },
    "pcb180"={
        stat2=stat_pcb180
        colIdPV=names(stat2)[grep("pv_",names(stat2))]
    },
    "pcb1260"={
        stat2=stat_pcb1260
        colIdPV=names(stat2)[grep("pv_",names(stat2))]
    }
    )
    colId=grep("pv_",names(stat2))
    i=c(i,which(stat2[,colId]<pThres[1]))
    colId=grep("qv_",names(stat2))
    i=c(i,which(stat2[,colId]<pThres[2]))
}
i=unique(i)
length(i)
cpgId=stat2$cpgId[i]
save(cpgId,file="cpgId_tmp.RData")

