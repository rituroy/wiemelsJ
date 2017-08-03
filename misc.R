dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))

##############################################
## Create list of probes in the same order as in the methylation data files

## Same order of probes in all files
datadir="docs/yuanyuan/Methylation/ALL-Guthrie/"
fName="guthrie_data_normtointernalprobes"
fName="guthrie_quantile_normalized"
fName="guthrie_quantile_arcsinesqrt"

datadir="docs/yuanyuan/Methylation/AML/"
fName="AML_methylation"

datadir="docs/all/LEU.data/"
fName="d.LEU"

datadir="docs/all/Guthrie.data/"
fName="d.GUTH"

blocksize=50000; N=485577
nblock=ceiling(N/blocksize)

print(nblock)
firstFlag=T
for(b in 1:nblock){
	print(b)
	nb=ifelse(b<nblock, blocksize, N%%blocksize)
	thisdata=read.delim(paste(datadir,fName,".txt",sep=""),header=FALSE,skip=blocksize*(b-1)+1, nrow=nb)
	probes=as.character(thisdata[,1])
	if (firstFlag) {
		write.table("probeId",file=paste("probeId.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)
		firstFlag=F
	}
	write.table(probes,file=paste("probeId.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=T)
}

f1=read.table("results/probeId.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
f2=read.table("probeId.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
table(f1$probeId==f2$probeId,exclude=NULL) ## All TRUE?
 
##############################################
## Create methylation data R object

fName="guthrie_quantile_normalized"
fName="guthrie_data_normtointernalprobes"
fName2=paste("docs/yuanyuan/Methylation/ALL-Guthrie/",fName,".txt",sep="")

fName="d.LEU"
fName2=paste("docs/all/LEU.data/",fName,".txt",sep="")

meth=read.delim(fName2,header=T)
probe=meth[,1]
meth=data.matrix(meth[,-1])
rownames(meth)=as.character(probe)
save(meth,file=paste(fName,".RData",sep=""))

##############################################
## Get methylation data summary

nProbe=100
nProbe=-1

cohortList=c("_allGuthSet1","_allGuthSet2","_allGuthSet1Set2_ctrlSubset")
cohortList=c("_allGuthSet1Set2_ctrlSubset")
tmp=rep(NA,length(cohortList))
out=data.frame(cohort=cohortList,minBeta=tmp,maxBeta=tmp,stringsAsFactors=F)
for (cohort in cohortList) {
    switch(cohort,
    "_allGuthSet1"={
        datadir="data/set1/"
        fName="beta_bmiq_allGuthSet1"
    },
    "_allGuthSet2"={
        datadir="data/set2/"
        fName="beta_bmiq_allGuthSet2"
    },
    "_allGuthSet1Set2_ctrlSubset"={
        datadir="data/set1set2/"
        fName="beta_bmiq_allGuthSet1Set2_ctrlSubset"
    }
    )
    meth=read.table(paste(datadir,fName,".txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=nProbe)
    rownames(meth)=meth$probeId
    meth=as.matrix(meth[,-1])
    save(meth,file=paste(fName,"_tmp.RData",sep=""))
    
    load(file=paste(fName,"_tmp.RData",sep=""))
    
    k=which(out$cohort==cohort)
    out$minBeta[k]=min(c(meth[which(meth!=0)]),na.rm=T)
    out$maxBeta[k]=max(c(meth[which(meth!=1)]),na.rm=T)
}
tbl=read.table(paste("data/summaryBeta.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=nProbe)
k=which(!out$cohort%in%tbl$cohort)
if (length(k)<nrow(out)) cat("Cohorts ",out$cohort[-k]," are already processed, so skipped!!!\n")
if (length(k)!=0) {
    tbl=rbind(tbl,out)
    write.table(tbl,file=paste("summaryBeta.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

##############################################
## Arcsine transform methylation data

fName="guthrie_quantile_normalized"
fNameOut="guthrie_quantile_arcsinesqrt"

fName="guthrie_data_normtointernalprobes"
fNameOut="guthrie_data_arcsinesqrt"

blocksize=50000; N=485577
nblock=ceiling(N/blocksize)

print(nblock)
firstFlag=T
for(b in 1:nblock){
	cat(b,":\n")
	nb=ifelse(b<nblock, blocksize, N%%blocksize)
	thisdata=read.delim(paste("docs/yuanyuan/Methylation/ALL-Guthrie/",fName,".txt",sep=""),header=FALSE,skip=blocksize*(b-1)+1, nrow=nb)
	probes=as.character(thisdata[,1])
	thisdata=data.matrix(thisdata[,-1])
	thisB=apply(thisdata,c(1,2),function(x) {asin(sqrt(x))})
	if (firstFlag) {
		write.table("probeId",file=paste("probeId.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)
		write.table(paste(c("TargetID",colnames(thisdata)),collapse="\t"),file=paste(fNameOut,".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)
		firstFlag=F
	}
	write.table(probes,file=paste("probeId.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=T)
	write.table(cbind(probes,as.data.frame(thisB)),file=paste(fNameOut,".txt",,sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=T)
}

########################################################################
## Cell mixture

datType="all"
datType="aml"

if (F) {
covdata=read.delim(paste("docs/yuanyuan/Samples/","ALLGuthrie-covdata-Feb232012.txt",sep=""),header=TRUE,sep="\t")
race=factor(covdata[,"int_ch_race"],levels=c(1:5))
ethn=factor(covdata[,"int_ch_ethnicity"],levels=c(1:3))
gestage=as.numeric(covdata[,"gestage"])
sex=factor(covdata[,"sex"])
covs=data.frame(sex=sex,race=race,ethn=ethn,gestage=gestage)

probes=read.table("results/probeId.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
ann=read.delim(paste("docs/yuanyuan/Annotation/","HumanMethylation450_15017482_v.1.1.csv",sep=""),header=TRUE, sep=",", skip=7)
ann=ann[match(probes$probeId,ann[,"IlmnID"]),]
rm(probes)
ann=as.matrix(ann)
#ann=ann[match(cgnames, ann[,1]),]
ann=ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq","Infinium_Design_Type", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
good=ann[,"CHR"]%in%c(1:22)
pos=as.integer(ann[,"MAPINFO"])

out=as.data.frame(t(sapply(as.character(ann[,"UCSC_RefGene_Group"]),function(x) {
	y=strsplit(x,";")[[1]]
	y=sort(unique(y))
	y2=y
	y2[which(y%in%c("TSS200","TSS1500"))]="TSS"
	y2=sort(unique(y2))
	c(paste(y,collapse=";"),length(y),paste(y2,collapse=";"),length(y2))
},USE.NAMES=F)),stringsAsFactors=F)
names(out)=c("group","numOfGroups","group2","numOfGroup2s")

ann2=cbind(pos,out)
}

datadir=paste("results/",datType,"/guthrie/norm/",sep="")
if (datType=="all") {
	covdata=read.delim(paste("docs/yuanyuan/Samples/","ALLGuthrie-covdata-Feb232012.txt",sep=""),header=TRUE,sep="\t")
	covdata$BeadChip=substr(covdata$Bead_Position,1,10)
	covdata$caco=covdata[,"Leukemia"]==1
	phen=covdata
	phen$caco=as.factor(phen$caco)
	phen$BeadChip=as.factor(phen$BeadChip)
	cellMix=read.table("results/all/guthrie/cellMixture/CellMixingCoefficients.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
	names(cellMix)=sub("X","",gsub(".","",names(cellMix),fixed=T))
} else {
	load(paste(datadir,"combatAdjustedData_aml.RData",sep=""))
	names(phen)[match("Beadchip",names(phen))]="BeadChip"
	rm(yMVals,annot)
	cellMix=read.table("results/aml/guthrie/cellMixture/CellMixCoef_aml.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
	j=match(phen$Subject_ID,cellMix$Subject_ID)
	print(table(is.na(j)))
	cellMix=cellMix[j,-1]
}


cellMix=as.matrix(cellMix)

cellMix1=cellMix

cellMix=cellMix1
for (k in 1:ncol(cellMix)) {
	cellMix[,k]=100*cellMix[,k]/sum(cellMix[,k])
}
grp=rep(colnames(cellMix),each=nrow(cellMix))

lim=c(0,85)

png(paste("cellMixBoxplot_",datType,".png",sep=""))
boxplot(c(cellMix1)~grp,ylim=lim,main=toupper(datType),ylab="Cell mixture coefficient")
dev.off()

if (F) {
png(paste("cellMixCoefAsSamplePercBoxplot_",datType,".png",sep=""))
boxplot(c(cellMix)~grp,ylim=lim,main=toupper(datType),ylab="Cell mixture coefficient (as sample %)")
dev.off()

#for (k in 1:ncol(cellMix)) {
#	boxplot(cellMix[,k],names=names(cellMix)[k],add=k>1,at=k,xlim=c(1,ncol(cellMix)))
#}
write.table(cellMix,file="CellMixingCoefAsSamplePercentage.txt", sep="\t", col.names=T, row.names=F, quote=F)
cellMix2=cellMix
}

cellMix=cellMix1
for (i in 1:nrow(cellMix)) {
	cellMix[i,]=100*cellMix[i,]/sum(cellMix[i,])
}
grp=rep(colnames(cellMix),each=nrow(cellMix))

png(paste("cellMixCoefSumTo100PercBoxplot_",datType,".png",sep=""))
boxplot(c(cellMix)~grp,ylim=lim,main=toupper(datType),ylab="Cell mixture coefficient (sum to 100%)")
dev.off()

grp2=apply(cbind(grp,rep(phen$caco,ncol(cellMix))),1,paste,collapse="_")
png(paste("cellMixCacoBoxplot_",datType,".png",sep=""))
boxplot(c(cellMix)~grp2,ylim=lim,main=toupper(datType),ylab="Cell mixture coefficient (sum to 100%)",las=2)
dev.off()

library(coin)
for (k in 1:ncol(cellMix)) {
	x=pvalue(wilcox_test(cellMix[,k]~as.factor(phen$caco),distribution="exact"))
	x2=summary(lm(cellMix[,k]~as.factor(phen$caco)))$coef[2,4]
#	cat(colnames(cellMix)[k],": pv ",signif(x,6),"\n",sep="")
	cat(colnames(cellMix)[k],": pv ",signif(x2,6),"\n",sep="")
}


AML: T-test
CD8T: pv 0.133021
CD4T: pv 0.160216
NK: pv 0.0137926
Bcell: pv 0.0656058
Mono: pv 0.0264283
Gran: pv 0.689565


AML: Wilcox
CD8T: pv 0.456419
CD4T: pv 0.174119
NK: pv 0.128589
Bcell: pv 0.0992841
Mono: pv 0.0483366
Gran: pv 0.884002


#for (k in 1:ncol(cellMix)) {
#	boxplot(cellMix[,k],names=names(cellMix)[k],add=k>1,at=k,xlim=c(1,ncol(cellMix)))
#}
write.table(cbind(Subject_ID=phen$Subject_ID,cellMix),file=paste("cellMixCoefSumTo100Perc_",datType,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
cellMix3=cellMix

########################################################################
## Create txt file from SAS file

dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))

library(sas7bdat)

## -----------------------------
tbl=read.sas7bdat("docs/aml/final.sas7bdat")
tbl$VAR7=sub("\r","",tbl$VAR7)
write.table(tbl, file="clinInfo_aml.txt", sep="\t", col.names=T, row.names=F, quote=F)

## -----------------------------
fNameClin=paste("clin_aml_20151002_2.txt",sep="")
clin=read.table(paste("docs/aml/",fNameClin,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

clin2=read.table(paste("Set1Set2_Guthrie_MethIDs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(clin2)=c("subjectId")

subjectId=c(clin2$subjectId,clin$subjectId[!clin$subjectId%in%clin2$subjectId])
tbl=data.frame(subjectId=subjectId,new85=rep("no",length(subjectId)),stringsAsFactors=F)
tbl$new85[!tbl$subjectId%in%clin2$subjectId]="yes"
write.table(tbl, file="Set1Set2_Guthrie_MethIDs_20151019.txt", sep="\t", col.names=T, row.names=F, quote=F)

## -----------------------------
for (datType in c("_allGuthSet1","_allGuthSet2","_allGuthSet1Set2")) {
    switch(datType,
    "_allGuthSet1"={
        datadirClin="set1/"
        fNameClin=paste("clin_allGuthSet1_20160523.txt",sep="")
        fNameOut="_20160928"
    },
    "_allGuthSet2"={
        datadirClin="set2/"
        fNameClin=paste("clin_allGuthSet2_20160523.txt",sep="")
        fNameOut="_20160928"
    },
    "_allGuthSet1Set2"={
        datadirClin="set1set2/"
        fNameClin=paste("clin_guthrieSet1Set2_20140619.txt",sep="")
        fNameOut="_20151022"
    }
    )
    clin=read.table(paste("docs/all/",datadirClin,fNameClin,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    clin1=read.table(paste("docs/all/set1/final.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(clin1)[match(c("Subject_ID","TargetID"),names(clin1))]=c("subjectId","guthrieId")
    clin2=read.table(paste("docs/all/set2/clin_guthrieSet2_20140619.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(clin2)[match(c("subjectID"),names(clin2))]=c("subjectId")
    clin$subjectId=""
    j=match(clin$guthrieId,clin1$guthrieId); j1=which(!is.na(j)); j2=j[j1]
    clin$subjectId[j1]=clin1$subjectId[j2]
    j=match(clin$guthrieId,clin2$guthrieId); j1=which(!is.na(j)); j2=j[j1]
    clin$subjectId[j1]=clin2$subjectId[j2]

    ##tbl=read.sas7bdat("docs/aml/jw_req_final.sas7bdat")
    #tbl=read.sas7bdat("docs/aml/jw_req_final_10192015.sas7bdat")
    tbl=read.sas7bdat("docs/all/set1set2/jw_req_final_10192015.sas7bdat")
    names(tbl)[match(c("ph","subjectID","income","motherht","motherwt","mothwtgn","ear_1yr","ear_1yr_freq","parity"),names(tbl))]=c("ph","subjectId","income","motherHt","motherWt","mothWtGn","ear1yr","ear1yrFreq","parity")
    #write.table(tbl, file="jw_req_final.txt", sep="\t", col.names=T, row.names=F, quote=F)
    #tbl=tbl[tbl$subjectId%in%tbl$subjectId,]
    for (k in 1:ncol(tbl)) {
        tbl[is.na(tbl[,k]),k]=NA
        print(names(tbl)[k])
        if (sum(!duplicated(tbl[,k]))<20) {
            print(table(tbl[,k],exclude=NULL))
        } else {
            print(summary(tbl[,k]))
        }
    }

    j=match(clin$subjectId,tbl$subjectId)
    if (any(is.na(j))) {
        tbl2=tbl[1:sum(is.na(j)),]
        for (k in 1:ncol(tbl2)) {
            tbl2[,k]=NA
        }
        tbl2$subjectId=clin$subjectId[is.na(j)]
        tbl=rbind(tbl,tbl2)
    }
    tbl=tbl[match(clin$subjectId,tbl$subjectId),,]
    summary(clin$income-tbl$income)

    tbl=cbind(clin,tbl[match(clin$subjectId,tbl$subjectId),!names(tbl)%in%names(clin)])
    write.table(tbl, file=paste("clin",datType,fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

## -----------------------------
fName="_20160523"

normFlag="_bmiq"

load("/Users/royr/UCSF/JoeWiemels/leukMeth/docs/SemiraGonsethNussle/ID_CS_Semira.RData")
names(sem)[match(c("TargetID","Subject_ID","bwt","GSage","birth"),names(sem))]=c("guthrieId","subjectId","birthWt","gestage","birth")
for (k in 1:ncol(sem)) {
    if (class(sem[,k])=="factor") sem[,k]=as.character(sem[,k])
}
sem[sem[,1]%in%sem[,1][duplicated(sem[,1])],]


for (datType in c("_allGuthSet1","_allGuthSet2","_allGuthSet1Set2")) {
    cat("\n\n===================",datType,"\n")
    switch(datType,
    "_allGuthSet2"={
        dirMeth=dirClin="docs/all/set2/"
        fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
        fNameClin="clin_guthrieSet2_20140619"
    },
    "_allGuthSet1"={
        dirMeth=dirClin="docs/all/set1/"
        fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),sep="")
        fNameClin="final"
        dirMethLeuk="docs/all/set1/LEU.data/"
        fNameMethLeuk="i.LEU.v2"
    },
    "_allGuthSet1Set2"={
        dirMeth=dirClin="docs/all/set1set2/"
        fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
        fNameClin="clin_guthrieSet1Set2_20151022"
    }
    )
    clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    k=which(names(clin)%in%c("TargetID","guthrieId"))
    names(clin)[k]="guthrieId"
    if ("Subject_ID"%in%names(clin)) clin=clin[,which(!names(clin)%in%c("Subject_ID"))]
    k=which(names(clin)%in%c("SubjectID","subjectID","subjectId"))
    names(clin)[k]="subjectId"
    k=which(names(clin)%in%c("birth_weight","birth_wt"))
    if (length(k)==1) {
        names(clin)[k]="birthWt"
        clin$birthWt=as.numeric(clin$birthWt)
    }
    #print(sort(names(clin)))
    #print(names(clin)[1:5])
    print(clin[clin$subjectId%in%clin$subjectId[duplicated(clin$subjectId)],1:5])
    j=match(clin$subjectId,sem$subjectId); j1=which(!is.na(j)); j2=j[j1]; j12=which(is.na(j))
    tbl=sem[j2,]
    if (length(j12)!=0) {
        tbl2=tbl[1:length(j12),]
        for (k in 1:ncol(tbl2)) {
            tbl2[,k]=NA
        }
        tbl2$subjectId=clin$subjectId[j12]
        tbl=rbind(tbl,tbl2)
    }
    k=which(!names(tbl)%in%names(clin))
    if (length(k)==0) {
        tbl=clin
    } else {
        nm=c(names(clin),names(tbl)[k])
        tbl=cbind(clin,tbl[match(clin$subjectId,tbl$subjectId),k])
        names(tbl)=nm
    }
    if (F) {
        ## Compare clin $ sem variables
        varList=data.frame(var1=c("gestage","birthWt"),var2=c("GSage","bwt"),stringsAsFactors=F)
        for (k in 1:nrow(varList)) {
            if (varList$var1[k]%in%names(tbl)) {
                print(varList[k,])
                if (varList$var1[k]=="birthWt") {
                    print(table(round(tbl[,varList$var1[k]])==tbl[,varList$var2[k]]))
                    j=which(round(tbl[,varList$var1[k]])!=tbl[,varList$var2[k]])
                    print(tbl[j,c("guthrieId",varList$var1[k],varList$var2[k])])
                } else {
                    print(table(tbl[,varList$var1[k]]==tbl[,varList$var2[k]]))
                }
                print(table(is.na(tbl[,varList$var1[k]])==is.na(tbl[,varList$var2[k]])))
            }
        }
    }
    write.table(tbl, file=paste("clin",datType,fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}
"
## Compare clin & sem variables:
## All gestages match
## Birth weights which do not match are like xxxx.5 & xxxx+1
    guthrieId birthWt  bwt
211     1979G  3118.5 3119
259     1713G  4252.5 4253
"

##############################################
## Txt file from RData

datType="aml"
datType="all"

datadir3=paste("results/",datType,"/guthrie/norm/",sep="")
datadir="docs/all/Guthrie.data/"
datadir2=datadir
fName="d.GUTH"
fName2="i.GUTH"

probes1=read.table("results/probeId.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
covdata=read.delim(paste(datadir2,fName2,".txt",sep=""),header=TRUE,sep="\t")
load("results/all/guthrie/norm/lmeAdjustedMValues_all_chipPos.RData")
j=match(colnames(yAdjM),paste("X",covdata$TargetID,sep=""))
j1=which(!is.na(j)); j2=j[j1]
yAdjM=yAdjM[,j1]
covdata=covdata[j2,]

tbl=cbind(IlmnID=probes,yAdjM)
write.table(cbind(IlmnID=probes,yAdjM), file="lmeAdjustedMValues_all_chipPos.txt", sep="\t", col.names=T, row.names=F, quote=F)

save(yAdjM,probes,file="mValues_chipPosLmeAdjusted_all.RData")

########################################################################
########################################################################
########################################################################
## Cluster clinical data to get meta variable for smoking

datType="allGuth"
datType="_allGuthSet2"
datType="_allGuthSet1"

for (datType in c("_allGuthSet2","_allGuthSet1")) {


if (datType=="allGuth") {
	datadir2="docs/all/"
	fInName2="final"

	varFlag="_leukBatchChipPos"
	datadir3=paste("results/all/guthrie/norm/",sep="")
	load(paste(datadir3,"wbcLmeAdjustedMDat_",datType,varFlag,".RData",sep=""))
	thisdata=mDat
	rm(mDat)

	clin=read.table(paste(datadir2,fInName2,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	clin$id=paste("X",clin$TargetID,sep="")
	j=match(colnames(thisdata),clin$id)
	min(diff(j)) ## >0 & not NA
	clin=clin[j,]
} else {
	switch(datType,
		   "_allGuthSet2"={
		   datadir=datadir2="docs/all/set2/"
		   fName="beta_funNorm_set2"
		   fName2="clin_guthrieSet2_20140619"
		   },
		   "_allGuthSet1"={
		   datadir=datadir2="docs/all/set1/"
		   fName="beta_funNorm_set1"
		   fName2="final"
		   datadir22="docs/all/set1/LEU.data/"
		   fName22="i.LEU.v2"
		   }
		   )

	clin=read.table(paste("docs/all/set2/clin_guthrieSet2_20140619.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	varSmoke=data.frame(varIn=c("passive_sm_home_post","smoke_mo_ever","smoke_mo_preg","M_SM_BF_DI","smoke_mo_3months_N","smoke_mo_preg_N","m_cignum_bf","smoke_mo_after_N","smoke_fa_ever","smoke_fa_3months","smoke_fa_3months_N","smoke_mo_bf_N"),
						varOut=c("paSmHmPost","moEver","moPreg","moBf","mo3mN",
							  "moPregN","moBfN","moAfterN",
							  "faEver","fa3m","fa3mN","moBf"),stringsAsFactors=F)
	k=match(varSmoke$varIn,names(clin)); k1=which(!is.na(k)); k2=k[k1]
	clin3=clin[,k2]
	
	switch(datType,
		"_allGuthSet2"={
		   clin=read.table(paste(datadir2,fName2,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		   clin$id=paste("X",clin$guthrieId,sep="")
		   clin$Leukemia=clin$caco
		   
		   clin$subtype=rep("",nrow(clin))
		   clin$subtype[which(clin$smhyper==1)]="hyperdiploid"
		   clin$subtype[which(clin$smtelaml==1)]="telaml" ## include subject with smhyper=1
		   clin$subtype[which(clin$smhyper==0 & clin$smtelaml==0)]="nonHypTelaml"
		   
		   beta=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
		   probeId=beta$probeId
		   beta=as.matrix(beta[,-1])
		   rownames(beta)=probeId
		   yMVals=beta
		   rm(probeId,beta)
		   
		   phen=clin
		   phen=phen[match(colnames(yMVals),phen$id),]
		},
		"_allGuthSet1"={
		   clin=read.table(paste(datadir2,fName2,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		   clin$id=paste("X",clin$TargetID,sep="")
		   clin$Beadchip=substr(clin$Bead_Position,1,10)
		   if ("Position1"%in%names(clin)) clin$Position=clin$Position1 else clin$Position=clin$Position.1
		   clin$caco=clin$Leukemia
		   clin$int_ch_ethnicity=clin$ch_ethnicity
		   
		   clin2=read.table(paste(datadir22,fName22,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		   j=match(clin$Subject_ID,clin2$Subject_ID)
		   j1=which(!is.na(j)); j2=j[j1]
		   clin$subtype=rep("",nrow(clin))
		   clin$subtype[j1][which(clin2$Subtype[j2]=="hyperdiploid")]="hyperdiploid"
		   clin$subtype[j1][which(clin2$Subtype[j2]=="t1221")]="telaml"
		   clin$subtype[j1][which(clin2$Subtype[j2]%in%c("mll","others","t119"))]="nonHypTelaml"
		   
		   beta=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
		   probeId=beta$probeId
		   beta=as.matrix(beta[,-1])
		   rownames(beta)=probeId
		   yMVals=beta
		   rm(probeId,beta)
		}
	)
}

fName1="_smokingVars"

varSetName=""
varSetName=c("","_noAfter","_noAfter_noNum")
varSetId=1

for (varSetId in 1:length(varSetName)) {
	
	fName=paste(fName1,varSetName[varSetId],datType,sep="")

	if (datType=="allGuth") {
		varPred=data.frame("paSmHmPost"=factor(clin[,"passive_sm_home_post"]), #96
		"moEver"=factor(clin[,"smoke_mo_ever"]), #150
		"moPreg"=factor(clin[,"smoke_mo_preg"]), #41
		"moBf"=factor(clin[,"M_SM_BF_DI"]), #14
		"mo3m"=factor(as.integer(varPred$mo3mN>0)),
		"mo3mN"=as.numeric(clin[,"smoke_mo_3months_N"]),
		"moPregN"=as.numeric(clin[,"smoke_mo_preg_N"]),
		"mo3mPregN"=as.numeric(clin[,"smoke_mo_3months_N"])+as.numeric(clin[,"smoke_mo_preg_N"]),
		"moBfN"=as.numeric(clin[,"m_cignum_bf"]), #14
		"moAfterN"=as.numeric(clin[,"smoke_mo_after_N"]), #62
		"faEver"=factor(clin[,"smoke_fa_ever"]),#213
		"fa3m"=factor(clin[,"smoke_fa_3months"]),#117
		"fa3mN"=as.numeric(clin[,"smoke_fa_3months_N"])) #107
	} else {
		x=names(clin)
		x=x[c(which(!x%in%varSmoke$varIn),which(x%in%varSmoke$varIn & x%in%names(clin3)))]
		clin=clin[,which(names(clin)%in%x)]
		k=match(varSmoke$varIn,names(clin)); k1=which(!is.na(k)); k2=k[k1]
		varPred=clin[,k2]
		names(varPred)=varSmoke$varOut[k1]
		varPred$mo3m=as.integer(varPred$mo3mN>0)
		varPred$mo3mPregN=varPred$mo3mN+varPred$moPregN
	}

	for (k in 1:ncol(varPred)) {
		varPred[,k]=as.integer(as.character(varPred[,k]))
	}

		if (varSetName[varSetId]=="_noAfter") {
			varPred=varPred[,which(names(varPred)%in%c("moPreg","mo3m","mo3mN","moPregN","mo3mPregN","fa3m","fa3mN"))]
		} else if (varSetName[varSetId]=="_noAfter_noNum") {
			varPred=varPred[,which(names(varPred)%in%c("moPreg","mo3m","fa3m"))]
		}

	#varPred=t(as.matrix(varPred))
	#plot(hclust(dist(t(varPred)),method="ward"))


	#######################################################
	library(marray)
	#source(paste(dirSrc,"functions/heatmap.4.3.R",sep=""))
	#source(paste(dirSrc,"functions/heatmapAcgh.6.R",sep=""))
	source(paste(dirSrc,"functions/heatmap.5.R",sep=""))
	source(paste(dirSrc,"functions/heatmapAcgh.7.R",sep=""))

	colHmap=c("blue", "red", "yellow")
	colList=c("skyblue","blue","yellow","purple")
	varList=colnames(varPred)
	varName=paste(varList," ",sep="")

	prId=1:ncol(varPred)
	samId=1:nrow(varPred)
	samId=which(apply(as.matrix(varPred[samId,prId]),1,function(x) mean(is.na(x)))<1)

	distMethod="pearson"
	distMethod="euclidean"
	#linkMethod="ward"
	linkMethod="ward.D2"

	phen2=cbind(clin[,c("id","Leukemia")],varPred)[samId,prId]
	arrayData=t(as.matrix(varPred[samId,prId]))

	cloneName=rep("",length(prId))
	cloneName=rownames(arrayData)
	#cloneCol=matrix(rep("white",nrow(arrayData)),nrow=1)
	cloneCol=NULL
	samName=clin$id[samId]
	samName=rep("",length(samId))


	lineSam=NULL
	phenList=colnames(varPred)
	phenName=paste(phenList," ",sep="")
	phenName2=phenName
	k=which(phenList%in%names(varPred))
	phenList=phenList[k]
	phenName=phenName[k]
	phenName2=phenName2[k]
	nm=cutoff=c()
	
	if (F) {
		samCol=matrix(rep(NA,length(samName)*length(phenName)),ncol=length(samName))
		lineSam=NULL
		k1=0
		for (phenId in 1:length(phenList)) {
			k1=k1+1
			kk=which(names(phen2)==phenList[phenId])
			x=phen2[,kk]
			grpUniq=sort(unique(x[!is.na(x)]))
			if (phenList[phenId]%in%c("type")) {
				samColUniq2=samColUniq
			} else if (phenList[phenId]%in%c("site")) {
				samColUniq2=rainbow(length(grpUniq))
			} else {
				samColUniq2=rainbow(length(grpUniq))
				samColUniq2=gray(0:(length(grpUniq)-1)/length(grpUniq))
			}
			for (k in 1:length(grpUniq)) {samCol[k1,which(phen2[,kk]==grpUniq[k])]=samColUniq2[k]}	
		}
		rownames(samCol)=phenName
	}
	samCol=NULL

	switch(distMethod,
	"pearson" = {distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
		clustC=hclust(distMat, method=linkMethod)
		distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
		clustR=hclust(distMat, method=linkMethod)},
	"spearman" = {distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
		distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))},
	"euclidean" = {distMat=dist(t(arrayData), method=distMethod)
		clustC=hclust(distMat, method=linkMethod)
		distMat=dist(arrayData, method=distMethod)
		clustR=hclust(distMat, method=linkMethod)}
	)

	summary(range(c(arrayData),na.rm=T))
	limit=c(-40,40)
	limit=c(-54,54)
	dat2=arrayData
	dat2[dat2==0]=limit[1]
	margins=c(.2,6)
	main=NULL

	nClust=3

	#png(paste("heatmap",fName,".png",sep=""),width=480*2,height=480*2); cexRow=cexCol=2
	pdf(paste("heatmap",fName,".pdf",sep="")); cexRow=cexCol=1
	#hcc=heatmap3(x=dat2, Rowv=as.dendrogram(clustR), Colv=as.dendrogram(clustC), distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit, cexCol=2, cexRow=2, high=colHmap[1], low=colHmap[2], mid=colHmap[3])
	hcc=heatmap3(x=dat2, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName, ncr=NA, ncc=nClust, scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit, cexCol=cexCol, cexRow=cexRow, high=colHmap[1], low=colHmap[2], mid=colHmap[3])
	dev.off()
	if (!is.null(samCol)) {
		for (varId in 1:length(varList)) {
			png(paste("heatmapSampleColorBarLegend_",varList[varId],".png",sep=""))
			x=as.character(pData(eset)[,varList[varId]]); x[x==""]=NA
			sampleColorLegend(tls=sort(unique(x)),col=colList,legendTitle=varName[varId])
			dev.off()
		}
	}

	if (!is.null(nClust)) {
		if (F) {
	#	png(paste("clusterSamples",fName,".png",sep=""))
		pdf(paste("clusterSamples",fName,".pdf",sep=""))
		plot(clustC,main=paste("Sample clusters with ",nClust," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust)
		dev.off()
		}
		
		clustId=cutree(clustC,k=nClust)[clustC$order]
		k1=which(!duplicated(clustId))
		for (k in 1:length(k1)) {
			clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
		}
		
		tbl=cbind(phen2[clustC$order,],clustId,sampleOrder=1:nrow(phen2))
		write.table(tbl, paste("clusterInfo",fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
	}
}
}
png("heatmapColorRange.png")
#heatmapColorBar(limit=limit,cols=colHmap)
maColorBar(c(0,1:limit[2]),col=maPalette(high=colHmap[1], low=colHmap[2], mid=colHmap[3],k=limit[2]+1)[c(1,(ceiling(limit[2]/2)+1):limit[2])])
dev.off()


## Smoking categories based on heatmaps
x=varPred
x$sm=NA
x$sm[which(x$moPreg==0 & x$mo3m==0 & x$fa3m==0)]=0
x$sm[which(x$moPreg==1)]=1
x$sm[which(x$mo3m==1 & x$moPreg==0)]=2
x$sm[which(x$fa3m==1 & x$mo3m==0 & x$moPreg==0)]=3

table(moPreg=varPred$moPreg,mo3m=varPred$mo3m,fa3m=varPred$fa3m,sm=x$sm,exclude=NULL)

########################################################################
########################################################################
########################################################################
## Check sample sizes of covariate categories

race=factor(clin[,"int_ch_race"],levels=c(1:5))
#ethn=factor(clin[,"int_ch_ethnicity"],levels=c(1:3))
ethn=factor(clin[,"ch_ethnicity"],levels=c(1:3))
gestage=as.numeric(clin[,"gestage"])
sex=factor(clin[,"sex"])
covs=data.frame(sex=sex,race=race,ethn=ethn,gestage=gestage)
covs$race4=as.integer(as.character(covs$race))
covs$race4[which(covs$race4==3)]=5
covs$race4=factor(covs$race4)
covs$race3=as.integer(as.character(covs$race))
covs$race3[which(covs$race3%in%c(2,3))]=5
covs$race3=factor(covs$race3)
table(covs$race,covs$ethn)
for (k in 1:ncol(covs)) {
	cat("\n\n================",names(covs)[k],"\n")
	print(table(covs[,k]))
	if (class(covs[,k])=="factor") {
		cat("pv ",fisher.test(covs[,k],clin$Leukemia)$p.value)
	} else {
		cat("pv ",wilcox.test(covs[,k]~clin$Leukemia)$p.value)
	}
}

########################################################################
########################################################################
########################################################################
## Birth weight

library(coin)

datadir2="docs/all/"
fInName2="final"
clin=read.table(paste(datadir2,fInName2,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clin$id=paste("X",clin$TargetID,sep="")
names(clin)[match(c("birth_weight"),names(clin))]=c("birthWt")
datadir22="docs/all/LEU.data/"
fName22="i.LEU.v2"
clin2=read.delim(paste(datadir22,fName22,".txt",sep=""),header=TRUE,sep="\t")
j=match(clin$Subject_ID,clin2$Subject_ID)
j1=which(!is.na(j)); j2=j[j1]
clin$Subtype=rep(NA,nrow(clin))
clin$Subtype[j1]=clin2$Subtype[j2]
#clin=clin[j1,]
#clin2=clin2[j2,]
#clin$Subtype=clin2$Subtype
#rm(datadir22,fName22,clin2)

clin2=read.table(paste(datadir2,"birthWeight/Ritu_birthweight_dataset - covariate data.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(clin2)[match(c("birth_wt"),names(clin2))]=c("birthWt")

j=match(clin$id,clin2$targetID)
table(is.na(j))
clin=clin[!is.na(j),]
clin2=clin2[j[!is.na(j)],]

tmp=rep(NA,nrow(clin))
varPred=data.frame(junk=tmp,hyperdip=tmp)
varPred$hyperdip=NA
varPred$hyperdip[which(clin$Subtype%in%c("mll","others","t119","t1221"))]=0
varPred$hyperdip[which(clin$Subtype=="hyperdiploid")]=1

race=factor(clin[,"int_ch_race"],levels=c(1:5))
#ethn=factor(clin[,"int_ch_ethnicity"],levels=c(1:3))
ethn=factor(clin[,"ch_ethnicity"],levels=c(1:3))
gestage=as.numeric(clin[,"gestage"])
sex=factor(clin[,"sex"])
covs=data.frame(sex=sex,race=race,ethn=ethn,gestage=gestage)
covs$race3=as.integer(as.character(covs$race))
covs$race3[which(covs$race3%in%c(2,3))]=5
covs$race3=factor(covs$race3)


j=which(clin[,"Leukemia"]==1)
boxplot(clin2$Percent[j]~varPred$hyperdip[j],main=signif(pvalue(wilcox_test(clin2$Percent[j]~as.factor(varPred$hyperdip[j]))),2))
boxplot(clin$gestage[j]~varPred$hyperdip[j],main=signif(pvalue(wilcox_test(clin$gestage[j]~as.factor(varPred$hyperdip[j]))),2))
fisher.test(varPred$hyperdip[j],clin$sex[j])
fisher.test(varPred$hyperdip[j],covs$ethn[j])
## pv=0.05258
fisher.test(varPred$hyperdip[j],covs$race[j])
fisher.test(varPred$hyperdip[j],covs$race3[j])

cor(clin2$Percent,clin$gestage,use="complete.obs",method="kendall")
boxplot(clin2$Percent~clin$sex,main=signif(pvalue(wilcox_test(clin2$Percent~as.factor(clin$sex))),2))

fisher.test(clin$Leukemia,covs$ethn)
## pv=0.002988
fisher.test(clin$Leukemia,covs$race)
boxplot(clin2$Percent~clin$Leukemia,main=signif(pvalue(wilcox_test(clin2$Percent~as.factor(clin$Leukemia))),2))


max(abs(clin$birthWt-clin2$birthWt),na.rm=T)/min(clin$birthWt,na.rm=T)
max(abs(clin$gestage-clin2$gestage),na.rm=T)/min(clin$gestage,na.rm=T)

table(clin2$gestage)
#26  28  32  33  34  35  36  37  38  39  40  41  42  43  44  45 
#1   1   4   4   3  10  13  31  60  79 164  60  22   2   2   3 
grp=table(clin2$gestage); grp=names(grp)[grp>9]; grp=as.integer(grp)
pdf("tmp.pdf")
for (k in 1:length(grp)) {
	par(mfrow=c(2,2))
	j=which(clin2$gestage==grp[k])
	pv=ks.test(clin2$birthWt[j], "pnorm", mean(clin2$birthWt[j],na.rm=T), sd(clin2$birthWt[j],na.rm=T))$p.value
	heading=paste("Gestage ",grp[k],": pv ",signif(pv,2),sep="")
	qqnorm(clin2$birthWt[j],main=heading); qqline(clin2$birthWt[j])
	plot(density(clin2$birthWt[j],na.rm=T),main=heading)
}
dev.off()

########################################################################
########################################################################
########################################################################
## Preprocess replication data

datadir="docs/all/set2/"
fIn="450K_replication_data_results_all columns"
fOut="methGuthrieRepl"

blocksize=10; N=485577
nblock=2

blocksize=50000; N=485577
nblock=ceiling(N/blocksize)

clin1=read.table(paste(datadir,"Final_dataset_june2012_replication methylation.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
clin2=read.table(paste(datadir,"I226_BeadChip_Sample_Layout.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
clin3=read.table(paste(datadir,"ID_labo_replication_CCLS_5Feb2014.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
clin4=read.table(paste(datadir,"covariates_052714.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
clin5=read.table(paste(datadir,"wiemels_income_061914.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

k=match(names(clin4),names(clin5)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
	if (any(clin4[,k1[k]]!=clin5[,k2[k]],na.rm=T) | any(is.na(clin4[,k1[k]])!=is.na(clin5[,k2[k]]))) {
		print(k)
	}
}

dat=clin4
dat=clin5
id=unique(dat$subjectID[duplicated(dat$subjectID)])
for (k in 1:ncol(dat)) {
for (jj in 1:length(id)) {
	j=which(dat$subjectID==id[jj])
	if (sum(!duplicated(dat[j,k]))!=1) {
		print(names(dat)[k])
		print(id[jj])
	}
}
}

names(clin1)[match(c("Guthrie.ID"),names(clin1))]=c("guthrieId")
names(clin5)[match(c("Guthrie.ID"),names(clin5))]=c("guthrieId")
names(clin2)[match(c("Sample.ID"),names(clin2))]=c("labId")
names(clin3)[match(c("labo_ID","GuthrieID"),names(clin3))]=c("labId","guthrieId")

clin5=clin5[match(clin1$guthrieId,clin5$guthrieId),]

for (datType in c("clin1","clin5")) {
	switch(datType,
		   "clin1"={dat=clin1},
		   "clin5"={dat=clin5}
		   )
	for (k in 1:ncol(dat)) {
		if (class(dat[,k])=="character") {
			dat[,k]=gsub("\"","",dat[,k])
		}
		if (names(dat)[k]%in%c("other_trans")) {
			dat[,k]=gsub(" ","",dat[,k])
		}
	}
	switch(datType,
		   "clin1"={clin1=dat},
		   "clin5"={clin5=dat}
		   )
}

x=as.integer(clin5$newid)
i=which(!is.na(x))
y=formatC(x[i],width=5,flag="0")
clin5$newid[i]=y

k=match(names(clin1),names(clin5)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
	if (any(clin1[,k1[k]]!=clin5[,k2[k]],na.rm=T) | any(is.na(clin1[,k1[k]])!=is.na(clin5[,k2[k]]))) {
		print(k)
		print(names(clin1)[k])
	}
}

clin1=cbind(clin1[,which(!names(clin1)%in%names(clin5))],clin5)
names(clin1)=sub("_new","",names(clin1))

clin2$arrayId=paste("X",clin2$Beadchip,"_",clin2$Position,sep="")

colNames=read.table(paste(datadir,fIn,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
colNames=names(colNames)
colIdBeta=grep("AVG_Beta",colNames)
colIdDetPV=grep("Detection.Pval",colNames)
arrayId=sub(".AVG_Beta","",colNames[colIdBeta],fixed=T)
table(arrayId%in%clin2$arrayId) ## All TRUE
clin2=clin2[match(arrayId,clin2$arrayId),]
table(clin3$labId%in%clin2$labId) ## All TRUE
clin2$guthrieId=NA
j=match(clin3$labId,clin2$labId); j3=which(!is.na(j)); j2=j[j3]
clin2$guthrieId[j2]=clin3$guthrieId[j3]
clin2[!clin2$guthrieId%in%clin1$guthrieId,]
#Plate Well  labId   Beadchip Position            arrayId guthrieId
#455     5  F11 BM1660 9702496164   R05C02 X9702496164_R05C02      <NA>
#456     5  F12 BM1683 9702496164   R06C02 X9702496164_R06C02      <NA>
j=match(clin2$guthrieId,clin1$guthrieId); j22=which(is.na(j)); j2=which(!is.na(j)); j1=j[j2]
colId2=c("Plate","Well","labId","Beadchip","Position")
clin11=cbind(clin1[j1,],clin2[j2,colId2])
clin12=clin1[1:length(j22),]
for (k in 1:ncol(clin12)) clin12[,k]=NA
clin12=cbind(clin12,clin2[j22,colId2])
clin=rbind(clin11,clin12)
clin=clin[match(arrayId,paste("X",clin$Beadchip,"_",clin$Position,sep="")),]
clin$caco=2-clin$caco
clin$int_ch_race=as.integer(clin$int_ch_race)

clin=clin[,which(!names(clin)%in%"hyperdip")]  ## Discrepancy is hyperdip status

#write.table(clin,file="clinGuthrieReplJune2012.txt", sep="\t", col.names=T, row.names=F, quote=F)
write.table(clin,file="clin_guthrieSet2_20140619.txt", sep="\t", col.names=T, row.names=F, quote=F)

samId=paste("X",clin$guthrieId,sep="")
j=which(is.na(clin$guthrieId))
samId[j]=paste("X",clin2$Beadchip,"_",clin2$Position,sep="")[j]

fName=paste(fOut,"_tmp",sep="")

write.table(paste(c("cpgId",samId),collapse="\t"),file=paste(fName,".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)
write.table("probeId",file=paste("probeId.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)
write.table(paste(c("cpgId",samId),collapse="\t"),file=paste(fOut,"_detP.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=F)

print(nblock)
for(b in 1:nblock) {
	print(b)
	nb=ifelse(b<nblock, blocksize, N%%blocksize)
	thisdata=read.table(paste(datadir,fIn,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,skip=blocksize*(b-1)+1, nrow=nb)
	dat=as.matrix(thisdata[,colIdBeta])
	dat2=as.matrix(thisdata[,colIdDetPV])
	write.table(cbind(thisdata[,2],dat2),file=paste(fOut,"_detP.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=T)
	dat[dat2>0.0001]=NA
	dat[apply(dat,1,function(x) mean(is.na(x)))>0.15,]=NA
	dat=cbind(thisdata[,2],dat)
	write.table(dat,file=paste(fName,".txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=T)
	probes=as.character(thisdata[,2])
	write.table(probes,file=paste("probeId.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F,append=T)
}

thisdata=read.table(paste(fName,".txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
for (j in 2:ncol(thisdata)) {
	thisdata[sapply(thisdata[,j],function(x) mean(is.na(x)),USE.NAMES=F)>0.15,j]=NA
}
write.table(thisdata,file=paste(fOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

if (F) {
samId=colnames(f1)[-1]
f3=as.matrix(f3[,-ncol(f3)])
f3=cbind(cpgId=f1[,1],f3)
write.table(f3,file=paste(fOut,"_detP.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)


f01=read.table("docs/all/replication/450K_replication_data_results_all columns.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=100)
f11=read.table("methGuthrieRepl.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=100)
f21=read.table("methGuthrieRepl_tmp.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=100)
f31=read.table("methGuthrieRepl_detP.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=100)
f11=as.matrix(f11[,-1])
f21=as.matrix(f21[,-1])
f31=as.matrix(f31[,-1])
	
	f01[,colIdBeta][1:5,1:5]
}


f3=read.table("methGuthrieRepl_detP.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
f3=as.matrix(f3[,-1])
f1=read.table("docs/all/Guthrie.data/d.GUTH.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
f1=as.matrix(f1[,-1])
x1=apply(f1,1,function(x) mean(!is.na(x)))
x3=apply(f3,1,function(x) mean(x<=0.0001,na.rm=T))
png("hist_loci.png")
par(mfrow=c(2,1))
hist(x1,main="Set1: Histogram of loci",xlab="proportion of samples with non-missing data")
hist(x3,main="Set2: Histogram of loci",xlab="proportion of samples with detection p-value >= 0.0001")
dev.off()

set.seed(645645)
samId=sample(1:ncol(f3),8,replace=F)
png("hist_sample_%02d.png")
par(mfrow=c(2,2))
for (j in samId) {
	hist(f3[,j],xlim=c(0,0.05),main=paste("Histogram: ",colnames(f3)[j],sep=""),xlab="Detection p-value")
}
dev.off()

f1=read.table("results/probeId.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
f2=read.table("probeId.txt", sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
table(f1$probeId==f2$probeId,exclude=NULL) ## All TRUE?

########################################################################
########################################################################
########################################################################
## Create AML clinical data table

dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))

datadir="docs/aml/"

clin1=read.table(paste(datadir,"11072011 Sample_Sheet (AML).csv",sep=""), sep=",", h=T, quote="", comment.char="",as.is=T,fill=T,skip=14)
clin2=read.table(paste(datadir,"i.AML.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
clin3=read.table(paste(datadir,"clinInfo_aml.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

names(clin1)[match(c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position"),names(clin1))]=c("guthrieId","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Beadchip","Position")
names(clin2)[match(c("CpG","CaseControl","Sex","Subject_ID"),names(clin2))]=c("guthrieId","caco","sex","subjectId")
names(clin3)[match(c("Subject_ID","CD8T","CD4T","NK","Bcell","Mono","VAR7","sex","ch_ethnicity","gestage"),names(clin3))]=c("subjectId","cd8t","cd4t","nk","bcell","mono","var7","sex","ch_ethnicity","gestage")

j=match(clin1$guthrieId,clin2$guthrieId); j1=which(!is.na(j)); j2=j[j1]
clin1=clin1[j1,]
clin2=clin2[j2,]

j=match(clin2$subjectId,clin3$subjectId); j1=which(!is.na(j)); j2=j[j1]
clin1=clin1[j1,]
clin2=clin2[j1,]
clin3=clin3[j2,]

clin11=clin1; clin12=clin2
for (k in which(names(clin11)%in%names(clin12))) {
	varName=names(clin11)[k]
	if (any(clin11[,varName]!=clin12[,varName],na.rm=T) | any(is.na(clin11[,varName])!=is.na(clin12[,varName]))) print(varName)
}
clin11=clin1; clin12=clin3
for (k in which(names(clin11)%in%names(clin12))) {
	varName=names(clin11)[k]
	if (any(clin11[,varName]!=clin12[,varName],na.rm=T) | any(is.na(clin11[,varName])!=is.na(clin12[,varName]))) print(varName)
}
clin11=clin2; clin12=clin3
for (k in which(names(clin11)%in%names(clin12))) {
	varName=names(clin11)[k]
	if (any(clin11[,varName]!=clin12[,varName],na.rm=T) | any(is.na(clin11[,varName])!=is.na(clin12[,varName]))) print(varName)
}

clin=cbind(clin1[,c("guthrieId","Beadchip","Position")],clin2[,c("caco","subjectId")],clin3[,c("cd8t","cd4t","nk","bcell","mono","var7","sex","ch_ethnicity","gestage")])
clin$caco=as.integer(clin$caco==1)

write.table(clin,file="clin_aml_20150114.txt", sep="\t", col.names=T, row.names=F, quote=F)

########################################################################
########################################################################
########################################################################
## Create birth defect clinical data table

dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))

datadir="docs/20150924_450K_Data/"
clin1=read.table(paste(datadir,"Cleft NTD 450K Samples 060716 - samples.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
names(clin1)[match(c("CBDMP.ID","Plate.ID","Source.Study","Status","Race.Ethnic","Sex","gDNA..ng.","Plate","Row","Column","Chip","Conversion.Batch","Conversion.Date..2015.","Duplicate.Replicate.Wells","Sentrix.ID","Sentrix.Posn"),names(clin1))]=c("cbdmpId","plateId","sourceStudy","status","raceEthn","sex","gDNA","plate","row","column","chip","conversionBatch","conversionDate2015","replicateWell","Beadchip","Position")
id=sapply(clin1$plateId,function(x) {
    y=gsub(" |-","_",x)
    if (!is.na(as.integer(substr(y,1,1)))) y=paste("X",y,sep="")
    y
},USE.NAMES=F)
clin=cbind(id,clin1)
clin$status=tolower(sub(" case","",clin$status))
write.table(clin,file="clin_birthDefect_20160607.txt", sep="\t", col.names=T, row.names=F, quote=F)

########################################################################
########################################################################
########################################################################
## Create data tables for Ivorra cohort

datadir="docs/misc/ivorra2014/"

## ------------------------------------

fName=paste(datadir,"GSE64316_series_matrix.txt",sep="")


clin1=read.table(paste(fName,sep=""), sep="\n", h=F, quote="", comment.char="",as.is=T,fill=T,nrow=400)
clin1=cbind(clin1,junk=rep("",nrow(clin1)))
clin1=clin1[which(gsub("\t","",clin1[,1])!=""),]
offset=grep("!series_matrix_table_begin",clin1[,1])+1
f1=read.table(paste(fName,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=offset)
names(f1)=sub("X.","",names(f1),fixed=T)
names(f1)=sub(".","",names(f1),fixed=T)
names(f1)[match("ID_REF",names(f1))]=c("probeId")
f1$probeId=gsub("\"","",f1$probeId)
f1=f1[1:(grep("!series_matrix_table_end",f1$probeId)-1),]
write.table(f1,file="beta_ivorra.txt", sep="\t", col.names=T, row.names=F, quote=F)


offset=3
fileList=dir(datadir,pattern="txt")
fileList=fileList[grep("GSM",fileList)]

id="GSM1568475"
idList=sapply(fileList,function(x) strsplit(x,"-")[[1]][1], USE.NAMES=F)

for (id in idList) {

	fName=paste(datadir,"GSM1568475-3018.txt",sep="")
	fName=paste(datadir,fileList[grep(id,fileList)],sep="")
	f2=read.table(paste(fName,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=offset)
	names(f2)[match(c("ID_REF","VALUE","Detection.Pval"),names(f2))]=c("probeId","beta","detP")
	f2=f2[match(f1$probeId,f2$probeId),]

	j=which(names(f1)==id)

	cat("\n\n================",id,"=============\n")
	print("table(f1[,j]==f2$beta,exclude=NULL)")
	print(table(f1[,j]==f2$beta,exclude=NULL))
	print("summary(f2$detP)")
	print(summary(f2$detP))
}


fName=paste(datadir,"GSM1568486-3029.txt",sep="")
f22=read.table(paste(fName,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=offset)

### --------------------------------
fName=paste(datadir,"GSE64316_series_matrix.txt",sep="")

clin1=read.table(paste(fName,sep=""), sep="\n", h=F, quote="", comment.char="",as.is=T,fill=T,nrow=400)
clin1=cbind(clin1,junk=rep("",nrow(clin1)))
clin1=clin1[which(gsub("\t","",clin1[,1])!=""),]
offset=grep("!series_matrix_table_begin",clin1[,1])+1
id=read.table(paste(fName,sep=""), sep="\t", h=F, quote="", comment.char="",as.is=T,fill=T,skip=offset,nrow=1)
id=as.vector(as.matrix(id)[1,-1])
id=gsub("\"","",id)
k=grep("!Sample_characteristics",clin1[,1]); offset=min(k)
clin=read.table(paste(fName,sep=""), sep="\t", h=F, quote="", comment.char="",as.is=T,fill=T,skip=offset,nrow=length(k))
clin=t(as.matrix(clin))
clin=gsub("\"","",clin)
nm=clin[1,]
clin=clin[-1,]
for (k in 1:ncol(clin)) {
	nm[k]=strsplit(clin[1,k],": ")[[1]][1]
	clin[,k]=sub(paste(nm[k],": ",sep=""),"",clin[,k])
}
colnames(clin)[match(c("prenatal tobacco exposure","newborn gender","tissue"),nm)]=c("smoke","sex","tissue")
clin=cbind(id,as.data.frame(clin,stringsAsFactors=F),stringsAsFactors=F)
clin=clin[!is.na(clin$id),]
clin$smoke=as.integer(as.factor(clin$smoke))-1
write.table(clin,file="clin_ivorra.txt", sep="\t", col.names=T, row.names=F, quote=F)

########################################################################
########################################################################
########################################################################
## Compare set1 & set2 clinical information & beta values

dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
dirSrc3=paste("code/",sep="")
dirSrc3=paste("/Users/royr/Downloads/wiemelsJ-all/",sep="")
setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))

##############################################

computerFlag="cluster"
computerFlag=""

##############################################

setFlag="set2"
setFlag=""

subsetFlag="propOfCacoEthnAsInSet2"
subsetFlag="case"
subsetFlag="ctrl"
subsetFlag=""

##############################################

source(paste(dirSrc3,"funcs.R",sep=""))
res=getClinData(setFlag=setFlag,subsetFlag=subsetFlag)
clin1=res$clin1
clin2=res$clin2
datadir11=res$dirClin1
datadir21=res$dirClin2
datadir12=res$dirMeth1
datadir22=res$dirMeth2
rm(res)

##############################################
## Section 1

## -------------------------------------------
## Compare set1 & set2 clinical variables

varName=c("caco","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
for (k1 in 1:length(varName)) {
	cat("\n\n===============",varName[k1],"===============\n")
	if (setFlag=="") {
		cat("set1")
		clin=clin1
		print(table(clin[,varName[k1]]))
		cat("\nset2")
	}
	clin=clin2
	print(table(clin[,varName[k1]]))
}

varName=c("caco","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
for (k1 in 1:length(varName)) {
	if (setFlag=="") {
		cat("\n\n=============== Set1 ===============\n\n")
		clin=clin1
		print(table(clin[,varName[k1]],dnn=varName[k1]))
		cat("\n\n=============== Set2 ===============\n\n")
	}
	clin=clin2
	print(table(clin[,varName[k1]],dnn=varName[k1]))
}

tbl=NULL
varName=c("caco","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
for (k1 in 1:(length(varName)-1)) {
	for (k2 in (k1+1):length(varName)) {
		if (setFlag=="") {
	#		cat("\n\n================",varName[k1]," and ",varName[k2],"================\n")
			clin=clin1
			j=!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]])
			x=table(clin[j,varName[k1]],clin[j,varName[k2]],dnn=varName[c(k1,k2)])
			if (nrow(x)>1 & ncol(x)>1) {
				cat("\n\n=============== Set1 ===============\n\n")
				print(table(clin[,varName[k1]],clin[,varName[k2]],dnn=varName[c(k1,k2)]))
				tbl1=c("set1",varName[k1],varName[k2],"Chi square test",chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,"")
				tbl=rbind(tbl,tbl1)
				#cat("set1 ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n")
				cat("\nChi-square test p-value: ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n",sep="")
				cat("\n\n=============== Set2 ===============\n\n")
			}
		}
		clin=clin2
		j=!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]])
		x=table(clin[j,varName[k1]],clin[j,varName[k2]],dnn=varName[c(k1,k2)])
		print(x)
		if (nrow(x)>1 & ncol(x)>1) {
			tbl1=c("set2",varName[k1],varName[k2],"Chi square test",chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,"")
			tbl=rbind(tbl,tbl1)
			cat("\nChi-square test p-value: ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n",sep="")
		}
		cat("\n\n==============================\n\n")
	}
}

varName=c("caco","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
for (k1 in 1:(length(varName)-1)) {
	for (k2 in (k1+1):length(varName)) {
		if (setFlag=="") {
	#		cat("\n\n================",varName[k1]," and ",varName[k2],"================\n")
			clin=clin1
			j=!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]])
			x=table(clin[j,varName[k1]],clin[j,varName[k2]],dnn=varName[c(k1,k2)])
			if (nrow(x)>1 & ncol(x)>1) {
				cat("\n\n=============== Set1 ===============\n\n")
				print(round(table(clin[,varName[k1]],clin[,varName[k2]],dnn=varName[c(k1,k2)])/(sum(!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]]))),2))
		#		cat("set1 ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n")
				cat("\nChi-square test p-value: ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n",sep="")
				cat("\n\n=============== Set2 ===============\n\n")
			}
		}
		clin=clin2
		j=!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]])
		x=table(clin[j,varName[k1]],clin[j,varName[k2]],dnn=varName[c(k1,k2)])
		if (nrow(x)>1 & ncol(x)>1) {
			print(round(x/(sum(!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]]))),2))
			cat("\nChi-square test p-value: ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n",sep="")
		}
	}
}

varName=c("caco","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
for (k1 in 1:(length(varName)-1)) {
	for (k2 in (k1+1):length(varName)) {
		cat("\n\n================",varName[k1]," and ",varName[k2],"================\n")
		if (setFlag=="") {
			clin=clin1
			j=!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]])
			x=table(clin[j,varName[k1]],clin[j,varName[k2]],dnn=varName[c(k1,k2)])
			if (nrow(x)>1 & ncol(x)>1) {
				cat("set1 ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n")
			}
		}
		clin=clin2
		j=!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]])
		x=table(clin[j,varName[k1]],clin[j,varName[k2]],dnn=varName[c(k1,k2)])
		if (nrow(x)>1 & ncol(x)>1) {
			cat("set2 ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n")
		}
	}
}

varName1=c("caco","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
varName2=c("gestage","DFE_nat","DFE_tot","DFE_Food","DFE_fort","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after_N","smoke_fa_3months_N","smoke_mo_bf_N","birthWt","dbirwt","pred_btw","pobw")
for (k1 in 1:length(varName1)) {
	for (k2 in 1:length(varName2)) {
		cat("\n\n================",varName1[k1]," and ",varName2[k2],"================\n")
		clin=clin2
		x=table(clin[,varName1[k1]])
		if (length(x)==2) {
			cat("================ Wilcoxon test \n")
			if (setFlag=="") {
				clin=clin1
				x=table(clin[,varName1[k1]])
				if (length(x)>1) {
					tbl1=c("set1",varName1[k1],varName2[k2],"Wilcoxon test",wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value,"")
					tbl=rbind(tbl,tbl1)
					cat("set1 ",signif(wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value,5),"\n")
				}
			}
			clin=clin2
			x=table(clin[,varName1[k1]])
			if (length(x)>1) {
				tbl1=c("set2",varName1[k1],varName2[k2],"Wilcoxon test",wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value,"")
				tbl=rbind(tbl,tbl1)
				if (setFlag=="") {
					ttl="set2 "
				} else {
					ttl="P-value "
				}
				cat(ttl,signif(wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value,5),"\n")
			}
		} else {
			cat("================ Kruskal-Wallis test \n")
			if (setFlag=="") {
				clin=clin1
				x=table(clin[,varName1[k1]])
				if (length(x)>1) {
					tbl1=c("set1",varName1[k1],varName2[k2],"Kruskal-Wallis test",kruskal.test(clin[,varName2[k2]]~as.factor(clin[,varName1[k1]]))$p.value,"")
					tbl=rbind(tbl,tbl1)
					cat("set1 ",signif(kruskal.test(clin[,varName2[k2]]~as.factor(clin[,varName1[k1]]))$p.value,5),"\n")
				}
			}
			clin=clin2
			x=table(clin[,varName1[k1]])
			if (length(x)>1) {
				tbl1=c("set2",varName1[k1],varName2[k2],"Kruskal-Wallis test",kruskal.test(clin[,varName2[k2]]~as.factor(clin[,varName1[k1]]))$p.value,"")
				tbl=rbind(tbl,tbl1)
				if (setFlag=="") {
					ttl="set2 "
				} else {
					ttl="P-value "
				}
				cat(ttl,signif(kruskal.test(clin[,varName2[k2]]~as.factor(clin[,varName1[k1]]))$p.value,5),"\n")
			}
		}
	}
}

varName=c("DFE_nat","DFE_tot","DFE_Food","DFE_fort","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after_N","smoke_fa_3months_N","smoke_mo_bf_N","birthWt","dbirwt","pred_btw","pobw","gestage")
for (k1 in 1:(length(varName)-1)) {
	for (k2 in (k1+1):length(varName)) {
		cat("\n\n================",varName[k1]," and ",varName[k2],"================\n")
		cat("================ Kendall's correlation test p-value:\n")
		if (setFlag=="") {
			clin=clin1
			res=cor.test(clin[,varName[k1]],clin[,varName[k2]],method="kendall")
			tbl1=c("set1",varName[k1],varName[k2],"Kendall's correlation test",res$p.value,res$estimate)
			tbl=rbind(tbl,tbl1)
			cat("set1 ",signif(res$p.value,5),"\n")
		}
		if (setFlag=="") {
			ttl="set2 "
		} else {
			ttl="P-value "
		}
		clin=clin2
		res=cor.test(clin[,varName[k1]],clin[,varName[k2]],method="kendall")
		tbl1=c("set2",varName[k1],varName[k2],"Kendall's correlation test",res$p.value,res$estimate)
		tbl=rbind(tbl,tbl1)
		cat(ttl,signif(res$p.value,5),"\n")
	}
}
rownames(tbl)=NULL
colnames(tbl)=c("set","variable1","variable2","test","pValue","corrCoef")
tbl=as.data.frame(tbl,stringsAsFactors=F)
for (colId in c(c("pValue","corrCoef"))) {
	tbl[,colId]=as.numeric(tbl[,colId])
}

for (varName in c("caco","sex","ch_hispanic_bc","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","mo_race","fa_race","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")) {
	cat("\n\n===============",varName,"===============\n")
	if (setFlag=="") {
		print(table(set1=clin1[,varName]))
		cat("\n")
	}
	print(table(set2=clin2[,varName]))
}

for (varName in c("caco","sex","ch_hispanic_bc","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","mo_race","fa_race","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")) {
	cat("\n\n===============",varName,"===============\n")
	if (setFlag=="") {
	#	print(table(set1=clin1[,varName]))
		print(round(table(set1=clin1[,varName])/sum(!is.na(clin1[,varName])),2))
		cat("\n")
	}
#	print(table(set2=clin2[,varName]))
	print(round(table(set2=clin2[,varName])/sum(!is.na(clin2[,varName])),2))
}

png("densityPlots_%02d.png")
for (varName in c("gestage","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after_N","smoke_fa_3months_N","DFE_nat","DFE_tot","DFE_Food","DFE_fort","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after_N","smoke_fa_3months_N","smoke_mo_bf_N","birthWt","dbirwt","pred_btw","pobw")) {
	par(mfrow=c(2,1))
	lim=range(c(clin1[,varName],clin2[,varName]),na.rm=T)
	if (setFlag=="") {
		plot(density(clin1[,varName],na.rm=T),xlim=lim,main="Set1",xlab=varName)
	}
	plot(density(clin2[,varName],na.rm=T),xlim=lim,main="Set2",xlab=varName)
}
dev.off()

png("varPairPlots_%02d.png",width=2*240,height=2*240)
par(mfcol=c(2,2))
for (var1 in c("birthWt")) {
	for (var2 in c("pred_btw","pobw")) {
		xlim=range(c(clin1[,var1],clin2[,var1]),na.rm=T)
		ylim=range(c(clin1[,var2],clin2[,var2]),na.rm=T)
		if (setFlag=="") {
			plot(clin1[,var1],clin1[,var2],xlim=xlim,ylim=ylim,main="Set1",xlab=var1,ylab=var2)
		}
		plot(clin2[,var1],clin2[,var2],xlim=xlim,ylim=ylim,main="Set2",xlab=var1,ylab=var2)
	}
}
dev.off()

tbl$variable=tbl$variable2
k=which(tbl$variable1%in%varSmoke$varIn)
tbl$variable[k]=tbl$variable1[k]
if (subsetFlag=="") {
	k=which(tbl$variable1=="caco" & !tbl$variable2%in%varSmoke$varIn)
	tbl$variable[k]=tbl$variable1[k]
}

set="set1"
set="set2"
var1=c("smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after_N","smoke_fa_3months_N","smoke_mo_bf_N","birthWt","dbirwt","pred_btw","pobw")
if (subsetFlag=="") {
	var1=c("caco",var1)
	var2=c("caco",var2)
}
var2=c("sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","gestage")
k=which(tbl$set==set & ((tbl$variable1%in%var1 & tbl$variable2%in%var2) | (tbl$variable1%in%var2 & tbl$variable2%in%var1)))
k=which(tbl$set==set & ((tbl$variable1%in%var1 & tbl$variable2%in%var2) | (tbl$variable1%in%var2 & tbl$variable2%in%var1)) & tbl$pValue<.07)
k=k[order(tbl$variable[k])]
tbl[k,]

save(tbl,clin1,clin2,file=paste("tmp_",subsetFlag,".RData",sep=""))


## ----------------------------------

## ----------------------------------
## Compare cell mixture estimates
## NLR = Granulocytes/(CD4T + CD8T + NK + Bcell)
## Refactor

library(coin)

datadir="results/cellMixture/"
cellMix1=read.table(paste(datadir,"cellMixCoef_allGuthSet1_leukChip.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
cellMix2=read.table(paste(datadir,"cellMixCoef_allGuthSet2_leukChip.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
cellMix1=cellMix1[match(paste("X",clin1$guthrieId,sep=""),cellMix1$guthrieId),]
cellMix2=cellMix2[match(paste("X",clin2$guthrieId,sep=""),cellMix2$guthrieId),]

cellMix1$nlr=cellMix1$Gran/(cellMix1$CD4T+cellMix1$CD8T+cellMix1$NK+cellMix1$Bcell)
cellMix2$nlr=cellMix2$Gran/(cellMix2$CD4T+cellMix2$CD8T+cellMix2$NK+cellMix2$Bcell)

#wilcox_test(cell, data, subset = NULL, weights = NULL, ...)

clin1=cbind(clin1,cellMix1[,which(!names(cellMix1)%in%names(clin1))])
clin2=cbind(clin2,cellMix2[,which(!names(cellMix2)%in%names(clin2))])


load("docs/SemiraGonsethNussle/refactor/Refactor_set1_and_2_K_6.RData")
rownames(Refactor_set1)=paste("X",rownames(Refactor_set1),sep="")
rownames(Refactor_set2)=paste("X",rownames(Refactor_set2),sep="")
colnames(Refactor_set1)=paste("refactor_",colnames(Refactor_set1),sep="")
colnames(Refactor_set2)=paste("refactor_",colnames(Refactor_set2),sep="")
Refactor_set1=Refactor_set1[match(clin1$id,rownames(Refactor_set1)),]
Refactor_set2=Refactor_set2[match(clin2$id,rownames(Refactor_set2)),]

clin1=cbind(clin1,Refactor_set1)
clin2=cbind(clin2,Refactor_set2)

for (varName in c("nlr","wbc",colnames(Refactor_set1))) {
	if (any(!is.na(clin1[,varName])) | any(!is.na(clin2[,varName]))) {
		cat("\n================",varName,"================\n")
		png(paste("densityPlots_",varName,".png",sep=""))
		par(mfrow=c(2,1))
		lim=range(c(clin1[,varName],clin2[,varName]),na.rm=T)
		ylim=c(0,1)
		if (varName=="wbc") {
			lim=c(0,100)
			ylim=c(0,.05)
		}
		if (setFlag=="") {
			j=which(clin1[,varName]<=lim[2])
			plot(density(clin1[j,varName],na.rm=T),xlim=lim,ylim=ylim,main=paste("Set1: Median ",round(median(clin1[j,varName],na.rm=T),2),sep=""),xlab=varName)
			axis(side=3)
			print(summary(clin1[,varName]))
		}
		j=which(clin2[,varName]<=lim[2])
		plot(density(clin2[j,varName],na.rm=T),xlim=lim,ylim=ylim,main=paste("Set2: Median ",round(median(clin2[j,varName],na.rm=T),2),sep=""),xlab=varName)
		axis(side=3)
		print(summary(clin2[,varName]))
		dev.off()
	}
}

tbl=NULL
#varName1=c("Subtype")
varName1=c("caco","hyperdipCtrl","telamlCtrl","nonHypTelamlCtrl","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","race3","DFE_sup","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
varName2=c("nlr",colnames(Refactor_set1))
for (k2 in 1:length(varName2)) {
	cat("================ Wilcoxon or Kruskal-Wallis test \n")
	for (k1 in 1:length(varName1)) {
		cat("\n================",varName1[k1]," and ",varName2[k2],"================\n")
		clin=clin2
		x=table(clin[,varName1[k1]])
		if (length(x)>1) {
			png(paste("boxPlots_",varName2[k2],"_groupBy",varName1[k1],".png",sep=""))
			par(mfrow=c(1,2))
			lim=range(c(clin1[,varName2[k2]],clin2[,varName2[k2]]),na.rm=T)
			ttl=varName1[k1]
			ttl=sub("Ctrl"," vs. ctrl",ttl)
			ttl=sub("hyperdipTelaml","hyperdip vs. telaml1",ttl)
			ttl=sub("hyperdipNonHypTelaml","hyperdip vs. non-hyperdip-telaml1",ttl)
			ttl=sub("telamlNonHypTelaml","telaml1 vs. non-hyperdip-telaml1",ttl)
			boxplot(clin1[,varName2[k2]]~clin1[,varName1[k1]],ylim=lim,main="Set1",xlab=ttl,ylab=varName2[k2])
			boxplot(clin2[,varName2[k2]]~clin2[,varName1[k1]],ylim=lim,main="Set2",xlab=ttl,ylab=varName2[k2])
			dev.off()
			clin=clin2
			x=table(clin[,varName1[k1]])
			if (length(x)==2) {
		#		cat("================ Wilcoxon test \n")
				if (setFlag=="") {
					clin=clin1
					x=table(clin[,varName1[k1]])
					if (length(x)>1) {
			#			pv=pValue(wilcox_test(clin[,varName2[k2]]~as.factor(clin[,varName1[k1]]),distribution="exact"))
						pv=wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value
						tbl1=c("set1",varName1[k1],varName2[k2],"Wilcoxon test",pv,"")
						tbl=rbind(tbl,tbl1)
						suf=""
						if (pv<0.05) {
							suf=" ****"
						} else if (pv<0.1) {
							suf=" **"
						}
						cat("set1 ",signif(pv,5),suf,"\n")
					}
				}
				clin=clin2
				x=table(clin[,varName1[k1]])
				if (length(x)>1) {
					pv=wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value
					tbl1=c("set2",varName1[k1],varName2[k2],"Wilcoxon test",pv,"")
					tbl=rbind(tbl,tbl1)
					if (setFlag=="") {
						ttl="set2 "
					} else {
						ttl="P-value "
					}
					suf=""
					if (pv<0.05) {
						suf=" ****"
					} else if (pv<0.1) {
						suf=" **"
					}
					cat(ttl,signif(pv,5),suf,"\n")
				}
			} else {
		#		cat("================ Kruskal-Wallis test \n")
				if (setFlag=="") {
					clin=clin1
					x=table(clin[,varName1[k1]])
					if (length(x)>1) {
						pv=kruskal.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value
						tbl1=c("set1",varName1[k1],varName2[k2],"Kruskal-Wallis test",pv,"")
						tbl=rbind(tbl,tbl1)
						suf=""
						if (pv<0.05) {
							suf=" ****"
						} else if (pv<0.1) {
							suf=" **"
						}
						cat("set1 ",signif(pv,5),suf,"\n")
					}
				}
				clin=clin2
				x=table(clin[,varName1[k1]])
				if (length(x)>1) {
					pv=kruskal.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value
					tbl1=c("set2",varName1[k1],varName2[k2],"Kruskal-Wallis test",pv,"")
					tbl=rbind(tbl,tbl1)
					if (setFlag=="") {
						ttl="set2 "
					} else {
						ttl="P-value "
					}
					suf=""
					if (pv<0.05) {
						suf=" ****"
					} else if (pv<0.1) {
						suf=" **"
					}
					cat(ttl,signif(pv,5),suf,"\n")
				}
			}
		}
	}
}

#varName1=c("dbirwt","pred_btw")
varName1=c("gestage","wbc","DFE_nat","DFE_tot","DFE_Food","DFE_fort","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after_N","smoke_fa_3months_N","smoke_mo_bf_N","birthWt","pobw")
varName2=c("nlr",colnames(Refactor_set1))
for (k2 in 1:length(varName2)) {
	cat("================ Kendall's correlation test p-value:\n")
	for (k1 in 1:length(varName1)) {
		cat("\n================",varName1[k1]," and ",varName2[k2],"================\n")
	#	cat("================ Kendall's correlation test p-value:\n")
		if (any(!is.na(clin1[,varName1[k1]]) & any(!is.na(clin1[,varName2[k2]])))) {
			png(paste("plot_",varName1[k1],"_vs_",varName2[k2],".png",sep=""))
			par(mfrow=c(1,2))
			xlim=range(c(clin1[,varName1[k1]],clin2[,varName1[k1]]),na.rm=T)
			ylim=range(c(clin1[,varName2[k2]],clin2[,varName2[k2]]),na.rm=T)
			plot(clin1[,varName1[k1]],clin1[,varName2[k2]],xlim=xlim,ylim=ylim,main="Set1",xlab=varName1[k1],ylab=varName2[k2])
			plot(clin2[,varName1[k1]],clin2[,varName2[k2]],xlim=xlim,ylim=ylim,main="Set2",xlab=varName1[k1],ylab=varName2[k2])
			dev.off()
			png(paste("rankPlot_",varName1[k1],"_vs_",varName2[k2],".png",sep=""))
			par(mfrow=c(1,2))
			x11=as.integer(as.factor(clin1[,varName1[k1]]))
			x21=as.integer(as.factor(clin2[,varName1[k1]]))
			x12=as.integer(as.factor(clin1[,varName2[k2]]))
			x22=as.integer(as.factor(clin2[,varName2[k2]]))
			x11=order(clin1[,varName1[k1]])
			x21=order(clin2[,varName1[k1]])
			x12=order(clin1[,varName2[k2]])
			x22=order(clin2[,varName2[k2]])
			xlim=range(c(x11,x21),na.rm=T)
			ylim=range(c(x12,x22),na.rm=T)
			plot(x11,x12,xlim=xlim,ylim=ylim,main="Set1",xlab=paste("Ranked ",varName1[k1],sep=""),ylab=paste("Ranked ",varName2[k2],sep=""))
			plot(x21,x22,xlim=xlim,ylim=ylim,main="Set2",xlab=paste("Ranked ",varName1[k1],sep=""),ylab=paste("Ranked ",varName2[k2],sep=""))
			dev.off()
			if (setFlag=="") {
				clin=clin1
				res=cor.test(clin[,varName1[k1]],clin[,varName2[k2]],method="kendall")
				pv=res$p.value
				tbl1=c("set1",varName1[k1],varName2[k2],"Kendall's correlation test",res$p.value,res$estimate)
				tbl=rbind(tbl,tbl1)
				suf=""
				if (pv<0.05) {
					suf=" ****"
				} else if (pv<0.1) {
					suf=" **"
				}
		#		cat("set1 ",signif(pv,5),suf,"\n")
				cat("set1 : Corr coef ",round(res$estimate,2)," (pv ",signif(pv,5),suf,")","\n",sep="")
			}
			if (setFlag=="") {
				ttl="set2 "
			} else {
				ttl="P-value "
			}
			clin=clin2
			res=cor.test(clin[,varName1[k1]],clin[,varName2[k2]],method="kendall")
			pv=res$p.value
			tbl1=c("set2",varName1[k1],varName2[k2],"Kendall's correlation test",res$p.value,res$estimate)
			tbl=rbind(tbl,tbl1)
			suf=""
			if (pv<0.05) {
				suf=" ****"
			} else if (pv<0.1) {
				suf=" **"
			}
		#	cat(ttl,signif(pv,5),suf,"\n")
			cat(ttl,": Corr coef ",round(res$estimate,2)," (pv ",signif(pv,5),suf,")","\n",sep="")
		}
	}
}

rownames(tbl)=NULL
colnames(tbl)=c("set","variable1","variable2","test","pValue","corrCoef")
tbl=as.data.frame(tbl,stringsAsFactors=F)
for (colId in c(c("pValue","corrCoef"))) {
	tbl[,colId]=as.numeric(tbl[,colId])
}

tbl[(substr(tbl$variable1,1,nchar("refactor"))=="refactor" | substr(tbl$variable2,1,nchar("refactor"))=="refactor") & tbl$pValue<0.05,]

## ----------------------------------
## SEM

if (computerFlag=="cluster") {
    ann=read.delim(paste("data/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
    snpVec=read.table(paste("data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
} else {
    ann=read.delim(paste("docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
    snpVec=read.table(paste("docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    dmr=read.table(paste("docs/SeungTae/leukemia.DMRs/leukemia.DMRs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
}
ann[which(ann[,"CHR"]=="X"),"CHR"]="23"
ann[which(ann[,"CHR"]=="Y"),"CHR"]="24"
ann[,"CHR"]=as.integer(ann[,"CHR"])
ann=ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
for (k in 1:ncol(ann)) if (class(ann[,k])=="factor") ann[,k]=as.character(ann[,k])

snpVec=snpVec[,1]
ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
ann$keep=as.integer(ann$snp==0 & ann$CHR%in%1:22)

load(paste(datadir21,"sem_",datType2,".RData",sep=""))
i=match(rownames(sem),ann[,"IlmnID"])
table(is.na(i))
ann=ann[i,]

x=apply(sem,2,function(x) {
    sum(!is.na(x))
})
clin2$sem=NA
clin2$sem[match(names(x),clin$id)]=x

## ----------------------------------
## Compare set1 & set2 methylation data

if (computerFlag=="cluster") {
	ann=read.delim(paste("data/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
	snpVec=read.table(paste("data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
} else {
	ann=read.delim(paste("docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
	snpVec=read.table(paste("docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	dmr=read.table(paste("docs/SeungTae/leukemia.DMRs/leukemia.DMRs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
}
ann[which(ann[,"CHR"]=="X"),"CHR"]="23"
ann[which(ann[,"CHR"]=="Y"),"CHR"]="24"
ann[,"CHR"]=as.integer(ann[,"CHR"])
ann=ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
for (k in 1:ncol(ann)) if (class(ann[,k])=="factor") ann[,k]=as.character(ann[,k])

snpVec=snpVec[,1]
ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
ann$keep=as.integer(ann$snp==0 & ann$CHR%in%1:22)

write.table(c("cpgId",ann$IlmnID[!ann$keep]),file="CpGs_to_exclude_Sept_24_fromSemira_plus_chrXY.txt", sep="\t", col.names=F, row.names=F, quote=F, append=F)

id=set=step=c()
id2=c("X1339G","X1298G","X0508G")
id=c(id,id2)
set=c(set,rep("set1",length(id2)))
step=c(step,rep("funNorm",length(id2)))
id2=c("X1762G","X0635G","X1588G","BM1660","BM1683")
id=c(id,id2)
set=c(set,rep("set2",length(id2)))
step=c(step,rep("funNorm",length(id2)))
id2=c("X0827G","X2220G","X1386G","X0319G","X0126G","X2222G","X0066G","X0133G")
id=c(id,id2)
set=c(set,rep("set2",length(id2)))
step=c(step,rep("bmiq",length(id2)))
write.table(data.frame(id,set,step),file="samples_to_exclude.txt", sep="\t", col.names=T, row.names=F, quote=F, append=F)

if (length(fileList1)>=2) {
    fName=fileList1[2]
    tmp=read.table(paste(datadir12,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
    beta1=read.table(paste(datadir12,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,col.names=names(tmp),nrow=nRow)
    fName=fileList2[2]
    tmp=read.table(paste(datadir22,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
    beta2=read.table(paste(datadir22,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,col.names=names(tmp),nrow=nRow)
}

if (F) {
    fName=fileList1[5]
    load(file=paste(fName,".RData",sep=""))
    beta12=beta
    fName=fileList2[5]
    load(file=paste(fName,".RData",sep=""))
    beta22=beta
    rm(beta)
}


fName=fileList1[betaFileId]
tmp=read.table(paste(datadir12,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
beta12=read.table(paste(datadir12,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,col.names=names(tmp),nrow=nRow)
fName=fileList2[betaFileId]
tmp=read.table(paste(datadir22,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
beta22=read.table(paste(datadir22,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,col.names=names(tmp),nrow=nRow)
rownames(beta12)=beta12$probeId
rownames(beta22)=beta22$probeId
if (length(fileList1)<2) {
    beta1=beta12
    beta2=beta22
}
beta12=as.matrix(beta12[-1])
beta22=as.matrix(beta22[-1])
save(beta1,beta2,clin1,clin2,beta12,beta22,ann,file="beta_tmp.RData")

rownames(beta1)=beta1$probeId
rownames(beta2)=beta2$probeId
beta1=as.matrix(beta1[-1])
beta2=as.matrix(beta2[-1])
clin1=clin1[match(colnames(beta1),clin1$id),]
clin2=clin2[match(colnames(beta2),clin2$id),]

save(beta1,beta2,clin1,clin2,beta12,beta22,ann,file="beta_tmp.RData")

load("beta_tmp.RData")

x1=apply(beta1,2,median,na.rm=T)
x2=apply(beta2,2,median,na.rm=T)
x12=apply(beta12,2,median,na.rm=T)
x22=apply(beta22,2,median,na.rm=T)

tbl=rbind(cbind(id=colnames(beta1),medianBeta_funNorm=round(x1,4),medianBeta_bmiq=round(x12,4),set="set1"),
	cbind(id=colnames(beta2),medianBeta_funNorm=round(x2,4),medianBeta_bmiq=round(x22,4),set="set2"))
write.table(tbl,file=paste("summaryBeta_forSample.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

for (dataset in c("_set1","_set2")) {
	switch(dataset,
		   "_set1"={
		   betaF=beta1
		   betaA=beta12
		   clin=clin1
		   },
		   "_set2"={
		   betaF=beta2
		   betaA=beta22
		   clin=clin2
		   }
		   )
	for (procFlag in c("_funNorm","_wbcAdj")) {
		switch(procFlag,
			   "_funNorm"={
			   beta=betaF
			   },
			   "_wbcAdj"={
			   beta=betaA
			   }
			   )
		iA=match(rownames(beta),ann$IlmnID)
		for (geneFlag in c("_noSnpAutosome","_snpAutosome","_autosome")) {
			switch(geneFlag,
				   "_noSnpAutosome"={
				   i=ann$keep[iA]==1
				   },
				   "_snpAutosome"={
				   i=ann$snp[iA]==1 & !ann$CHR[iA]%in%c(23,24)
				   },
				   "_autosome"={
				   i=!ann$CHR[iA]%in%c(23,24)
				   }
				   )
			for (subsetFlag in c("","_ctrlSubset","_caseSubset")) {
				j=1:nrow(clin)
				switch(dataset,
					   "_ctrlSubset"={
					   j=which(clin$caco==0)
					   },
					   "_caseSubset"={
					   j=which(clin$caco==1)
					   }
					   )
				x=apply(beta[i,j],1,mad,na.rm=T)
				png(paste("hist_cpgMad",geneFlag,subsetFlag,procFlag,dataset,".png",sep=""))
				hist(x,main="",xlab="median absolute deviation")
				dev.off()
				png(paste("cpgMad",geneFlag,subsetFlag,procFlag,dataset,".png",sep=""))
				plot(sort(x),ylab="median absolute deviation")
				dev.off()
				write.table(cbind(ann[iA,][i,],mad=x),file=paste("cpgMad",geneFlag,subsetFlag,procFlag,dataset,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
			}
		}
	}
}

i=match(rownames(beta1),rownames(beta2))
i1=which(!is.na(i)); i2=i[i1]
iA=match(rownames(beta1)[i1],ann$IlmnID)

j1=which(clin1$caco==0)
j2=which(clin2$caco==0)
clust1=hclust(dist(t(as.matrix(beta1[i1,j1][which(ann$keep[iA]==1),]))),method="ward.D2")
clust2=hclust(dist(t(as.matrix(beta2[i2,j2][which(ann$keep[iA]==1),]))),method="ward.D2")
save(clust1,clust2,file="clust_ctrlSubset_tmp.RData")

clust1=hclust(dist(t(as.matrix(beta1[i1,][which(ann$keep[iA]==1),]))),method="ward.D2")
clust2=hclust(dist(t(as.matrix(beta2[i2,][which(ann$keep[iA]==1),]))),method="ward.D2")
save(clust1,clust2,file="clust_tmp.RData")

## ---------------------------------
## Calculate no. of CpGs & samples excluded

library(minfi)
load("rgsetRaw_set1.RData")
rgsetRaw1=rgsetRaw
load("rgsetRaw_set2.RData")
rgsetRaw2=rgsetRaw
load("rgsetFN_set1.RData")
rgsetFN1=rgsetFN
load("rgsetFN_set2.RData")
rgsetFN2=rgsetFN

class(rgsetRaw1$assayData)
table(featureNames(rgsetRaw1)==featureNames(rgsetRaw2),exclude=NULL)
table(featureNames(rgsetFN1)==featureNames(rgsetFN2),exclude=NULL)
table(ann$IlmnID%in%featureNames(rgsetFN1))
i1=which(!ann$IlmnID%in%featureNames(rgsetFN1))
i2=which(ann$IlmnID%in%featureNames(rgsetFN1))
table(ann$CHR[i1],exclude=NULL)
table(ann$CHR[i2],exclude=NULL)
table(ann$snp[i1],exclude=NULL)
table(ann$snp[i2],exclude=NULL)
ann[i1,][which(ann$snp[i1]!=1),][10,]
table(ann$IlmnID[i1]=="[Controls]",exclude=NULL)
table(ann$IlmnID[i2]=="[Controls]",exclude=NULL)
table(duplicated(ann$IlmnID[i1]))
table(duplicated(ann$Name[i1]))

iA=match(featureNames(rgsetFN1),ann$IlmnID)
i1=match(beta1)

i11=apply(beta1,1,function(x) mean(is.na(x)))
i21=apply(beta2,1,function(x) mean(is.na(x)))

nrow(rgsetFN1)

table(qc=i11!=1)
table(qc=i21!=1)

table(qc=i11!=1,keep=ann$keep[iA])
table(qc=i21!=1,keep=ann$keep[iA])

ncol(rgsetRaw1)
ncol(rgsetFN1)
ncol(rgsetRaw2)
ncol(rgsetFN2)

## ---------------------------------
i=which(ann$keep[iA]==1)

subsetList=c(" ","_ctrlSubset","_natCatLow_ctrlSubset","_natCatHigh_ctrlSubset")
lim=c(0,1)
for (subsetId in 1:length(subsetList)) {
	switch(subsetList[subsetId],
		" "={
		   j1=1:nrow(clin1)
		   j2=1:nrow(clin2)
		},
		"_ctrlSubset"={
		   j1=which(clin1$caco==0)
		   j2=which(clin2$caco==0)
		},
		"_natCatLow_ctrlSubset"={
		   j1=which(clin1$caco==0 & clin1$DFE_natCat==0)
		   j2=which(clin2$caco==0 & clin2$DFE_natCat==0)
		},
		"_natCatHigh_ctrlSubset"={
		   j1=which(clin1$caco==0 & clin1$DFE_natCat==1)
		   j2=which(clin2$caco==0 & clin2$DFE_natCat==1)
		}
	)
	y1=apply(beta1[i1,j1][i,],1,median,na.rm=T)
	y2=apply(beta2[i2,j2][i,],1,median,na.rm=T)
	heading=paste("Median beta-value (Pearson Corr ",round(cor(y1,y2,use="complete.obs",method="pearson"),2),")",sep="")
	png(paste("plotMedianBeta",sub(" ","",subsetList[subsetId]),"_set2VsSet1.png",sep=""))
	plot(y1,y2,xlim=lim,ylim=lim,main=heading,xlab="Set 1",ylab="Set 2")
	abline(c(0,1),col="red")
	dev.off()
	y1=apply(beta1[i1,j1][i,],1,mean,na.rm=T)
	y2=apply(beta2[i2,j2][i,],1,mean,na.rm=T)
	heading=paste("Mean beta-value (Pearson Corr ",round(cor(y1,y2,use="complete.obs",method="pearson"),2),")",sep="")
	png(paste("plotMeanBeta",sub(" ","",subsetList[subsetId]),"_set2VsSet1.png",sep=""))
	plot(y1,y2,xlim=lim,ylim=lim,main=heading,xlab="Set 1",ylab="Set 2")
	abline(c(0,1),col="red")
	dev.off()
}

subsetList=c(" ","_ctrlSubset")
varList="DFE_natCat"
varId=1
lim=NULL
for (subsetId in 1:length(subsetList)) {
	switch(subsetList[subsetId],
		   " "={
		   j1=1:nrow(clin1)
		   j2=1:nrow(clin2)
		   },
		   "_ctrlSubset"={
		   j1=which(clin1$caco==0)
		   j2=which(clin2$caco==0)
		   },
		   "_natCatLow_ctrlSubset"={
		   j1=which(clin1$caco==0 & clin1$DFE_natCat==0)
		   j2=which(clin2$caco==0 & clin2$DFE_natCat==0)
		   },
		   "_natCatHigh_ctrlSubset"={
		   j1=which(clin1$caco==0 & clin1$DFE_natCat==1)
		   j2=which(clin2$caco==0 & clin2$DFE_natCat==1)
		   }
		   )
	j11=j1[which(clin1[j1,varList[varId]]==0)]
	j12=j1[which(clin1[j1,varList[varId]]==1)]
	j21=j2[which(clin2[j2,varList[varId]]==0)]
	j22=j2[which(clin2[j2,varList[varId]]==1)]
	
	y1=apply(beta1[i1,j12][i,],1,mean,na.rm=T)-apply(beta1[i1,j11][i,],1,mean,na.rm=T)
	y2=apply(beta2[i2,j22][i,],1,mean,na.rm=T)-apply(beta2[i2,j21][i,],1,mean,na.rm=T)
	heading=paste(varList[varId],": Mean Beta Difference (Pearson Corr ",round(cor(y1,y2,use="complete.obs",method="pearson"),2),")",sep="")
	png(paste("plotMeanBetaDiff",sub(" ","",subsetList[subsetId]),"_set2VsSet1.png",sep=""))
	plot(y1,y2,xlim=lim,ylim=lim,main=heading,xlab="Set 1",ylab="Set 2")
	abline(c(0,1),col="red")
	dev.off()
	
	y11=apply(beta1[i1,j11][i,],1,mean,na.rm=T)
	y12=apply(beta1[i1,j12][i,],1,mean,na.rm=T)
	y21=apply(beta2[i2,j21][i,],1,mean,na.rm=T)
	y22=apply(beta2[i2,j22][i,],1,mean,na.rm=T)
	cor(y11,y21,use="complete.obs",method="pearson")
	cor(y12,y22,use="complete.obs",method="pearson")
	cor(y12-y11,y22-y21,use="complete.obs",method="pearson")
	cor(log(y12)-log(y11),log(y22)-log(y21),use="complete.obs",method="pearson")
	summary(y11)
	summary(y12)
	summary(y21)
	summary(y22)
	summary(y12-y11)
	summary(y22-y21)
	
	x1=beta1; x1=log2(x1/(1-x1))
	x2=beta2; x2=log2(x2/(1-x2))
	y11=apply(x1[i1,j11][i,],1,function(x) mean(x[!is.infinite(x)],na.rm=T))
	y12=apply(x1[i1,j12][i,],1,function(x) mean(x[!is.infinite(x)],na.rm=T))
	y21=apply(x2[i2,j21][i,],1,function(x) mean(x[!is.infinite(x)],na.rm=T))
	y22=apply(x2[i2,j22][i,],1,function(x) mean(x[!is.infinite(x)],na.rm=T))
	ii=which(!is.infinite(y11) & !is.infinite(y12) & !is.infinite(y21) & !is.infinite(y22))
	cor(y11[ii],y21[ii],use="complete.obs",method="pearson")
	cor(y12[ii],y22[ii],use="complete.obs",method="pearson")
	cor(y12[ii]-y11[ii],y22[ii]-y21[ii],use="complete.obs",method="pearson")
	
}

cor(y11,y21,use="complete.obs",method="pearson")
[1] 0.9983
cor(y12,y22,use="complete.obs",method="pearson")
[1] 0.9981654
cor(y12-y11,y22-y21,use="complete.obs",method="pearson")
[1] -0.1182888

summary(y11)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.007999 0.083150 0.525800 0.480900 0.868300 0.988100 
summary(y12)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.008217 0.084650 0.526000 0.480800 0.866700 0.987700 
summary(y21)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NAs 
0.008276 0.077670 0.485500 0.464900 0.844700 0.987000        1 
summary(y22)
Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NAs 
0.007872 0.077480 0.487800 0.465900 0.846400 0.987000        1 
summary(y12-y11)
summary(y22-y21)


## ---------------------------------
library(marray)
source(paste(dirSrc,"functions/heatmap.5.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.R",sep=""))

#load("beta_tmp.RData")

if (F) {
i1=i2=iA=prId=1:nrow(beta1)
beta1=as.matrix(beta1[,-1])
beta2=as.matrix(beta2[,-1])
}

subsetFlag=""
subsetFlag="_ctrlSubset"

load(paste("clust",subsetFlag,"_tmp.RData",sep=""))

if (F) {
prId=which(ann$keep[iA]==1)
prId=prId[1:10]

beta1=beta1[i1,][prId,]
beta2=beta2[i2,][prId,]
ann=ann[iA,][prId,]
i1=i2=iA=prId=1:nrow(ann)
clin1=clin1[match(colnames(beta1),clin1$id),]
clin2=clin2[match(colnames(beta2),clin2$id),]
}

colList=c("skyblue","blue","yellow","purple")
samColUniq=c("skyblue","blue","brown","red","orange","yellow","green","cyan","pink","magenta","purple")
colHM=c("red","blue","grey")

varList=c("type")
varName=c("type ")

distMethod="pearson"
distMethod=""
linkMethod="ward.D2"

varCatList=c("Batch","caco","sex","int_ch_ethnicity","ch_hispanic_bc","race3","Subtype","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","DFE_natCat","DFE_foodCat","DFE_totCat","DFE_supCat","smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke3","smoke2")
varContList=c("gestage","Beadchip")
varList=varName=c(varCatList,varContList)

for (datasetFlag in c("_allGuthrieSet1","_allGuthrieSet2")) {
	fName2=paste(subsetFlag,datasetFlag,sep="")
	switch(datasetFlag,
		"_allGuthrieSet1"={
		   heading="Set 1"
		   arrayData=beta1[i1,][prId,]
		   phen=clin1
		   clustC=clust1
		   },
		"_allGuthrieSet2"={
		   heading="Set 2"
		   arrayData=beta2[i2,][prId,]
		   phen=clin2
		   clustC=clust2
		}
	)
	if (subsetFlag=="_ctrlSubset") {
		j=which(phen$caco==0)
		arrayData=arrayData[,j]
		phen=phen[j,]
	}
#	centr=apply(arrayData,1,median,na.rm=T)
#	arrayData=arrayData-centr
	cloneName=rep("",nrow(arrayData))
	if (F) {
		cloneName=geneSym[id]
		clonƒphneName=featureNames(eset)[id]
		cloneName=rep("",length(id))
		#cloneCol=matrix(rep("white",nrow(arrayData)),nrow=1)
	}
	cloneCol=NULL
	samName=phen$id
	samName=rep("",nrow(phen))
	samCol=NULL
	if (F) {
		#samCol=matrix(sapply(as.character(phen$type),function(x) {if (x=="ACF") {y="skyblue"} else if(x=="ECF") {y="blue"} else {y="yellow"}; y}, USE.NAMES=F),nrow=1)
		samCol=matrix(nrow=length(varList),ncol=nrow(phen))
		for (varId in 1:length(varList)) {
			x=as.character(phen[,varList[varId]]); x[x==""]=NA; x=as.integer(as.factor(x))
			grpUniq=sort(unique(x))
			if (length(grpUniq)<=length(colList)) {
				samCol[varId,]=colList[x]
			} else {
				samCol[varId,]=rainbow(length(grpUniq))[x]
			}
		}
		rownames(samCol)=varName
	}
	samCol=matrix(rep(NA,length(samName)*length(varName)),ncol=length(samName))
	lineSam=NULL
	k1=0
	for (varId in 1:length(varList)) {
		k1=k1+1
		kk=which(names(phen)==varList[varId])
		x=phen[,kk]
		grpUniq=sort(unique(x[!is.na(x)]))
		if (varList[varId]%in%varCatList) {
			samColUniq2=samColUniq
		} else if (varList[varId]%in%c("site")) {
			samColUniq2=rainbow(length(grpUniq))
		} else {
			samColUniq2=rainbow(length(grpUniq))
			samColUniq2=gray(0:(length(grpUniq)-1)/length(grpUniq))
		}
		for (k in 1:length(grpUniq)) {samCol[k1,which(phen[,kk]==grpUniq[k])]=samColUniq2[k]}	
	}
	rownames(samCol)=varName
	
	summary(range(c(arrayData),na.rm=T))
	limit=c(0,1)
	margins=c(1,1)
	main=NULL

	switch(distMethod,
		   "pearson" = {distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
		   clustC=hclust(distMat, method=linkMethod)
		   distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
		   clustR=hclust(distMat, method=linkMethod)},
		   "spearman" = {distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
		   distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))},
		   "euclidean" = {distMat=dist(t(arrayData), method=distMethod)
		   distMat=dist(arrayData, method=distMethod)}
		   )
	png(paste("heatmap",fName2,".png",sep=""),width=480*2,height=480*2)
	hcc=heatmap3(x=arrayData, Rowv=NA, Colv=as.dendrogram(clustC), distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=samCol, RowSideColors=cloneCol, labCol=samName, labRow=cloneName,scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3])
	dev.off()
	png("heatmapColorRange.png",width=480,height=140)
	heatmapColorBar(limit=limit,cols=colHM)
	dev.off()
	if (!is.null(samCol)) {
		for (varId in 1:length(varList)) {
			png(paste("heatmapSampleColorBarLegend_",varList[varId],".png",sep=""))
			x=as.character(phen[,varList[varId]]); x[x==""]=NA
			grpUniq=table(x)
			grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
			if (varList[varId]%in%varCatList) {
				samColUniq2=samColUniq
			} else if (varList[varId]%in%c("site")) {
				samColUniq2=rainbow(length(grpUniq))
			} else {
				samColUniq2=gray(0:(length(grpUniq)-1)/length(grpUniq))
			}
			sampleColorLegend(tls=grpUniq,col=samColUniq2,legendTitle=varName[varId])
		dev.off()
		}
	}
}

## ---------------------------------
if (computerFlag=="cluster") {
	datadir=""
} else {
	datadir="results/comparison/"
}

stat11=read.table(paste(datadir,"stat_refFreeEWAS_dfeNat_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat12=read.table(paste(datadir,"stat_refFreeEWAS_dfeTot_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat13=read.table(paste(datadir,"stat_refFreeEWAS_dfeFood_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat14=read.table(paste(datadir,"stat_refFreeEWAS_dfeFort_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat15=read.table(paste(datadir,"stat_refFreeEWAS_dfeNatCat_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat16=read.table(paste(datadir,"stat_refFreeEWAS_dfeTotCat_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat17=read.table(paste(datadir,"stat_refFreeEWAS_dfeFoodCat_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat18=read.table(paste(datadir,"stat_refFreeEWAS_dfeSupCat_ctrlSubset_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat21=read.table(paste(datadir,"stat_refFreeEWAS_dfeNat_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat22=read.table(paste(datadir,"stat_refFreeEWAS_dfeTot_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat23=read.table(paste(datadir,"stat_refFreeEWAS_dfeFood_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat24=read.table(paste(datadir,"stat_refFreeEWAS_dfeFort_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat25=read.table(paste(datadir,"stat_refFreeEWAS_dfeNatCat_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat26=read.table(paste(datadir,"stat_refFreeEWAS_dfeTotCat_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat27=read.table(paste(datadir,"stat_refFreeEWAS_dfeFoodCat_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat28=read.table(paste(datadir,"stat_refFreeEWAS_dfeSupCat_ctrlSubset_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## -------------------
stat1=stat15; stat2=stat25; fName="refFreeEWAS_dfeNatCat"; compName="Natural folate (categorical)"; xlab="Set1: Estimates"; ylab="Set2: Estimates"

ii=order(stat2$pvBeta)
ii=ii[stat2$cpgId%in%ann$IlmnID[iA][which(ann$keep[iA]==1)]]
		
ii1=ii[1]
i=which(ann$IlmnID[iA]==stat2$cpgId[ii1])

summary(lm(beta1[i1[i],]~clin1$DFE_natCat))$coef
summary(lm(beta2[i2[i],]~clin2$DFE_natCat))$coef

########################################################################
########################################################################
########################################################################
## ALL Guthrie Set1 + Set2
## Process

cohort="_allGuthSet1Set2_ctrlSubset"
datadir="data/set1set2"
fName=paste("beta_funNorm",cohort,".txt",sep="")
tmp=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
beta11=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,col.names=names(tmp),nrow=nRow)
fName=paste("beta_bmiq",cohort,".txt",sep="")
tmp=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
beta12=read.table(paste(datadir,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,col.names=names(tmp),nrow=nRow)

x=apply(beta1,2,median,na.rm=T)

tbl=rbind(cbind(id=colnames(beta1),medianBeta_funNorm=round(x1,4),medianBeta_bmiq=round(x12,4),set="set1"),
cbind(id=colnames(beta2),medianBeta_funNorm=round(x2,4),medianBeta_bmiq=round(x22,4),set="set2"))
write.table(tbl,file=paste("summaryBeta_forSample",cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)


########################################################################
########################################################################
########################################################################
## PCA

computerFlag=""
computerFlag="cluster"

## ---------------------------------
subsetName2=""
datType="allGuthSet2"
datType="allGuthSet1"
datType="allGuthSet1Set2"

#for (subsetName2 in c("_set1","_set2")) {

#for (datType in c("allGuthSet2","allGuthSet1")) {


## ---------------------------------

if (computerFlag=="cluster") {
	switch(datType,
		   "allGuthSet2"={
		   datName="ALL Guthrie, set2"
		   datadir1="data/set2/"
		   datadir2=datadir1
		   fileList=c("mDat_funNorm_set2","beta_funNorm_set2")
		   fName2="clin_guthrieSet2_20140619"
		   },
		   "allGuthSet1"={
		   datName="ALL Guthrie, set1"
		   datadir1="data/set1/"
		   datadir2=datadir1
		   fileList=c("mDat_funNorm_set1","beta_funNorm_set1")
		   fName2="final"
		   },
		   "allGuthSet1Set2"={
		   datName="ALL Guthrie, set1 & set2"
		   datadir1="data/set1set2/"
		   datadir2=datadir1
		   fileList=c("mDat_funNorm_set1set2","beta_funNorm_set1set2")
		   fileList=c("beta_funNorm_set1set2")
		   fName2="clin_guthrieSet1Set2_20140619"
		   }
		   )	
	nRow=-1
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
	setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
	switch(datType,
		   "allGuthSet2"={
		   datName="ALL Guthrie, set2"
		   datadir1="docs/all/set2/"
		   datadir2=datadir1
		   fileList=c("beta_funNorm_set2")
		   fName2="clin_guthrieSet2_20140619"
		   },
		   "allGuthSet1"={
		   datName="ALL Guthrie, set1"
		   datadir1="docs/all/set1/"
		   datadir2=datadir1
		   fileList=c("beta_funNorm_set1")
		   fName2="final"
		   },
		   "allGuthSet1Set2"={
		   datName="ALL Guthrie, set1 & set2"
		   datadir1="docs/all/set1set2/"
		   datadir2=datadir1
		   fileList=c("beta_funNorm_set1set2")
		   fName2="clin_guthrieSet1Set2_20140619"
		   }
		   )
	nRow=1000
}

for (fName in fileList) {
	clin=read.table(paste(datadir1,fName2,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	tmp=read.table(paste(datadir2,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=2)
	dat=read.table(paste(datadir2,fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,col.names=names(tmp),nrow=nRow)
	
	probeId=dat$probeId
	dat=dat[,-1]
	rownames(dat)=probeId
	rm(probeId)
	
	switch(datType,
		   "allGuthSet2"={
		   if (F) {
		   clin$caco2=sapply(clin$newid,function(x) {
							 y=substr(x,5,5)
							 y
							 },USE.NAMES=F)
		   table(clin$caco2,clin$caco,exclude=NULL)
		   }
		   },
		   "allGuthSet1"={
		   clin$guthrieId=clin$TargetID
		   clin$Beadchip=substr(clin$Bead_Position,1,10)
		   clin$caco=clin$Leukemia
		   }
		   )
	clin$id=paste("X",clin$guthrieId,sep="")
	clin=clin[match(colnames(dat),clin$id),]
	
	if (subsetName2!="") {
		clin=clin[which(clin$set==sub("_","",subsetName2)),]
		j=match(colnames(dat),clin$id); j1=which(!is.na(j)); j2=j[j1]
		dat=dat[,j1]
		clin=clin[j2,]
	}
	
	for (centerFlag in c("centered","notCentered")) {
		if (centerFlag=="centered") {
			x=as.data.frame(t(apply(dat,1,function(x) {x-median(x,na.rm=T)})))
		} else {
			x=dat
		}
		y=apply(x,1,function(x) mean(is.na(x)))
		cat("\n\n----------------------------------------------------\n")
		cat(fName,", ",centerFlag,"\n")
		
		print("dim(dat)")
		print(dim(dat))
		
		y=apply(x,1,function(x) mean(is.na(x)))
		
		fmla=as.formula(paste(" ~ ", paste(colnames(x), collapse= "+")))
		fit=prcomp(fmla, center=F, scale=F, data = x, scores=T)
#fit=prcomp(x, center=F, scale=F)
		fit1=prcomp(t(x[y==0,]), center=F, scale=F)
		fit1=prcomp(t(x[y==0,]), center=F, scale=F)
		for (varId in which(names(clin)%in%c("Beadchip","caco","set"))) {
			png(paste("pc_",names(clin)[varId],"_",centerFlag,fName,".png",sep=""))
			par(mfrow=c(2,2))
			cat("-----------------",names(clin)[varId],", ",centerFlag,", ",fName,"------------------\n")
			for (k in 1:4) {
				fit2=lm(fit$x[k,]~as.factor(clin[,varId]))
				pv=anova(fit2)[1,5]
				x=table(clin[,varId])
				if (length(x)==2) {
					fit2=wilcox.test(fit$x[k,]~as.factor(clin[,varId]))
					testName="Wilcoxon"
				} else {
					fit2=kruskal.test(fit$x[k,]~as.factor(clin[,varId]))
					testName="Kruskal-Wallis"
				}
				pv=fit2$p.value
				cat("PC",k,": ",pv,"\n",sep="")
#print(fit2)
#plot(sort(fit$x[k,]))
				if (names(clin)[varId]=="caco") {
					j=which(fit$x[k,]>min(fit$x[k,],na.rm=T))
					j=1:ncol(fit$x)
				} else {
					j=1:ncol(fit$x)
				}
				j=which(fit$x[k,]>min(fit$x[k,],na.rm=T))
#				fit2=kruskal.test(fit$x[k,]~as.factor(clin$caco))
#				pv=fit2$p.value
#				if (sum(!duplicated(clin[j,varId])))<11) {
				boxplot(fit$x[k,j]~as.factor(clin[j,varId]),main=paste(datName,"\n",testName," test p-value: ",signif(pv,2),sep=""),xlab=names(clin)[varId],ylab=paste("PC",k,sep=""))
#				}
			}
			dev.off()
		}
	}
}

#}

########################################################################
########################################################################
########################################################################

clin=read.table(paste("docs/all/Guthrie.data/i.GUTH.v2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clin$id=paste("X",clin$TargetID,sep="")
clin1=clin
clin=read.table(paste("docs/all/final.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clin$id=paste("X",clin$TargetID,sep="")

clin1=clin1[match(clin$TargetID, clin1$TargetID),]

clin2=read.table("docs/all/LEU.data/i.LEU.v2.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

j=match(clin$Subject_ID,clin2$Subject_ID)
j1=which(!is.na(j)); j2=j[j1]

k=match(names(clin),names(clin2))
k1=which(!is.na(k)); k2=k[k1]

clin=clin[j1,k1]
clin2=clin2[j2,k2]

for (k in grep("hyperdiploid_",names(clin2))) {
	cat("\n\n==============",k,names(clin2)[k],"===============\n\n")
	tbl=table(clin2[,k],Subtype=clin2$Subtype,exclude=NULL)
	
}

plot(clin2$)

########################################################################
########################################################################
########################################################################

clin=read.table(paste("docs/all/final.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clin2=read.table(paste("docs/all/folate.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(clin2)=sub("_new","",names(clin3))

j=match(clin$Subject_ID,clin2$Subject_ID)
j1=which(!is.na(j)); j2=j[j1]

k=match(names(clin),names(clin2))
k1=which(!is.na(k)); k2=k[k1]

clin=clin[j1,k1]
clin2=clin2[j2,k2]


for (k in grep("DFE_",names(clin2))) {
	if (any(clin[,k]!=clin2[,k],na.rm=T) | any(is.na(clin[,k])!=is.na(clin2[,k]))) cat(k,names(clin2)[k],"\n")
	if (F) {
	if (any(clin[,k]!=clin2[,k],na.rm=T)) cat("Different value:",k,names(clin2)[k],"\n")
	if (any(is.na(clin[,k])!=is.na(clin2[,k]))) {
		cat("Different missing:",k,names(clin2)[k],"\n")
		print(table(is.na(clin[,k]),is.na(clin2[,k])))
	}
	}
}

########################################################################
########################################################################
########################################################################

clin=read.table(paste("docs/all/LEU.data/i.LEU.v2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clin21=read.table(paste("docs/all/normalBCell/BCELL.methylation.data.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=10)
clin22=read.table(paste("docs/all/normalBCell/manifest 45K set one.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
clin22=clin22[1:16,7:ncol(clin22)]

out=data.frame(names(clin2)[which(substr(names(clin21),4,6)=="_S2")],function(x) {
	id=x
	
}, USE.NAMES=F)


########################################################################
########################################################################
########################################################################

datadir="docs/all/set2/"
datadir="data/set2/"

load(paste(datadir,"lmeAdjustedMValues_allGuth_set2_leukChip.RData",sep=""))
dat2=yAdjM
rm(yAdjM)

#dat1=read.table(paste(datadir,"mDat_funNorm_set2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

tmp=read.table(paste(datadir,"mDat_funNorm_set2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=50)
colName=names(tmp)
colClass=rep("",ncol(tmp))
for (k in 1:length(colClass)) {
	colClass[k]=class(tmp[,k])
}
dat1=NULL
block=1:50
block=1:5
nProbe=1000
set.seed(476767)
offset=sample(seq(1,(nrow(dat2)-nProbe),by=nProbe),length(block),replace=F)
cat("------------- Read fun norm data -------------\n\n")
timeStamp=Sys.time()
for (k in 1:length(block)) {
	if (k%%10==0) print(k)
	tmp=read.table(paste(datadir,"mDat_funNorm_set2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=offset[k],nrow=nProbe,col.names=colName,colClasses=colClass)
	dat1=rbind(dat1,tmp)
}
timeStamp=c(timeStamp,Sys.time())
print(diff(timeStamp))

i=match(dat1$probeId,rownames(dat2)); i1=which(!is.na(i)); i2=i[i1]
j=match(names(dat1),colnames(dat2)); j1=which(!is.na(j)); j2=j[j1]
dat1=as.matrix(dat1[i1,j1])
dat2=dat2[i2,j2]
rownames(dat1)=rownames(dat2)
colnames(dat1)=colnames(dat2)

nProbe=1000
nProbe=10
set.seed(6489560)
i=sample(1:nrow(dat1),2*nProbe,replace=F)
iSame=i[1:nProbe]
iDiff=i[(nProbe+1):(2*nProbe)]
tmp=rep(NA,2*nProbe)
tmpC=rep("",2*nProbe)
out=data.frame(id1=tmpC,id2=tmpC,corP=tmp,stringsAsFactors=F)
nrow(out)
x=nrow(out)/10
k=1
cat("------------- Get correlation -------------\n\n")
timeStamp=Sys.time()
for (i in 1:length(iSame)) {
	i1=i2=iSame[i]
	if (k%%x==0) print(k)
	out$id1[k]=rownames(dat1)[i1]
	out$id2[k]=rownames(dat2)[i2]
	out$corP[k]=cor(dat1[i1,],dat2[i2,],method="pearson",use="complete.obs")
	k=k+1
}
for (i in 1:length(iSame)) {
	i1=iSame[i]
	i1=iDiff[i]
	if (k%%x==0) print(k)
	out$id1[k]=rownames(dat1)[i1]
	out$id2[k]=rownames(dat2)[i2]
	out$corP[k]=cor(dat1[i1,],dat2[i2,],method="pearson",use="complete.obs")
	k=k+1
}
timeStamp=c(timeStamp,Sys.time())
print(diff(timeStamp))
k1=which(out$id1==out$id2)
k2=(1:nrow(out))[-k1]
print("summary(out$corP[k1]) - same probe")
print(summary(out$corP[k1]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.5543  0.9638  0.9857  0.9644  0.9954  1.0000 
print("summary(out$corP[k2]) - diff probes")
print(summary(out$corP[k2]))
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.265300  0.008542  0.134700  0.139000  0.260700  0.554200 
ttl=rep("same probe",nrow(out))
ttl[k2]="diff probes"
png("boxplot_funNorm_combat_set2.png")
boxplot(out$corP~ttl,main="Set 2 beta values\nCorrelation: FunNorm vs. FunNorm+ComBat",ylab="Pearson correlation coefficient")
dev.off()

save.image("tmp_funNorm_combat_set2.RData")

########################################################################
########################################################################
########################################################################

# Methylation - SNP association

load("beta_tmp.RData")

subsetFlag="_CCLS_GWAS"
subsetFlag="_CCLS_hispanic_GWAS"

for (subsetFlag in c("_CCLS_hispanic_GWAS","_CCLS_GWAS")) {
		
	datadir1="docs/misc/"; datadir2=""
	datadir1="data/"; datadir2=""

	if (subsetFlag=="_CCLS_hispanic_GWAS") {
		snp=read.table(paste(datadir1,"ARID5B_CCLS_Hispanic_GWAS.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		names(snp)[match(c("ID","CaCo","rs7089424"),names(snp))]=c("subjectId","caco","rs7089424")
	} else {
		snp=read.table(paste(datadir1,"ARID5B_CCLS_Genotyped.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		names(snp)[match(c("Subject.ID","rs7089424"),names(snp))]=c("subjectId","rs7089424")
	}
	stat=read.table(paste(datadir2,"stat_refFree_pobw_ctrl_set1set2_bmiq_pv0.05.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)

	j=match(snp$subjectId,clin1$subjectId); j11=which(!is.na(j)); j12=j[j11]
	j=match(snp$subjectId,clin2$subjectId); j21=which(!is.na(j)); j22=j[j21]

	library("exactRankTests")

	colId="rs7089424"
	ii=grep("ARID5B",stat$UCSC_RefGene_Name)
	ii=grep("ARID5B",ann$UCSC_RefGene_Name)
	out=matrix(nrow=length(ii),ncol=4,dimnames=list(ann$IlmnID[ii],c("pv_ctrl_set1","pv_case_set1","pv_ctrl_set2","pv_case_set2")))
	for (iii in 1:length(ii)) {
		i=ii[iii]
		j=which(clin1$caco[j12]==0)
		y=beta12[i,j12[j]]; x=snp[j11[j],colId]
		if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
			pv=kruskal.test(y~x)$p.value
			out[iii,"pv_ctrl_set1"]=pv
		}

		j=which(clin1$caco[j12]==1)
		x=snp[j11[j],colId]; y=beta12[i,j12[j]]
		if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
			pv=kruskal.test(y~x)$p.value
			out[iii,"pv_case_set1"]=pv
		}

		j=which(clin2$caco[j22]==0)
		x=snp[j21[j],colId]; y=beta22[i,j22[j]]
		if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
			pv=kruskal.test(y~x)$p.value
			out[iii,"pv_ctrl_set2"]=pv
		}

		j=which(clin2$caco[j22]==1)
		x=snp[j21[j],colId]; y=beta22[i,j22[j]]
		if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
			pv=kruskal.test(y~x)$p.value
			out[iii,"pv_case_set2"]=pv
		}
	}
	summary(c(out))

	snpUniq=sort(unique(snp[,colId]))
	compList=c("0v1","0v2","1v2")
	out=matrix(nrow=length(ii),ncol=length(compList)*4,dimnames=list(ann$IlmnID[ii],paste("pv_snp",rep(compList,each=4),"_",c("ctrl_set1","case_set1","ctrl_set2","case_set2"),sep="")))
	for (iii in 1:length(ii)) {
		i=ii[iii]
		for (k1 in 1:(length(snpUniq)-1)) {
			for (k2 in (k1+1):length(snpUniq)) {
				colIdOut=paste(snpUniq[c(k1,k2)],collapse="v")
				j=which(clin1$caco[j12]==0)
				x=snp[j11[j],colId]; y=beta12[i,j12[j]]
				if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
					pv=wilcox.exact(y[which(x==snpUniq[k1])],y[which(x==snpUniq[k2])])$p.value
					out[iii,paste("pv_snp",colIdOut,"_ctrl_set1",sep="")]=pv
				}
				
				j=which(clin1$caco[j12]==1)
				x=snp[j11[j],colId]; y=beta12[i,j12[j]]
				if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
					pv=wilcox.exact(y[which(x==snpUniq[k1])],y[which(x==snpUniq[k2])])$p.value
					out[iii,paste("pv_snp",colIdOut,"_case_set1",sep="")]=pv
				}
				
				j=which(clin2$caco[j22]==0)
				x=snp[j21[j],colId]; y=beta22[i,j22[j]]
				if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
					pv=wilcox.exact(y[which(x==snpUniq[k1])],y[which(x==snpUniq[k2])])$p.value
					out[iii,paste("pv_snp",colIdOut,"_ctrl_set2",sep="")]=pv
				}
				
				j=which(clin2$caco[j22]==1)
				x=snp[j21[j],colId]; y=beta22[i,j22[j]]
				if (sum(!duplicated(x[!is.na(x) & !is.na(y)]))!=0) {
					pv=wilcox.exact(y[which(x==snpUniq[k1])],y[which(x==snpUniq[k2])])$p.value
					out[iii,paste("pv_snp",colIdOut,"_case_set2",sep="")]=pv
				}
			}
		}
	}
	summary(c(out))
	pThres=0.05
	pThres=0.001
	x=apply(out,2,function(x,pThres) {
			mean(x<pThres,na.rm=T)
			},pThres=pThres)
	names(x)=NULL
	sort(x)
	colnames(out)[which(x>0)]
	
	k=which(x>0)
	i2=which(out[,k]<pThres)
	i1=match(rownames(out)[i2],ann$IlmnID)
	for (iii in 1:length(i1)) {
		i=i1[iii]
#		if (i==1) ylim=c(.9,1) else ylim=NULL
		ylim=NULL
		png(paste("methSnp_",ann$IlmnID[i],subsetFlag,"_ARID5B.png",sep=""))
		par(mfrow=c(2,2))
		j=which(clin1$caco[j12]==0)
		pv=out[i2[iii],grep("_ctrl_set1",colnames(out))]
		boxplot(beta12[i,j12[j]]~snp[j11[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",paste(signif(pv,2)," (",compList,")",sep="",collapse=", "),sep=""),xlab="set1 ctrl: Genotype",ylab="methylation",notch=T)
		j=which(clin1$caco[j12]==1)
		pv=out[i2[iii],grep("_case_set1",colnames(out))]
		boxplot(beta12[i,j12[j]]~snp[j11[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",paste(signif(pv,2)," (",compList,")",sep="",collapse=", "),sep=""),xlab="set1 case: Genotype",ylab="methylation",notch=T)
		j=which(clin2$caco[j22]==0)
		pv=out[i2[iii],grep("_ctrl_set2",colnames(out))]
		boxplot(beta22[i,j22[j]]~snp[j21[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",paste(signif(pv,2)," (",compList,")",sep="",collapse=", "),sep=""),xlab="set2 ctrl: Genotype",ylab="methylation",notch=T)
		j=which(clin2$caco[j22]==1)
		pv=out[i2[iii],grep("_case_set2",colnames(out))]
		boxplot(beta22[i,j22[j]]~snp[j21[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",paste(signif(pv,2)," (",compList,")",sep="",collapse=", "),sep=""),xlab="set2 case: Genotype",ylab="methylation",notch=T)
		dev.off()
	}
	
	tbl=cbind(ann[ii,c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","Relation_to_UCSC_CpG_Island")],out,stringsAsFactors=F)
	write.table(tbl,file=paste("methSnp",subsetFlag,"_ARID5B.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
	switch(subsetFlag,
		   "_CCLS_GWAS"={
		   out1=out
		   },
		   "_CCLS_hispanic_GWAS"={
		   out2=out
		   }
	)

	ii=match(c("cg25953130","cg02863179"),ann$IlmnID)
	for (i in ii) {
		if (i==1) ylim=c(.9,1) else ylim=NULL
		png(paste("methSnp_",ann$IlmnID[i],"_ARID5B.png",sep=""))
		par(mfrow=c(2,2))
		j=which(clin1$caco[j12]==0)
		pv=kruskal.test(beta12[i,j12[j]]~snp[j11[j],colId])$p.value
		boxplot(beta12[i,j12[j]]~snp[j11[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",signif(pv,2),sep=""),xlab="set1 ctrl: Genotype",ylab="methylation",notch=T)
		j=which(clin1$caco[j12]==1)
		pv=kruskal.test(beta12[i,j12[j]]~snp[j11[j],colId])$p.value
		boxplot(beta12[i,j12[j]]~snp[j11[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",signif(pv,2),sep=""),xlab="set1 case: Genotype",ylab="methylation",notch=T)
		j=which(clin2$caco[j22]==0)
		pv=kruskal.test(beta22[i,j22[j]]~snp[j21[j],colId])$p.value
		boxplot(beta22[i,j22[j]]~snp[j21[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",signif(pv,2),sep=""),xlab="set2 ctrl: Genotype",ylab="methylation",notch=T)
		j=which(clin2$caco[j22]==1)
		pv=kruskal.test(beta22[i,j22[j]]~snp[j21[j],colId])$p.value
		boxplot(beta22[i,j22[j]]~snp[j21[j],colId],ylim=ylim,main=paste(ann$IlmnID[i],"\nPV: ",signif(pv,2),sep=""),xlab="set2 case: Genotype",ylab="methylation",notch=T)
		dev.off()
	}

	table(clin1$caco[j12])
	table(clin2$caco[j22])
}

########################################################################
########################################################################
########################################################################
## AML: Clinical table

nProbe=-1
nProbe=101

idNoRefactor=c("X1288G","X1466G","X0153G","X0077G","X0201G","X1191G","X1264G","X1219G","X0873G","X0927G","X1075G","X1766G","X0691G")

dirClin=dirMeth="data/aml/"
fNameMeth=paste("beta_bmiq_aml.txt",sep="")
fNameClin=paste("clin_aml_20150114.txt",sep="")
fNameClin=paste("clin_aml_20151002_2.txt",sep="")
clin=read.table(paste(dirClin,fNameClin,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
meth=read.table(paste(dirMeth,fNameMeth,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
clin$id=paste("X",clin$guthrieId,sep="")
j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
if (fNameClin=="clin_aml_20150114.txt") {
	clin$include="yes"
	clin$include[-j2]="excludeByFunNorm"
	clin$include[which(clin$id%in%idNoRefactor)]="excludeByRefactor"
	load("R_est_init_aml.RData")
	dat=matrix(nrow=nrow(clin),ncol=ncol(R_est$ind$coord))
	colnames(dat)=paste("prinComp",1:ncol(dat),sep="")
	dat[match(rownames(R_est$ind$coord),clin$id),]=R_est$ind$coord
	clin=cbind(clin,dat)
	write.table(clin,file=paste("clin_aml_20151002.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F,append=F)
}
meth=meth[j1,]
clin=clin[j2,]
clin0=clin

for (subsetFlag in c("","excludedByRefactor")) {
	cat("\n\n#####################",subsetFlag,"################\n")
	if (subsetFlag=="excludedByRefactor") {
		## Excluded by refactor
		j=which(!clin0$id%in%idNoRefactor)
		clin2=clin0[j,]
	} else {
		clin2=clin0
	}

	tbl=NULL
	varName=c("caco","sex","ch_ethnicity","include","downSyndrome")
	for (k1 in 1:(length(varName)-1)) {
		for (k2 in (k1+1):length(varName)) {
			clin=clin2
			j=!is.na(clin[,varName[k1]]) & !is.na(clin[,varName[k2]])
			x=table(clin[j,varName[k1]],clin[j,varName[k2]],dnn=varName[c(k1,k2)])
			print(x)
			if (nrow(x)>1 & ncol(x)>1) {
				tbl1=c("set2",varName[k1],varName[k2],"Chi square test",chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,"")
				tbl=rbind(tbl,tbl1)
				cat("\nChi-square test p-value: ",signif(chisq.test(clin[,varName[k1]],clin[,varName[k2]])$p.value,5),"\n",sep="")
			}
			cat("\n\n==============================\n\n")
		}
	}


	varName1=c("caco","sex","ch_ethnicity","include","downSyndrome")
	varName2=c("gestage",paste("prinComp",1:6,sep=""))
	for (k1 in 1:length(varName1)) {
		for (k2 in 1:length(varName2)) {
			cat("\n\n================",varName1[k1]," and ",varName2[k2],"================\n")
			clin=clin2
			x=table(clin[,varName1[k1]])
			if (length(x)==2) {
				cat("================ Wilcoxon test \n")
				clin=clin2
				x=table(clin[,varName1[k1]])
				if (length(x)>1) {
					tbl1=c("set2",varName1[k1],varName2[k2],"Wilcoxon test",wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value,"")
					tbl=rbind(tbl,tbl1)
					ttl="P-value "
					cat(ttl,signif(wilcox.test(clin[,varName2[k2]]~clin[,varName1[k1]])$p.value,5),"\n")
				}
			} else if (length(x)>2) {
				cat("================ Kruskal-Wallis test \n")
				clin=clin2
				x=table(clin[,varName1[k1]])
				if (length(x)>1) {
					tbl1=c("set2",varName1[k1],varName2[k2],"Kruskal-Wallis test",kruskal.test(clin[,varName2[k2]]~as.factor(clin[,varName1[k1]]))$p.value,"")
					tbl=rbind(tbl,tbl1)
					ttl="P-value "
					cat(ttl,signif(kruskal.test(clin[,varName2[k2]]~as.factor(clin[,varName1[k1]]))$p.value,5),"\n")
				}
			}
		}
	}
}

clin=clin0
grp=as.character(clin$downSyndrome)
grp[which(grp=="0")]="noDownSyndrome"
grp[which(grp=="1")]="downSyndrome"
grp=paste(grp,clin$include,sep="/")
grp=sub("yes","include",grp)
grpUniq=sort(unique(grp))
for (gId in 1:length(grpUniq)) {
	cat("\n===========",grpUniq[gId],"===========\n")
	for (k in which(names(clin)%in%paste("prinComp",1:2,sep=""))) {
		cat("summary(abs(",names(clin)[k],"))\n",sep="")
		print(summary(abs(clin[which(grp==grpUniq[gId]),k])))
	}
}

########################################################################
########################################################################
########################################################################

load("data.RData")
phen=pdata
phen$subjectId=phen$subjectID

# 661724 - T-ALL, 328606 - T-ALL/B-ALL
x=c("661724","328606")
table(x%in%phen$subjectId)

########################################################################
########################################################################
########################################################################
## QC analysis

library(FactoMineR)

## ---------------

dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
dirSrc3=paste("code/",sep="")
dirSrc3=paste("/Users/royr/Downloads/wiemelsJ-all/",sep="")
setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))

##############################################

computerFlag="cluster"
computerFlag=""

##############################################

setFlag="set2"
setFlag=""

subsetFlag="propOfCacoEthnAsInSet2"
subsetFlag="case"
subsetFlag="ctrl"
subsetFlag=""


colList=c("black","red")
colList=c("skyblue","blue","yellow","purple","black","red","orange","green","cyan","darkgreen")

##############################################

source(paste(dirSrc3,"funcs.R",sep=""))
res=getClinData(setFlag=setFlag,subsetFlag=subsetFlag)
clin1=res$clin1
clin2=res$clin2
datadir11=res$dirClin1
datadir21=res$dirClin2
datadir12=res$dirMeth1
datadir22=res$dirMeth2
rm(res)

clin1$id=paste("X",clin1$guthrieId,sep="")
clin2$id=paste("X",clin2$guthrieId,sep="")
clin1$beadPos=paste(clin1$Beadchip,"_",clin1$Position,sep="")
clin2$beadPos=paste(clin2$Beadchip,"_",clin2$Position,sep="")

## ---------------
datadir="docs/HelenHansen/"
qcX=read.table(paste(datadir,"SET 1& 2 BAC control probes_GC samples_xo.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
qcC=read.table(paste(datadir,"SET 1& 2 BAC control probes_GC samples_countOfFlags.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
qcR=read.table(paste(datadir,"SET 1& 2 BAC control probes_GC samples_rawBacScores.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
qcThres=read.table(paste(datadir,"SET 1& 2 BAC control probes_GC samples_bacThresholds.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

nameInfo=data.frame(name1=c(names(qcX),names(qcC),names(qcR)),stringsAsFactors=F)

names(qcX)[match(c("Sample.Name","Sentrix.Barcode","Sentrix.Position","Any.non.background.bilufite.conversion.flag","Any.Bisulfite.Conversion..Flag","Any.Flag"),names(qcX))]=
    c("sampleId","Beadchip","Position","convNBFlag","convFlag","anyFlag")
names(qcC)[match(c("Sample.Name","Sentrix.Barcode","Sentrix.Position","Staining.Green","Staining.Red","Extension.Green","Extension.Red","Hybridization.High.Medium","HybridizationMedium.Low","Target.Removal1","Target.Removal2","Specificity.1.Green","Specificity.1.Red","Specificity.2","Specificity.2.Background","NonPolymorphic.Green","NonPolymorphic.Red","Bisulfite.Conversion.1.Green","Bisulfite.Conversion.1.Red","Bisulfite.Conversion.2","BisulfiteConversion.1.Background.Green","Bisulfite.Conversion.1.Background.Red","Bisulfite.Conversion.2.Background","Any.non.background.bilufite.conversion.flag","Any.Bisulfite.Conversion..Flag","Any.Flag"),names(qcC))]=
    c("sampleId","Beadchip","Position","stainG","stainR","extG","extR","hybHiMed","hybMedLo","tarRem1","tarRem2","spec1G","spec1R","spec2","spec2Bgd","nonPolG","nonPolR","conv1G","conv1R","conv2","conv1BgdG","conv1BgdR","conv2Bgd","convNBFlag","anyConvFlag","anyFlag")
names(qcR)[match(c("Sample.Name","Sentrix.Barcode","Sentrix.Position","Restoration","Staining.Green","Staining.Red","Extension.Green","Extension.Red","Hybridization.High.Medium","HybridizationMedium.Low","Target.Removal1","Target.Removal2","Bisulfite.Conversion.1.Green","BisulfiteConversion.1.Background.Green","Bisulfite.Conversion.1.Red","Bisulfite.Conversion.1.Background.Red","Bisulfite.Conversion.2","Bisulfite.Conversion.2.Background","Specificity.1.Green","Specificity.1.Red","Specificity.2","Specificity.2.Background","NonPolymorphicGreen","NonPolymorphic.Red","sample.group","X450k.EPIC","DNA.source","Restored..Y.N.","conversion.lab","Illumina.Core","Approx.date.of.project","Flagged.by.MethylLight..only.USC.","P..0.05","Excluded.for.QC.by.RR","Illumina.known.hyb.probe..labeling.error","Total.ng.DNA.for.bis.conv","Concentrated"),names(qcR))]=
    c("sampleId","Beadchip","Position","restore","stainG","stainR","extG","extR","hybHiMed","hybMedLo","tarRem1","tarRem2","conv1G","conv1BgdG","conv1R","conv1BgdR","conv2","conv2Bgd","spec1G","spec1R","spec2","spec2Bgd","nonPolG","nonPolR","group","chipType","dnaSrc","restored","convLab","ilmnCore","projDate","methylLightFlag","pv.05","qcExclRR","IlmnPrLabelErr","totDna","conc")
names(qcThres)[match(c("X","BAC.cut.off.for.passing.QC"),names(qcThres))]=c("name","thres")
qcThres$id[match(c("Staining Green","Staining Red","Extension Green","Extension Red","Hybridization High Medium","HybridizationMedium Low","Target Removal1","Target Removal2","Bisulfite Conversion 1 Green","BisulfiteConversion 1 Background Green","Bisulfite Conversion 1 Red","Bisulfite Conversion 1 Background Red","Bisulfite Conversion 2","Bisulfite Conversion 2 Background","Specificity 1 Green","Specificity 1 Red","Specificity 2","Specificity 2 Background","NonPolymorphicGreen","NonPolymorphic Red"),qcThres$name)]=
    c("stainG","stainR","extG","extR","hybHiMed","hybMedLo","tarRem1","tarRem2","conv1G","conv1BgdG","conv1R","conv1BgdR","conv2","conv2Bgd","spec1G","spec1R","spec2","spec2Bgd","nonPolG","nonPolR")
qcThres$thres=as.numeric(sub(">","",qcThres$thres))

nameInfo$name2=c(names(qcX),names(qcC),names(qcR))
nameInfo=nameInfo[!duplicated(nameInfo$name1),]

qcC=qcC[match(qcX$sampleId,qcC$sampleId),]
qcR=qcR[match(qcX$sampleId,qcR$sampleId),]

if (F) {
    qcComb=data.frame(sampleId=qcX$sampleId,stringsAsFactors=F)
    tbl=qcR[,which(names(qcR)%in%qcThres$id)]
    for (k2 in 1:ncol(tbl)) {
        k1=match(names(tbl)[k2],qcThres$id)
        tbl[,k2]=as.integer(tbl[,k2]>qcThres$thres[k1])
    }
    tbl=qcC[,which(names(qcC)%in%qcThres$id)]
    for (k2 in names(tbl)) {
        cat("\n\n==============",k2,"\n")
        k1=match(k2,qcThres$id)
        tbl[,k2]=as.integer(qcR[,k2]>qcThres$thres[k1])
        print(table(qcC=qcC[,k2],qcR=tbl[,k2],exclude=NULL))
    }
    qcComb=cbind(sampleId=qcX$sampleId,tbl)
}

clin1=clin1[which(clin1$beadPos%in%qcX$sampleId),]
clin2=clin2[which(clin2$beadPos%in%qcX$sampleId),]

## ---------------

nProbe=10001
nProbe=101
nProbe=-1

## ---------------
## Run this if necessary

cohortList=c("_set1","_set2")
for (cohortFlag in cohortList) {
    if (nProbe==(-1)) {
        load(paste("beta",cohortFlag,".RData",sep=""))
    } else {
        cohortFlag="_set1"
        datadir=paste("docs/all/",sub("_","",cohortFlag),"/",sep="")
        fName=paste("beta_funNorm",cohortFlag,".txt",sep="")
        meth1=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        
        cohortFlag="_set2"
        datadir=paste("docs/all/",sub("_","",cohortFlag),"/",sep="")
        fName=paste("beta_funNorm",cohortFlag,".txt",sep="")
        meth2=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
    }
    switch(cohortFlag,
    "_set1"={
        rownames(meth1)=meth1$probeId
        j=match(clin1$id,colnames(meth1)); j1=which(!is.na(j)); j2=j[j1]
        betaThis=t(meth1[,j2])
        rm(meth1)
    },
    "_set2"={
        rownames(meth2)=meth2$probeId
        j=match(clin2$id,colnames(meth2)); j1=which(!is.na(j)); j2=j[j1]
        betaThis=t(meth2[,j2])
        rm(meth2)
    }
    )
    #cpgId=1:nrow(betaThis)
    #samId=1:ncol(betaThis)
    #betaThis=t(betaThis)
    cpgId=apply(betaThis,2,function(x) mean(!is.na(x))!=0)
    samId=apply(betaThis,1,function(x) mean(!is.na(x))!=0)
    betaThis=betaThis[samId,cpgId]
    
    cpgInfo=as.data.frame(t(apply(betaThis,2,function(x) {
        meanThis=mean(x,na.rm=T)
        sdThis=sd(x,na.rm=T)
        y=quantile(x,na.rm=T)
        zeroThis=mean(x==0,na.rm=T)
        oneThis=mean(x==1,na.rm=T)
        c(meanThis,sdThis,zeroThis,oneThis,y)
    })),stringsAsFactors=F)
    names(cpgInfo)=c("mean","sd","mean0","mean1",paste("perc",c(0,25,50,75,100),sep=""))
    cpgInfo=cbind(cpgId=colnames(betaThis),cpgInfo)
    save(cpgInfo,file=paste("cpgInfo",cohortFlag,".RData",sep=""))
    
    fitPCA=PCA(betaThis,graph=F,ncp=6)
    save(fitPCA,file=paste("fitPCA",cohortFlag,".RData",sep=""))
}


## ---------------
## Plots of some methylated & demethyated CpGs

qcFlag="_xo"
varList=paste(c("anyFlag"),"_xo",sep="")
varId=1
varList=paste(c("convNBFlag","convFlag","anyFlag"),"_xo",sep="")

qcFlag="_rpart"
varList=c("outlier","outlierPred")

cohortList=c("_set1","_set2")
for (cohortFlag in cohortList) {
    if (nProbe==(-1)) {
        load(paste("beta",cohortFlag,".RData",sep=""))
    } else {
        cohortFlag="_set1"
        datadir=paste("docs/all/",sub("_","",cohortFlag),"/",sep="")
        fName=paste("beta_funNorm",cohortFlag,".txt",sep="")
        meth1=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        
        cohortFlag="_set2"
        datadir=paste("docs/all/",sub("_","",cohortFlag),"/",sep="")
        fName=paste("beta_funNorm",cohortFlag,".txt",sep="")
        meth2=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
    }
    switch(cohortFlag,
    "_set1"={
        rownames(meth1)=meth1$probeId
        j=match(clin1$id,colnames(meth1)); j1=which(!is.na(j)); j2=j[j1]
        betaThis=t(meth1[,j2])
        rm(meth1)
    },
    "_set2"={
        rownames(meth2)=meth2$probeId
        j=match(clin2$id,colnames(meth2)); j1=which(!is.na(j)); j2=j[j1]
        betaThis=t(meth2[,j2])
        rm(meth2)
    }
    )
    cpgId=apply(betaThis,2,function(x) mean(!is.na(x))!=0)
    samId=apply(betaThis,1,function(x) mean(!is.na(x))!=0)
    betaThis=betaThis[samId,cpgId]
    load(file=paste("cpgInfo",cohortFlag,".RData",sep=""))
    switch(cohortFlag,
    "_set1"={
        clin=clin1
        dat=dat12
    },
    "_set2"={
        clin=clin2
        dat=dat22
    }
    )
    cohortName=capWords(sub("_","",cohortFlag))
    tbl2=qcX
    names(tbl2)=paste(names(tbl2),"_xo",sep="")
    tbl2=cbind(tbl2,qcC)
    tbl1=qcR
    names(tbl1)=paste(names(tbl1),"_raw",sep="")
    tbl2=cbind(tbl2,tbl1)
    if (qcFlag=="_rpart") {
        if (is.numeric(dat$classObs)) tmp=rep(NA,nrow(tbl2)) else tmp=rep("",nrow(tbl2))
        tbl1=data.frame(classObs=tmp,classPred=tmp,stringsAsFactors=F)
        j=match(dat$beadPos,tbl2$sampleId); j=j[!is.na(j)]
        for (k in names(tbl1)) {
            tbl1[j,k]=dat[,k]
        }
        names(tbl1)=c("outlier","outlierPred")
        tbl2=cbind(tbl2,tbl1)
    }
    nm=c(names(clin),varList)
    clin=cbind(clin,tbl2[match(clin$beadPos,tbl2$sampleId),varList])
    names(clin)=nm
    for (varId in 1:length(varList)) {
        if (is.character(clin[,varList[varId]])) {
            if (all(clin[,varList[varId]]!="0",na.rm=T)) {clin[is.na(clin[,varList[varId]]),varList[varId]]="0"}
            clin[which(clin[,varList[varId]]=="X"),varList[varId]]="1"
            clin[,varList[varId]]=as.integer(clin[,varList[varId]])
        } else {
            if (all(clin[,varList[varId]]!=0,na.rm=T)) {clin[is.na(clin[,varList[varId]]),varList[varId]]=0}
        }
    }
    clin=clin[match(rownames(betaThis),clin$id),]
    i1=order(cpgInfo$perc50,decreasing=F)[1:10]
    i2=order(cpgInfo$perc50,decreasing=T)[1:10]
    i1=which(cpgInfo$perc100<.1)[1:10]
    i2=which(cpgInfo$perc0>.9)[1:10]
    i=1
    offset=.1
    for (varId in 1:length(varList)) {
        parInfo=list(pch=20,cex=1,cex.axis=.8,las=3,cex.leg=.6)
        png(paste("scatterPlot",cohortFlag,"_",varList[varId],qcFlag,"_%1d.png",sep=""),width=3*240,height=2*240)
        par(mfrow=c(2,1))
        for (pId in 1:2) {
            if (pId==1) {
                header=paste(cohortName,": Demethylated CpGs",sep="")
                ii=i1
            } else {
                header=paste(cohortName,": Methylated CpGs",sep="")
                ii=i2
            }
            xlim=c(1,length(ii)+1)
            ylim=range(c(betaThis[,ii]),na.rm=T)
            i=1
            j1=which(!is.na(betaThis[,ii[i]]))
            y=betaThis[j1,ii[i]]
            x=rep(i,length(y))
            plot(x,y,xlim=xlim,ylim=ylim,main=header,xlab="",ylab="beta-value",pch=parInfo$pch,cex=parInfo$cex,cex.axis=parInfo$cex.axis,las=parInfo$las,type="n",xaxt="n")
            axis(side=1,at=xlim[1]:xlim[2],labels=c(colnames(betaThis)[ii],""),tick=F,cex.axis=parInfo$cex.axis,las=parInfo$las)
            grpUniq=sort(unique(clin[j1,varList[varId]]))
            if (pId==1) legend(x=length(ii)+.5,y=ylim[2],fill=colList[1:length(grpUniq)],legend=grpUniq,title=varList[varId],cex=parInfo$cex.leg)
            for (i in 1:length(ii)) {
                j1=which(!is.na(betaThis[,ii[i]]))
                y=betaThis[j1,ii[i]]
                x=rep(i,length(y))
                y0=round(y*100000)
                y1=unique(y0[duplicated(y0)])
                if (length(y1)!=0) {
                    for (j1 in 1:length(y1)) {
                        j=which(y0==y1[j1])
                        x2=(0:floor(length(j)/2))*offset
                        j2=j[1:length(x2)]
                        x[j2]=x[j2]+x2
                        if (length(j2)!=length(j)) {
                            j2=j[(length(j2)+1):length(j)]
                            x2=-(1:length(j2))*offset
                            x[j2]=x[j2]+x2
                        }
                    }
                }
                colVec=rep("black",nrow(clin))
                j1=which(!is.na(betaThis[,ii[i]]))
                grpUniq=sort(unique(clin[j1,varList[varId]]))
                for (gId in 1:length(grpUniq)) {
                    colVec[which(clin[j1,varList[varId]]==grpUniq[gId])]=colList[gId]
                }
                points(x,y,pch=parInfo$pch,cex=parInfo$cex,col=colVec)
                j=which(clin[j1,varList[varId]]!=grpUniq[1])
                points(x[j],y[j],pch=parInfo$pch,cex=parInfo$cex,col=colVec[j])
            }
        }
        dev.off()
    }
}

## ---------------

library(coin)

datadir="results/qc/"

qcFlag="_cntOfFlg"
varList=c()
for (k in 1:ncol(qcC)) {
    x=sum(!duplicated(qcC[,k]))
    if (x==2) varList=c(varList,names(qcC)[k])
    if (x>2 & x<21) {
        cat("\n\n==============",k,names(qcC)[k],x,"\n")
        print(table(qcC[,k],exclude=NULL))
    }
}
#varList=varList[1:2]

qcFlag="_xo"
varList=paste(c("convNBFlag","convFlag","anyFlag"),"_xo",sep="")

qcFlag="_rpart"
#varList=paste(c("classObs","classPred"),"_rpart",sep="")
varList=c("outlier","outlierPred")

cohortList=c("_set1","_set2")

xlim=ylim=c(-2000,1000)
xlim=ylim=NULL
load(file=paste(datadir,"fitPCA_set1.RData",sep=""))
xlim=range(fitPCA$ind$coord[,1],na.rm=T)
ylim=range(fitPCA$ind$coord[,2],na.rm=T)
load(file=paste(datadir,"fitPCA_set2.RData",sep=""))
xlim=range(c(xlim,fitPCA$ind$coord[,1]),na.rm=T)
ylim=range(c(ylim,fitPCA$ind$coord[,2]),na.rm=T)


png(paste("pcaScreePlot.png",sep=""),width=3*240,height=2*240)
par(mfrow=c(1,2))
for (cohortFlag in cohortList) {
    cohortName=capWords(sub("_","",cohortFlag))
    load(file=paste(datadir,"fitPCA",cohortFlag,".RData",sep=""))
    x=fitPCA$eig[1:20,"percentage of variance"]
    names(x)=paste("PC",1:length(x),sep="")
    barplot(x,ylim=c(0,17),main=cohortName,xlab="Principal component",ylab="Variance explained (%)",las=3)
}
dev.off()

tbl=NULL
for (cohortFlag in cohortList) {
    cohortName=capWords(sub("_","",cohortFlag))
    switch(cohortFlag,
        "_set1"={
            clin=clin1
            dat=dat12
        },
        "_set2"={
            clin=clin2
            dat=dat22
        }
    )
    tbl2=qcX
    names(tbl2)=paste(names(tbl2),"_xo",sep="")
    tbl2=cbind(tbl2,qcC)
    tbl1=qcR
    names(tbl1)=paste(names(tbl1),"_raw",sep="")
    tbl2=cbind(tbl2,tbl1)
    if (qcFlag=="_rpart") {
        if (is.numeric(dat$classObs)) tmp=rep(NA,nrow(tbl2)) else tmp=rep("",nrow(tbl2))
        tbl1=data.frame(classObs=tmp,classPred=tmp,stringsAsFactors=F)
        j=match(dat$beadPos,tbl2$sampleId); j=j[!is.na(j)]
        for (k in names(tbl1)) {
            tbl1[j,k]=dat[,k]
            #tbl1[-j,k]=NA
        }
        #names(tbl1)=paste(names(tbl1),"_rpart",sep="")
        names(tbl1)=c("outlier","outlierPred")
        tbl2=cbind(tbl2,tbl1)
    }
    nm=c(names(clin),varList)
    clin=cbind(clin,tbl2[match(clin$beadPos,tbl2$sampleId),varList])
    names(clin)=nm

    for (varId in 1:length(varList)) {
        if (is.character(clin[,varList[varId]])) {
            if (all(clin[,varList[varId]]!="0",na.rm=T)) {clin[is.na(clin[,varList[varId]]),varList[varId]]="0"}
            clin[which(clin[,varList[varId]]=="X"),varList[varId]]="1"
            clin[,varList[varId]]=as.integer(clin[,varList[varId]])
        } else {
            if (all(clin[,varList[varId]]!=0,na.rm=T)) {clin[is.na(clin[,varList[varId]]),varList[varId]]=0}
        }
    }
    load(file=paste(datadir,"fitPCA",cohortFlag,".RData",sep=""))
    
    #Refactor_dat=fitPCA$ind$coord[,c(1:6)]
    if (nrow(fitPCA$ind$coord)==nrow(clin)) {
        png(paste("pca",cohortFlag,qcFlag,"_%1d.png",sep=""),width=3*240,height=2*240)
        #png(paste("pca",cohortFlag,qcFlag,"_%1d.png",sep=""))
        par(mfrow=c(2,3))
        for (varId in 1:length(varList)) {
            y=table(clin[,varList[varId]])
            
            if (F) {
                par(mfrow=c(1,1))
                pv1=pv2=NA
            }

            if (T) {
                if (length(y)==2) {
                    testType="wilcox"
                    res=wilcox_test(fitPCA$ind$coord[,1]~as.factor(clin[,varList[varId]]),distribution="exact")
                    pv1=pvalue(res)
                    res=wilcox_test(fitPCA$ind$coord[,2]~as.factor(clin[,varList[varId]]),distribution="exact")
                    pv2=pvalue(res)
                } else if (length(y)>2) {
                    testType="kruskal"
                    res=kruskal_test(fitPCA$ind$coord[,1]~as.factor(clin[,varList[varId]]),distribution="exact")
                    pv1=pvalue(res)
                    res=kruskal_test(fitPCA$ind$coord[,2]~as.factor(clin[,varList[varId]]),distribution="exact")
                    pv2=pvalue(res)
                } else {
                    pv1=pv2=NA
                }
            }
            tbl2=c(sub("_","",cohortFlag),sub("_","",subsetFlag),varList[varId],sum(clin[,varList[varId]]==1,na.rm=T),pv1,pv2)
            tbl=rbind(tbl,tbl2)
            colVec=rep("black",nrow(clin))
            grpUniq=sort(unique(clin[,varList[varId]]))
            for (gId in 1:length(grpUniq)) {
                colVec[which(clin[,varList[varId]]==grpUniq[gId])]=colList[gId]
            }
            plot(fitPCA,xlim=xlim,ylim=ylim,label="none",col.ind=colVec,title=paste("PCA: ",cohortName,"\nPV: ",varList[varId]," vs PC1 ",signif(pv1,2),", vs PC2 ",signif(pv2,2),sep=""))
            abline(c(0,1),lty="dotted")
            abline(v=c(-500,500),lty="dotted")
            abline(h=c(-500,500),lty="dotted")
            legend(x=-2000,y=-500,fill=colList[1:length(grpUniq)],legend=grpUniq,title=varList[varId])
        }
        dev.off()
    } else {
        cat("No. of samples differ after PCA!!!\n")
    }
}
rownames(tbl)=NULL
tbl=as.data.frame(tbl,stringsAsFactors=F)
names(tbl)=c("dataset","subset","variable","numFlagged","pvPC1","pvPC2")
for (k in c("numFlagged","pvPC1","pvPC2")) tbl[,k]=as.numeric(tbl[,k])
tbl2=tbl
for (k in c("pvPC1","pvPC2")) {
    sig=rep("",nrow(tbl))
    sig[which(tbl[,k]<0.1)]="*"
    sig[which(tbl[,k]<0.05)]="**"
    sig[which(tbl[,k]<0.001)]="***"
    tbl2[,k]=paste(signif(tbl[,k],2)," ",sig,sep="")
}
colId=c("dataset","variable","varDescription","numFlagged","pvPC1","pvPC2")
tbl2$varDescription=nameInfo$name1[match(tbl2$variable,nameInfo$name2)]
tbl2=tbl2[,colId]
tbl2=cbind(tbl2[which(tbl2$dataset=="set1"),c("variable","varDescription","numFlagged","pvPC1","pvPC2")],tbl2[which(tbl2$dataset=="set2"),c("numFlagged","pvPC1","pvPC2")])
names(tbl2)=c("variable","varDescription",paste(c("numSampleFlagged","pvPC1","pvPC2"),rep(cohortList,each=3),sep=""))
write.table(tbl2,file=paste("pca_qcFlag",qcFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

## ---------------

########################################################################
########################################################################
########################################################################

if (F) {
    nProbe=-1

    cohortFlag="_set1"
    datadir=paste("data/",sub("_","",cohortFlag),"/",sep="")
    fName=paste("beta_funNorm",cohortFlag,".txt",sep="")
    meth1=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
    save(meth1,file="beta_set1.RData")
    rm(meth1)

    cohortFlag="_set2"
    datadir=paste("data/",sub("_","",cohortFlag),"/",sep="")
    fName=paste("beta_funNorm",cohortFlag,".txt",sep="")
    meth2=read.table(paste(datadir,fName,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
    save(meth2,file="beta_set2.RData")
    rm(meth2)
}

########################################################################
########################################################################
########################################################################
## rpart
## Prepare input data for rpart then run rpart.R

library(rpart)

datadir="results/qc/"

varResp="classObs"
varSamId="beadPos"

qcFlag="_cntOfFlg"
varList=c()
for (k in 1:ncol(qcC)) {
    x=sum(!duplicated(qcC[,k]))
    if (x==2) varList=c(varList,names(qcC)[k])
    if (x>2 & x<21) {
        cat("\n\n==============",k,names(qcC)[k],x,"\n")
        print(table(qcC[,k],exclude=NULL))
    }
}

qcFlag="_xo"
varList=paste(c("convNBFlag","convFlag","anyFlag"),"_xo",sep="")

qcFlag="_raw"
varList=c()
for (k in 1:ncol(qcR)) {
    x=sum(!duplicated(qcR[,k]))
    if (x>1 & is.numeric(qcR[,k])) varList=c(varList,names(qcR)[k])
    if (x>1 & x<21 & is.character(qcR[,k])) {
        cat("\n\n==============",k,names(qcR)[k],x,"\n")
        print(table(qcR[,k],exclude=NULL))
    }
}
varList=c("stainG","stainR","extG","extR","hybHiMed","hybMedLo","tarRem1","tarRem2","conv1G","conv1BgdG","conv1R","conv1BgdR","conv2","conv2Bgd","spec1G","spec1R","spec2","spec2Bgd","nonPolG","nonPolR")
varList=paste(varList,"_raw",sep="")

thres=c(-500,500)

for (datType in c("_allGuthSet2","_allGuthSet1")) {
    cat("\n\n=================== ",datType," ==================\n\n",sep="")
    dat2Type=tolower(sub("allGuth","",datType))
    load(file=paste(datadir,"fitPCA",dat2Type,".RData",sep=""))
    switch(datType,
    "_allGuthSet2"={
        clin=clin2
    },
    "_allGuthSet1"={
        clin=clin1
    }
    )
    tbl2=qcX
    names(tbl2)=paste(names(tbl2),"_xo",sep="")
    tbl2=cbind(tbl2,qcC)
    tbl1=qcR
    names(tbl1)=paste(names(tbl1),"_raw",sep="")
    tbl2=cbind(tbl2,tbl1)
    nm=c(names(clin),varList)
    clin=cbind(clin,tbl2[match(clin$beadPos,tbl2$sampleId),varList])
    names(clin)=nm
    for (varId in 1:length(varList)) {
        if (is.character(clin[,varList[varId]])) {
            if (all(clin[,varList[varId]]!="0",na.rm=T)) {clin[is.na(clin[,varList[varId]]),varList[varId]]="0"}
            clin[which(clin[,varList[varId]]=="X"),varList[varId]]="1"
            clin[,varList[varId]]=as.integer(clin[,varList[varId]])
        } else {
            if (all(clin[,varList[varId]]!=0,na.rm=T)) {clin[is.na(clin[,varList[varId]]),varList[varId]]=0}
        }
    }
    clin$classObs=rep("0",nrow(clin))
    clin$classObs[which(fitPCA$ind$coord[,1]<thres[1] | fitPCA$ind$coord[,1]>thres[2] | fitPCA$ind$coord[,2]<thres[1] | fitPCA$ind$coord[,2]>thres[2])]="X"
    
    if (F) {
        clin$flagLev=rep("",nrow(clin))
        clin$flagLev[which(fitPCA$ind$coord[,1]<thres[1] | fitPCA$ind$coord[,1]>thres[2] | fitPCA$ind$coord[,2]<thres[1] | fitPCA$ind$coord[,2]>thres[2])]="99"
        clin$flagLev[which(fitPCA$ind$coord[,1]<thres[1] & fitPCA$ind$coord[,2]<thres[1])]="ll"
        clin$flagLev[which(fitPCA$ind$coord[,1]<thres[1] & fitPCA$ind$coord[,2]>=thres[1])]="ln"
    }
    
    j=match(clin$beadPos,tbl2$sampleId)
    classObs=clin$classObs
    classPred=tbl2$anyFlag_xo[j]
    classPred[which(classPred=="0")]="0"
    classPred[which(classPred=="X")]="X"
    cat("\nN: ",length(classObs),"\n",sep="")
    if (is.numeric(classPred)) {
        classPred=as.numeric(classPred)
        cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
    } else {
        print(table(classObs,classPred))
        cat("\nMisclassification rate: ",round(mean(classObs!=classPred,na.rm=T),2),sep="")
    }

    if (is.numeric(clin[,varResp])) tmp=rep(NA,nrow(clin)) else tmp=rep("",nrow(clin))
    dat=data.frame(clin[,c(varSamId,varList,varResp)],tmp,stringsAsFactors=F)
    names(dat)=c(varSamId,varList,"classObs","classPred")
    switch(datType,
    "_allGuthSet2"={
        dat22=dat
    },
    "_allGuthSet1"={
        dat12=dat
    }
    )
    j=1:nrow(clin)
    for (k in c(varResp,varList)) {
        j=j[which(!is.na(clin[j,k]))]
    }
    clin=clin[j,]
    for (k in c(varResp)) {
        clin[,k]=as.factor(clin[,k])
    }
    switch(datType,
    "_allGuthSet2"={
        dat2=clin[,c(varSamId,varResp,varList)]
    },
    "_allGuthSet1"={
        dat1=clin[,c(varSamId,varResp,varList)]
    }
    )
}

"
=================== _allGuthSet2 ==================


N: 451
classPred
classObs   0   X
0 363  63
X   9  16

Misclassification rate: 0.16

=================== _allGuthSet1 ==================


N: 479
classPred
classObs   0   X
0 363  97
X   6  13

Misclassification rate: 0.22

"
## ---------------------------------------------
## NOT USED

clin=dat2
j=1:nrow(clin)
for (k in c(varResp,varList)) {
    j=j[which(!is.na(clin[j,k]))]
}
for (k in c(varResp,varList)) {
    #clin[,k]=as.factor(clin[,k])
}
modelThis=as.formula(paste(varResp,"~",paste(varList,collapse="+"),sep=""))
ctrlThis=rpart.control(minsplit = 2, minbucket = 1, cp = 0.01,maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10,surrogatestyle = 0, maxdepth = 30)
ctrlThis=rpart.control(minsplit = 20, minbucket = round(20/3), cp = 0.01,maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10,surrogatestyle = 0, maxdepth = 30)
ctrlThis=rpart.control(minsplit = 5, minbucket = 2, cp = 0.01,maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10,surrogatestyle = 0, maxdepth = 30)
fit <- rpart(modelThis, data = clin[j,],method = "class",control=ctrlThis)
png("rpart_set2.png")
plot(fit)
par(xpd = TRUE)
text(fit, use.n = TRUE)
dev.off()


x="2) nonPolG< 0.5 434 15 0 (0.96543779 0.03456221)
4) spec1R< 0.5 433 14 0 (0.96766744 0.03233256)
8) conv2< 0.5 431 13 0 (0.96983759 0.03016241)
16) conv1G< 0.5 429 12 0 (0.97202797 0.02797203) *
17) conv1G>=0.5 2  1 0 (0.50000000 0.50000000)
34) stainG< 0.5 1  0 0 (1.00000000 0.00000000) *
35) stainG>=0.5 1  0 1 (0.00000000 1.00000000) *
9) conv2>=0.5 2  1 0 (0.50000000 0.50000000)
18) conv1G>=0.5 1  0 0 (1.00000000 0.00000000) *
19) conv1G< 0.5 1  0 1 (0.00000000 1.00000000) *
5) spec1R>=0.5 1  0 1 (0.00000000 1.00000000) *
3) nonPolG>=0.5 17  7 1 (0.41176471 0.58823529)
6) stainG>=0.5 2  0 0 (1.00000000 0.00000000) *
7) stainG< 0.5 15  5 1 (0.33333333 0.66666667)
14) conv1G>=0.5 1  0 0 (1.00000000 0.00000000) *
15) conv1G< 0.5 14  4 1 (0.28571429 0.71428571) *"
y=sapply(strsplit(x,"\n")[[1]],function(x) {
    y=sub("< ","<",gsub(" +"," ",x))
    y=strsplit(y," ")[[1]]
},USE.NAMES=F)


library(MASS)
fit <- glm(modelThis, family="binomial", data = clin[j,])
fit2 <- stepAIC(fit, k=log(nrow(clin[j,])), trace = FALSE)
fit2$anova
summary(glm(classObs ~ nonPolG + nonPolR, family="binomial", data = clin[j,]))



########################################################################
########################################################################
########################################################################
## NOT USED

kk=c()
for (k in 1:ncol(qcR)) if (!is.numeric(qcR[,k])) {
    kk=c(kk,k)
    cat(k,names(qcR)[k],"\n")
}
kk=c()
for (k in 1:ncol(qcR)) {
    if (is.numeric(qcR[,k])) {
        kk=c(kk,k)
        cat(k,names(qcR)[k],sum(is.na(qcR[,k])),"\n")
    }
}
names(qcR)[!names(qcR)%in%names(qcR)]
