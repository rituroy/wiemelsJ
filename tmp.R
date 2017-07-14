computerFlag="cluster"
computerFlag=""

if (computerFlag=="cluster") {
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
	setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

####################################################################
####################################################################

cohort="_aml"
cohort="_allGuthSet2"
cohort="_leuk"
cohort="_allGuthSet1Set2"
#cohort="_tcgaGbm"
cohort="_birthDefect"
cohort="_uscEpic"
cohort="_allGuthSet1"

subsetFlag="_ctrl"
subsetFlag="_case"
subsetFlag=""

subsetName=ifelse(subsetFlag=="",subsetFlag,paste(subsetFlag,"Subset",sep=""))

### R code from vignette source 'vignettes/minfi/inst/doc/minfi.Rnw'

options(width=70)
require(minfi)
#require(minfiData)
if (cohort=="_uscEpic") {
    library(IlluminaHumanMethylationEPICmanifest)
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
} else {
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
}

if (computerFlag=="cluster") {
	if (cohort=="_allGuthSet1") {
		baseDir="data/set1/idat/"
		targets=read.table("data/set1/final.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		clin2=read.table("data/set1/matches.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		
		names(clin2)[match(c("Array.ID","Old.match","New.match"),names(clin2))]=c("arrayId","guthrieIdO","guthrieIdN")
		clin2$Beadchip=substr(clin2$arrayId,1,10)
		clin2$Position=substr(clin2$arrayId,12,17)
		#clin2$guthrieIdO=paste("X",clin2$guthrieIdO,sep="")
		#clin2$guthrieIdN=paste("X",clin2$guthrieIdN,sep="")
		#clin2$type=substr(clin2$guthrieIdN,6,6)
		clin2$type=substr(clin2$guthrieIdN,5,5)
		clin2[!clin2$type%in%c("G","C"),]
		clin2$type[!clin2$type%in%c("G","C")]=NA
		
		targets$arrayId=paste(targets$Beadchip,targets$Position,sep="_")
		#targets$guthrieId=paste("X",targets$TargetID,sep="")
		targets$guthrieId=targets$TargetID
		
		targets$Beadchip=targets$Position=NA
		j=match(clin2$guthrieIdN,targets$guthrieId); j2=which(!is.na(j)); j1=j[j2]
		targets$Beadchip[j1]=clin2$Beadchip[j2]
		targets$Position[j1]=clin2$Position[j2]
        targets$caco=targets$Leukemia
    } else if (cohort=="_allGuthSet2") {
		baseDir="data/set2/idat/"
		#targets=read.table("data/set2/clinGuthrieReplJune2012.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		targets=read.table("data/set2/clin_guthrieSet2_20140619.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	} else if (cohort=="_allGuthSet1Set2") {
		baseDir="data/idat/"
		targets=read.table("data/set1set2/clin_guthrieSet1Set2_20140619.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		#targets$guthrieId=paste("X",targets$guthrieId,sep="")
	} else if (cohort=="_leuk") {
		baseDir="data/set1/idat/"

		clin=read.delim("data/set1/i.LEU.v2.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		names(clin)[match(c("Plate","Sentrix_ID","Sentrix_Position"),names(clin))]=c("Batch","Beadchip","Position")
		clin$id=sapply(clin$Sample,function(x) {if (is.na(as.integer(substr(x,1,1)))) x else paste("X",x,sep="")},USE.NAMES=F)
		clin$group="Leukemia"
		clin$group2="Leuk"
		
		tbl1=read.table("data/0708011 Sample_Sheet (Fetal blood).csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=14)
		names(tbl1)[match(c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position"),names(tbl1))]=c("id1","Sample_Well","Sample_Plate","group","Pool_ID","Beadchip","Position")
		tbl1$id=tbl1$id1
		j=which(!is.na(as.integer(tbl1$id1)))
		tbl1$group[j]=paste("S",tbl1$group[j],sep="")
		tbl1$group2=tbl1$group
		tbl1$group2[which(tbl1$group=="Allcell-B")]="B"
		tbl1$group2[which(tbl1$group=="Allcell-non-B")]="nB"
		tbl1$id=paste(tbl1$group2,"_",tbl1$id1,sep="")
		tbl1$sex=NA
		tbl1[tbl1$id%in%tbl1$id[duplicated(tbl1$id)],]
		
		names(clin)[names(clin)%in%names(tbl1)]
		#names(clin)[!names(clin)%in%names(tbl1)]
		#names(tbl1)[!names(tbl1)%in%names(clin)]
		
		k=match(names(clin),names(tbl1)); k1=which(!is.na(k)); k2=k[k1]
		targets=rbind(clin[,k1],tbl1[!duplicated(tbl1$id),k2])
		x=gsub("_+","_",gsub("(","_",gsub(" |-|)","_",targets$id),fixed=T))
		j=which(substr(x,nchar(x),nchar(x))=="_")
		if (length(j)!=0) x[j]=substr(x[j],1,nchar(x[j])-1)
		x[!(x%in%colnames(yMVals))]
		targets$id=x
		targets$Basename=paste("/",targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,sep="")

		if (F) {
			targets$tmp=paste(targets$Beadchip,"_",targets$Position,sep="")
			targets$tmp=paste(targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,"_Grn.idat",sep="")
			x=dir(path="data/set1/idat",pattern="_Grn.idat",recursive=T)
			table(targets$tmp%in%x)
			
			tbl1=read.table("data/0708011 Sample_Sheet (Fetal blood).csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=14)
			names(tbl1)[match(c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position"),names(tbl1))]=c("id1","Sample_Well","Sample_Plate","group","Pool_ID","Beadchip","Position")
			tbl1$id=tbl1$id1
			j=which(!is.na(as.integer(tbl1$id1)))
			tbl1$id[j]=paste("X",tbl1$id1[j],"_S",tbl1$group[j],sep="")
			
			targets=tbl1
			targets$Basename=paste("/",targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,sep="")

			targets$tmp=paste(targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,"_Grn.idat",sep="")
			x=dir(path="data/set1/idat",pattern="_Grn.idat",recursive=T)
			table(targets$tmp%in%x)
		}
    } else if (cohort=="_birthDefect") {
        baseDir="data/birthDefect/idat/"
        targets=read.table("data/birthDefect/clin_birthDefect_20160607.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    } else if (cohort=="_tcgaGbm") {
        baseDir="data/tcgaGbm/idat/"
        x=dir(baseDir)
        x=x[grep("_Grn",x)]
        targets=as.data.frame(t(sapply(x,function(x) {
            y=sub("_Grn.idat","",x,fixed=T)
            y=strsplit(y,"_")[[1]]
            c(paste("X",y[1],"_",y[2],sep=""),y)
        },USE.NAMES=F)),stringsAsFactors=F)
        names(targets)=c("id","Beadchip","Position")
        rm(x)
	} else {
		baseDir="data/aml/I169_DataFiles_part4_110311/"
		datadir="data/"
		#tbl1=read.table("data/List of 110 pilot ANB specimens to Joe 12-15-2011.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		#tbl2=read.table(paste(datadir,"AML-Sample-Layout.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
		#tbl3=read.table(paste(datadir,"aml/i.AML.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		tbl4=read.table(paste(datadir,"aml/11072011 Sample_Sheet (AML).csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=14)
		names(tbl4)[match(c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position"),names(tbl4))]=c("guthrieID","Sample_Well","Sample_Plate","sampleGroup","Pool_ID","Beadchip","Position")
		targets=read.table(paste(datadir,"aml/i.AML.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		names(targets)=c("guthrieId","caco","sex","subjectID")
		j=match(targets$guthrieID,tbl4$guthrieID)
		targets=cbind(targets,tbl4[,c("sampleGroup","Beadchip","Position")])
	}
} else {
	if (cohort=="_allGuthSet1") {
        baseDir="docs/all/set1/idat/"
        targets=read.table("docs/all/set1/final.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        clin2=read.table("docs/all/set1/matches.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        
        names(clin2)[match(c("Array.ID","Old.match","New.match"),names(clin2))]=c("arrayId","guthrieIdO","guthrieIdN")
        clin2$Beadchip=substr(clin2$arrayId,1,10)
        clin2$Position=substr(clin2$arrayId,12,17)
        #clin2$guthrieIdO=paste("X",clin2$guthrieIdO,sep="")
        #clin2$guthrieIdN=paste("X",clin2$guthrieIdN,sep="")
        #clin2$type=substr(clin2$guthrieIdN,6,6)
        clin2$type=substr(clin2$guthrieIdN,5,5)
        clin2[!clin2$type%in%c("G","C"),]
        clin2$type[!clin2$type%in%c("G","C")]=NA
        
        targets$arrayId=paste(targets$Beadchip,targets$Position,sep="_")
        #targets$guthrieId=paste("X",targets$TargetID,sep="")
        targets$guthrieId=targets$TargetID
        
        targets$Beadchip=targets$Position=NA
        j=match(clin2$guthrieIdN,targets$guthrieId); j2=which(!is.na(j)); j1=j[j2]
        targets$Beadchip[j1]=clin2$Beadchip[j2]
        targets$Position[j1]=clin2$Position[j2]
        targets$caco=targets$Leukemia
	} else if (cohort=="_allGuthSet2") {
		baseDir="docs/all/set2/idat/"
		targets=read.table("docs/all/set2/clinGuthrieReplJune2012.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    } else if (cohort=="_uscEpic") {
        baseDir="docs/uscEpic/idat/"
        fileList=list.files(baseDir)
        targets=read.table("docs/uscEpic/1436 (UCSF-16).txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(targets)[match(c("Plate","Bead.chip..","Well.Position","Well.Position.1","Sentrix.ID","Terminus","Complete.Barcode","Tube.Label","X","Volume.For.Bis","Volume.water.to.45","ALU.C4..HB.313..Ct.value","CONV.100...HB.365..Ct.value","CONV.50...HB.382..Ct.value","CONV.0...HB.368..Ct.value","Detected.CpG..0.05.","Percent.of.Probes..P.0.05."),names(targets))]=
            c("plate","group","wellPosition","wellPosition2","Beadchip","Position","id","tubeLabel","X","volBis","volWater45","aluC4","conv100","conv50","conv0","detCpG","probePerc")
        out=as.data.frame(t(sapply(fileList[grep("_Grn",fileList)],function(x) {
            y=strsplit(x,"_")[[1]]
            c(y[1:2])
        },USE.NAMES=F)),stringsAsFactors=F)
        names(out)=c("Beadchip","Position")
        targets$id=paste(targets$Beadchip,"_",targets$Position,sep="")
	} else {
	}
}
#tbl=targets[,which(!names(targets)%in%c("arrayId","TargetID","Leukemia","Position1","Bead_Position","Blinded_Study_ID"))]
#write.table(tbl,file=paste("clinGuthrieSet1Orig.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

#list.files(baseDir)
#list.files(file.path(baseDir, "6042308046"))


#targets=read.450k.sheet(baseDir)
if (!cohort%in%c("_tcgaGbm","_uscEpic")) {
    targets=targets[targets$Beadchip%in%list.files(baseDir),]
}

## -------------------------------------------------
if (F) {

	baseDir="data/set1/idat/"
	targets=read.table("data/set1/final.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	clin2=read.table("data/set1/matches.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
	meth=read.table("data/set1/d.GUTH.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=5)
	meth=meth[,-1]
	methL=read.table("data/set1/d.LEU.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=5)
	methL=methL[,-1]

	names(clin2)[match(c("Array.ID","Old.match","New.match"),names(clin2))]=c("arrayId","guthrieIdO","guthrieIdN")
	clin2$Beadchip=substr(clin2$arrayId,1,10)
	clin2$Position=substr(clin2$arrayId,12,17)
	clin2$guthrieIdO=paste("X",clin2$guthrieIdO,sep="")
	clin2$guthrieIdN=paste("X",clin2$guthrieIdN,sep="")
	clin2$type=substr(clin2$guthrieIdN,6,6)
	clin2[!clin2$type%in%c("G","C"),]
	clin2$type[!clin2$type%in%c("G","C")]=NA

	targets$arrayId=paste(targets$Beadchip,targets$Position,sep="_")
	targets$guthrieId=paste("X",targets$TargetID,sep="")

	targets$Beadchip=targets$Position=NA
	j=match(clin2$guthrieIdN,targets$guthrieId); j2=which(!is.na(j)); j1=j[j2]
	targets$Beadchip[j1]=clin2$Beadchip[j2]
	targets$Position[j1]=clin2$Position[j2]

	targets=targets[j1,]

	targets=targets[match(names(meth),targets$guthrieId),]

	x=sort(unique(as.character(targets$Beadchip)))
	y=sort(list.files(baseDir))
	x[!x%in%y]
	y[!y%in%x]

	j=which(clin2$type=="G")
	x=sort(unique(as.character(clin2$Beadchip[j])))
	y=sort(list.files(baseDir))
	x[!x%in%y]
	y[!y%in%x]
	
}

## -------------------------------------------------

#sub(baseDir, "", targets$Basename)
targets$Basename=paste("/",targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,sep="")
if (cohort=="_allGuthSet1") {
	targets$caco=targets$Leukemia
	targets$cacoDesc=""
	targets$cacoDesc[which(targets$caco==0)]="ctrl"
	targets$cacoDesc[which(targets$caco==1)]="case"
} else if (cohort=="_allGuthSet2") {
	#targets$cacoDesc=as.character(targets$caco)
	targets$cacoDesc="bm"
	targets$cacoDesc[which(targets$caco==1)]="ctrl"
	targets$cacoDesc[which(targets$caco==2)]="case"
} else if (cohort=="_allGuthSet1Set2") {
	targets$cacoDesc=""
	targets$cacoDesc[which(targets$caco==0)]="ctrl"
	targets$cacoDesc[which(targets$caco==1)]="case"
} else if (cohort=="_leuk") {
	targets$caco=targets$group2
	targets$cacoDesc=targets$group
} else if (cohort=="_birthDefect") {
    targets$caco=targets$status
    targets$cacoDesc=targets$status
    targets$cacoDesc[which(targets$caco=="")]="miscellaneous"
} else if (cohort=="_uscEpic") {
    targets$caco=targets$group
    targets$cacoDesc=paste("group",targets$group)
} else if (cohort%in%c("_tcgaGbm")) {
    targets$caco=0
    targets$cacoDesc=""
} else {
	targets$cacoDesc=""
	targets$cacoDesc[which(targets$caco==1)]="ctrl"
	targets$cacoDesc[which(targets$caco==2)]="case"
}

#rgsetRaw=read.450k.exp(base=baseDir, targets=targets)
#rgsetRaw=read.450k.exp(base=baseDir, recursive=T)
#save(rgsetRaw,file="rgsetRaw_read.450k.exp.RData")

rgsetRaw=read.metharray.exp(base=baseDir, recursive=T)
#save(rgsetRaw,file="rgsetRaw_read.metharray.exp.RData")

if (cohort%in%c("_leuk","_birthDefect","_tcgaGbm","_uscEpic")) {
	k=which(names(targets)%in%c("id","Beadchip","Position"))
} else {
	k=which(names(targets)%in%c("guthrieId","Beadchip","Position"))
}
targets[!paste(targets$Beadchip,"_",targets$Position,sep="")%in%sampleNames(rgsetRaw),k]
sampleNames(rgsetRaw)[grep("6057833044",sampleNames(rgsetRaw))]

j=match(sampleNames(rgsetRaw),paste(targets$Beadchip,"_",targets$Position,sep="")); j1=which(!is.na(j)); j2=j[j1]
rgsetRaw=rgsetRaw[,j1]
targets=targets[j2,]
save(rgsetRaw,file=paste("rgsetRaw",cohort,".RData",sep=""))

load(paste("rgsetRaw",cohort,".RData",sep=""))
#targets=targets[]
j=match(sampleNames(rgsetRaw),paste(targets$Beadchip,"_",targets$Position,sep="")); j1=which(!is.na(j)); j2=j[j1]
rgsetRaw=rgsetRaw[,j1]
targets=targets[j2,]

mani=getManifest(rgsetRaw)

rgsetRaw

if (!cohort%in%c("_leuk","_birthDefect","_tcgaGbm","_uscEpic")) {
	j=which(is.na(targets$guthrieId))
	targets$id=paste("X",targets$guthrieId,sep="")
	if (length(j)!=0) {targets$id[j]=paste(targets$Beadchip,"_",targets$Position,sep="")[j]}
}
pData(rgsetRaw)=targets[match(sampleNames(rgsetRaw),paste(targets$Beadchip,"_",targets$Position,sep="")),]
sampleNames(rgsetRaw)=pData(rgsetRaw)$id

pd=pData(rgsetRaw)
#pd[1:5,1:4]

#rgsetRaw2=read.450k.exp(file.path(baseDir, "5723646052"))
#rgsetRaw3=read.450k.exp(baseDir, recursive=TRUE)

#targets2=read.csv(file.path(baseDir, "SampleSheet.csv"), stringsAsFactors=FALSE, skip=7)
#targets2

#targets2$Basename=file.path(baseDir, targets2$Sentrix_ID, paste0(targets2$Sentrix_ID, targets2$Sentrix_Position))

qcReport(rgsetRaw, sampNames=pd$id, sampGroups=pd$cacoDesc, pdf=paste("qcReport_rawData",cohort,".pdf",sep=""))

###################################################

j=1:ncol(rgsetRaw)
switch(subsetFlag,
"_ctrl"={
    j=which(targets$caco==0)
},
"_case"={
    j=which(targets$caco==1)
}
)
rgsetRaw=rgsetRaw[,j]
targets=targets[j,]
save.image("tmp.RData")

rgsetFN=preprocessFunnorm(rgsetRaw)

save(rgsetFN,file=paste("rgsetFN",cohort,subsetName,".RData",sep=""))
#load(paste("rgsetFN",cohort,subsetName,".RData",sep=""))

if (cohort%in%c("_leuk","_birthDefect")) {
	colId=c("id","sex","predictedSex","Beadchip","Position")
	k="id"
} else if (cohort%in%c("_tcgaGbm")) {
    colId=c("id","predictedSex","Beadchip","Position")
    k="id"
} else if (cohort%in%c("_uscEpic")) {
    colId=c("id","group","predictedSex","Beadchip","Position")
    k="id"
} else {
	#colId=c("labId","sex","predictedSex","Beadchip","Position")]
	colId=c("guthrieId","sex","predictedSex","Beadchip","Position")
	k="guthrieId"
}

table(observedSex=colData(rgsetFN)$sex, predictedSex=colData(rgsetFN)$predictedSex,exclude=NULL)
if (cohort=="_birthDefect") {
    j=which(colData(rgsetFN)$sex!=colData(rgsetFN)$predictedSex)
} else if (cohort%in%c("_tcgaGbm","_uscEpic")) {
    j=c()
} else {
    table(colData(rgsetFN)$sex,2-as.integer(as.factor(colData(rgsetFN)$predictedSex)))
    colData(rgsetFN)[which(colData(rgsetFN)$sex!=(3-as.integer(as.factor(colData(rgsetFN)$predictedSex)))),colId]
    j=which(is.na(colData(rgsetFN)$sex))
}
if (length(j)!=0) {
    tbl=as.data.frame(colData(rgsetFN))[j,c("id","sex","predictedSex")]
    rownames(tbl)=NULL
    print(tbl)
}
if (length(j)!=0) colData(rgsetFN)[j,colId]

colData(rgsetFN)[which(colData(rgsetFN)[,k]%in%c("X0508G","X1298G")),colId]

if (F) {
    tmpC=rep("",nrow(rgsetFN))
    prInfo=data.frame(Name=featureNames(rgsetFN),Type=tmpC,stringsAsFactors=F)
    x=getProbeInfo(mani, type=c("I"))
    i=match(prInfo$Name,x$Name); i1=which(!is.na(i)); i2=i[i1]
    prInfo$Type[i1]="I"
    x=getProbeInfo(mani, type=c("II"))
    i=match(prInfo$Name,x$Name); i1=which(!is.na(i)); i2=i[i1]
    prInfo$Type[i1]="II"

    par(mfrow=c(1,2))
    plotBetasByType(msetRaw[,1], main="Raw")
    plotBetasByType(getBeta(rgsetFN)[,1], probeTypes=prInfo,main="Funnorm")

    png(paste("mdsPlot",cohort,subsetName,".png",sep=""))
    par(mfrow=c(1,2))
    mdsPlot(getBeta(rgsetRaw), numPositions=1000, sampGroups=pd$cacoDesc, sampNames=pd$id,main="Raw")
    mdsPlot(getBeta(rgsetFN), numPositions=1000, sampGroups=pd$cacoDesc, sampNames=pd$id,main="Funnorm")
    dev.off()
}

###################################################

save.image("tmp1.RData")

load("tmp1.RData")

if (cohort=="_allGuthSet1") {
	j=which(!sampleNames(rgsetRaw)%in%c("X1339G","X1298G","X0508G"))
} else if (cohort=="_allGuthSet2") {
	j=which(!sampleNames(rgsetRaw)%in%c("X1762G","X0635G","X1588G","9702496164_R05C02","9702496164_R06C02"))
} else if (cohort=="_allGuthSet1Set2") {
	j=which(!sampleNames(rgsetRaw)%in%c("X1339G","X1298G","X0508G","X1762G","X0635G","X1588G","9702496164_R05C02","9702496164_R06C02"))
} else if (cohort%in%c("_birthDefect")) {
    j=which(!sampleNames(rgsetRaw)%in%c("X8216","HB0004"))
} else if (cohort%in%c("_leuk","_tcgaGbm","_uscEpic")) {
	j=1:ncol(rgsetRaw)
} else {
	j=which(!sampleNames(rgsetRaw)%in%c("X1095G"))
}
rgsetRawTmp=rgsetRaw[,j]
save(rgsetRaw,rgsetRawTmp,cohort,subsetName,file=paste("rgsetRawEtc",cohort,subsetName,".RData",sep=""))
rgsetFN=preprocessFunnorm(rgsetRawTmp)
rm(rgsetRawTmp)
#save(rgsetFN,file="tmp2.RData")
#rgsetFN=preprocessFunnorm(rgsetRaw[,j])
save(rgsetFN,file=paste("rgsetFN",cohort,subsetName,".RData",sep=""))

#load(paste("rgsetFN",cohort,subsetName,".RData",sep=""))

###################################################
detP=detectionP(rgsetRaw, type="m+u")
rm(rgsetRaw)
detP=detP[match(rownames(rgsetFN),rownames(detP)),match(colnames(rgsetFN),colnames(detP))]
write.table(cbind(probeId=rownames(detP),detP),file=paste("detP",cohort,subsetName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

dat=getBeta(rgsetFN)
rm(rgsetFN)
dat[detP>0.01]=NA
rm(detP)
dat[apply(dat,1,function(x) mean(is.na(x)))>0.15,]=NA
write.table(cbind(probeId=rownames(dat),dat),file=paste("beta_funNorm",cohort,subsetName,"_tmp.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

dat[,apply(dat,2,function(x) mean(is.na(x)))>0.15]=NA

beta=dat
rm(dat)
write.table(cbind(probeId=rownames(beta),beta),file=paste("beta_funNorm",cohort,subsetName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

###################################################
mDat = beta
rm(beta)
#mDat = ifelse(mDat<=0, 0.00001, mDat)
#mDat = ifelse(mDat>=1, 0.99999, mDat)
mDat = ifelse(mDat<0.00001, 0.00001, mDat)
mDat = ifelse(mDat>0.99999, 0.99999, mDat)
mDat = log2(mDat)-log2(1-mDat)
save(mDat,file=paste("mDat_funNorm",cohort,subsetName,".RData",sep=""))
tbl=cbind(probeId=rownames(mDat),mDat)
write.table(tbl,file=paste("mDat_funNorm",cohort,subsetName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
rm(tbl)

###################################################
for (datType in c("mDat","beta")) {
	datadir=paste("data/",cohort,"/",sep="")
	dat=read.table(paste(datadir,datType,"_funNorm",cohort,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
	save.image(paste("tmp_",datType,cohort,".RData",sep=""))
	load(paste("tmp_",datType,cohort,".RData",sep=""))
	rownames(dat)=dat$probeId
	dat=as.matrix(dat[,-1])
	if (cohort%in%c("_leuk","_birthDefect","_tcgaGbm","_uscEpic",)) {
		j=match(colnames(dat),targets$id); j1=which(!is.na(j)); j2=j[j1]
	} else {
		j=match(colnames(dat),paste("X",targets$guthrieId,sep="")); j1=which(!is.na(j)); j2=j[j1]
	}
	dat=dat[,j1]
	targets=targets[j2,]

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

	ann$keep = ann$snp==0 & !ann$CHR%in%c(23,24)

	i=match(rownames(dat),ann$IlmnID); i1=which(!is.na(i)); i2=i[i1]
	dat=dat[i1,j1]
	ann=ann[i2,]

	if (cohort=="_allGuthSet1") {
		targets$caco=targets$Leukemia
	}

	## ----------------------------------------------
	colList=c("red","blue","orange","yellow","green","cyan","brown","olivedrab","skyblue","pink","magenta","purple")
		
	if (cohort%in%c("_leuk","_birthDefect")) {
		varList=c("caco","sex","Beadchip","Position")
		subset2List=c("",paste("_",unique(targets$group2),sep=""))
    } else if (cohort%in%c("_tcgaGbm")) {
            varList=c("Beadchip","Position")
            subset2List=c("",paste("_",unique(targets$group2),sep=""))
    } else if (cohort%in%c("_uscEpic")) {
        varList=c("group","Beadchip","Position")
        subset2List=""
	} else {
		varList=c("group2","sex","Beadchip","Position")
		subset2List=c("","_ctrl","_case")
	}

	i=ann$keep
	for (subset2Flag in subset2List) {
		if (subset2Flag=="") {
			j=1:nrow(targets)
		} else {
		if (cohort=="_leuk") {
			j=which(paste("_",targets$group2,sep="")==subset2Flag)
		} else {
			if (subset2Flag=="_ctrl") {
				j=which(targets$caco==0)
			} else if (subset2Flag=="_case") {
				j=which(targets$caco==1)
			}
		}
		for (varId in 1:length(varList)) {
			x=unique(targets[j,varList[varId]]); x=x[!is.na(x)]
			if (length(x)>1) {
				png(paste("mdsPlot_",datType,"_",varList[varId],subset2Flag,cohort,".png",sep=""))
				mdsPlot(dat[i,j], sampNames=targets$id[j], sampGroups=targets[j,varList[varId]],main=paste(datType,", ",varList[varId],": MDS plot\n1000 most variable CpG positions",sep=""),pal=colList,legendNCol=3)
				dev.off()
			}
		}
	}
}

###################################################
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

ann$keep = ann$snp==0 & !ann$CHR%in%c(23,24)

ann1=ann
dat=mDat

i=match(rownames(dat),ann$IlmnID); i1=which(!is.na(i)); i2=i[i1]
dat=dat[i1,j1]
ann=ann[i2,]

sdVec=apply(mDat,1,sd,na.rm=T)
save.image("tmp_10.RData")

nPr=10000
nPr=30000
#dat2=mDat[1:nPr,]
dat2=mDat[order(sdVec,decreasing=T)[1:nPr],]
#colnames(dat2)=paste(targets$group,"_",targets$wellPosition,sep="")
colnames(dat2)=targets$wellPosition

header=paste(nPr," most variable CpGs",sep="")

png(paste("hclust_epic_",nPr,"cpg.png",sep=""))
plot(hclust(dist(t(dat2)),method="ward.D2"),main=header,xlab="")
dev.off()

library(FactoMineR)
#fit <- PCA(t(dat), graph = F, ncp =6)
fit <- PCA(t(dat2), graph = F, ncp =6)
limAll=NULL
png(paste("pca_epic_",nPr,"cpg.png",sep=""))
plot(fit,xlim=limAll,ylim=limAll,sub=paste("PCA: ",header,sep=""))
if (F) {
    j=which(targets$group==1)
    points(fit$ind$coord[j,1],fit$ind$coord[j,1],col="red")
    j=which(targets$group==2)
    points(fit$ind$coord[j,1],fit$ind$coord[j,1],col="green")
}
abline(c(0,1),lty="dotted")
dev.off()

####################################################################
####################################################################
## Estimate cell mixture effects

#setwd("/Users/royr/UCSF/MargaretWrensch/gbm2014")

library(minfi)

cohortIn="_allGuthSet2"; cohort="_allGuthSet2"
cohortIn="_allGuthSet1Set2"; cohort="_allGuthSet1Set2"
cohortIn="_allGuthSet1"; cohort="_allGuthSet1"
cohortIn="_allGuthSet1Set2"; cohort="_allGuthSet1Set2_ctrlSubset"

compCellType="Blood"; compCellTypeName=""; cellTypeVec=c("CD8T","CD4T","NK","Bcell","Mono","Gran")
compCellType="CordBlood"; compCellTypeName="_cordBlood"; cellTypeVec=c("nRBC","CD8T","CD4T","NK","Bcell","Mono","Gran") ## default

load(file=paste("rgsetFN",cohortIn,".RData",sep=""))
phen=pData(rgsetFN)
rm(rgsetFN)
load(file=paste("rgsetRaw",cohortIn,".RData",sep=""))
phen$id2=paste(phen$Beadchip,"_",phen$Position,sep="")
#rgsetRaw=rgsetRaw[,1:5]
j=match(sampleNames(rgsetRaw),phen$id2); j1=which(!is.na(j)); j2=j[j1]
if (any(!phen$id2%in%sampleNames(rgsetRaw))) cat("There are non-matching IDs !!!\n")
rgsetRaw=rgsetRaw[,j1]
phen=phen[j2,]
sampleNames(rgsetRaw)=phen$id

if (length(grep("_ctrlSubset",cohort))==1) {
	j=which(phen$caco==0)
	rgsetRaw=rgsetRaw[,j]
	phen=phen[j,]
}

counts <- estimateCellCounts(rgsetRaw,compositeCellType=compCellType,cellTypes=cellTypeVec,meanPlot=FALSE)

save(counts,file=paste("counts_tmp",compCellTypeName,cohort,".RData",sep=""))

write.table(cbind(id=rownames(counts),counts),file=paste("cellCount_minfi",compCellTypeName,cohort,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

####################################################################
####################################################################
