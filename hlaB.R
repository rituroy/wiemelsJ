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
## Run regression_refactor_part2.R first

datadir="docs/all/hlaB/"
probG=read.table(paste(datadir,"HLA-B Alleles_Probabilities.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
callG=read.table(paste(datadir,"HLA-B Alleles_Best Guess Calls.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(probG)[match("Subject.ID",names(probG))]="subjectId"
names(callG)[match("Columns",names(callG))]="subjectId"
callG=callG[match(probG$subjectId,callG$subjectId),match(names(probG),names(callG))]
annG=data.frame(id=names(probG)[-1],chr=6,stringsAsFactors=F)
phenG=data.frame(id=probG$subjectId,subjectId=probG$subjectId,stringsAsFactors=F)
probG=t(as.matrix(probG[,-1]))
callG=t(as.matrix(callG[,-1]))
rownames(probG)=rownames(callG)=annG$id
colnames(probG)=colnames(callG)=phenG$id

annC=read.table(paste(datadir,"annotation_res_concordance.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

candGeneInfo=data.frame(geneSym=c("HLA-B","AS","HLAB","B-4901"),stringsAsFactors=F)
x=ann[grep("HLA-B|AS|HLAB|B-4901",toupper(ann$UCSC_RefGene_Name)),]
x=ann[grep("HLA-B|AS|HLAB|B-4901",toupper(ann$geneSym)),]
x=ann[which(toupper(ann$geneSym)%in%toupper(candGeneInfo$geneSym)),]


callG2=probG
callG2[probG<.5]=0
callG2[probG>=.5 & probG<=1.5]=1
callG2[probG>1.5]=2
callG2[is.na(probG)]=NA
table(callG=c(callG),callG2=c(callG2),exclude=NULL)
i=which(c(callG)!=c(callG2))
png("tmp.png")
boxplot(c(probG)[i]~c(callG)[i])
dev.off()


## ----------------------------------------------

stat0=stat_hlabCo12
stat1=stat_hlabCa1
stat2=stat_hlabCa2

pThres=0.05
lim=range(c(stat0$coef_hlaBallele1v0,stat1$coef_hlaBallele1v0,stat2$coef_hlaBallele1v0,stat0$coef_hlaBallele2v0,stat1$coef_hlaBallele2v0,stat2$coef_hlaBallele2v0),na.rm=T)
lim=c(-5,5)
varName="Coefficient"
cexMain=3; cexLab=3; cexAxis=2
colList=c("green","cyan","red","yellow")
#png("coefPlot_hlaB.png",width=3*480,height=2*480)
#par(mfrow=c(2,3))
#varList=c("hlaBallele1v0","hlaBallele2v0")
png("coefPlot_hlaB.png",width=3*480,height=1*480)
par(mfrow=c(1,3))
varList=c("hlaBallele1v0")
par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(5, 6, 4, 2) + 0.1)
for (vId in 1:length(varList)) {
    colId=paste("coef_",varList[vId],sep="")
    colIdPV=paste("holm_",varList[vId],sep="")
    header=paste(varList[vId])
    for (compFlag in c("co12Vca1","co12Vca2","ca1Vca2")) {
        switch(compFlag,
        "co12Vca1"={
            stat_1=stat0; stat_2=stat1
            ttl=paste(c("Set1+Set2 ctrl","Set1 case"),": ",varName,sep="")
        },
        "co12Vca2"={
            stat_1=stat0; stat_2=stat2
            ttl=paste(c("Set1+Set2 ctrl","Set2 case"),": ",varName,sep="")
        },
        "ca1Vca2"={
            stat_1=stat1; stat_2=stat2
            ttl=paste(c("Set1 case","Set2 case"),": ",varName,sep="")
        },
        )
        i=match(paste(stat_1$cpgId,stat_1$gene_genotype),paste(stat_2$cpgId,stat_2$gene_genotype)); i1=which(!is.na(i)); i2=i[i1]
        stat_1=stat_1[i1,]; stat_2=stat_2[i2,]
        x=table(stat_1[,colIdPV]<pThres,stat_2[,colIdPV]<pThres)
        plot(stat_1[,colId],stat_2[,colId],xlim=lim,ylim=lim,main=header,xlab=ttl[1],ylab=ttl[2],cex.main=cexMain,cex.lab=cexLab,cex.axis=cexAxis)
        i=which(stat_1[,colIdPV]<pThres & stat_2[,colIdPV]>=pThres)
        points(stat_1[i,colId],stat_2[i,colId],col=colList[1])
        i=which(stat_1[,colIdPV]>=pThres & stat_2[,colIdPV]<pThres)
        points(stat_1[i,colId],stat_2[i,colId],col=colList[2])
        i=which(stat_1[,colIdPV]<pThres & stat_2[,colIdPV]<pThres)
        points(stat_1[i,colId],stat_2[i,colId],col=colList[3])
        #i=match(stat0$cpgId,stat1$cpgId); i0=which(!is.na(i)); i1=i[i1]
        #i=match(stat0$cpgId[i0],stat2$cpgId); i0=i0[which(!is.na(i))]; i1=i1[which(!is.na(i))]; i2=i[which(!is.na(i))]
        #i=which(stat0[i0,colIdPV]<pThres & stat1[i1,colIdPV]<pThres & stat2[i2,colIdPV]<pThres)
        #points(stat_1[i,colId],stat_2[i,colId],col=colList[4])
        abline(h=0,lty="dotted"); abline(v=0,lty="dotted"); abline(c(0,1),lty="dotted");
        #ttl2=c("x-axis","y-axis","both","all 3")
        ttl2=c("x-axis","y-axis","both")
        if (nrow(x)==ncol(x)) ttl2=paste(ttl2," (",c(x)[2:4],")",sep="")
        legend(x=lim[1],y=lim[2],title=paste("Holm<=",pThres,sep=""),legend=ttl2,fill=colList,cex=2)
    }
}
dev.off()


#i=match(stat0$cpgId,stat1$cpgId); i0=which(!is.na(i)); i1=i[i1]
#i=match(stat0$cpgId[i0],stat2$cpgId); i0=i0[which(!is.na(i))]; i1=i1[which(!is.na(i))]; i2=i[which(!is.na(i))]
#i=which(stat0[i0,colIdPV]<pThres & stat1[i1,colIdPV]<pThres & stat2[i2,colIdPV]<pThres)

## ----------------------------------------------
loadDataFlag=T
loadDataFlag=F

transformFlag=""
transformFlag="_mVal"


cpgId=unique(c(stat0$cpgId,stat1$cpgId,stat2$cpgId))

if (loadDataFlag) {
    dirCom="docs/all/"
    datType="_allGuthSet2"
    load("methAll/data_allGuthSet2.RData")
    meth=meth[rownames(meth)%in%cpgId,]
    if (transformFlag=="_mVal") {
        methInfo=read.table(paste(dirCom,"summaryBeta.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
        k=which(methInfo$cohort==datType)
        meth[which(meth==0)]=methInfo$minBeta*0.5; meth[which(meth==1)]=methInfo$maxBeta+(1-methInfo$maxBeta)*0.5
        cat("Min meth:",min(meth),min(meth,na.rm=T),"\n")
        cat("Max meth:",max(meth),max(meth,na.rm=T),"\n")
        meth=log2(meth/(1-meth))
    }
}
iA2=match(rownames(meth),ann$IlmnID)
j=match(phen$subjectId,phenG$subjectId); j1=which(!is.na(j)); j2=j[j1]

i=match(paste(stat0$cpgId,stat0$gene_genotype),paste(stat1$cpgId,stat1$gene_genotype)); i0=which(!is.na(i)); i1=i[i0]
i=match(paste(stat0$cpgId,stat0$gene_genotype)[i0],paste(stat2$cpgId,stat2$gene_genotype)); i0=i0[which(!is.na(i))]; i1=i1[which(!is.na(i))]; i2=i[which(!is.na(i))]

pThres=0.05
adjPFlag="holm"
varList=c("hlaBallele1v0")
vId=1
colId=paste("coef_",varList[vId],sep="")
colIdPV=paste("holm_",varList[vId],sep="")
compList=c("signifCaCo","signifCoNotCa","signifCaNotCo","signifCoCa1","signifCoCa2","signifCoNotCa1","signifCoNotCa2","signifCa1NotCo","signifCa2NotCo")
compList=c("signifCaCo","signifCoNotCa","signifCaNotCo","signifCoCa1","signifCoCa2","signifCoNotCa1","signifCoNotCa2","signifCa1NotCo","signifCa2NotCo")

for (modelFlag in c("meth~caco+genotype","meth~caco*genotype")) {
    tbl=NULL
    for (compId in 1:length(compList)) {
        switch(compList[compId],
            "signifCaCo"={
                header="Significant in set1+set2 ctrl, set1 case, set2 case"
                i=which(stat0[i0,colIdPV]<pThres & stat1[i1,colIdPV]<pThres & stat2[i2,colIdPV]<pThres)
            },
            "signifCoNotCa"={
                header="Significant in set1+set2 ctrl but not in set1 case or set2 case"
                i=which(stat0[i0,colIdPV]<pThres & stat1[i1,colIdPV]>=pThres & stat2[i2,colIdPV]>=pThres)
            },
            "signifCaNotCo"={
                header="Significant set1 case & set2 case but not in in set1+set2 ctrl"
                i=which(stat0[i0,colIdPV]>=pThres & stat1[i1,colIdPV]<pThres & stat2[i2,colIdPV]<pThres)
            },
            "signifCoCa1"={
                header="Significant in set1+set2 ctrl & set1 case"
                i=which(stat0[i0,colIdPV]<pThres & stat1[i1,colIdPV]<pThres)
            },
            "signifCoCa2"={
                header="Significant in set1+set2 ctrl & set2 case"
                i=which(stat0[i0,colIdPV]<pThres & stat2[i2,colIdPV]<pThres)
            },
            "signifCoNotCa1"={
                header="Significant in set1+set2 ctrl but not in set1 case"
                i=which(stat0[i0,colIdPV]<pThres & stat1[i1,colIdPV]>=pThres)
            },
            "signifCoNotCa2"={
                header="Significant in set1+set2 ctrl but not in set2 case"
                i=which(stat0[i0,colIdPV]<pThres & stat2[i2,colIdPV]>=pThres)
            },
            "signifCa1NotCo"={
                header="Significant set1 case but not in in set1+set2 ctrl"
                i=which(stat0[i0,colIdPV]>=pThres & stat1[i1,colIdPV]<pThres)
            },
            "signifCa2NotCo"={
                header="Significant set2 case but not in in set1+set2 ctrl"
                i=which(stat0[i0,colIdPV]>=pThres & stat2[i2,colIdPV]<pThres)
            }
        )
        cat("\n\n",header,"\n",sep="")
        cat("No. of loci:",length(i),"\n",sep="")
        if (length(i)==0) {cat("No significant CpGs !!!\n"); next}
        tmp=rep(NA,length(i)); tmpC=rep("",length(i))
        switch(modelFlag,
            "meth~caco+genotype"={
                colId=c("caco1","genotype1","genotype2")
            },
            "meth~caco*genotype"={
                colId=c("caco1","genotype1","genotype2","caco1.genotype1","caco1.genotype2")
            }
        )
        out=matrix(nrow=length(i),ncol=3*length(colId))
        colnames(out)=paste(c("coef","pv",adjPFlag),"_",rep(colId,each=3),sep="")
        tbl2=data.frame(comprison=rep(header,length(i)),stat0[i0[i],c("cpgId","gene_genotype")],stringsAsFactors=F)
        iList=i
        for (i in 1:length(iList)) {
            ii=iList[i]
            ii1=which(rownames(meth)==stat0$cpgId[i0][ii])
            ii2=which(rownames(callG)==stat0$gene_genotype[i0][ii])
            dat=cbind(meth=meth[ii1,j1],phen[j1,],genotype=as.factor(callG[ii2,j2]))
            if (all(dat$genotype!=0,na.rm=T)) {cat("No allele0!!!  Exiting ...\n"); break}
            fit=try(lm(as.formula(modelFlag),data=dat))
            if (class(fit)=="try-error") {cat("Cannot fit !!!\n"); next}
            res=summary(fit)$coef
            rownames(res)=sub(":",".",rownames(res))
            res=res[2:nrow(res),]
            out[i,match(paste("coef_",rownames(res),sep=""),colnames(out))]=res[,"Estimate"]
            out[i,match(paste("pv_",rownames(res),sep=""),colnames(out))]=res[,"Pr(>|t|)"]
        }
        tbl=rbind(tbl,cbind(tbl2,out))
    }
    for (k in grep("pv_",names(tbl))) {
        tbl[,sub("pv",adjPFlag,names(tbl)[k])]=p.adjust(tbl[,k],method=adjPFlag)
    }
    for (k in grep("coef_",names(tbl))) tbl[,k]=round(tbl[,k],2)
    for (k in grep(paste("pv_|",adjPFlag,"_",sep=""),names(tbl))) tbl[,k]=signif(tbl[,k],2)
    names(tbl)=gsub("gene_genotype","genotype",gsub("genotype1","hlaBallele1v0",gsub("genotype2","hlaBallele2v0",names(tbl))))
    switch(modelFlag,
        "meth~caco+genotype"={
            tblA=tbl
        },
        "meth~caco*genotype"={
            tblM=tbl
        }
    )
    fName=sub("~","Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T))
    fName=paste("stat_",fName,"_hlaB_allGuthSet2",transformFlag,".txt",sep="")
    write.table(tbl,file=fName, append=F,col.names=T,row.names=F, sep="\t",quote=F)
}

for (modelFlag in c("meth~caco+genotype","meth~caco*genotype")) {
    switch(modelFlag,
    "meth~caco+genotype"={
        tbl=tblA
    },
    "meth~caco*genotype"={
        tbl=tblM
    }
    )
    cat("\n\nModel: ",modelFlag,"\n",sep="")
    cat("-----------------------------\n",sep="")
    for (k in grep(adjPFlag,names(tbl))) {
        i=which(tbl[,k]<pThres)
        cat(names(tbl)[k],": No. of loci significant ",length(i),"\n",sep="")
        if (length(i)==0) next
        #print(tbl[i,])
    }
}
"
Model: meth~caco+genotype
-----------------------------
holm_caco1: No. of loci significant 0
holm_hlaBallele1v0: No. of loci significant 91
holm_hlaBallele2v0: No. of loci significant 0


Model: meth~caco*genotype
-----------------------------
holm_caco1: No. of loci significant 3
holm_hlaBallele1v0: No. of loci significant 63
holm_hlaBallele2v0: No. of loci significant 0
holm_caco1.hlaBallele1v0: No. of loci significant 0
holm_caco1.hlaBallele2v0: No. of loci significant 0
"

modelFlag="meth~caco*genotype"
tbl1=tblM
ii=which(tbl1$holm_caco1<pThres & tbl1$holm_hlaBallele1v0<pThres)
header="Loci with both main effects significant in multiplicative model\nNo locus has both caco and genotype as significant in additive model\nThere are no significant interaction effects"
for (modelFlag in c("meth~caco+genotype","meth~caco*genotype")) {
    switch(modelFlag,
    "meth~caco+genotype"={
        tbl=tblA
    },
    "meth~caco*genotype"={
        tbl=tblM
    }
    )
    tbl=tbl[ii,]
    kk=c()
    for (k in 1:ncol(tbl)) if (any(!is.na(tbl[,k]))) kk=c(kk,k)
    tbl=tbl[,kk]
    fName=paste(sub("~","Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T)),"_",adjPFlag,pThres,sep="")
    fName=paste("stat_",fName,"_hlaB_allGuthSet2",transformFlag,".txt",sep="")
    write.table(header,file=fName,append=F,col.names=F,row.names=F, sep="\t",quote=F)
    write.table(paste("\nStatistics from model: ",modelFlag,sep=""),file=fName,append=T,col.names=F,row.names=F, sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F, sep="\t",quote=F)
}
ii=ii[1]
ii1=which(rownames(meth)==tbl1$cpgId[ii])
ii2=which(rownames(callG)==tbl1$genotype[ii])
grp=paste(phen$caco[j1],callG[ii2,j2])
j=which(callG[ii2,j2]!=2)
grp=grp[j]
grpUniqAll=paste(rep(0:1,each=3)," ",0:2,sep="")
grpUniq=sort(unique(grp))
ttl=paste(paste("hlaBallele",0:2,sep=""),"\n",rep(c("ctrl","case"),each=3),sep="")
ttl=ttl[match(grpUniq,grpUniqAll)]
x=table(grp)
ttl=paste(ttl," (",x,")",sep="")
png(paste("methGenotypeBoxplot_signif_hlaB_allGuthSet2",transformFlag,".png",sep=""))
boxplot(meth[ii1,j1[j]]~grp,names=ttl,main=paste("Set2: ",tbl1$cpgId[ii],", ",tbl1$genotype[ii],sep=""),ylab="M-value",las=0)
abline(h=1,lty="dotted")
dev.off()


####################################################################
####################################################################
