computerFlag="cluster"
computerFlag=""

if (computerFlag=="cluster") {
	setwd("/home/royr/project/JoeWiemels")
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
	setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

##############################################
if (computerFlag=="cluster") {
    dirSmirnov="data/IvanSmirnov/"
    dirHg="tmp/"
} else {
    dirSmirnov="docs/IvanSmirnov/"
    dirHg=""
    dirHg="annotation/"
}

##############################################
## Section 1

## Mouse Dec. 2011 (GRCm38/mm10) Assembly
nProbe=1000
nProbe=-1
datadir=dirSmirnov
tbl=read.table(paste(datadir,"EC_HH_4473_mincov10_CpG_meth.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=nProbe)
colId=3:ncol(tbl)
dge=list(counts=as.matrix(tbl[,colId])/100,ann=cbind(geneId=paste("gene",1:nrow(tbl),sep=""),tbl[,-colId]),samples=data.frame(id=names(tbl)[colId],group=rep(c("corn oil","pcb"),each=3),stringsAsFactors=F))
rownames(dge$counts)=dge$ann$geneId
dge$ann$chr=as.integer(sub("chr","",dge$ann$chr))
rm(tbl)


##############################################
## Map mouse - human

if (F) {
    
## -----------------------------------------
## Mouse to human
tbl=paste("chr",dge$ann$chr,":",dge$ann$pos,"-",dge$ann$pos,sep="")
write.table(tbl,paste("ann_forUCSCLiftover_mm10.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)
datadir=""
annC=read.table(paste(datadir,"ann_fromUCSCliftover_mmu10toHg19.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
annU=read.table(paste(datadir,"ann_fromUCSCliftover_unconverted_mmu10toHg19.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
i2=which(substr(annU[,1],1,nchar("chr"))=="chr")
x2=annU[i2,1]
x1=rep(NA,length(x))
x1[!x%in%x2]=annC[,1]
out=t(sapply(x1,function(x) {
    y=strsplit(x,":")[[1]]
    chr=y[1]
    if (!is.na(chr)) {
        if (chr=="chrX") {chr="chr23"} else if (chr=="chrY") {chr="chr24"}
    }
    chr=as.integer(sub("chr","",chr))
    pos=strsplit(y[2],"-")[[1]]
    c(chr,pos[1],pos[2])
},USE.NAMES=F))
tmp=rep("",nrow(ann))
ann2=matrix("",nrow=nrow(ann),ncol=3)
colnames(ann2)=c("chr","start","end")
ann2[i,]=out
ann2=as.data.frame(ann2,stringsAsFactors=F)
ann2$chr=sub("chr","",ann2$chr)
for (k in which(names(ann2)%in%c("chr","start","end"))) {
    ann2[,k]=as.integer(ann2[,k])
}
ann2=cbind(ann[,c("geneId","geneSym")],ann2)

## -----------------------------------------
## Human to mouse
load("annAll.RData")
tbl=paste(ann$CHR,":",ann$MAPINFO,"-",ann$MAPINFO,sep="")
write.table(tbl,paste("ann_forUCSCLiftover_hg19.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)
ann=dge$ann
datadir=""
annC=read.table(paste(datadir,"ann_bumphunter_fwer0.05_fromUCSCLiftover_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
annU=read.table(paste(datadir,"ann_bumphunter_fwer0.05_fromUCSCLiftover_unconverted_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
i2=which(substr(annU[,1],1,nchar("chr"))=="chr")
x2=annU[i2,1]
x1=rep(NA,length(x))
x1[!x%in%x2]=annC[,1]
out=t(sapply(x1,function(x) {
    y=strsplit(x,":")[[1]]
    chr=y[1]
    if (!is.na(chr)) {
        if (chr=="chrX") {chr="chr23"} else if (chr=="chrY") {chr="chr24"}
    }
    chr=as.integer(sub("chr","",chr))
    pos=strsplit(y[2],"-")[[1]]
    c(chr,pos[1],pos[2])
},USE.NAMES=F))
tmp=rep("",nrow(ann))
ann2=matrix("",nrow=nrow(ann),ncol=3)
colnames(ann2)=c("chr","start","end")
ann2[i,]=out
ann2=as.data.frame(ann2,stringsAsFactors=F)
ann2$chr=sub("chr","",ann2$chr)
for (k in which(names(ann2)%in%c("chr","start","end"))) {
    ann2[,k]=as.integer(ann2[,k])
}
ann2=cbind(ann[,c("geneId","geneSym")],ann2)

}

## -----------------------------------------
## Top genes: Human to mouse

genesetFlag="pcb"
genesetFlag="pcb_logged_PCB_138_SRS"
genesetFlag="pcb_logged_PCB_105_SRS"
genesetFlag="pcb_logged_PCB_aroclor1260"
genesetFlag="candGene"
genesetFlag="bumphunter"


#for (genesetFlag in c("pcb_logged_PCB_138_SRS","pcb_logged_PCB_105_SRS")) {

datadir=dirHg

if (genesetFlag=="bumphunter") {
    pvName="fwer"
    pThres=0.05
    pThres=0.1
    
    pvName="pv"
    pThres=0.05
    pThres=0.2
    
    fNameC=paste("_",genesetFlag,"_",pvName,pThres,sep="")
    
    if (F) {
        ann=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_forUCSCLiftover_hg19.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
        ann=as.data.frame(t(sapply(ann[,1],function(x) {
            y=strsplit(x,":")[[1]]
            chr=y[1]
            if (!is.na(chr)) {
                if (chr=="chrX") {chr="chr23"} else if (chr=="chrY") {chr="chr24"}
            }
            chr=as.integer(sub("chr","",chr))
            pos=strsplit(y[2],"-")[[1]]
            c(chr,pos[1],pos[2])
        },USE.NAMES=F)),stringsAsFactors=F)
        names(ann)=paste(c("chr","start","end"),"_hg",sep="")
        for (k in which(names(ann2)%in%paste(c("chr","start","end"),"_hg",sep=""))) {
            ann[,k]=as.integer(ann[,k])
        }
    }
    ann=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_forUCSCLiftover_hg19.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    x=ann[,1]
    ann=read.table(paste(datadir,"stat_all_",genesetFlag,"_",pvName,pThres,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(ann)[match(c("p.value"),names(ann))]=c("pvalue")
    annC=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_fromUCSCLiftover_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    annU=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_fromUCSCLiftover_unconverted_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
} else if (substr(genesetFlag,1,3)=="pcb") {
    pvName="qv"
    pThres=0.05
    
    pvName="pv"
    pThres=0.001
    pThres=0.01
    
    fNameC=paste("_",genesetFlag,"_",pvName,pThres,sep="")
    
    load("annAll.RData")
    ann1=ann
    if (F) {
        ann=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_forUCSCLiftover_hg19.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
        ann=as.data.frame(t(sapply(ann[,1],function(x) {
            y=strsplit(x,":")[[1]]
            chr=y[1]
            if (!is.na(chr)) {
                if (chr=="chrX") {chr="chr23"} else if (chr=="chrY") {chr="chr24"}
            }
            chr=as.integer(sub("chr","",chr))
            pos=strsplit(y[2],"-")[[1]]
            c(chr,pos[1],pos[2])
        },USE.NAMES=F)),stringsAsFactors=F)
        names(ann)=paste(c("chr","start","end"),"_hg",sep="")
        for (k in which(names(ann2)%in%paste(c("chr","start","end"),"_hg",sep=""))) {
            ann[,k]=as.integer(ann[,k])
        }
    }
    ann=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_forUCSCLiftover_hg19.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    x=ann[,1]
    ann=read.table(paste(datadir,"stat_",genesetFlag,"_",pvName,pThres,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(ann)[match(c("coef","pv","qv"),names(ann))]=c("value","pvalue","qvalue")
    i=match(ann$cpgId,ann1$IlmnID)
    ann$chr=ann1$CHR[i]
    ann$start=ann$end=ann1$MAPINFO[i]
    rm(ann1)
    ann=ann[,c("comparison","cpgId","chr","start","end","value","pvalue","qvalue")]
    annC=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_fromUCSCLiftover_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    annU=read.table(paste(datadir,"ann_",genesetFlag,"_",pvName,pThres,"_fromUCSCLiftover_unconverted_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
} else {
    
    fNameC=paste("_candGene",sep="")
    
    ann=read.table(paste(datadir,"ann",fNameC,"_forUCSCLiftover_hg19.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    x=ann[,1]
    #ann=read.table(paste(datadir,"stat",fNameC,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    tmp=rep(NA,length(x))
    out=t(sapply(x,function(x) {
        y=strsplit(sub("chr","",x),":")[[1]]
        chr=y[1]
        pos=strsplit(y[2],"-")[[1]]
        print(chr)
        print(pos)
        as.integer(c(chr,pos))
    },USE.NAMES=F))
    ann=data.frame(comparison=c("H19","IGF2"),chr=out[,1],start=out[,2],end=out[,3],value=tmp,pvalue=tmp,fwer=tmp,stringsAsFactors=F)
    annC=read.table(paste(datadir,"ann",fNameC,"_fromUCSCLiftover_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
    annU=read.table(paste(datadir,"ann",fNameC,"_fromUCSCLiftover_unconverted_hg19toMmu10.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
}

i=1:nrow(ann)
x0=as.character(ann$chr)
x0[which(x0=="22")]="X"
x0[which(x0=="23")]="Y"
x0=paste("chr",x0,":",ann$start,"-",ann$end,sep="")
i2=which(substr(annU[,1],1,nchar("chr"))=="chr")
x2=annU[i2,1]
x1=rep(NA,length(x))
x1[!x%in%x2]=annC[,1]

if (F) {
    out=t(sapply(x1,function(x) {
        y=strsplit(x,":")[[1]]
        chr=y[1]
        if (!is.na(chr)) {
            if (chr=="chrX") {chr="chr23"} else if (chr=="chrY") {chr="chr24"}
        }
        chr=as.integer(sub("chr","",chr))
        pos=strsplit(y[2],"-")[[1]]
        c(chr,pos[1],pos[2])
    },USE.NAMES=F))
}
out=t(sapply(x1[match(x0,x)],function(x) {
    if (!is.na(x)) {
        y=strsplit(x,":")[[1]]
        chr=y[1]
        if (chr=="chrX") {chr="chr23"} else if (chr=="chrY") {chr="chr24"}
        chr=as.integer(sub("chr","",chr))
        pos=strsplit(y[2],"-")[[1]]
    } else {
        chr=NA
        pos=c(NA,NA)
    }
    c(chr,pos[1],pos[2])
},USE.NAMES=F))

tmp=rep("",nrow(ann))
ann2=matrix("",nrow=nrow(ann),ncol=3)
colnames(ann2)=paste(c("chr","start","end"),"_mmu",sep="")
ann2[i,]=out
ann2=as.data.frame(ann2,stringsAsFactors=F)
for (k in which(names(ann2)%in%paste(c("chr","start","end"),"_mmu",sep=""))) {
    ann2[,k]=as.integer(ann2[,k])
}
i=which(!is.na(ann2$chr_mmu))
ann=ann[i,]
ann2=ann2[i,]
ann2=cbind(ann,ann2)
#ann2=ann2[1,]
#ann2=ann2[2,]
tbl=NULL
thres=10^6
thres=1000
thres=100
thres=10000
cat("No. of human CpGs: ",nrow(ann2),"\n",sep="")
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
for (i1 in 1:nrow(ann2)) {
    #for (i1 in 1:10) {
    if (i1%%100==0) print(i1)
    i=which(dge$ann$chr==ann2$chr_mmu[i1])
    if (length(i)!=0) {
        if (F) {
            i2=i[which(dge$ann$pos[i]>=ann2$start_mmu[i1] & dge$ann$pos[i]<=ann2$end_mmu[i1])]
            if (length(i2)==0) {
                i2=dist2gene=c()
                ii=max(i[dge$ann$pos[i]<ann2$start_mmu[i1]])
                if (!is.infinite(ii)) {
                    i2=c(i2,ii)
                    dist2gene=c(dist2gene,dge$ann$pos[ii]-ann2$start_mmu[i1])
                }
                ii=min(i[dge$ann$pos[i]>ann2$pos[i1]])
                if (!is.infinite(ii)) {
                    i2=c(i2,ii)
                    dist2gene=c(dist2gene,dge$ann$pos[ii]-ann2$end_mmu[i1])
                }
            } else {
                dist2gene=rep(0,length(i2))
            }
        }
        if (F) {
            i2=dist2gene=c()
            ii=i[dge$ann$pos[i]<ann2$start_mmu[i1]]
            ii=ii[order(dge$ann$pos[ii]-ann2$start_mmu[i1],decreasing=T)]
            ii=ii[!is.infinite(ii)]
            ii=ii[1:min(c(10,length(ii)))]
            if (length(ii)!=0) {
                i2=c(i2,ii)
                dist2gene=c(dist2gene,dge$ann$pos[ii]-ann2$start_mmu[i1])
            }
            ii=i[dge$ann$pos[i]>ann2$end_mmu[i1]]
            ii=ii[order(dge$ann$pos[ii]-ann2$end_mmu[i1],decreasing=T)]
            ii=ii[!is.infinite(ii)]
            ii=ii[1:min(c(10,length(ii)))]
            if (length(ii)!=0) {
                i2=c(i2,ii)
                dist2gene=c(dist2gene,dge$ann$pos[ii]-ann2$end_mmu[i1])
            }
        }
        if (T) {
            i2=dist2gene=c()
            ii=i[which(dge$ann$pos[i]>=(ann2$start_mmu[i1]) & dge$ann$pos[i]<=(ann2$end_mmu[i1]))]
            if (length(ii)!=0) {
                i2=c(i2,ii)
                dist2gene=c(dist2gene,rep(0,length(ii)))
            }
            ii=i[which(dge$ann$pos[i]>=(ann2$start_mmu[i1]-thres) & dge$ann$pos[i]<ann2$start_mmu[i1])]
            if (length(ii)!=0) {
                i2=c(i2,ii)
                dist2gene=c(dist2gene,dge$ann$pos[ii]-ann2$start_mmu[i1])
            }
            ii=i[which(dge$ann$pos[i]>ann2$end_mmu[i1] & dge$ann$pos[i]<=(ann2$end_mmu[i1]+thres))]
            if (length(ii)!=0) {
                i2=c(i2,ii)
                dist2gene=c(dist2gene,dge$ann$pos[ii]-ann2$end_mmu[i1])
            }
        }
        #print(length(i2))
        if (length(i2)!=0) {
            if (substr(genesetFlag,1,3)=="pcb") {
                tbl2=cbind(comparison=rep(ann2$comparison[i1],length(i2)),cpgId=ann2$cpgId[i1],matrix(rep(unlist(ann2[i1,!names(ann2)%in%c("comparison","cpgId")]),length(i2)),nrow=length(i2),byrow=T),dge$ann[i2,],dist2gene)
            } else {
                tbl2=cbind(comparison=rep(ann2$comparison[i1],length(i2)),matrix(rep(unlist(ann2[i1,!names(ann2)%in%c("comparison","cpgId")]),length(i2)),nrow=length(i2),byrow=T),dge$ann[i2,],dist2gene)
            }
            tbl=rbind(tbl,tbl2)
        }
    }
}
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))
#names(tbl)=c(paste(c("chr","start","end"),"_hg",sep=""),paste(c("geneId","chr","pos"),"_mmu",sep=""),"dist")
names(tbl)=c(paste(names(ann),"_hg",sep=""),paste(c("chr","start","end"),"_mmu_mapped",sep=""),paste(c("geneId","chr","pos"),"_mmu",sep=""),"dist2mapped_mmu")
rownames(tbl)=NULL
candGene=tbl

if (F) {
    grp=paste(tbl$chr_hg,tbl$start_hg,tbl$end_hg)
    grpUniq=unique(grp)
    for (gId in 1:length(grpUniq)) {
        i=which(grp==grpUniq[gId])
        #print(summary(tbl$dist2mapped_mmu[i]))
        print(length(i))
    }
}
rm(tbl,tbl2)
save(candGene,file=paste("candGene",fNameC,".RData",sep=""))
#}

## -----------------------------------------


## -----------------------------------------
## Section 2
##
## Voom
## Run section 1 first

library(limma)
library(edgeR)

fName3="_mmu"

pThres=0.05

pThres=0.1

dgeF=dge
cpm=cpm(dge$counts,log=F)
lcpm=cpm(dge$counts,log=TRUE)
maxCntVec=apply(cpm,1,max,na.rm=T)
maxCntVec2=apply(dge$counts,1,max,na.rm=T)

minCnt=0
modelThis=paste("~group",sep="")
modelThis=as.formula(modelThis)
design=model.matrix(modelThis,data=dgeF$samples)
if (F) {
    png("tmp.png")
    plot(fit)
    dev.off()
}
minCntList=c(0,1,2)
minCntList=c("_allCpG","_perc67","_perc80")
for (minCnt in minCntList) {
    #i=which(maxCntVec>=minCnt)
    if (substr(minCnt,1,5)=="_perc") {
        cutoff=quantile(maxCntVec2,probs=as.numeric(sub("_perc","",minCnt))/100,na.rm=T)
        i=which(maxCntVec2>=cutoff)
    } else if (minCnt=="_allCpG") {
        i=1:length(maxCntVec2)
    } else {
        i=which(maxCntVec>=minCnt)
    }
    dat <- voom(dgeF$counts[i,],design,save.plot=F)
    fit <- lmFit(dat,design)
    rm(dat)
    fit <- eBayes(fit)
    colId=2
    top=data.frame(geneId=rownames(fit$coef),log2FC=fit$coef[,colId],pvalue=fit$p.value[,colId],stringsAsFactors=F)
    write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

    fName=paste("_mmu_minCnt",minCnt,sep="")
    header=paste("Mouse:\nCpGs with at least one sample having counts per million (CPM) >= ",minCnt,sep="")
    fName=paste("_mmu",minCnt,sep="")
    if (minCnt=="_allCpG") {
        header="Mouse: All CpGs"
    } else {
        header=paste("Mouse: N = ",length(i),"\nCpGs with at least one sample having\nbeta-value >= ",cutoff," (",sub("_perc","",minCnt)," percentile)",sep="")
    }
    
    png(paste("hist",fName,".png",sep=""))
    hist(top$pvalue,main=header,xlab="P-value")
    dev.off()
    
    switch(as.character(minCnt),
        "0"={
            top0=top
        },
        "1"={
            top1=top
        },
        "2"={
            top2=top
        },
        "_allCpG"={
            top0=top
        },
        "_perc67"={
            top67=top
        },
        "_perc80"={
            top80=top
        }
    )
}

## ----------------------------------
## Candidate genes

if (genesetFlag%in%c("candGene") | (genesetFlag=="bumphunter" & fNameC=="_bumphunter_fwer0.1")) {
    top=top0
    #top=top1
    top=top67
    
    load(file=paste("annotation/candGene",fNameC,".RData",sep=""))
    i=match(candGene$geneId_mmu,top$geneId); i1=which(!is.na(i)); i2=i[i1]
    tbl=top[i2,c("log2FC","pvalue")]
    names(tbl)=paste(names(tbl),"_mmu",sep="")
    tbl=cbind(candGene[i1,],tbl)
    for (k in c("log2FC_mmu","pvalue_mmu")) tbl[,k]=round(tbl[,k],4)
    rownames(tbl)=NULL

    tbl$bonf_mmu=p.adjust(tbl$pvalue_mmu,method="bonferroni")
    tbl$holm_mmu=p.adjust(tbl$pvalue_mmu,method="holm")
    grp=tbl$comparison_hg
    grpUniq=unique(grp)
    for (gId in 1:length(grpUniq)) {
        i=which(grp==grpUniq[gId])
        tbl$bonf_mmu[i]=p.adjust(tbl$pvalue_mmu[i],method="bonferroni")
        tbl$holm_mmu[i]=p.adjust(tbl$pvalue_mmu[i],method="holm")
    }

    if (genesetFlag=="bumphunter") {
        colIdPV="pvalue_hg"; pThres=c(0.1,1); lfc=c(0,0); dist2Reg=10000
    } else {
        colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(0,0); dist2Reg=0; compFlag="H19"
        colIdPV="pvalue_hg"; pThres=c(0.05,1); lfc=c(0,0); dist2Reg=0; compFlag="H19"
    }
    i=abs(tbl$dist2mapped_mmu)<=dist2Reg
    #tbl[i,]
    i=tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=dist2Reg
    #tbl=tbl[i,]

    names(tbl)[match("comparison_hg",names(tbl))]="geneSym"
    tbl=tbl[,c("geneSym","chr_hg","start_hg","end_hg","chr_mmu_mapped","start_mmu_mapped","end_mmu_mapped","geneId_mmu","chr_mmu","pos_mmu","dist2mapped_mmu","log2FC_mmu","pvalue_mmu","bonf_mmu","holm_mmu")]
    fName=paste(fNameC,"_pv_mmu",pThres[2],"_dist2reg",dist2Reg,sep="")
    write.table(tbl[i,],paste("stat",fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

## ----------------------------------
## Quadrant plot of statistics of hg vs. mmu

top=top0
#top=top1
top=top67

if (F) {
    pvName="fwer"
    pThres=0.05
    pThres=0.1
    pvName="pv"
    pThres=0.05
    pThres=0.2
}

pvInfo=data.frame(colIdPV=c("pvalue_hg","fwer_hg","qvalue_hg"),pvName=c("pv","fwer","qv"),stringsAsFactors=F)

datadir=""
datadir="annotation/"
load(file=paste(datadir,"candGene",fNameC,".RData",sep=""))
for (k in c("comparison_hg")) candGene[,k]=as.character(candGene[,k])

i=match(candGene$geneId_mmu,top$geneId); i1=which(!is.na(i)); i2=i[i1]
tbl=top[i2,c("log2FC","pvalue")]
names(tbl)=paste(names(tbl),"_mmu",sep="")
tbl=cbind(candGene[i1,],tbl)
for (k in c("log2FC_mmu","pvalue_mmu")) tbl[,k]=round(tbl[,k],4)
rownames(tbl)=NULL
if (F) {
    pThres=c(0.1,0.1)
    pThres=c(0.05,0.1)
    #tbl[tbl$pvalue_mmu<pThres,]
    i=tbl$pvalue_mmu<pThres
    i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms" & tbl$fwer<pThres[1] & tbl$pvalue_mmu<pThres[2]
    #tbl=tbl[i,]
    #write.table(tbl[i,],paste("stat_candGene",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}

colIdPV="pvalue_hg"; pThres=c(0.05,0.1)
colIdPV="pvalue_hg"; pThres=c(0.2,0.2)
colIdPV="fwer_hg"; pThres=c(0.05,0.05)
colIdPV="pvalue_hg"; pThres=c(0.05,0.05)
colIdPV="pvalue_hg"; pThres=c(0.05,0.1)
colIdPV="pvalue_hg"; pThres=c(0.1,0.1)
colIdPV="pvalue_hg"; pThres=c(0.2,0.2); dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_105_SRS, 1000maxGap, 1000perms"
colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=.1; dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms"
colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(.1,.1); dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_105_SRS, 1000maxGap, 1000perms"
colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(.1,.1); dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_105_SRS, 1000maxGap, 1000perms"
colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(.1,.1); dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms"

colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(0,0); dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms"
colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(0,0); dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_105_SRS, 1000maxGap, 1000perms"
colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(0,0); dist2Reg=10; compFlag="bumphunter: allGuthSet1Set2, logged_PCB_138_SRS, 1000maxGap, 1000perms"

if (substr(genesetFlag,1,3)=="pcb") {
    varList=c("logged_PCB_105_SRS","logged_PCB_138_SRS","logged_PCB_aroclor1260")
}
varList=unique(as.character(tbl$comparison_hg))
if (genesetFlag=="bumphunter") {
    varList=varList[grep(" allGuthSet1Set2,",varList)]
    varList=varList[grep(", 1000maxGap",varList)]
}
for (varId in varList) {
    if (all(tbl$comparison_hg!=varId)) next
    if (substr(genesetFlag,1,3)=="pcb") {
        xlab="Locus-wise, allGuthSet1Set2: Effect size"
    } else {
        xlab="Bumphunter, allGuthSet1Set2: Effect size"
    }
    colIdPV="pvalue_hg"; pThres=c(0.1,0.1); lfc=c(0,0); dist2Reg=10; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(.1,0); dist2Reg=10; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(.1,0); dist2Reg=5; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(.1,.1); dist2Reg=5; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(.1,.1); dist2Reg=100; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(0,0); dist2Reg=5; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(0,0); dist2Reg=10; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(0,0); dist2Reg=100; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(0,0); dist2Reg=1000; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(0,0); dist2Reg=10; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.05,0.05); lfc=c(0,0); dist2Reg=1000; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="pvalue_hg"; pThres=c(0.1,0.1); lfc=c(0,0); dist2Reg=1000; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
    colIdPV="qvalue_hg"; pThres=c(0.05,0.05); lfc=c(0,0); dist2Reg=10000; compFlag=varId
    colIdPV="pvalue_hg"; pThres=c(0.001,0.001); lfc=c(0,0); dist2Reg=1000; compFlag=varId
    colIdPV="pvalue_hg"; pThres=c(0.001,0.001); lfc=c(0,0); dist2Reg=100; compFlag=varId
    colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(0,0); dist2Reg=max(tbl$dist2mapped_mmu); compFlag=varId

    if (substr(genesetFlag,1,3)=="pcb") {
        colIdPV="pvalue_hg"; pThres=c(0.01,0.01); lfc=c(0,0); dist2Reg=1000; compFlag=varId
    } else {
        colIdPV="pvalue_hg"; pThres=c(0.2,0.2); lfc=c(0,0); dist2Reg=1000; compFlag=varId
    }

    pThresInit=pThres
    
    sort(unique(tbl$comparison_hg))
    if (F) {
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_138_SRS, 1000maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=10
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_105_SRS, 1000maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=10
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=5
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 400maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=100
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 100maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=1000
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=10
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2]
        i=tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=10
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_105_SRS, 1000maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=10
        i=tbl$comparison_hg=="bumphunter: allGuthSet1Set2, logged_PCB_aroclor1260, 1000maxGap, 1000perms" & tbl[,colIdPV]<pThres[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$dist2mapped_mmu)<=10
    }
    i=tbl$comparison_hg==compFlag & tbl[,colIdPV]<pThres[1] & abs(tbl$value_hg)>=lfc[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$log2FC_mmu)>=lfc[2] & abs(tbl$dist2mapped_mmu)<=dist2Reg
    if (any(i)) {
        i=which(i)
        i=i[!duplicated(paste(tbl$chr_mmu[i],tbl$start_mmu[i],tbl$end_mmu[i]))]
        x=table(sign(tbl$value_hg[i]),sign(tbl$log2FC_mmu[i]))
        x
        res=try(mcnemar.test(x=x,correct=F))
        pv=ifelse(class(res)=="try-error",NA,res$p.value)
        cat("\n\n",varId,", No. of mmu CpGs: ",length(i),", prop in same dirn with hg: ",round(mean(sign(tbl$value_hg[i])==sign(tbl$log2FC_mmu[i])),2),", PV: ",signif(pv,2),"\n",sep="")
        cor(tbl$value_hg[i],tbl$log2FC_mmu[i],use="complete.obs",method="pearson")
        fName=paste("_",gsub(": |, ","_",compFlag),"_",pvInfo$pvName[match(colIdPV,pvInfo$colIdPV)],"_hg",pThres[1],"_pv_mmu",pThres[2],"_dist2reg",dist2Reg,sep="")
        write.table(tbl[i,],paste("stat_candGene",fName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

        header=paste(compFlag,", ",pvInfo$pvName[match(colIdPV,pvInfo$colIdPV)],"_hg<",pThres[1],"\npv_mmu<",pThres[2],", abs(dist2mapped_mmu)<=",dist2Reg,sep="")
        png(paste("effectSizePlot_hgVmmu",fName,".png",sep=""))
        plot(tbl$value_hg[i],tbl$log2FC_mmu[i],main=paste(header,"\nNo. of mouse CpGs: ",length(i),", proportion in same direction with human: ",round(mean(sign(tbl$value_hg[i])==sign(tbl$log2FC_mmu[i])),2),sep=""),xlab=xlab,ylab="Mouse: Effect size",cex.main=.8)
        abline(h=0,lty="dotted")
        abline(v=0,lty="dotted")
        dev.off()
    }
    
    if (colIdPV=="pvalue_hg" & all(lfc==0)) {
        for (dist2Reg in c(10,100,1000,10000)) {
            pThres=pThresInit
            fName=paste("_",gsub(": |, ","_",compFlag),"_",pvInfo$pvName[match(colIdPV,pvInfo$colIdPV)],"_hg",pThres[1],"_dist2reg",dist2Reg,sep="")
            #if (colIdPV=="pvalue_hg" & all(lfc==0) & dist2Reg==1000) {
            cat("\n\nNo. of mmu CpGs mapped:\n")
            if (substr(genesetFlag,1,3)=="pcb") {
                xVec=c(.001,.00075,0.0005,0.00025,0.0001)
                xVec=c(0.01,0.0075,0.0050,0.0025,0.0001)
            } else {
                xVec=c(.05,0.1,0.15,0.2,0.25)
            }
            yVec=rep(NA,length(xVec))
            ttl=rep(NA,length(xVec))
            for (vId in 1:length(xVec)) {
                if (substr(genesetFlag,1,3)=="pcb") {
                    pThres[2]=xVec[vId]
                    header=paste(compFlag,"\npvalue_hg<",pThres[1],", abs(dist2mapped_mmu)<=",dist2Reg,sep="")
                    xlab="Mouse: P-value cutoff"
                } else {
                    pThres=rep(xVec[vId],2)
                    header=paste(compFlag,"\nabs(dist2mapped_mmu)<=",dist2Reg,sep="")
                    xlab="P-value cutoff"
                }
                #colIdPV="pvalue_hg"; lfc=c(0,0); dist2Reg=1000; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
                i=tbl$comparison_hg==compFlag & tbl[,colIdPV]<pThres[1] & abs(tbl$value_hg)>=lfc[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$log2FC_mmu)>=lfc[2] & abs(tbl$dist2mapped_mmu)<=dist2Reg
                i=which(i)
                i=i[!duplicated(paste(tbl$chr_mmu[i],tbl$start_mmu[i],tbl$end_mmu[i]))]
                x=table(sign(tbl$value_hg[i]),sign(tbl$log2FC_mmu[i]))
                x
                cat("PV ",pThres[2],": ",length(i),"\n",sep="")
                yVec[vId]=mean(sign(tbl$value_hg[i])==sign(tbl$log2FC_mmu[i]))
                ttl[vId]=paste(xVec[vId],"\n(N=",length(i),")",sep="")
            }
            if (any(!is.na(yVec))) {
                png(paste("effectSizeVpv_propSameDir_hgVmmu",fName,".png",sep=""))
                plot(xVec,yVec,main=paste(header,sep=""),xlab=xlab,ylab="Proportion of mouse CpGs in same direction as in human",cex.main=.8,xaxt="n")
                axis(side=1,at=xVec,labels=ttl)
                dev.off()
            }
        }
    }
    
    if (colIdPV=="pvalue_hg" & all(lfc==0)) {
        pThres=pThresInit
        fName=paste("_",gsub(": |, ","_",compFlag),"_",pvInfo$pvName[match(colIdPV,pvInfo$colIdPV)],"_hg",pThres[1],"_pv_mmu",pThres[2],sep="")
        cat("\n\nNo. of mmu CpGs mapped:\n")
        xVec=c(1,10,100,1000,10000)
        
        yVec=rep(NA,length(xVec))
        ttl=rep(NA,length(xVec))
        for (vId in 1:length(xVec)) {
            dist2Reg=xVec[vId]
            if (substr(genesetFlag,1,3)=="pcb") {
                pThres=c(0.01,0.01)
            } else {
                pThres=c(0.05,0.05)
            }
            header=paste(compFlag,"\n",pvInfo$pvName[match(colIdPV,pvInfo$colIdPV)],"_hg<",pThres[1],", pv_mmu<",pThres[2],sep="")
            xlab="Mouse: dist2mappedRegion"
            #colIdPV="pvalue_hg"; lfc=c(0,0); dist2Reg=1000; compFlag=paste("bumphunter: allGuthSet1Set2, ",varId,", 1000maxGap, 1000perms",sep="")
            i=tbl$comparison_hg==compFlag & tbl[,colIdPV]<pThres[1] & abs(tbl$value_hg)>=lfc[1] & tbl$pvalue_mmu<pThres[2] & abs(tbl$log2FC_mmu)>=lfc[2] & abs(tbl$dist2mapped_mmu)<=dist2Reg
            i=which(i)
            i=i[!duplicated(paste(tbl$chr_mmu[i],tbl$start_mmu[i],tbl$end_mmu[i]))]
            x=table(sign(tbl$value_hg[i]),sign(tbl$log2FC_mmu[i]))
            x
            cat("dist2reg ",pThres[2],": ",length(i),"\n",sep="")
            yVec[vId]=mean(sign(tbl$value_hg[i])==sign(tbl$log2FC_mmu[i]))
            ttl[vId]=paste(xVec[vId],"\n(N=",length(i),")",sep="")
        }
        if (any(!is.na(yVec))) {
            png(paste("effectSizeVdist2reg_propSameDir_hgVmmu",fName,".png",sep=""))
            plot(log10(xVec),yVec,main=paste(header,sep=""),xlab=xlab,ylab="Proportion of mouse CpGs in same direction as in human",cex.main=.8,xaxt="n")
            axis(side=1,at=log10(xVec),labels=ttl)
            dev.off()
        }
    }
    
}

####################################################################################
####################################################################################

if (F) {
genesetFlag="pcb_logged_PCB_aroclor1260"
genesetFlag="pcb_logged_PCB_138_SRS"
genesetFlag="pcb_logged_PCB_105_SRS"

pvName="pv"
pThres=0.001
pThres=0.01

fNameC=paste("_",genesetFlag,"_",pvName,pThres,sep="")
}
