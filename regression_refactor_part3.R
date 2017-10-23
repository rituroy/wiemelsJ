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

## ----------------------------------------------


#############################################################
## take a look at significant loci enrichment
#############################################################

datFlag="meth"
datFlag=""

ann$region <- sapply(as.character(ann[, "UCSC_RefGene_Group"]), function(z){
    if(z=="") {
        return("Intergenic")
    } else {
        zz <- strsplit(z,";")[[1]]
        if(all(zz%in%c("TSS1500","TSS200","5'UTR","1stExon"))) {
            return("Promoter")
        } else if(all(zz%in%c("Body","3'UTR"))) {
            return("Body")
        } else {
            return("Uncertain") #uncertain if the locus mapped to different regions in different transcripts
        }
    }
})
ann$relationToIsland <- ann[,"Relation_to_UCSC_CpG_Island"]; ann$relationToIsland <- gsub("S_","",gsub("N_","",ann$relationToIsland)); ann$relationToIsland[ann$relationToIsland==""] <- "notCGI"

## -------------------
if (datFlag=="meth") {
    load("data_tmp_allGuthSet1Set2_ctrl.RData")
}

colIdPV="qv"; pThres=0.05

stat2=stat_pcb170
names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
if (datFlag=="meth") {
    i=match(stat2$cpgId,rownames(meth)); i1=which(!is.na(i)); i2=i[i1]
    stat2=stat2[i1,]
    meth=meth[i2,]
}
i=match(stat2$cpgId,ann$IlmnID); i2=which(ann$keep[i]==1)
stat2=stat2[i2,]
if (datFlag=="meth") {
    meth=meth[i2,]
}
iA2=match(stat2$cpgId,ann$IlmnID)
i2=1:nrow(stat2)

id1=as.integer(stat2$cpgId%in%cpgId & stat2[,colIdPV]<pThres)

if (F) {
    id1=as.character(id1)
    id1[which(id1=="0")]="Not in our gene list"
    id1[which(id1=="1")]="In our gene list"
}
id1[is.na(ann$keep[iA2][i2]) | !ann$keep[iA2][i2]]=NA

## -------------------
id2=ann$region[iA2][i2]
id2[ann$region[iA2][i2]=="Uncertain"]=NA

if (F) {
    i=which(ann$Infinium_Design_Type[iA2]=="I")
    fisher.test(table(id1[i], id2[i]=="Promoter"))
    i=which(ann$Infinium_Design_Type[iA2]=="II")
    fisher.test(table(id1[i], id2[i]=="Promoter"))
    
    i=which(ann$Infinium_Design_Type[iA2]=="I")
    fisher.test(table(id1[i], id2[i]))
    i=which(ann$Infinium_Design_Type[iA2]=="II")
    fisher.test(table(id1[i], id2[i]))
}

for (designFlag in c("I","II")) {
    cat("\n=========== Infinium_Design_Type ",designFlag," ====================\n",sep="")
    ii=which(ann$Infinium_Design_Type[iA2][i2]==designFlag)
    x=table(id1[ii],id2[ii])
    rownames(x)=c("Not in our gene list","In our gene list")
    print(x)
    cat("\n")
    for (k1 in 1:(ncol(x)-1)) {
        for (k2 in (1+k1):ncol(x)) {
            i=ii[which(id2[ii]%in%c(colnames(x)[c(k1,k2)]))]
            pv=fisher.test(id1[i],id2[i])$p.value
            suf=""
            if (pv<0.05) {
                suf=" ****"
            } else if (pv<0.1) {
                suf=" **"
            }
            cat(colnames(x)[k1]," vs. ",colnames(x)[k2],": pv ",signif(pv,4),suf,"\n",sep="")
        }
    }
}

## -------------------
id2=ann$relationToIsland[iA2][i2]

for (designFlag in c("I","II")) {
    cat("\n=========== Infinium_Design_Type ",designFlag," ====================\n",sep="")
    ii=which(ann$Infinium_Design_Type[iA2][i2]==designFlag)
    x=table(id1[ii],id2[ii])
    rownames(x)=c("Not in our gene list","In our gene list")
    print(x)
    cat("\n")
    for (k1 in 1:(ncol(x)-1)) {
        for (k2 in (1+k1):ncol(x)) {
            i=ii[which(id2[ii]%in%c(colnames(x)[c(k1,k2)]))]
            suf=""
            x2=table(id1[i],id2[i])
            res=try(fisher.test(x2))
            if (class(res)=="try-error") {
                pv=NA
            } else {
                pv=res$p.value
                if (pv<0.05) {
                    suf=" ****"
                } else if (pv<0.1) {
                    suf=" **"
                }
            }
            cat(colnames(x)[k1]," vs. ",colnames(x)[k2],": pv ",signif(pv,4),suf,"\n",sep="")
            
        }
    }
}

#############################################################
## Is there a global demethylation associated with caco-pobw interaction
#############################################################

dirn=sign(stat2$coef[i2])
dirn[dirn==0]=NA
dirn[is.na(ann$keep[iA2][i2]) | !ann$keep[iA2][i2]]=NA

#Is there a global demethylation associated with caco-pobw interaction?
for (designFlag in c("I","II")) {
    cat("\n=========== Infinium_Design_Type ",designFlag," ====================\n",sep="")
    ii=which(ann$Infinium_Design_Type[iA2][i2]==designFlag)
    
    ## global: binomial test against 0.5
    x=table(dirn[ii])
    #print(x)
    res=binom.test(x=sum(dirn[ii]==-1,na.rm=TRUE), n=sum(!is.na(dirn[ii])), p=0.5, alternative="greater")
    #print(res)
    pv=res$p.value
    suf=""
    if (pv<0.05) {
        suf=" ****"
    } else if (pv<0.1) {
        suf=" **"
    }
    cat("Number total: ",length(ii),"\n",sep="")
    cat("Proportion down globally: ",round(res$estimate,2),"\npv (vs. 0.5) ",signif(pv,4),suf,"\n",sep="")
}

#Is there greater demethylation for significant caco-pobw interaction loci than globally?
for (designFlag in c("I","II")) {
    cat("\n=========== Infinium_Design_Type ",designFlag," ====================\n",sep="")
    ii=which(ann$Infinium_Design_Type[iA2][i2]==designFlag)
    i=ii[which(id1[ii]==1)]
    i_2=ii[which(id1[ii]==0)]
    
    ## among P<0.05
    x=table(dirn[ii])
    #print(x)
    res=prop.test(c(sum(dirn[i]==-1,na.rm=TRUE), sum(dirn[ii]==-1,na.rm=TRUE)),c(sum(!is.na(dirn[i])), sum(!is.na(dirn[ii]),na.rm=TRUE)), alternative="greater")
    #print(res)
    pv=res$p.value
    suf=""
    if (pv<0.05) {
        suf=" ****"
    } else if (pv<0.1) {
        suf=" **"
    }
    cat("Number total = ",length(ii),", number signif = ",length(i),"\n",sep="")
    cat("Proportion down: For signif loci (",round(res$estimate[1],2),") vs. globally (",round(res$estimate[2],2),"):\npv ",signif(pv,4),suf,"\n",sep="")
    
    ## among P>=0.05
    res=prop.test(c(sum(dirn[i]==-1,na.rm=TRUE), sum(dirn[i_2]==-1,na.rm=TRUE)),c(sum(!is.na(dirn[i])), sum(!is.na(dirn[i_2]),na.rm=TRUE)), alternative="greater")
    #print(res)
}




## -----------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------
