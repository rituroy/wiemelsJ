## ---------------------------------

computerFlag=""
computerFlag="cluster"

## ---------------------------------

nProbe=101
nProbe=-1

## ---------------------------------
datType="_aml"
datType="_allGuthSet2"
datType="_allGuthSet1"
datType="_leuk"
datType="_allGuthSet1Set2"

subsetFNFlag="ctrl"
subsetFNFlag=""
subsetFNFlag="case"

subsetFlag="ctrl"
subsetFlag=""
subsetFlag="case"

## ---------------------------------

##############################################

if (computerFlag=="cluster") {
	setwd("/home/royr/project/JoeWiemels")
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
	setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

subsetName=subsetFlag
if (subsetFlag!="") {
	subsetName=paste("_",subsetFlag,"Subset",sep="")
}

subsetFNName=subsetFNFlag
if (subsetFNFlag!="") {
    subsetFNName=paste("_",subsetFNFlag,"SubsetFunNorm",sep="")
}

heading=paste(c(", ",subsetFlag,", ",", ",datType,", ",subsetFNName),collapse="")
cat("\n\n============================",subsetFlag,", ",datType,", ",subsetFNName,"===========================\n\n")

##############################################

timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))

dirMethRef=""
if (computerFlag=="cluster") {
	dirClin=dirMeth="data/"
	switch(datType,
		"_allGuthSet2"={
		   dirClin=dirMeth="data/set2/"
		   fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
		   fNameClin=paste("clin_guthrieSet2_20140619.txt",sep="")
		},
		"_allGuthSet1"={
		   dirClin=dirMeth="data/set1/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
		   fNameClin=paste("final.txt",sep="")
		},
		"_allGuthSet1Set2"={
		   dirClin=dirMeth="data/set1set2/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
           #fNameClin=paste("clin_guthrieSet1Set2_20140619.txt",sep="")
           fNameClin=paste("clin_allGuthSet1Set2_20160523.txt",sep="")
        },
        "_leuk"={
            dirClin=dirMeth=dirClin2="data/set1/"
            fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
            fNameClin="i.LEU.v2.txt"
            fNameClin2="0708011 Sample_Sheet (Fetal blood)"
        },
		"_aml"={
		   dirClin=dirMeth="data/aml/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
		   fNameClin=paste("clin_aml_20150114.txt",sep="")
		}
	)
} else {
	switch(datType,
		"_allGuthSet2"={
		   dirClin=dirMeth="docs/all/set2/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
		   fNameClin=paste("clin_guthrieSet2_20140619.txt",sep="")
		},
		"_allGuthSet1"={
		   dirClin=dirMeth="docs/all/set1/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
		   fNameClin=paste("final.txt",sep="")
	   },
	   "_allGuthSet1Set2"={
		   dirClin=dirMeth="docs/all/set1set2/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
           #fNameClin=paste("clin_guthrieSet1Set2_20140619.txt",sep="")
           fNameClin=paste("clin_allGuthSet1Set2_20160523.txt",sep="")
       },
       "_leuk"={
           dirMeth="docs/all/set1/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
           dirClin="docs/all/set1/LEU.data/"
           fNameClin="i.LEU.v2.txt"
           dirClin2="docs/all/set1/preBcell/"
           fNameClin2="0708011 Sample_Sheet (Fetal blood)"
       },
		"_aml"={
		   dirClin=dirMeth="docs/aml/"
           fNameMeth=paste("beta_funNorm",datType,ifelse(subsetFNFlag=="","",paste("_",subsetFNFlag,"Subset",sep="")),".txt",sep="")
		   fNameClin=paste("clin_aml_20150114.txt",sep="")
		}
	)
}

phen=read.table(paste(dirClin,fNameClin,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
meth=read.table(paste(dirMeth,fNameMeth,sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
rownames(meth)=meth$probeId
meth=as.matrix(meth[,-1])
switch(datType,
	"_allGuthSet2"={
		names(phen)[match(c("subjectID","birth_wt"),names(phen))]=c("subjectId","birthWt")
		phen$id=paste("X",phen$guthrieId,sep="")
	},
	"_allGuthSet1"={
		names(phen)[match(c("Subject_ID","birth_weight"),names(phen))]=c("subjectId","birthWt")
		phen$id=paste("X",phen$TargetID,sep="")
        phen$caco=phen$Leukemia
	},
	"_allGuthSet1Set2"={
    },
    "_leuk"={
        names(phen)[match(c("Plate","Sentrix_ID","Sentrix_Position"),names(phen))]=c("Batch","Beadchip","Position")
        phen$id=sapply(phen$Sample,function(x) {if (is.na(as.integer(substr(x,1,1)))) x else paste("X",x,sep="")},USE.NAMES=F)
        phen$group="Leukemia"
        phen$group2="Leuk"
        
        tbl1=read.table(paste(dirClin2,fNameClin2,".csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=14)
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
        k=match(names(phen),names(tbl1)); k1=which(!is.na(k)); k2=k[k1]
        phen=rbind(phen[,k1],tbl1[!duplicated(tbl1$id),k2])
        
        x=gsub("_+","_",gsub("(","_",gsub(" |-|)","_",phen$id),fixed=T))
        j=which(substr(x,nchar(x),nchar(x))=="_")
        if (length(j)!=0) x[j]=substr(x[j],1,nchar(x[j])-1)
        phen$id=x
        
        colnames(meth)=sub("_$","",gsub("_+","_",gsub(".","_",colnames(meth),fixed=T)))
    },
"_aml"={
    phen$id=paste("X",phen$guthrieId,sep="")
	}
)
j=match(colnames(meth),phen$id); j1=which(!is.na(j)); j2=j[j1]
meth=meth[,j1]
phen=phen[j2,]
rownames(phen)=phen$id

if (subsetFlag!="") {
	switch(subsetFlag,
		   "case"={samId=which(phen$caco==1)},
		   "ctrl"={samId=which(phen$caco==0)}
	)
	meth=meth[,samId]
	phen=phen[samId,]
	rm(samId)
}

## ----------------------------------------------
if (computerFlag=="") {
    #load(file="ann.RData")
    load(file="annAll.RData")
} else {
    if (computerFlag=="cluster") {
        ann=read.delim(paste(dirMethRef,"data/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
        snpVec=read.table(paste(dirMethRef,"data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    } else {
        ann=read.delim(paste(dirMethRef,"docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
        snpVec=read.table(paste(dirMethRef,"docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    }
}
designType=as.character(ann$Infinium_Design_Type)
designType[which(designType=="I")]="1"
designType[which(designType=="II")]="2"
designType=as.integer(designType)
i=match(rownames(meth),ann[,"IlmnID"])
designType=designType[i]

save.image(paste("tmp1",datType,subsetFNName,subsetName,".RData",sep=""))

################################################
# Run BMIQ
################################################

beta=matrix(nrow=nrow(meth),ncol=ncol(meth),dimnames=list(rownames(meth),colnames(meth)))

timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))

if (computerFlag=="cluster") {
	source(paste(dirMethRef,"docs/bmiq/BMIQ_1.3.my.R",sep=""))
} else {
	source(paste(dirMethRef,"docs/bmiq/BMIQ_1.3.my.R",sep=""))
}

for (samId in 1:ncol(meth)) {
	i=!is.na(meth[,samId])
	res=BMIQ(beta.v=meth[i,samId],design.v=designType[i],nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=colnames(meth)[samId],fNameSuf=paste(datType,subsetFNName,subsetName,sep=""))
	beta[i,samId]=res$nbeta
}

timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))

print(diff(timeStamp))

write.table(cbind(probeId=rownames(beta),beta),file=paste("beta_bmiq",datType,subsetFNName,subsetName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

save.image(paste("tmp2",datType,subsetFNName,subsetName,".RData",sep=""))

x=apply(meth,2,median,na.rm=T)
print("summary(beta_funNorm)")
print(summary(x))
x1=x

x=apply(beta,2,median,na.rm=T)
print("summary(beta_bmiq)")
print(summary(x))
x2=x

tbl=cbind(id=colnames(beta),medianBeta_funNorm=round(x1,4),medianBeta_bmiq=round(x2,4))
write.table(tbl,file=paste("summaryBeta_forSample",datType,subsetFNName,subsetName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

######################################################
######################################################


######################################################
######################################################
