## Wiemel project

argv=commandArgs(TRUE)
varThis=argv[1]
subsetFlag=argv[2]

if (F) {
varThis="caco_covSexPlateCelltype"
varThis="caco_covSexPlate"
subsetFlag="dsal"

varThis="cacoXdsalperi_covSexPlate"
subsetFlag=""
}

x=strsplit(varThis,"_")[[1]]
varThis=x[1]
k=which(substr(x,1,3)=="cov")
if (length(k)==1) covThis=paste("_",x[k],sep="") else covThis=""

cat("\n\n============================ dmrcate.R ===========================\n\n",sep="")

## ---------------------------------

computerFlag=""
computerFlag="cluster"

## ---------------------------------

nProbe=101
nProbe=10001
nProbe=-1

## ---------------------------------
subsetName2=""

normFlag=""
normFlag="_funNorm"
normFlag="_bmiq"

transformFlag=""
transformFlag="_mVal"

datType="_periDsal"; subsetName2=""

computeFlag="covariance"
computeFlag="regression"


## ---------------------------------
subsetFFlag="_noEpiStr"
subsetFFlag=""

setFlag=ifelse(subsetName2=="",tolower(sub("allGuth","",datType)),subsetName2)

if (F) {
subsetFlag="hispCtrl"
subsetFlag="noHispWtCtrl"

subsetFlag="female"
subsetFlag="male"

#subsetFlag="noHisp"
subsetFlag="hisp"
subsetFlag="noHispWt"

subsetFlag="dsal"
subsetFlag="peri"

subsetFlag=""

subsetFlag="case"
subsetFlag="ctrl"
}

## ---------------------------------

varFlag="_season"; covFlag=covThis; varName=""; termName=""; modelFlag=paste("meth~",varThis,sep=""); computeFlag[2]="linear"

## ---------------------------------

varFlag="_cacoXdsalperi"; covFlag=""; varName=""; termName=c("caco","dsalPeri","interaction")

varFlag="_dsalPeri"; covFlag=covThis; varName=""; termName=""; modelFlag=paste("meth~",varThis,sep=""); computeFlag[2]="linear"

## ---------------------------------
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag=paste("meth~",varThis,sep=""); computeFlag[2]="linear"
varFlag="_caco"; covFlag="_covSexPlate"; varName=""; termName=""; modelFlag=paste("meth~",varThis,sep=""); computeFlag[2]="linear"
varFlag="_caco"; covFlag="_covSexPlateCelltype"; varName=""; termName=""; modelFlag=paste("meth~",varThis,sep=""); computeFlag[2]="linear"
## ---------------------------------

varFlag=paste("_",varThis,sep=""); covFlag=covThis; varName=""; termName=""; modelFlag=paste("meth~",varThis,sep=""); computeFlag[2]="linear"

## ---------------------------------

covPCFlag="_covPrinComp123456"
covPCFlag="_covPrinComp12"
covPCFlag=""
covPCFlag="_covPrinComp1234"

covESFlag=""
covESFlag="_covEpStr"

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
subsetFNName=""
if (subsetFlag!="") {
	subsetName=paste("_",subsetFlag,"Subset",sep="")
	switch(subsetFlag,
       "case"={varName=sub(")"," cases)",varName)},
       "ctrl"={varName=sub(")"," controls)",varName)},
       "female"={varName=sub(")"," females)",varName)},
       "male"={varName=sub(")"," males)",varName)},
       "hisp"={varName=sub(")"," hispanics)",varName)},
       "noHisp"={varName=sub(")"," non-hispanics)",varName)},
       "hyperdip"={varName=sub(")"," hyperdiploids)",varName)},
       "telaml"={varName=sub(")"," tel/aml1s)",varName)},
       "noHypTelaml"={varName=sub(")"," non-hyperdiploids/non-tel/aml1s)",varName)},
       "hispCtrl"={varName=sub(")"," hispanic controls)",varName)},
       "peri"={varName=sub(")"," perinatal)",varName)},
       "dsal"={varName=sub(")"," downsyndrome)",varName)}
    )
}

heading=paste(c(varFlag,", ",subsetFlag,", ",covFlag,", ",covPCFlag,", ",covESFlag,", ",subsetFFlag,", ",datType,subsetName2,", ",normFlag,", ",transformFlag),collapse="")
cat("\n\n============================ ",paste(computeFlag,collapse=", "),", Refactor ===========================\n\n",sep="")
cat("\n\n============================",modelFlag,", ",varFlag,", ",subsetFlag,", ",covFlag,", ",covPCFlag,", ",covESFlag,", ",subsetFFlag,", ",datType,subsetName2,", ",normFlag,", ",transformFlag,"===========================\n\n")

##############################################
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))

if (computerFlag=="cluster") {
	dirMeth="data/"
	dirClin=dirMeth
	dirBW="data/"
	dirCom=dirMeth
    dirRefactor=dirMeth
	switch(datType,
		"_allGuthSet2"={
		   dirMeth=dirClin="data/set2/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
           fNameClin="clin_allGuthSet2_20180409"
        },
		"_allGuthSet1"={
		   dirMeth=dirClin="data/set1/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),sep="")
           fNameClin="clin_allGuthSet1_20180409"
           dirMethLeuk=dirMeth
           fNameMethLeuk="clin_leuk_20180409"
        },
        "_allGuthSet1Set2"={
            dirMeth=dirClin="data/set1set2/"
            dirMeth=dirClin="data/set1set2/new/"
            fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1set2",datType),subsetFNName,subsetName,sep="")
            fNameClin="clin_allGuthSet1Set2_20181003"
        },
		"_leuk"={
		   dirMeth=dirClin=dirClin2="data/set1/"
		   fNameMeth=paste("beta",normFlag,datType,sep="")
           fNameClin="clin_leuk_20180409"
        },
		"_allGuthSet1Set2Combat"={
		   dirMeth=dirClin="data/set1set2/"
		   fNameMeth="combatAdjustedBeta_allGuthSet1Set2_set"
		   fNameClin="clin_guthrieSet1Set2_20151022"
		},
		"_aml"={
		   dirClin=dirMeth=dirRefactor="data/aml/"
		   fNameMeth=paste("beta",normFlag,"_aml",sep="")
		   fNameClin=paste("clin_aml_20150114",sep="")
        },
        "_periDsal"={
            dirMeth="/home/royr/project/JoeWiemels/normalize/results/"
            dirClin=dirRefactor="/home/royr/project/JoeWiemels/data/periDsal/"
            if (transformFlag=="_mVal") {
                fNameMeth="mValBmiq.RData"
            } else {
                fNameMeth="betaBmiq.RData"
            }
            fNameClin="sampleInfo_periDsal_20181119"
        },
		"_ivorra"={
		   dirMeth=dirClin="data/ivorra2014/"
		   fNameMeth="beta_ivorra"
		   fNameClin="clin_ivorra"
		}
	)
} else {
	dirBW="docs/birthWeight/"
	dirCom="docs/all/"
    switch(datType,
		"_allGuthSet2"={
		   dirMeth=dirClin="docs/all/set2/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
           fNameClin="clin_allGuthSet2_20180409"
           dirRefactor=dirClin
	   },
	   "_allGuthSet1"={
			dirMeth=dirClin="docs/all/set1/"
			fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),sep="")
            fNameClin="clin_allGuthSet1_20180409"
            dirMethLeuk="docs/all/set1/LEU.data/"
            fNameMethLeuk="clin_leuk_20180409"
            dirRefactor=dirClin
       },
       "_allGuthSet1Set2"={
           dirMeth=dirClin="docs/all/set1set2/"
           dirMeth=dirClin="docs/all/set1set2/new/"
           fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),subsetFNName,subsetName,sep="")
           fNameClin="clin_allGuthSet1Set2_20181003"
           dirRefactor=dirClin
       },
		"_leuk"={
		   dirMeth="docs/all/set1/"
		   fNameMeth=paste("beta",normFlag,datType,sep="")
		   dirClin="docs/all/set1/LEU.data/"
           fNameClin="clin_leuk_20180409"
           dirRefactor=dirClin
        },
		"_aml"={
		   dirClin=dirMeth="docs/aml/"
		   fNameMeth=paste("beta",normFlag,"_aml",sep="")
		   fNameClin=paste("clin_aml_20150114",sep="")
           dirRefactor="docs/aml/refactor/"
        },
        "_periDsal"={
            dirMeth="/Users/royr/UCSF/JoeWiemels/leukMeth/epic/results/"
            dirClin=dirRefactor="/Users/royr/UCSF/JoeWiemels/leukMeth/docs/periDsal/"
            if (transformFlag=="_mVal") {
                fNameMeth="mValBmiq.RData"
            } else {
                fNameMeth="betaBmiq.RData"
            }
            fNameClin="sampleInfo_periDsal_20181119"
        },
	   "_ivorra"={
		   dirMeth=dirClin="docs/misc/ivorra2014/"
		   fNameMeth="beta_ivorra"
		   fNameClin="clin_ivorra"
           dirRefactor=dirClin
       }
	)
}

if (datType%in%c("_allGuthSet1Set2") & subsetFlag!="") {
    samInfo=read.table(paste(dirClin,"summaryBeta_forSample",datType,subsetFNName,subsetName,".txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
} else if (datType%in%c("_allGuthSet2","_allGuthSet1","_allGuthSet1Set2","_allGuthSet1Set2Combat")) {
    samInfo=read.table(paste(dirCom,"summaryBeta_forSample.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
} else {
    if (exists("samInfo")) rm(samInfo)
}

switch(datType,
	"_allGuthSet2"={
		clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        
        beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
		probeId=beta$probeId
		beta=as.matrix(beta[,-1])
		rownames(beta)=probeId
		meth=beta
		rm(probeId,beta)
		
		if (normFlag=="_bmiq") {
			j=which(colnames(meth)%in%samInfo$id[which(samInfo$set=="set2" & samInfo$medianBeta_bmiq>=0.3)])
			meth=meth[,j]
		}
		
        clin2=read.table(paste(dirClin,"epistructure",datType,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        id=clin$id[!clin$id%in%clin2$id]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$id=id
        clin2=rbind(clin2,tmp)
        clin=cbind(clin,clin2[match(clin$id,clin2$id),grep("epistr",names(clin2))])
        rm(clin2,tmp,id)

		phen=clin
		phen=phen[match(colnames(meth),phen$id),]
        
        load(paste(dirRefactor,"Refactor_dat",datType,".RData",sep=""))
        colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
        j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=cbind(phen[j1,],Refactor_dat[j2,])
	},
	"_allGuthSet1"={
        clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

        beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        probeId=beta$probeId
        beta=as.matrix(beta[,-1])
        rownames(beta)=probeId
        meth=beta
        rm(probeId,beta)
        if (subsetName2=="_set2EthnProp") {
            clin2=read.table(paste(dirClin,"guthrieId_propOfCacoEthnAsInSet2_set1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            clin=clin[match(clin2$guthrieId,clin$id),]
            j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
            meth=meth[,j1]
            clin=clin[j2,]
        } else if (subsetName2=="_noNonRndChip") {
           x=table(clin$Beadchip,clin$caco)
           y=apply(x,1,function(x) {all(x>2)})
           chisq.test(x[y,])
           clin=clin[which(clin$Beadchip%in%names(y)[y]),]
           j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
           meth=meth[,j1]
           clin=clin[j2,]
           rm(x,y)
        }

        clin2=read.table(paste(dirClin,"epistructure",datType,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        id=clin$id[!clin$id%in%clin2$id]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$id=id
        clin2=rbind(clin2,tmp)
        clin=cbind(clin,clin2[match(clin$id,clin2$id),grep("epistr",names(clin2))])
        rm(clin2,tmp,id)

        phen=clin
        phen=phen[match(colnames(meth),phen$id),]

        load(paste(dirRefactor,"Refactor_dat",datType,".RData",sep=""))
        colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
        j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=cbind(phen[j1,],Refactor_dat[j2,])
	},
	"_leuk"={
		clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        
		beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        colnames(beta)=sub("_$","",gsub("_+","_",gsub(".","_",colnames(beta),fixed=T)))
		probeId=beta$probeId
		beta=as.matrix(beta[,-1])
		rownames(beta)=probeId
		meth=beta
		rm(probeId,beta)

		phen=clin

		phen=phen[match(colnames(meth),phen$id),]
	},
	"_allGuthSet1Set2"={
		clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        clin$Batch=paste(clin$set,clin$Batch)
        
        if (length(dir(dirClin,pattern=paste("epistructure",datType,subsetFNName,subsetName,".txt",sep="")))!=0) {
            clin2=read.table(paste(dirClin,"epistructure",datType,subsetFNName,subsetName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            id=clin$id[!clin$id%in%clin2$id]
            tmp=clin2[1:length(id),]
            for (k in 1:ncol(tmp)) tmp[,k]=NA
            tmp$id=id
            clin2=rbind(clin2,tmp)
            clin=cbind(clin,clin2[match(clin$id,clin2$id),grep("epistr",names(clin2))])
            rm(clin2,tmp,id)
        }
        
        if (length(dir(dirClin,pattern=paste("cellCount_minfi_cordBlood",datType,subsetFNName,subsetName,".txt",sep="")))!=0) {
            clin2=read.table(paste(dirClin,"cellCount_minfi_cordBlood",datType,subsetFNName,subsetName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            id=clin$id[!clin$id%in%clin2$id]
            tmp=clin2[1:length(id),]
            for (k in 1:ncol(tmp)) tmp[,k]=NA
            tmp$id=id
            clin2=rbind(clin2,tmp)
            clin=cbind(clin,clin2[match(clin$id,clin2$id),which(!names(clin2)%in%names(clin))])
            rm(clin2,tmp,id)
        }
        beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        probeId=beta$probeId
		beta=as.matrix(beta[,-1])
		rownames(beta)=probeId
		meth=beta
		rm(probeId,beta)
        
        if (normFlag=="_bmiq") {
            j=which(colnames(meth)%in%samInfo$id[which(samInfo$medianBeta_bmiq>=0.3)])
            meth=meth[,j]
        }
        
		j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
		meth=meth[,j1]
		clin=clin[j2,]
		if (subsetName2!="") {
			clin=clin[which(clin$set==sub("_","",subsetName2)),]
			j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
			meth=meth[,j1]
			clin=clin[j2,]
		}
		phen=clin
		phen=phen[match(colnames(meth),phen$id),]
        
        if (length(dir(dirRefactor,pattern=paste("Refactor_dat",datType,subsetFNName,subsetName,".RData",sep="")))!=0) {
            load(paste(dirRefactor,"Refactor_dat",datType,subsetFNName,subsetName,".RData",sep=""))
            colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
            j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
            meth=meth[,j1]
            phen=cbind(phen[j1,],Refactor_dat[j2,])
        }
	},
	"_allGuthSet1Set2Combat"={
		clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
		probeId=beta$cpgId
		beta=as.matrix(beta[,-1])
		rownames(beta)=probeId
		meth=beta
		rm(probeId,beta)
		j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
		meth=meth[,j1]
		clin=clin[j2,]
		if (subsetName2!="") {
			clin=clin[which(clin$set==sub("_","",subsetName2)),]
			j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
			meth=meth[,j1]
			clin=clin[j2,]			
		}
		phen=clin
		phen=phen[match(colnames(meth),phen$id),]
	},
	"_aml"={
		clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		clin$id=paste("X",clin$guthrieId,sep="")
		clin$int_ch_ethnicity=clin$ch_ethnicity
		beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
		probeId=beta$probeId
		beta=as.matrix(beta[,-1])
		rownames(beta)=probeId
		meth=beta
		rm(probeId,beta)
		j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
		meth=meth[,j1]
		clin=clin[j2,]
		if (subsetName2!="") {
			clin=clin[which(clin$set==sub("_","",subsetName2)),]
			j=match(colnames(meth),clin$id); j1=which(!is.na(j)); j2=j[j1]
			meth=meth[,j1]
			clin=clin[j2,]
		}
		
		clin2=read.table(paste(dirClin,"clin_aml_20151002+SNPeigenvalues.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		clin=cbind(clin,clin2[match(clin$guthrieId,clin2$guthrieId),grep("SNP_PC",names(clin2))])
		rm(clin2)
        
        clin2=read.table(paste(dirClin,"epistructure",datType,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        id=clin$id[!clin$id%in%clin2$id]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$id=id
        clin2=rbind(clin2,tmp)
        clin=cbind(clin,clin2[match(clin$id,clin2$id),grep("epistr",names(clin2))])
        rm(clin2,tmp,id)
        
		phen=clin
		phen=phen[match(colnames(meth),phen$id),]
		
		load(paste(dirRefactor,"Refactor_dat_aml.RData",sep=""))
		colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
		j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
		meth=meth[,j1]
		phen=cbind(phen[j1,],Refactor_dat[j2,])
    },
    "_periDsal"={
        phen=read.table(paste(dirClin,fNameClin,".txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
        names(phen)[match(c("plate"),names(phen))]=c("Batch")
        phen$caco=as.integer(phen$caco=="case")
        
        load(paste(dirMeth,fNameMeth,sep=""))
        if (transformFlag=="_mVal") {
            meth=mVal
            rm(mVal)
        } else {
            meth=betaBmiq
            rm(betaBmiq)
        }
        if (nProbe!=-1) meth=meth[1:nProbe,]
        
        j=match(colnames(meth),phen$id); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=phen[j2,]
        
        samId=which(phen$quality!="3. low from bmiq" & !is.na(phen$periDsal) & !is.na(phen$caco))
        meth=meth[,samId]
        phen=phen[samId,]
        
        if (subsetFlag%in%c("peri","dsal")) {
            clin2=read.table(paste(dirClin,"epistructure",datType,subsetName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        } else {
            clin2=read.table(paste(dirClin,"epistructure",datType,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        }
        id=phen$id[!phen$id%in%clin2$id]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$id=id
        clin2=rbind(clin2,tmp)
        phen=cbind(phen,clin2[match(phen$id,clin2$id),grep("epistr",names(clin2))])
        rm(clin2,tmp,id)
        
        if (subsetFlag%in%c("peri","dsal")) {
            load(paste(dirRefactor,"Refactor_dat",datType,subsetName,".RData",sep=""))
        } else {
            load(paste(dirRefactor,"Refactor_dat",datType,".RData",sep=""))
        }
         colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
        j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=cbind(phen[j1,],Refactor_dat[j2,])
    },
    "_ivorra"={
    clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		
		beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
		probeId=beta$probeId
		beta=as.matrix(beta[,-1])
		rownames(beta)=probeId
		meth=beta
		rm(probeId,beta)
		
		phen=clin
		phen=phen[match(colnames(meth),phen$id),]
	}
)
rownames(phen)=phen$id
if (length(grep("_allGuth",datType)!=0)) {
	phen$smoke_mo_3months=as.integer(phen$smoke_mo_3months_N>0)
	phen$smoke_mo_after=as.integer(phen$smoke_mo_after_N>0)
	phen$smoke_mo_bf=as.integer(phen$smoke_mo_bf_N>0)
}
if ("birthWt"%in%names(phen)) {
    cutoff=median(phen$birthWt,na.rm=T)
    cutoff=3500 ## The median between 3400 & 3500 for set1 ctrl, set2 ctrl & set 2 caco datasets
    cat("Birth weight cutoff: ",cutoff,"\n\n",sep="")
    print("summary(phen$birthWt)")
    print(summary(phen$birthWt))
    phen$birthWtBi=as.integer(phen$birthWt>cutoff)
}
if ("pobw"%in%names(phen)) {
    cutoff=median(phen$pobw,na.rm=T)
    cutoff=1 ## The median is ~ 1 for set1 ctrl, set2 ctrl & set 2 caco datasets
    cat("POBW cutoff: ",cutoff,"\n\n",sep="")
    print("summary(phen$pobw)")
    print(summary(phen$pobw))
    phen$pobwBi=as.integer(phen$pobw>cutoff)
}
if ("motherWt"%in%names(phen)) {
    phen$moBmi=phen$motherWt/((phen$motherHt/100)^2)
    phen$moBmiCat=NA
    phen$moBmiCat[which(round(phen$moBmi,1)>=30)]="obese"
    phen$moBmiCat[which(round(phen$moBmi,1)>=25 & round(phen$moBmi,1)<30)]="overweight"
    phen$moBmiCat[which(round(phen$moBmi,1)>=18.5 & round(phen$moBmi,1)<25)]="normal"
    phen$moBmiCat[which(round(phen$moBmi,1)<18.5)]="underweight"
    phen$moBmiOverOrObese=phen$moBmiOver=phen$moBmiObese=phen$moBmiUnder=phen$moBmiCat
    phen$moBmiOverOrObese=as.integer(phen$moBmiCat%in%c("overweight","obese"))
    phen$moBmiOverOrObese[!phen$moBmiCat%in%c("normal","overweight","obese")]=NA
    phen$moBmiOver=as.integer(phen$moBmiCat%in%c("overweight"))
    phen$moBmiOver[!phen$moBmiCat%in%c("normal","overweight")]=NA
    phen$moBmiObese=as.integer(phen$moBmiCat%in%c("obese"))
    phen$moBmiObese[!phen$moBmiCat%in%c("normal","obese")]=NA
    phen$moBmiUnder=as.integer(phen$moBmiCat%in%c("underweight"))
    phen$moBmiUnder[!phen$moBmiCat%in%c("normal","underweight")]=NA
}
if (length(grep("PCB",names(phen)))!=0) {
    pcb=read.table(paste(dirCom,"pcb.txt",sep=""), sep=" ", h=F, quote="", comment.char="",as.is=T,fill=T,skip=5)
    names(pcb)=c("pcbNo","clPos","aroclor1016","aroclor1242","aroclor1248A3.5","aroclor1248G3.5","aroclor1254Late","aroclor1254","aroclor1260")
    for (k in 1:ncol(pcb)) if (is.character(pcb[,k])) pcb[,k]=gsub("\"","",pcb[,k])
    for (k in c("pcbNo")) pcb[,k]=as.integer(pcb[,k])
    for (k in c("aroclor1016","aroclor1242","aroclor1248A3.5","aroclor1248G3.5","aroclor1254Late","aroclor1254","aroclor1260")) pcb[,k]=as.numeric(pcb[,k])
    k=grep("logged_PCB",names(phen))
    dat=as.matrix(phen[,k])
    x=as.integer(gsub("logged_PCB_|_SRS","",colnames(dat)))
    x=x[!is.na(x)]
    phen$logged_PCB_aroclor1260=dat%*%pcb$aroclor1260[match(x,pcb$pcbNo)]
}

if (subsetFlag!="") {
	switch(subsetFlag,
		   "case"={samId=which(phen$caco==1)},
		   "ctrl"={samId=which(phen$caco==0)},
           "female"={samId=which(phen$sex==2)},
           "male"={samId=which(phen$sex==1)},
           "hisp"={samId=which(phen$int_ch_ethnicity==1)},
		   "noHispWt"={samId=which(phen$int_ch_ethnicity==2)},
		   "noHisp"={
			   samId=which(phen$subtype==2)
			   samId=NULL
		   },
		   "hyperdip"={samId=which(phen$subtype=="hyperdiploid")},
		   "telaml"={samId=which(phen$subtype=="telaml")},
		   "noHypTelaml"={samId=which(phen$subtype=="nonHypTelaml")},
           "hispCtrl"={samId=which(phen$int_ch_ethnicity==1 & phen$caco==0)},
           "noHispWtCtrl"={samId=which(phen$int_ch_ethnicity==2 & phen$caco==0)},
           "peri"={samId=which(phen$periDsal=="peri")},
           "dsal"={samId=which(phen$periDsal=="dsal")}
           )
	meth=meth[,samId]
	phen=phen[samId,]
	rm(samId)
}
if (subsetFFlag!="") {
    switch(subsetFlag,
    "_noEpiStr"={varName=sub(")"," cases)",varName)}
    )
    meth=meth[,samId]
    phen=phen[samId,]
    rm(samId)
}

if (length(grep("_allGuth",datType)!=0)) {
	phen$subtype[which(phen$caco==0)]="control"

	phen$caco=as.factor(phen$caco)
	phen$Beadchip=as.factor(phen$Beadchip)

	phen$race=as.factor(phen$int_ch_race)
	phen$ethn=as.factor(phen$int_ch_ethnicity)
	phen$sex=as.factor(phen$sex)

	phen$race3=as.integer(as.character(phen$race))
	phen$race3[which(phen$race3%in%c(2,3))]=5
	phen$race3=as.factor(phen$race3)

    if (length(grep("Race3",covFlag))==1) {
		phen$race2=phen$race3
    } else if (length(grep("Hisp",covFlag))==1) {
		phen$race2=phen$ch_hispanic_bc
	} else {
		phen$race2=phen$ethn
	}
}
if (datType%in%c("_leuk","_aml")) {
    if ("int_ch_race"%in%names(phen)) phen$race=as.factor(phen$int_ch_race)
    if ("int_ch_ethnicity"%in%names(phen)) phen$ethn=as.factor(phen$int_ch_ethnicity)
    phen$race2=phen$ethn
}

for (k in which(names(phen)%in%c("caco","sex","Beadchip","race","ethn","periDsal"))) {
	phen[,k]=as.factor(phen[,k])
}

if (exists("clin")) rm(clin)
switch(varFlag,
	"_rasChip"={
	   phen$caco=phen$RAS
	},
	"_hyperdipChip"={
		phen$caco=NA
		phen$caco[which(phen$Subtype%in%c("mll","others","t119","t1221"))]=0
		phen$caco[which(phen$Subtype=="hyperdiploid")]=1
	},
	"_hyperdipTelamlChip"={
		phen$caco=NA
		phen$caco[which(phen$Subtype%in%c("t1221"))]=0
		phen$caco[which(phen$Subtype=="hyperdiploid")]=1
	},
	"_hyperdipCtrl"={
		phen$caco=NA
		phen$caco[which(phen$subtype=="control")]=0
		phen$caco[which(phen$subtype=="hyperdiploid")]=1
	},
	"_telamlCtrl"={
	   phen$caco=NA
	   phen$caco[which(phen$subtype=="control")]=0
	   phen$caco[which(phen$subtype=="telaml")]=1
	},
	"_noHypTelamlCtrl"={
	   phen$caco=NA
	   phen$caco[which(phen$subtype=="control")]=0
	   phen$caco[which(phen$subtype=="nonHypTelaml")]=1
    },
    "_telamlOther"={
           phen$caco=NA
           phen$caco[which(phen$group=="Leukemia")]=0
           phen$caco[which(phen$subtype=="telaml")]=1
    },
	"_leukPreb"={
		phen$caco=NA
		phen$caco[which(phen$group%in%c("S2","S3","S4"))]=0
		phen$caco[which(phen$group=="Leukemia")]=1
	},
	"_leukSmokeChip"={
		phen$smoke=NA
		phen$smoke[which(phen$smoke_mo_preg==0 & phen$smoke_mo_3months_N==0 & phen$smoke_fa_3months==0)]=0
		phen$smoke[which(phen$smoke_mo_preg==1)]=1
		phen$smoke[which(phen$smoke_mo_3months_N>0 & phen$smoke_mo_preg==0)]=2
		phen$smoke[which(phen$smoke_fa_3months==1 & phen$smoke_mo_3months_N==0 & phen$smoke_mo_preg==0)]=3
		phen$smoke[which(phen$smoke==2)]=1
		phen$smoke=factor(phen$smoke)
		phen$caco=phen$smoke
		
		clin$smoke=NA
		clin$smoke[which(clin$moPreg==0 & clin$mo3m==0 & clin$fa3m==0)]=0
		clin$smoke[which(clin$moPreg==1)]=1
		clin$smoke[which(clin$mo3m==1 & clin$moPreg==0)]=2
		clin$smoke[which(clin$fa3m==1 & clin$mo3m==0 & clin$moPreg==0)]=3
		clin$smoke[which(clin$smoke==2)]=1
		clin$smoke=factor(clin$smoke)

	},
	"_hyperdipTelaml"={
		phen$caco=phen$hyperdipTelaml
	},
	"_hyperdipNonHypTelaml"={
		phen$caco=phen$hyperdipNonHypTelaml
	},
	"_telamlNonHypTelaml"={
		phen$caco=phen$telamlNonHypTelaml
	},
	"_dfeFood"={
		phen$caco=phen$DFE_Food
	},
	"_dfeFort"={
		phen$caco=phen$DFE_fort
	},
	"_dfeNat"={
		phen$caco=phen$DFE_nat
	},
	"_dfeTot"={
		phen$caco=phen$DFE_tot
	},
	"_dfeFoodCat"={
		phen$caco=as.integer(phen$DFE_Food>=400)
	},
	"_dfeNatCat"={
		phen$caco=as.integer(phen$DFE_nat>=200)
	},
	"_dfeSupCat"={
		phen$caco=as.integer(phen$DFE_sup>0)
	},
	"_dfeTotCat"={
		phen$caco=as.integer(phen$DFE_tot>=400)
	},
	"_smoke3"={
		## Smoking 3 categories
		phen$caco=NA
		phen$caco[which(phen$smoke_mo_preg==0 & phen$smoke_mo_3months==0 & phen$smoke_fa_3months==0)]=0
		phen$caco[which(phen$smoke_mo_preg==1)]=1
		phen$caco[which(phen$smoke_mo_3months>0 & phen$smoke_mo_preg==0)]=2
		phen$caco[which(phen$smoke_fa_3months==1 & phen$smoke_mo_3months==0 & phen$smoke_mo_preg==0)]=3
		phen$caco[which(phen$caco==2)]=1
	},
	"_smoke2"={
		## Smoking 2 categories
		phen$caco=NA
		phen$caco[which(phen$smoke_mo_preg==0 & phen$smoke_mo_3months==0 & phen$smoke_fa_3months==0)]=0
		phen$caco[which(phen$smoke_mo_preg==1)]=1
		phen$caco[which(phen$smoke_mo_3months>0 & phen$smoke_mo_preg==0)]=2
		phen$caco[which(phen$smoke_fa_3months==1 & phen$smoke_mo_3months==0 & phen$smoke_mo_preg==0)]=3
		phen$caco[which(phen$caco==2)]=1
		phen$caco[which(phen$caco==3)]=1
	},
	"_smoke"={
		phen$caco=phen$smoke
	},
	"_birthWt"={
		phen$caco=phen$birthWt
	},
	"_pobw"={
		phen$caco=phen$pobw
	},
	"_cacoXbirthwt"={
		phen$predVar=phen$birthWt
	},
	"_cacoXpobw"={
		phen$predVar=phen$pobw
    },
    "_chipPos"={phen$caco=phen$Position
    },
	"_season"={
        phen$caco=phen$season
        phen$caco[which(phen$season=="autumn")]="1_autumn"
        phen$caco[which(phen$season=="winter")]="2_winter"
        phen$caco[which(phen$season=="spring")]="3_spring"
        phen$caco[which(phen$season=="summer")]="4_summer"
    },
    "_dsalPeri"={phen$caco=as.integer(phen$periDsal=="dsal")
    },
    "_cacoXdsalperi"={
        phen$predVar=as.integer(phen$periDsal=="dsal")
    }
)

switch(varThis,
    "telamlCtrl"={
           phen$telamlCtrl=NA
           phen$telamlCtrl[which(phen$subtype=="control")]=0
           phen$telamlCtrl[which(phen$subtype=="telaml")]=1
    },
    "telamlOther"={
        phen$telamlOther=NA
        phen$telamlOther[which(phen$group=="Leukemia")]=0
        phen$telamlOther[which(phen$subtype=="telaml")]=1
    },
    "telamlNonTelaml"={
        phen$telamlNonTelaml=as.integer(phen$subtype=="telaml")
        phen$telamlNonTelaml[phen$caco!=1 | phen$subtype==""]=NA
    },
    "hyperdipTelaml"={
        phen$hyperdipTelaml=as.integer(phen$subtype=="hyperdiploid")
        phen$hyperdipTelaml[!phen$subtype%in%c("hyperdiploid","telaml")]=NA
    },
    "dsalPeri"={phen$dsalPeri=as.integer(phen$periDsal=="dsal")
    }
)

if ("caco"%in%names(phen) & !varFlag%in%c("_dfeFood","_dfeFort","_dfeNat","_dfeTot",paste("_smoke_",c("mo3mN","moPregN","moAfterN","fa3mN","moBfN"),sep=""),"_smoke","_birthWt","_pobw","_cacoVbirthWt","_cacoXpobw","_cacoXdsalperi")) {
	phen$caco=as.factor(phen$caco)
}

## ----------------------------------------------
if (datType=="_periDsal") {
    load(paste(dirMeth,"ann850k.RData",sep=""))
    ann=as.data.frame(ann850k)
    rm(ann850k)
    rownames(ann)=NULL
    names(ann)[match(c("Name","chr","pos"),names(ann))]=c("IlmnID","CHR","MAPINFO")
    ann$CHR[which(ann$CHR=="chrX")]="chr23"
    ann$CHR[which(ann$CHR=="chrY")]="chr24"
    ann$CHR=as.integer(sub("chr","",ann$CHR))
    ann$snp=as.integer(!is.na(ann$Probe_rs))
    
    i=match(rownames(meth),ann[,"IlmnID"])
    table(is.na(i))
    ann=ann[i,]

    ann$keep=TRUE
    
    ann$geneSym=sapply(toupper(ann$UCSC_RefGene_Name),function(x) {
        strsplit(x,";")[[1]][1]
    },USE.NAMES=F)
    ann$geneSym[is.na(ann$geneSym)]=""
} else {
    if (computerFlag=="") {
        #load(file="ann.RData")
        load(file="annAll.RData")
    } else {
        if (computerFlag=="cluster") {
            ann=read.delim(paste("data/","HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
            snpVec=read.table(paste("data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        } else {
            ann=read.delim(paste("docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
            snpVec=read.table(paste("docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        }
        ann[which(ann[,"CHR"]=="X"),"CHR"]="23"
        ann[which(ann[,"CHR"]=="Y"),"CHR"]="24"
        ann[,"CHR"]=as.integer(ann[,"CHR"])
        ann=ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
        for (k in 1:ncol(ann)) if (class(ann[,k])=="factor") ann[,k]=as.character(ann[,k])
        
        snpVec=snpVec[,1]
        ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
    }
    ann$geneSym=sapply(toupper(ann$UCSC_RefGene_Name),function(x) {
        strsplit(x,";")[[1]][1]
    },USE.NAMES=F)
    ann$geneSym[is.na(ann$geneSym)]=""
    
    i=match(rownames(meth),ann[,"IlmnID"])
    table(is.na(i))
    ann=ann[i,]
    
    #ann$keep=ann$CHR%in%1:22
    ann$keep=ann$CHR%in%1:22 & apply(meth,1,function(x) {any(!is.na(x))})
    ann$keep=ann$snp==0 & ann$CHR%in%1:22 & apply(meth,1,function(x) {any(!is.na(x))})
    ann$keep=apply(meth,1,function(x) {any(!is.na(x))})
}

## ----------------------------------------------

if ("sem"%in%names(phen) & datType=="_allGuthSet2") {
    phen$logSemNoSnpNoSexChr=phen$semNoSnpNoSexChr=logSemNoSexChr=phen$semNoSexChr=phen$logSem=phen$sem=NA
    x=apply(sem,2,function(x) {sum(!is.na(x))})
    j=match(phen$id,names(x)); j1=which(!is.na(j)); j2=j[j1]
    phen$sem[j1]=x[j2]
    phen$logSem[j1]=log(x[j2])
    x=apply(sem[which(ann$CHR%in%1:22),],2,function(x) {sum(!is.na(x))})
    j=match(phen$id,names(x)); j1=which(!is.na(j)); j2=j[j1]
    phen$semNoSexChr[j1]=x[j2]
    phen$logSemNoSexChr[j1]=log(x[j2])
    x=apply(sem[which(ann$snp==0 & ann$CHR%in%1:22),],2,function(x) {sum(!is.na(x))})
    j=match(phen$id,names(x)); j1=which(!is.na(j)); j2=j[j1]
    phen$semNoSnpNoSexChr[j1]=x[j2]
    phen$logSemNoSnpNoSexChr[j1]=log(x[j2])
}

save.image(paste(paste("tmp_dmrcate_",sub("~","Resp_",gsub("*","X",gsub("+","_",modelFlag,fixed=T),fixed=T)),subsetName,covFlag,covPCFlag,covESFlag,subsetFFlag,datType,subsetName2,normFlag,transformFlag,".RData",sep=""),sep=""))

datObj=list(meth=meth,phen=phen)

################################################
# Regression
################################################
library(DMRcate)

colId=2
if (varFlag%in%c("_cacoXbirthwt","_cacoXpobw","_cacoXdsalperi")) {
    ## CHECK !!!
    model1="~caco*predVar"
} else if (varFlag%in%c("_cacoXhlaB")) {
    ## CHECK !!!
    model1="~caco*hlaB"
    #colId=c("caco1","hlaBallele1v0","hlaBallele2v0","hlaBallele2v1","caco1:hlaBallele1v0","caco1:hlaBallele2v0","caco1:hlaBallele2v1"),colnames(Xunadj))
    
} else {
    ## CHECK !!!
    model1="~caco"
    model1=paste("~",sub("~","+",modelFlag),sep="")
    k=grep("Sex",covFlag)
    if (length(k)==1) {
        model1=paste(model1,"+sex",sep="")
    }
    k=grep("Race3|Ethn|Hisp",covFlag)
    if (length(k)==1) {
        model1=paste(model1,"+race2",sep="")
    }
    k=grep("Gestage",covFlag)
    if (length(k)==1) {
        model1=paste(model1,"+gestage",sep="")
    }
    k=grep("Plate",covFlag)
    if (length(k)==1) {
        model1=paste(model1,"+Batch",sep="")
    }
    k=grep("Celltype",covFlag)
    if (length(k)==1) {
        #model1=paste(model1,"+nRBC+CD8T+CD4T+NK+Bcell+Mono+Gran",sep="")
        model1=paste(model1,"+nRBC+CD8T+CD4T+NK+Bcell+Mono",sep="")
    }
    k=grep("SnpPC",covFlag)
    if (length(k)==1) {
        x=strsplit(covFlag,"SnpPC")[[1]][2]
        for (k2 in 1:nchar(x)) {
            x2=substr(x,k2,k2)
            if (is.na(as.integer(x2))) {
                break
            } else {
                model1=paste(model1,"+SNP_PC",x2,sep="")
            }
        }
    }
    k=grep("Set",covFlag)
    if (length(k)==1) {
        model1=paste(model1,"+set",sep="")
    }
}
if (covPCFlag!="") {
    model1=paste(model1,"+",paste("prinComp",strsplit(sub("_covPrinComp","",covPCFlag),"")[[1]],collapse="+",sep=""),sep="")
}
if (covESFlag!="") {
    model1=paste(model1,"+",paste("epistr",1:5,collapse="+",sep=""),sep="")
}
model2=sub("meth*","",sub("meth+","",model1,fixed=T),fixed=T)
model2=as.formula(model2)
Xunadj=model.matrix(model2,data=phen)
if (varFlag%in%c("_cacoXbirthwt","_cacoXpobw","_cacoXdsalperi")) {
    #colId=match(c("caco1","predVar","caco1:predVar"),colnames(Xunadj))
    colId=match(c("caco1:predVar"),colnames(Xunadj))
}
samId=match(rownames(Xunadj),phen$id)
meth=meth[ann$keep,samId]
phen=phen[samId,]

if (transformFlag=="_mVal") {
    if (datType%in%c("_allGuthSet1","_allGuthSet2","_allGuthSet1Set2_ctrlSubset")) {
        methInfo=read.table(paste(dirCom,"summaryBeta.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
        k=which(methInfo$cohort==datType)
        meth[which(meth==0)]=methInfo$minBeta*0.5; meth[which(meth==1)]=methInfo$maxBeta+(1-methInfo$maxBeta)*0.5
        cat("Min meth:",min(meth),min(meth,na.rm=T),"\n")
        cat("Max meth:",max(meth),max(meth,na.rm=T),"\n")
        meth=log2(meth/(1-meth))
        
    } else if (datType%in%c("_periDsal")) {
    } else {
        cat("Require m-value data !!! \n\n")
        meth=NULL
    }
}

## ---------------------------
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))

meth <- rmSNPandCH(meth, dist=2, mafcut=0.05, rmXY=T)
myannotation <- cpg.annotate("array", meth, what="M", arraytype = "EPIC",analysis.type="differential", design=Xunadj, coef=colId)

save.image(paste(paste("tmp1_dmrcate_",sub("~","Resp_",gsub("*","X",gsub("+","_",modelFlag,fixed=T),fixed=T)),subsetName,covFlag,covPCFlag,covESFlag,subsetFFlag,datType,subsetName2,normFlag,transformFlag,".RData",sep=""),sep=""))

dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)
save.image(paste(paste("tmp2_dmrcate_",sub("~","Resp_",gsub("*","X",gsub("+","_",modelFlag,fixed=T),fixed=T)),subsetName,covFlag,covPCFlag,covESFlag,subsetFFlag,datType,subsetName2,normFlag,transformFlag,".RData",sep=""),sep=""))

results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges

timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))

save.image(paste(paste("tmp3_dmrcate_",sub("~","Resp_",gsub("*","X",gsub("+","_",modelFlag,fixed=T),fixed=T)),subsetName,covFlag,covPCFlag,covESFlag,subsetFFlag,datType,subsetName2,normFlag,transformFlag,".RData",sep=""),sep=""))

tbl=as.data.frame(results.ranges)
write.table(tbl,file=paste("dmrcate_",sub("~","Resp_",gsub("*","X",gsub("+","_",modelFlag,fixed=T),fixed=T)),subsetName,covFlag,covPCFlag,covESFlag,subsetFFlag,datType,subsetName2,normFlag,transformFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

