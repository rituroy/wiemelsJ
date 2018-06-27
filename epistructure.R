## Use subsetted clinical & refactor data instead of all samples
## Before 04/03/18, used files *_allGuthSet1Set2_ctrlSubset* instead of *_allGuthSet1Set2_ctrlSubsetFunNorm_ctrlSubset* files
## and summaryBeta_forSample for allGuthSet1Set2/ctrlSubset
## ---------------------------------

computerFlag=""
computerFlag="cluster"

## ---------------------------------

nProbe=101
nProbe=-1
nProbe=10001

## ---------------------------------
subsetName2=""

normFlag=""
normFlag="_funNorm"
normFlag="_bmiq"

datType="_leuk"; subsetName2=""
datType="_aml"; subsetName2=""
datType="_allGuthSet1"; subsetName2=""
datType="_allGuthSet2"; subsetName2=""
datType="_allGuthSet1Set2"; subsetName2=""

setFlag=ifelse(subsetName2=="",tolower(sub("allGuth","",datType)),subsetName2)

#subsetFlag="noHisp"
subsetFlag="noHispWt"
subsetFlag="hisp"

subsetFlag=""

subsetFlag="case"
subsetFlag="ctrl"

covFlag=""
mediationFlag=F
varFlag=""
covPCFlag="_covPrinComp1234"
winsorTail=0.05

library(limma)
library(FactoMineR)

#for (subsetName2 in c("_set1","_set2")) {
#for (datType in c("_allGuthSet1","_allGuthSet2")) {
		
## ---------------------------------

## All smoke variables from set1 & set2
## All smoke variables in set2 are in set1

## ---------------------------------

## ---------------------------------

#for (subsetFlag in c("noHisp","hisp")) {

#	for (varFlag in c("_hyperdipCtrl","_telamlCtrl","_noHypTelamlCtrl")) {
#for (varFlag in c("_dfeFoodCat","_dfeSupCat","_dfeNat","_dfeTot","_dfeFood","_dfeFort","_dfeNatCat","_dfeTotCat")) {
#for (varFlag in paste("_smoke_",varSmoke$varOut,sep="")) {

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
	#varFlag=paste(varFlag,"_",subsetFlag,"Subset",sep="")
	subsetName=paste("_",subsetFlag,"Subset",sep="")
    if (F) {
        switch(subsetFlag,
		   "case"={varName=sub(")"," cases)",varName)},
		   "ctrl"={varName=sub(")"," controls)",varName)},
		   "hisp"={varName=sub(")"," hispanics)",varName)},
		   "noHisp"={varName=sub(")"," non-hispanics)",varName)},
		   "hyperdip"={varName=sub(")"," hyperdiploids)",varName)},
		   "telaml"={varName=sub(")"," tel/aml1s)",varName)},
		   "noHypTelaml"={varName=sub(")"," non-hyperdiploids/non-tel/aml1s)",varName)}
        )
    }
    if (datType=="_allGuthSet1Set2") subsetFNName=paste("_",subsetFlag,"SubsetFunNorm",sep="")
}

heading=paste(c(varFlag,", ",subsetFlag,", ",covFlag,", ",covPCFlag,", ",datType,subsetName2,", ",normFlag),collapse="")
cat("\n\n============================ Epistructure ===========================\n\n")
cat("\n\n============================",varFlag,", ",subsetFlag,", ",covFlag,", ",covPCFlag,", ",datType,subsetName2,", ",normFlag,"===========================\n\n")

##############################################

timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))

if (computerFlag=="cluster") {
	dirMeth="data/"
	dirClin=dirMeth
	dirBW="data/"
	dirCom=dirMeth
	dirRefactor=dirEpistructure=dirMeth
	switch(datType,
		"_allGuthSet2"={
		   dirMeth=dirClin="data/set2/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),subsetName,sep="")
		   fNameClin="clin_guthrieSet2_20140619"
           fNameClin="clin_allGuthSet2_20160928"
		},
		"_allGuthSet1"={
		   dirMeth=dirClin="data/set1/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),subsetName,sep="")
		   fNameClin="final"
           fNameClin="clin_allGuthSet1_20160928"
		   dirMethLeuk=dirMeth
		   fNameMethLeuk="i.LEU.v2"
        },
        "_allGuthSet1Set2"={
            dirMeth=dirClin="data/set1set2/"
            fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1set2",datType),subsetFNName,subsetName,sep="")
            fNameClin="clin_guthrieSet1Set2_20140619"
            fNameClin="clin_guthrieSet1Set2_20151022"
            fNameClin="clin_allGuthSet1Set2_20160523"
        },
		"_leuk"={
		   dirMeth=dirClin=dirClin2="data/set1/"
		   fNameMeth=paste("beta",normFlag,datFlag,subsetName,sep="")
		   fNameClin="i.LEU.v2"
		   fNameClin2="0708011 Sample_Sheet (Fetal blood)"
		},
		"_allGuthSet1Set2Combat"={
		   dirMeth=dirClin="data/set1set2/"
		   fNameMeth=paste("combatAdjustedBeta_allGuthSet1Set2",subsetName,"_set",sep="")
		   fNameClin="clin_guthrieSet1Set2_20140619"
		},
		"_aml"={
		   dirClin=dirMeth=dirRefactor="data/aml/"
		   fNameMeth=paste("beta",normFlag,"_aml",subsetName,sep="")
		   fNameClin=paste("clin_aml_20150114",sep="")
		},
		"_ivorra"={
		   dirMeth=dirClin="data/ivorra2014/"
		   fNameMeth=paste("beta_ivorra",subsetName,sep="")
		   fNameClin="clin_ivorra"
		}
	)
} else {
	dirBW="docs/birthWeight/"
	dirCom="docs/all/"
    #dirRefactor="docs/SemiraGonsethNussle/refactor/"
    dirEpistructure="docs/SemiraGonsethNussle/epistructure/"
    switch(datType,
		"_allGuthSet2"={
		   dirMeth=dirClin="docs/all/set2/"
           dirRefactor=dirClin
           fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),subsetName,sep="")
		   fNameClin="clin_guthrieSet2_20140619"
           fNameClin="clin_allGuthSet2_20160928"
	   },
	   "_allGuthSet1"={
			dirMeth=dirClin="docs/all/set1/"
            dirRefactor=dirClin
            fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),subsetName,sep="")
			fNameClin="final"
            fNameClin="clin_allGuthSet1_20160928"
			dirMethLeuk="docs/all/set1/LEU.data/"
			fNameMethLeuk="i.LEU.v2"
       },
       "_allGuthSet1Set2"={
           dirMeth=dirClin="docs/all/set1set2/"
           dirRefactor=dirClin
           fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),subsetFNName,subsetName,sep="")
           fNameClin="clin_allGuthSet1Set2_20160523"
       },
       "_leuk"={
		   dirMeth="docs/all/set1/"
           dirRefactor=dirClin
           fNameMeth=paste("beta",normFlag,datFlag,subsetName,sep="")
		   dirClin="docs/all/set1/LEU.data/"
		   fNameClin="i.LEU.v2"
		   dirClin2="docs/all/set1/preBcell/"
		   fNameClin2="0708011 Sample_Sheet (Fetal blood)"
		},
		"_aml"={
		   dirClin=dirMeth="docs/aml/"
		   dirRefactor="docs/aml/refactor/"
		   fNameMeth=paste("beta",normFlag,"_aml",subsetName,sep="")
		   fNameClin=paste("clin_aml_20150114",sep="")
	   },
	   "_ivorra"={
		   dirMeth=dirClin="docs/misc/ivorra2014/"
           dirRefactor=dirClin
           fNameMeth=paste("beta_ivorra",subsetName,sep="")
		   fNameClin="clin_ivorra"
		}
	)
}

if (datType%in%c("_allGuthSet1Set2") & subsetFlag!="") {
    samInfo=read.table(paste(dirClin,"summaryBeta_forSample",datType,subsetFNName,subsetName,".txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
} else if (datType%in%c("_allGuthSet2","_allGuthSet1","_allGuthSet1Set2","_allGuthSet1Set2Combat")) {
    samInfo=read.table(paste(dirCom,"summaryBeta_forSample.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
}

switch(datType,
	"_allGuthSet2"={
        clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        #names(clin)[match(c("subjectID","birth_wt"),names(clin))]=c("subjectId","birthWt")
        clin$id=paste("X",clin$guthrieId,sep="")
        #clin$Leukemia=clin$caco
        
        clin$subtype=rep("",nrow(clin))
        clin$subtype[which(clin$smhyper==1)]="hyperdiploid"
        clin$subtype[which(clin$smtelaml==1)]="telaml" ## include subject with smhyper=1
        clin$subtype[which(clin$smhyper==0 & clin$smtelaml==0)]="nonHypTelaml"
        clin$hyperdipTelaml=as.integer(clin$subtype=="hyperdiploid")
        clin$hyperdipTelaml[!clin$subtype%in%c("hyperdiploid","telaml")]=NA
        clin$hyperdipNonHypTelaml=as.integer(clin$subtype=="hyperdiploid")
        clin$hyperdipNonHypTelaml[!clin$subtype%in%c("hyperdiploid","nonHypTelaml")]=NA
        clin$telamlNonHypTelaml=as.integer(clin$subtype=="telaml")
        clin$telamlNonHypTelaml[!clin$subtype%in%c("telaml","nonHypTelaml")]=NA
        
        clin$birthWt=as.numeric(clin$birthWt)
        bw=read.table(paste(dirBW,"pobw-2014-12-26.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(bw)[match("SubjectID",names(bw))]=c("subjectId")
        j=match(clin$subjectId,bw$subjectId)
        j1=which(!is.na(j)); j2=j[j1]
        clin$pobw=clin$pred_btw=clin$dbirwt=rep(NA,nrow(clin))
        clin$dbirwt[j1]=bw$dbirwt[j2]
        clin$pred_btw[j1]=bw$pred_btw[j2]
        clin$pobw[j1]=bw$pobw[j2]
        
        if (F) {
            load(paste(dirRefactor,"Refactor_set1_and_2_K_6.RData",sep=""))
            Refactor_set=Refactor_set2
            rownames(Refactor_set)=paste("X",rownames(Refactor_set),sep="")
            colnames(Refactor_set)=paste("refactor_",colnames(Refactor_set),sep="")
            j=match(clin$id,rownames(Refactor_set))
            j1=which(!is.na(j)); j2=j[j1]
            j11=which(is.na(j))
            tmp=matrix(nrow=nrow(clin),ncol=ncol(Refactor_set))
            colnames(tmp)=colnames(Refactor_set)
            tmp[j1,]=Refactor_set[j2,]
            clin=cbind(clin,tmp)
            rm(Refactor_set1,Refactor_set2,Refactor_set,tmp)
        }
        
        
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
        
        phen=clin
        phen=phen[match(colnames(meth),phen$id),]
        
        load(paste(dirRefactor,"Refactor_dat",datType,subsetName,".RData",sep=""))
        colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
        j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=cbind(phen[j1,],Refactor_dat[j2,])
    },
    "_allGuthSet1"={
        clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        #names(clin)[match(c("Subject_ID","birth_weight"),names(clin))]=c("subjectId","birthWt")
        #clin$id=paste("X",clin$TargetID,sep="")
        clin$id=paste("X",clin$guthrieId,sep="")
        clin$Beadchip=substr(clin$Bead_Position,1,10)
        if ("Position1"%in%names(clin)) clin$Position=clin$Position1 else clin$Position=clin$Position.1
        clin$caco=clin$Leukemia
        clin$int_ch_ethnicity=clin$ch_ethnicity
        
        clin2=read.table(paste(dirMethLeuk,fNameMethLeuk,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        j=match(clin$subjectId,clin2$Subject_ID)
        j1=which(!is.na(j)); j2=j[j1]
        
        clin$subtype=rep("",nrow(clin))
        clin$subtype[j1][which(clin2$Subtype[j2]=="hyperdiploid")]="hyperdiploid"
        clin$subtype[j1][which(clin2$Subtype[j2]=="t1221")]="telaml"
        clin$subtype[j1][which(clin2$Subtype[j2]%in%c("mll","others","t119"))]="nonHypTelaml"
        clin$hyperdipTelaml=as.integer(clin$subtype=="hyperdiploid")
        clin$hyperdipTelaml[!clin$subtype%in%c("hyperdiploid","telaml")]=NA
        clin$hyperdipNonHypTelaml=as.integer(clin$subtype=="hyperdiploid")
        clin$hyperdipNonHypTelaml[!clin$subtype%in%c("hyperdiploid","nonHypTelaml")]=NA
        clin$telamlNonHypTelaml=as.integer(clin$subtype=="telaml")
        clin$telamlNonHypTelaml[!clin$subtype%in%c("telaml","nonHypTelaml")]=NA
        
        bw=read.table(paste(dirBW,"pobw-2014-12-26.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(bw)[match("SubjectID",names(bw))]=c("subjectId")
        j=match(clin$subjectId,bw$subjectId)
        j1=which(!is.na(j)); j2=j[j1]
        clin$pobw=clin$pred_btw=clin$dbirwt=rep(NA,nrow(clin))
        clin$dbirwt[j1]=bw$dbirwt[j2]
        clin$pred_btw[j1]=bw$pred_btw[j2]
        clin$pobw[j1]=bw$pobw[j2]
        
        if (F) {
            load(paste(dirRefactor,"Refactor_set1_and_2_K_6.RData",sep=""))
            Refactor_set=Refactor_set1
            rownames(Refactor_set)=paste("X",rownames(Refactor_set),sep="")
            colnames(Refactor_set)=paste("refactor_",colnames(Refactor_set),sep="")
            j=match(clin$id,rownames(Refactor_set))
            j1=which(!is.na(j)); j2=j[j1]
            j11=which(is.na(j))
            tmp=matrix(nrow=nrow(clin),ncol=ncol(Refactor_set))
            colnames(tmp)=colnames(Refactor_set)
            tmp[j1,]=Refactor_set[j2,]
            clin=cbind(clin,tmp)
            rm(Refactor_set1,Refactor_set2,Refactor_set,tmp)
        }
        
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
        
        phen=clin
        phen=phen[match(colnames(meth),phen$id),]
        
        load(paste(dirRefactor,"Refactor_dat",datType,subsetName,".RData",sep=""))
        colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
        j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=cbind(phen[j1,],Refactor_dat[j2,])
    },
    "_leuk"={
        clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
		names(clin)[match(c("Plate","Sentrix_ID","Sentrix_Position"),names(clin))]=c("Batch","Beadchip","Position")
		clin$id=sapply(clin$Sample,function(x) {if (is.na(as.integer(substr(x,1,1)))) x else paste("X",x,sep="")},USE.NAMES=F)
		clin$group="Leukemia"
		clin$group2="Leuk"
		
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
		k=match(names(clin),names(tbl1)); k1=which(!is.na(k)); k2=k[k1]
		clin=rbind(clin[,k1],tbl1[!duplicated(tbl1$id),k2])
		
		
        #x=gsub("(",".",gsub(" |-|)",".",clin$id),fixed=T)
		x=gsub("_+","_",gsub("(","_",gsub(" |-|)","_",clin$id),fixed=T))
		j=which(substr(x,nchar(x),nchar(x))=="_")
		if (length(j)!=0) x[j]=substr(x[j],1,nchar(x[j])-1)
		clin$id=x
		
		beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
		probeId=beta$probeId
		beta=as.matrix(beta[,-1])
		rownames(beta)=probeId
		meth=beta
		rm(probeId,beta)

		phen=clin
        #phen$id[!(phen$id%in%colnames(meth))]
		phen=phen[match(colnames(meth),phen$id),]
	},
	"_allGuthSet1Set2"={
		clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
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
		phen=clin
		phen=phen[match(colnames(meth),phen$id),]
        
        load(paste(dirRefactor,"Refactor_dat",datType,subsetFNName,subsetName,".RData",sep=""))
        colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
        j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=cbind(phen[j1,],Refactor_dat[j2,])
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
		
		phen=clin
		phen=phen[match(colnames(meth),phen$id),]
		
		load(paste(dirRefactor,"Refactor_dat_aml.RData",sep=""))
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

if (subsetFlag!="") {
	switch(subsetFlag,
		   "case"={samId=which(phen$caco==1)},
		   "ctrl"={samId=which(phen$caco==0)},
		   "hisp"={
			   if ("ch_hispanic_bc"%in%names(phen)) {
					samId=which(phen$ch_hispanic_bc==1)
			   } else {
					samId=which(phen$int_ch_ethnicity==1)
			   }
		   },
		   "noHispWt"={samId=which(phen$int_ch_ethnicity==2)},
		   "noHisp"={
			   samId=which(phen$subtype==2)
			   samId=NULL
		   },
		   "hyperdip"={samId=which(phen$subtype=="hyperdiploid")},
		   "telaml"={samId=which(phen$subtype=="telaml")},
		   "noHypTelaml"={samId=which(phen$subtype=="nonHypTelaml")}
	)
	meth=meth[,samId]
	phen=phen[samId,]
	rm(samId)
}

if (mediationFlag) {
    samId=which(!is.na(phen$caco) &!is.na(phen$pobw) & !is.na(phen$sex) & !is.na(phen$ch_hispanic_bc))
    meth=meth[,samId]
    phen=phen[samId,]
    rm(samId)
}

if (length(grep("_allGuth",datType)!=0)) {
	if ("subtype"%in%names(phen)) phen$subtype[which(phen$caco==0)]="control"

	phen$caco=as.factor(phen$caco)
	phen$Beadchip=as.factor(phen$Beadchip)

	phen$race=as.factor(phen$int_ch_race)
	phen$ethn=as.factor(phen$int_ch_ethnicity)
	phen$sex=as.factor(phen$sex)
	#phen$DFE_sup=as.factor(phen$DFE_sup>0)

	phen$race3=as.integer(as.character(phen$race))
	phen$race3[which(phen$race3%in%c(2,3))]=5
	phen$race3=as.factor(phen$race3)

    #if (covFlag%in%c("_covSexRace3Gestage","_covRace3")) {
    if (length(grep("Race3",covFlag))==1) {
		phen$race2=phen$race3
    #} else if (covFlag%in%c("_covHisp","_covSexHisp")) {
    } else if (length(grep("Hisp",covFlag))==1) {
		phen$race2=phen$ch_hispanic_bc
	} else {
		phen$race2=phen$ethn
	}
}
if (datType=="_aml") {
	phen$ethn=as.factor(phen$int_ch_ethnicity)
	phen$race2=phen$ethn
}

for (k in which(names(phen)%in%c("caco","sex","Beadchip"))) {
	phen[,k]=as.factor(phen[,k])
}

rm(clin)

## ----------------------------------------------
if (F) {
    if (computerFlag=="") {
        load(file="ann.RData")
    } else {
        if (computerFlag=="cluster") {
            ann=read.delim(paste("data/","HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
            #snpVec=read.table(paste("data/CpGs to exclude_FINAL.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
            snpVec=read.table(paste("data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        } else {
            ann=read.delim(paste("docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
            #snpVec=read.table(paste("docs/ShwetaChoudhry/CpGs to exclude_FINAL.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
            snpVec=read.table(paste("docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        }
        ann[which(ann[,"CHR"]=="X"),"CHR"]="23"
        ann[which(ann[,"CHR"]=="Y"),"CHR"]="24"
        ann[,"CHR"]=as.integer(ann[,"CHR"])
        ann=ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
        for (k in 1:ncol(ann)) if (class(ann[,k])=="factor") ann[,k]=as.character(ann[,k])
        i=match(rownames(meth),ann[,"IlmnID"])
        table(is.na(i))
        ann=ann[i,]

        snpVec=snpVec[,1]
        ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
    }


    ann$geneSym=sapply(toupper(ann$UCSC_RefGene_Name),function(x) {
        strsplit(x,";")[[1]][1]
    },USE.NAMES=F)
    ann$geneSym[is.na(ann$geneSym)]=""

    #keep=!ann$CHR%in%c(23,24)
    keep=!ann$CHR%in%c(23,24) & apply(meth,1,function(x) {any(!is.na(x))})
    keep=ann$snp==0 & !ann$CHR%in%c(23,24) & apply(meth,1,function(x) {any(!is.na(x))})
    keep=ann$snp==0 & ann$CHR%in%1:22 & apply(meth,1,function(x) {any(!is.na(x))})
}

################################################
# Regression
################################################

save.image(paste("tmp",datType,subsetFNName,subsetName,".RData",sep=""))
load(paste("tmp",datType,subsetFNName,subsetName,".RData",sep=""))
ann=data.frame(IlmnID=rownames(meth),id=rownames(meth),stringsAsFactors=T)

prId=read.table(paste(dirEpistructure,"epistructure_reference_sites.txt",sep=""), sep="\t", h=F, quote="", comment.char="",as.is=T,fill=T)
prId=prId[,1]

keep=which(ann$IlmnID%in%prId)

model1=""
model1=paste(model1,"+",paste("prinComp",strsplit(sub("_covPrinComp","",covPCFlag),"")[[1]],collapse="+",sep=""),sep="")
model1=sub("+","~",model1,fixed=T)
model1=as.formula(model1)
Xunadj=model.matrix(model1,data=phen)
samId=match(rownames(Xunadj),phen$id)

lm.unadj=eBayes(lmFit(meth[keep,samId], Xunadj),robust=ifelse(is.null(winsorTail),F,T), winsor.tail.p=winsorTail)
res=residuals(lm.unadj,y=meth[keep,samId])

head(lm.unadj$coef)
head(lm.unadj$p.value)


table(is.na(res))
res_nona=na.omit(res)
res_PCA=PCA(res_nona, graph = F)

epistr1=res_PCA $var$coord[, c(1)]
epistr2=res_PCA $var$coord[, c(2)]
epistr3=res_PCA $var$coord[, c(3)]
epistr4=res_PCA $var$coord[, c(4)]
epistr5=res_PCA $var$coord[, c(5)]
epistr=cbind(epistr1, epistr2, epistr3, epistr4, epistr5)

tbl=cbind(id=rownames(epistr),as.data.frame(epistr))
write.table(tbl,file=paste("epistructure",datType,subsetFNName,subsetName,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
