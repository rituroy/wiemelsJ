## Linear or logistic regression
## Additive or multiplicative model

## ---------------------------------

computerFlag=""
computerFlag="cluster"

## ---------------------------------

nProbe=10001
nProbe=101
nProbe=-1

## ---------------------------------
subsetName2=""

normFlag=""
normFlag="_funNorm"
normFlag="_bmiq"

transformFlag="_mVal"
transformFlag=""

datType="_leuk"; subsetName2=""
datType="_allGuthSet1"; subsetName2="_noNonRndChip"
datType="_allGuthSet1Set2"; subsetName2=""
datType="_aml"; subsetName2=""
datType="_allGuthSet1"; subsetName2=""
datType="_allGuthSet2"; subsetName2=""

mediationFlag=T
mediationFlag=F

computeFlag="covariance"
computeFlag="regression"


## ---------------------------------
winsorTail=NULL
winsorTail=0.3
winsorTail=0.2
winsorTail=0.1
winsorTail=0.05

subsetFFlag="_noEpiStr"
subsetFFlag=""

setFlag=ifelse(subsetName2=="",tolower(sub("allGuth","",datType)),subsetName2)

subsetFlag="hispCtrl"
subsetFlag="noHispWtCtrl"

subsetFlag="case"
subsetFlag="ctrl"

subsetFlag="female"
subsetFlag="male"

#subsetFlag="noHisp"
subsetFlag="hisp"
subsetFlag="noHispWt"

subsetFlag=""

#for (subsetName2 in c("_set1","_set2")) {
#for (datType in c("_allGuthSet1","_allGuthSet2")) {
		
## ---------------------------------

## All smoke variables from set1 & set2
## All smoke variables in set2 are in set1
varSmoke=data.frame(varIn=c("passive_sm_home_post","smoke_mo_ever","smoke_mo_preg","M_SM_BF_DI","smoke_mo_3months","smoke_mo_3months_N","smoke_mo_preg_N","m_cignum_bf","smoke_mo_after","smoke_mo_after_N","smoke_fa_ever","smoke_fa_3months","smoke_fa_3months_N","smoke_mo_bf","smoke_mo_bf_N","smoke3","smoke2"),
varOut=c("paSmHmPost","moEver","moPreg","moBf","mo3m","mo3mN","moPregN","moBfN","moAfter","moAfterN","faEver","fa3m","fa3mN","moBf","moBfN","smoke3","smoke2"),stringsAsFactors=F)

## Smoke variables in common between set1 & set2
varSmoke=data.frame(varIn=c("smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after","smoke_mo_after_N","smoke_fa_ever","smoke_fa_3months","smoke_fa_3months_N","smoke_mo_bf","smoke_mo_bf_N","smoke3","smoke2"),
varOut=c("moEver","moPreg","mo3m","mo3mN","moPregN","moAfter","moAfterN","faEver","fa3m","fa3mN","moBf","moBfN","smoke3","smoke2"),stringsAsFactors=F)

## Ctrl subset: Covariates associated with a smoking variable in set2 are removed for both set1 & set2 analysis. None of those variables are associated with set1 smoking variables
## Ethnicity kept for smoke_fa_3months_N (p=0.070)
## Sex excluded for smoke_mo_3months_N (p=0.0630)
k=nrow(varSmoke)
setList=c("_set1","_set2")
covInfo=data.frame(set=rep(c("_set1","_set2"),each=k),subset=rep("ctrl",k*2),variable=rep(c("smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_fa_ever","smoke_fa_3months","smoke_mo_bf","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after_N","smoke_fa_3months_N","smoke_mo_bf_N","smoke2","smoke3"),2),covariates=rep("_covSexEthnGestage",k*2),stringsAsFactors=F)
covInfo$covariates[which(covInfo$set%in%setList & covInfo$subset=="ctrl" & covInfo$variable%in%c("smoke_fa_3months"))]="_covSexRace3Gestage"
covInfo$covariates[which(covInfo$set%in%setList & covInfo$subset=="ctrl" & covInfo$variable%in%c("smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_after","smoke_mo_preg_N","smoke_mo_after_N","smoke3"))]="_covSexEthn"
covInfo$covariates[which(covInfo$set%in%setList & covInfo$subset=="ctrl" & covInfo$variable%in%c("smoke_mo_3months_N"))]="_covEthn"

varFlag="_dfeNat"; covFlag="_covSexEthnGestage"; varName=""; termName=""
varFlag="_dfeNat"; covFlag="_covSexRace3Gestage"; varName=""; termName=""
varFlag="_dfeNatCat"; covFlag=""; varName=""; termName=""
varFlag="_dfeNatCat"; covFlag=""; varName=""; termName=""

varFlag="_smoke_moPreg"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_mo3m"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_mo3mN"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_moPregN"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_moAfter"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_moAfterN"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_faEver"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_fa3m"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_fa3mN"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_moBf"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_moBfN"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke3"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke2"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""
varFlag="_smoke_moEver"; covFlag=covInfo$covariates[which(covInfo$set==setFlag & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_smoke_",varSmoke$varOut,sep="")==varFlag)])]; varName=""; termName=""


if (datType%in%c("_allGuthSet1Set2","_allGuthSet1Set2Combat")) {
    covFlag=covInfo$covariates[which(covInfo$set=="_set2" & covInfo$subset==subsetFlag & covInfo$variable==varSmoke$varIn[which(paste("_",varSmoke$varOut,sep="")==varFlag)])]
}

varFlag="_smoke"; covFlag=""; varName=""; termName=""

if (subsetFlag=="ctrl") {varFlag="_birthWt"; covFlag="_covHisp"; varName=""; termName=""}
if (subsetFlag=="case") {varFlag="_birthWt"; covFlag="_covSexHisp"; varName=""; termName=""}
varFlag="_cacoXbirthwt"; covFlag=""; varName=""; termName=c("caco","birthWt","interaction")
varFlag="_pobw"; covFlag=""; varName=""; termName=""
varFlag="_cacoXpobw"; covFlag=""; varName=""; termName=c("caco","pobw","interaction")

varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="meth~pobw"; computeFlag[2]="linear"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="meth~pobwBi"; computeFlag[2]="linear"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*pobw"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*pobwBi"; computeFlag[2]="logistic"

varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*birthWt"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*birthWtBi"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="meth~birthWt"; computeFlag[2]="linear"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="meth~birthWtBi"; computeFlag[2]="linear"


varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*moBmi"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*moBmiOverOrObese"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*moBmiOver"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*moBmiObese"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag=""; varName=""; termName=""; modelFlag="caco~meth*moBmiUnder"; computeFlag[2]="logistic"

varFlag="_hyperdipTelaml"; covFlag="_covEthnGestage"; varName=""; termName=""
varFlag="_hyperdipNonHypTelaml"; covFlag="_covSexRace3Gestage"; varName=""; termName=""
varFlag="_telamlNonHypTelaml"; covFlag="_covRace3"; varName=""; termName=""

varFlag="_leukPreb"; covFlag=""; varName=""; termName=""
varFlag="_caco"; covFlag="_covSexRace3Gestage"; varName=""; termName=""
varFlag="_caco"; covFlag=""; varName=""; termName=""
varFlag="_caco"; covFlag="_covSexEthnGestage"; varName=""; termName=""

varFlag="_caco"; covFlag="_covSexGestageSnpPC1"; varName=""; termName=""
varFlag="_caco"; covFlag="_covSexGestageSnpPC12"; varName=""; termName=""
varFlag="_caco"; covFlag="_covSexGestageSnpPC123"; varName=""; termName=""
varFlag="_caco"; covFlag="_covSexGestageSnpPC1234"; varName=""; termName=""

varFlag="_caco"; covFlag="_covSex"; varName=""; termName=""
varFlag="_caco"; covFlag="_covEthn"; varName=""; termName=""
varFlag="_caco"; covFlag="_covGestage"; varName=""; termName=""
varFlag="_caco"; covFlag="_covSexEthn"; varName=""; termName=""
varFlag="_caco"; covFlag="_covEthnGestage"; varName=""; termName=""
varFlag="_caco"; covFlag="_covSexEthnGestage"; varName=""; termName=""; modelFlag="meth~caco"; computeFlag[2]="linear"
varFlag="_caco"; covFlag="_covSexEthnGestage"; varName=""; termName=""; modelFlag="caco~meth"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="caco~meth"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag="_covGestage"; varName=""; termName=""; modelFlag="meth~caco"; computeFlag[2]="linear"
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="meth~caco"; computeFlag[2]="linear"


if (mediationFlag) {
    computeFlag=c("regression","")
    modelFlag="caco~pobw"; computeFlag[2]="logistic" # Model 1
    modelFlag="meth~pobw"; computeFlag[2]="linear" # Model 2
    modelFlag="caco~meth+pobw"; computeFlag[2]="logistic" # Model 3
    
    computeFlag="covariance"
    modelFlag="caco~meth+pobw" # Compute variance/covarance
    
    computeFlag=c("regression","")
    modelFlag="pobw~meth"; computeFlag[2]="linear" # Model 2
    modelFlag="caco~meth+pobw"; computeFlag[2]="logistic" # Model 3
    modelFlag="caco~meth"; computeFlag[2]="logistic" # Model 1
    
    varFlag="_caco"; covFlag="_covSexHisp"; varName=""
}

covPCFlag="_covPrinComp123456"
covPCFlag="_covPrinComp12"
covPCFlag=""
covPCFlag="_covPrinComp1234"

covESFlag="_covEpStr" # do not adjust for cell mixture since refactor components are regressed out when estimating epistructures
covESFlag=""

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
if (subsetFlag!="") {
	#varFlag=paste(varFlag,"_",subsetFlag,"Subset",sep="")
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
           "noHispWtCtrl"={varName=sub(")"," non-hispanic white controls)",varName)}
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
		   fNameClin="clin_guthrieSet2_20140619"
           fNameClin="clin_allGuthSet2_20160928"
        },
		"_allGuthSet1"={
		   dirMeth=dirClin="data/set1/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),sep="")
		   fNameClin="final"
           fNameClin="clin_allGuthSet1_20160928"
		   dirMethLeuk=dirMeth
		   fNameMethLeuk="i.LEU.v2"
        },
        "_allGuthSet1Set2"={
            dirMeth=dirClin="data/set1set2/"
            fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1set2",datType),sep="")
            fNameClin="clin_guthrieSet1Set2_20151022"
        },
		"_leuk"={
		   dirMeth=dirClin=dirClin2="data/set1/"
		   fNameMeth=paste("beta",normFlag,datFlag,sep="")
		   fNameClin="i.LEU.v2"
		   fNameClin2="0708011 Sample_Sheet (Fetal blood)"
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
		"_ivorra"={
		   dirMeth=dirClin="data/ivorra2014/"
		   fNameMeth="beta_ivorra"
		   fNameClin="clin_ivorra"
		}
	)
} else {
	dirBW="docs/birthWeight/"
	dirCom="docs/all/"
    #dirRefactor="docs/SemiraGonsethNussle/refactor/"
    switch(datType,
		"_allGuthSet2"={
		   dirMeth=dirClin="docs/all/set2/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
		   fNameClin="clin_guthrieSet2_20140619"
           fNameClin="clin_allGuthSet2_20160928"
           dirRefactor=dirClin
	   },
	   "_allGuthSet1"={
			dirMeth=dirClin="docs/all/set1/"
			fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),sep="")
			fNameClin="final"
            fNameClin="clin_allGuthSet1_20160928"
			dirMethLeuk="docs/all/set1/LEU.data/"
			fNameMethLeuk="i.LEU.v2"
            dirRefactor=dirClin
       },
       "_allGuthSet1Set2"={
           dirMeth=dirClin="docs/all/set1set2/"
           fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
           fNameClin="clin_guthrieSet1Set2_20151022"
           dirRefactor=dirClin
       },
		"_leuk"={
		   dirMeth="docs/all/set1/"
		   fNameMeth=paste("beta",normFlag,datFlag,sep="")
		   dirClin="docs/all/set1/LEU.data/"
		   fNameClin="i.LEU.v2"
		   dirClin2="docs/all/set1/preBcell/"
		   fNameClin2="0708011 Sample_Sheet (Fetal blood)"
           dirRefactor=dirClin
        },
		"_aml"={
		   dirClin=dirMeth="docs/aml/"
		   fNameMeth=paste("beta",normFlag,"_aml",sep="")
		   fNameClin=paste("clin_aml_20150114",sep="")
           dirRefactor="docs/aml/refactor/"
        },
	   "_ivorra"={
		   dirMeth=dirClin="docs/misc/ivorra2014/"
		   fNameMeth="beta_ivorra"
		   fNameClin="clin_ivorra"
           dirRefactor=dirClin
       }
	)
}

if (datType%in%c("_allGuthSet1Set2") & length(grep("ctrl",tolower(subsetFlag)))==1) {
    samInfo=read.table(paste(dirClin,"summaryBeta_forSample_allGuthSet1Set2_ctrlSubset.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
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
		
		
#		x=gsub("(",".",gsub(" |-|)",".",clin$id),fixed=T)
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
#		phen$id[!(phen$id%in%colnames(meth))]

		phen=phen[match(colnames(meth),phen$id),]
	},
	"_allGuthSet1Set2"={
		clin=read.table(paste(dirClin,fNameClin,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        
        bw=read.table(paste(dirBW,"pobw-2014-12-26.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(bw)[match("SubjectID",names(bw))]=c("subjectId")
        j=match(clin$subjectId,bw$subjectId)
        j1=which(!is.na(j)); j2=j[j1]
        clin$pobw=clin$pred_btw=clin$dbirwt=rep(NA,nrow(clin))
        clin$dbirwt[j1]=bw$dbirwt[j2]
        clin$pred_btw[j1]=bw$pred_btw[j2]
        clin$pobw[j1]=bw$pobw[j2]
        
        if (length(grep("ctrl",tolower(subsetFlag)))==1) {
            beta=read.table(paste(dirMeth,fNameMeth,"_ctrlSubset.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        } else {
            beta=read.table(paste(dirMeth,fNameMeth,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        }
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
        
        if (length(grep("ctrl",tolower(subsetFlag)))==1) {
            load(paste(dirRefactor,"Refactor_dat",datType,"_ctrlSubset.RData",sep=""))
        } else {
            load(paste(dirRefactor,"Refactor_dat",datType,".RData",sep=""))
        }
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

####################################################################
####################################################################


##########################################
############## e.g. of choice of indexes for exclusion
# toremove<- which(phenotypes$sex!=phenotypes$Gender)



#################### Regression Models


library(data.table)# to process results
library(MASS) # rlm function for robust linear regression
library(sandwich) #Huberís estimation of the standard error
library(lmtest) # to use coeftest
#library(parallel) # to use multicore approach - part of base R


removeOutliers <- function(probes) {
    require(matrixStats)
    if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
    rowIQR <- rowIQRs(probes, na.rm = T)
    row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
    maskL <- probes < row2575[,1] - 3 * rowIQR
    maskU <- probes > row2575[,2] + 3 * rowIQR
    initial_NAs <- rowSums(is.na(probes))
    probes[maskL] <- NA
    removed_lower <- rowSums(is.na(probes))-initial_NAs
    probes[maskU] <- NA
    removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
    N_for_probe <- rowSums(!is.na(probes))
    Log <- data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
    return(list(probes, Log))
}


require(matrixStats)
getSEM <- function(probes) {
    rowIQR <- rowIQRs(probes, na.rm = T)
    row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
    out=apply(probes,2,function(x,y) {
        z1=x-y[,1]
        z2=x-y[,2]
        z=z1
        z[which(z1>=0)]=NA
        i=which(z2>0)
        z[i]=z2[i]
        z
    },y=cbind(row2575[,1] - 3 * rowIQR,row2575[,2] + 3 * rowIQR))
    return(out)
}
system.time(sem <- getSEM(meth))

save(sem,file=paste("sem",datType,".RData",sep=""))


####################################################################
####################################################################
