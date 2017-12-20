## Linear or logistic regression
## Additive or multiplicative model

argv=commandArgs(TRUE)
varThis=argv[1]

if (F) {
    '
    cd /Users/royr/UCSF/JoeWiemels/leukMeth
    
    qRscript regression_refactor_1.R "logged_PCB_105_SRS"
    qRscript regression_refactor_1.R "logged_PCB_118_SRS"
    qRscript regression_refactor_1.R "logged_PCB_138_SRS"
    qRscript regression_refactor_1.R "logged_PCB_153_SRS"
    qRscript regression_refactor_1.R "logged_PCB_170_SRS"
    qRscript regression_refactor_1.R "logged_PCB_180_SRS"
    qRscript regression_refactor_1.R "logged_PCB_aroclor1260"
    '
    
    varThis="logged_PCB_105_SRS"
    varThis="logged_PCB_118_SRS"
    varThis="logged_PCB_138_SRS"
    varThis="logged_PCB_153_SRS"
    varThis="logged_PCB_170_SRS"
    varThis="logged_PCB_180_SRS"
    varThis="logged_PCB_aroclor1260"
}

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

transformFlag=""
transformFlag="_mVal"

datType="_leuk"; subsetName2=""
datType="_allGuthSet1"; subsetName2="_noNonRndChip"
datType="_aml"; subsetName2=""
datType="_allGuthSet1Set2"; subsetName2=""
datType="_allGuthSet2"; subsetName2=""
datType="_allGuthSet1"; subsetName2=""

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

subsetFlag="female"
subsetFlag="male"

subsetFlag=""

#subsetFlag="noHisp"
subsetFlag="hisp"
subsetFlag="noHispWt"

subsetFlag="case"
subsetFlag="ctrl"

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

varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="caco~meth*sem"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="caco~meth*logSem"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="caco~meth*semNoSexChr"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="caco~meth*logSemNoSexChr"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="caco~meth*semNoSnpNoSexChr"; computeFlag[2]="logistic"
varFlag="_caco"; covFlag="_covSexGestage"; varName=""; termName=""; modelFlag="caco~meth*logSemNoSnpNoSexChr"; computeFlag[2]="logistic"

varFlag="_caco"; covFlag="_covSet"; varName=""; termName=""; modelFlag=paste("meth~",varThis,sep=""); computeFlag[2]="linear"


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

covESFlag=""
covESFlag="_covEpStr" # do not adjust for cell mixture since refactor components are regressed out when estimating epistructures

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
        
        clin2=read.table(paste(dirCom,"chemicals.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(clin2)[match(c("subjectID"),names(clin2))]="subjectId"
        id=clin$subjectId[!clin$subjectId%in%clin2$subjectId]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$subjectId=id
        clin2=rbind(clin2,tmp)
        clin=cbind(clin,clin2[match(clin$subjectId,clin2$subjectId),grep("_PCB",names(clin2))])
        rm(clin2,tmp,id)
        
        load(paste(dirClin,"sem",datType,".RData",sep=""))

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

        clin2=read.table(paste(dirCom,"chemicals.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(clin2)[match(c("subjectID"),names(clin2))]="subjectId"
        id=clin$subjectId[!clin$subjectId%in%clin2$subjectId]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$subjectId=id
        clin2=rbind(clin2,tmp)
        clin=cbind(clin,clin2[match(clin$subjectId,clin2$subjectId),grep("_PCB",names(clin2))])
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
        
        bw=read.table(paste(dirBW,"pobw-2014-12-26.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(bw)[match("SubjectID",names(bw))]=c("subjectId")
        j=match(clin$subjectId,bw$subjectId)
        j1=which(!is.na(j)); j2=j[j1]
        clin$pobw=clin$pred_btw=clin$dbirwt=rep(NA,nrow(clin))
        clin$dbirwt[j1]=bw$dbirwt[j2]
        clin$pred_btw[j1]=bw$pred_btw[j2]
        clin$pobw[j1]=bw$pobw[j2]
        
        clin2=read.table(paste(dirClin,"epistructure",datType,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        id=clin$id[!clin$id%in%clin2$id]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$id=id
        clin2=rbind(clin2,tmp)
        clin=cbind(clin,clin2[match(clin$id,clin2$id),grep("epistr",names(clin2))])
        rm(clin2,tmp,id)
        
        clin2=read.table(paste(dirCom,"chemicals.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(clin2)[match(c("subjectID"),names(clin2))]="subjectId"
        id=clin$subjectId[!clin$subjectId%in%clin2$subjectId]
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$subjectId=id
        clin2=rbind(clin2,tmp)
        clin=cbind(clin,clin2[match(clin$subjectId,clin2$subjectId),grep("_PCB",names(clin2))])
        rm(clin2,tmp,id)
        
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
if (length(grep("PCB",names(phen)))!=0) {
    pcb=read.table(paste(dirCom,"pcb.txt",sep=""), sep=" ", h=F, quote="", comment.char="",as.is=T,fill=T,skip=5)
    names(pcb)=c("pcbNo","clPos","aroclor1016","aroclor1242","aroclor1248A3.5","aroclor1248G3.5","aroclor1254Late","aroclor1254","aroclor1260")
    for (k in 1:ncol(pcb)) if (is.character(pcb[,k])) pcb[,k]=gsub("\"","",pcb[,k])
    #pcb=pcb[!is.na(as.integer(pcb$pcbNo)) & !is.na(as.integer(pcb$aroclor1016)),]
    for (k in c("pcbNo")) pcb[,k]=as.integer(pcb[,k])
    for (k in c("aroclor1016","aroclor1242","aroclor1248A3.5","aroclor1248G3.5","aroclor1254Late","aroclor1254","aroclor1260")) pcb[,k]=as.numeric(pcb[,k])
    k=grep("logged_PCB",names(phen))
    dat=as.matrix(phen[,k])
    x=as.integer(gsub("logged_PCB_|_SRS","",colnames(dat)))
    x=x[!is.na(x)]
    phen$logged_PCB_aroclor1260=dat%*%pcb$aroclor1260[match(x,pcb$pcbNo)]
    #j=2
    #sum(dat[j,]*pcb$aroclor1260[match(x,pcb$pcbNo)])
    #phen$logged_PCB_aroclor1260[j]
}

if (subsetFlag!="") {
	switch(subsetFlag,
		   "case"={samId=which(phen$caco==1)},
		   "ctrl"={samId=which(phen$caco==0)},
           "female"={samId=which(phen$sex==2)},
           "male"={samId=which(phen$sex==1)},
           "hisp"={
               if (F) {
                   if ("ch_hispanic_bc"%in%names(phen)) {
                        samId=which(phen$ch_hispanic_bc==1)
                   } else {
                        samId=which(phen$int_ch_ethnicity==1)
                   }
               }
               samId=which(phen$int_ch_ethnicity==1)
		   },
		   "noHispWt"={samId=which(phen$int_ch_ethnicity==2)},
		   "noHisp"={
			   samId=which(phen$subtype==2)
			   samId=NULL
		   },
		   "hyperdip"={samId=which(phen$subtype=="hyperdiploid")},
		   "telaml"={samId=which(phen$subtype=="telaml")},
		   "noHypTelaml"={samId=which(phen$subtype=="nonHypTelaml")},
           "hispCtrl"={samId=which(phen$int_ch_ethnicity==1 & phen$caco==0)},
           "noHispWtCtrl"={samId=which(phen$int_ch_ethnicity==2 & phen$caco==0)}
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

if (mediationFlag) {
    samId=which(!is.na(phen$caco) &!is.na(phen$pobw) & !is.na(phen$sex) & !is.na(phen$ch_hispanic_bc))
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
#colnames(meth)=phen$Subject_ID
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
	"_chipPos"={phen$caco=phen$Position}
)
if (varFlag%in%paste("_smoke_",varSmoke$varOut,sep="")) {
	k=match(varSmoke$varIn,names(phen)); k2=which(!is.na(k)); k1=k[k2]
	names(phen)[k1]=varSmoke$varOut[k2]
	phen$caco=phen[,sub("_smoke_","",varFlag)]
}

#phen=phen[,c("caco","Leukemia","sex","Beadchip")]
if (!varFlag%in%c("_dfeFood","_dfeFort","_dfeNat","_dfeTot",paste("_smoke_",c("mo3mN","moPregN","moAfterN","fa3mN","moBfN"),sep=""),"_smoke","_birthWt","_pobw","_cacoVbirthWt","_cacoXpobw")) {
	phen$caco=as.factor(phen$caco)
}

if (F) {
	phen$sex=as.factor(phen$sex)
	phen$Beadchip=as.factor(phen$Beadchip)
}

## ----------------------------------------------
if (computerFlag=="") {
    #load(file="ann.RData")
    load(file="annAll.RData")
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

	snpVec=snpVec[,1]
	ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
}

i=match(rownames(meth),ann[,"IlmnID"])
table(is.na(i))
ann=ann[i,]

#keep=ann$CHR%in%1:22
keep=ann$CHR%in%1:22 & apply(meth,1,function(x) {any(!is.na(x))})
keep=ann$snp==0 & ann$CHR%in%1:22 & apply(meth,1,function(x) {any(!is.na(x))})
keep=apply(meth,1,function(x) {any(!is.na(x))})

## ----------------------------------------------

if (datType=="_allGuthSet2") {
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

################################################
# Regression
################################################

colId=2
if (varFlag%in%c("_cacoXbirthwt","_cacoXpobw")) {
    ## CHECK !!!
    model1="~caco*predVar"
    colId=match(c("caco1","predVar","caco1:predVar"),colnames(Xunadj))
} else {
    if (mediationFlag) {
        #model1=paste("~",strsplit(modelFlag,"~")[[1]][2],sep="")
        model1=paste("~",sub("~","+",modelFlag),sep="")
    } else {
        ## CHECK !!!
        model1="~caco"
        model1=paste("~",sub("~","+",modelFlag),sep="")
    }
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
samId=match(rownames(Xunadj),phen$id)
meth=meth[keep,samId]
phen=phen[samId,]

if (transformFlag=="_mVal") {
    meth1=meth
    methInfo=read.table(paste(dirCom,"summaryBeta.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
    k=which(methInfo$cohort==datType)
    meth[which(meth==0)]=methInfo$minBeta*0.5; meth[which(meth==1)]=methInfo$maxBeta+(1-methInfo$maxBeta)*0.5
    cat("Min meth:",min(meth),min(meth,na.rm=T),"\n")
    cat("Max meth:",max(meth),max(meth,na.rm=T),"\n")
    meth=log2(meth/(1-meth))
}


timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
#if (mediationFlag) {
    if (computeFlag[1]=="covariance") {
        varId=gsub("+","",gsub("~|caco|meth","",modelFlag),fixed=T)
        nm=c(varId,"meth")
        out=matrix(nrow=nrow(meth),ncol=length(nm)+1,dimnames=list(rownames(meth),c(paste("var_",nm,sep=""),"cov")))
        for (i in 1:nrow(meth)) {
            j=!is.na(meth[i,])
            out[i,1]=var(phen[j,varId])
            out[i,2]=var(meth[i,j])
            out[i,3]=cov(phen[j,varId],meth[i,j])
        }
        tbl=data.frame(cpgId=ann$IlmnID[keep],out,stringsAsFactors=F)
        write.table(tbl,file=paste("stat_mediation_covariance_",sub("~","Resp_",gsub("+","_",modelFlag,fixed=T)),subsetName,covFlag,covPCFlag,covESFlag,subsetFFlag,datType,subsetName2,normFlag,transfomrFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    } else {
        #nm=c("intercept",strsplit(strsplit(modelFlag,"~")[[1]][2],"+",fixed=T)[[1]])
        nm=c("intercept",strsplit(strsplit(sub("*","+",modelFlag,fixed=T),"~")[[1]][2],"+",fixed=T)[[1]])
        if (length(grep("*",modelFlag,fixed=T))) nm=c(nm,sub("*",":",strsplit(strsplit(modelFlag,"~")[[1]][2],"+",fixed=T)[[1]],fixed=T))
        tmp=matrix(nrow=nrow(meth),ncol=length(nm),dimnames=list(rownames(meth),nm))
        out=list(coef=tmp,se=tmp,pValue=tmp)
        x=strsplit(modelFlag,"~")[[1]][1]
        model1=paste(x,sub(paste(x,"+",sep=""),"",model1,fixed=T),sep="")
        model1=sub("meth","meth[i,j]",model1)
        colId2=nm
        colId=sub("caco","caco1",sub("intercept","(Intercept)",sub("meth","meth[i, j]",nm)))
        model1=as.formula(model1)
        if (computeFlag[2]=="linear") {
            # Run linear regresssion model
            for (i in 1:nrow(meth)) {
                j=!is.na(meth[i,])
                fit=tryCatch(lm(model1, data=phen[j,]),error = function(e) e)
                if (!inherits(fit,"error")) {
                    res=summary(fit)$coef
                    out$coef[i,]=res[colId,1]
                    out$se[i,]=res[colId,2]
                    out$pValue[i,]=res[colId,4]
                }
            }
        } else {
            # Run logistic regresssion model
            for (i in 1:nrow(meth)) {
                j=!is.na(meth[i,])
                fit=tryCatch(glm(model1, family="binomial",data=phen[j,]),error = function(e) e)
                if (!inherits(fit,"error")) {
                    res=summary(fit)$coef
                    out$coef[i,]=res[colId,1]
                    out$se[i,]=res[colId,2]
                    out$pValue[i,]=res[colId,4]
                }
            }
        }
        for (kk in 2:ncol(out$coef)) {
            k=colId2[kk]
            tbl2=cbind(out$coef[,k],out$se[,k],out$pValue[,k])
            if (kk==2) {
                tbl=tbl2
            } else {
                tbl=cbind(tbl,tbl2)
            }
        }
        colnames(tbl)=paste(c("coef","se","pv"),"_",rep(colnames(out$coef)[2:ncol(out$coef)],each=3),sep="")
        tbl=data.frame(cpgId=ann$IlmnID[keep],tbl,stringsAsFactors=F)
        write.table(tbl,file=paste("stat_",ifelse(mediationFlag,"mediation_",""),sub("~","Resp_",gsub("*","X",gsub("+","_",modelFlag,fixed=T),fixed=T)),subsetName,covFlag,covPCFlag,covESFlag,subsetFFlag,datType,subsetName2,normFlag,transformFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    }
#}
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))
