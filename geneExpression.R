## Some sections have to be run on cluster. Those with save() statements for files tmp_expr.RData, tmp_expr_apobec.RData
## ---------------------------------

computerFlag="cluster"
computerFlag=""

## ---------------------------------

nProbe=10001
nProbe=-1
nProbe=101

## ---------------------------------
subsetName2=""

normFlag=""
normFlag="_bmiq"
normFlag="_funNorm"

transformFlag=""
transformFlag="_mVal"

datType="_allGuthSet1"; subsetName2="_noNonRndChip"
datType="_allGuthSet1Set2"; subsetName2=""
datType="_aml"; subsetName2=""
datType="_allGuthSet1"; subsetName2=""
datType="_allGuthSet2"; subsetName2=""
datType="_leuk"; subsetName2=""

for (datType in c("_leuk","_allGuthSet1")) {
    ## ---------------------------------
    
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
    
    ##############################################
    
    if (computerFlag=="cluster") {
        setwd("/home/royr/project/JoeWiemels")
    } else {
        dirSrc="/Users/royr/UCSF/"
        dirSrc2=dirSrc
        setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
    }
    
    varFlag=""
    
    subsetName=subsetFlag
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
        "noHispWtCtrl"={varName=sub(")"," non-hispanic white controls)",varName)}
        )
    }
    
    heading=paste(c(varFlag,", ",subsetFlag,", ",subsetFFlag,", ",datType,subsetName2,", ",normFlag,", ",transformFlag),collapse="")
    cat("\n\n============================ ",", ===========================\n\n",sep="")
    cat("\n\n============================",varFlag,", ",subsetFlag,", ",subsetFFlag,", ",datType,subsetName2,", ",normFlag,", ",transformFlag,"===========================\n\n")
    
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
            fNameMeth=paste("beta",normFlag,datType,sep="")
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
            fNameMeth=paste("beta",normFlag,datType,sep="")
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
        
        phen=clin
        phen=phen[match(colnames(meth),phen$id),]
        
        load(paste(dirRefactor,"Refactor_dat",datType,".RData",sep=""))
        colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
        j=match(colnames(meth),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
        meth=meth[,j1]
        phen=cbind(phen[j1,],Refactor_dat[j2,])
        
        phenG=phen
        methG=meth
        rm(phen,meth)
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
        #phen$id[!(phen$id%in%colnames(meth))]
        
        phen=phen[match(colnames(meth),phen$id),]
        
        phenL=phen
        methL=meth
        rm(phen,meth)
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
}
#save(phenG,methG,phenL,methL,file="tmp_expr.RData")

## -------------------------------------------
## -------------------------------------------
load("tmp_expr.RData")

nProbe=-1
phenE=read.table(paste(dirMethLeuk,"i.LEU.exp.v2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
expr=read.table(paste(dirMethLeuk,"d.LEU.exp.v2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
phen=as.data.frame(t(sapply(names(expr)[2:118],function(x) strsplit(x,"_")[[1]],USE.NAMES=F)),stringsAsFactors=F)
names(phen)=c("id","subtype")
phenE$id=paste("X",phenE$Sample,sep="")

exprP=read.table(paste(dirMethLeuk,"BCELL.expression.data.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
exprP=exprP[,-1]
phenP=as.data.frame(t(sapply(names(exprP),function(x) strsplit(x,".",fixed=T)[[1]][1:2],USE.NAMES=F)),stringsAsFactors=F)
names(phenP)=c("id","subtype")
names(exprP)=paste(phenP$id,"_",phenP$subtype,sep="")
expr=cbind(expr,exprP)
phen=rbind(phen,phenP)
rm(exprP,phenP)

if (F) {
    datadir="/Users/royr/Downloads/HuGene-1_0-st-v1-na36-hg19-probeset-csv/"
    tmp=read.table(paste(datadir,"HuGene-1_0-st-v1.na36.hg19.probeset.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=100)
    offset=grep("probeset_id",tmp[,1])-1
    annE=read.table(paste(datadir,"HuGene-1_0-st-v1.na36.hg19.probeset.csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=offset,nrow=nProbe)
    names(annE)=substr(names(annE),3,nchar(names(annE))-1)
    annE=annE[,c("probeset_id","transcript_cluster_id","seqname","start","stop","gene_assignment")]
    for (k in 1:ncol(annE)) {
        if (is.character(annE[,k])) annE[,k]=gsub("\"","",annE[,k])
    }
    for (k in which(names(annE)%in%c("probeset_id","transcript_cluster_id","start","stop"))) {
        annE[,k]=as.integer(annE[,k])
    }
    annE$geneSym=sapply(annE$gene_assignment,function(x) {
        if (x=="---") {
            gene=NA
        } else {
            gene=strsplit(x," /// ")[[1]]
            gene=paste(unique(sapply(gene,function(x) {strsplit(x," // ")[[1]][2]},USE.NAMES=F)),collapse=" ")
        }
        gene
    },USE.NAMES=F)
    save(annE,file="annExprP.RData")
}

datadir="/Users/royr/Downloads/HuGene-1_0-st-v1-na36-hg19-transcript-csv/"
if (F) {
    tmp=read.table(paste(datadir,"HuGene-1_0-st-v1.na36.hg19.transcript.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=100)
    offset=grep("transcript_cluster_id",tmp[,1])-1
    offset=24
    tmp=read.table(paste(datadir,"HuGene-1_0-st-v1.na36.hg19.transcript.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,skip=offset-1,nrow=1)
    #annE=read.table(paste(datadir,"HuGene-1_0-st-v1.na36.hg19.transcript.csv",sep=""),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,skip=offset,nrow=nProbe)
    annE=read.table(paste(datadir,"HuGene-1_0-st-v1.na36.hg19.transcript.csv",sep=""),sep="\n",h=F,quote="",comment.char="",as.is=T,fill=T,skip=offset,nrow=nProbe)
    names(annE)=gsub("\"","",unlist(tmp))
}
offset=23
annE=read.csv(paste(datadir,"HuGene-1_0-st-v1.na36.hg19.transcript.csv",sep=""), header = TRUE, sep = ",", quote = "\"",dec = ".", as.is=T,fill=T,comment.char = "",skip=offset,nrow=nProbe)
annE=annE[,c("probeset_id","transcript_cluster_id","seqname","start","stop","gene_assignment")]
for (k in 1:ncol(annE)) {
    if (is.character(annE[,k])) annE[,k]=gsub("\"","",annE[,k])
}
for (k in which(names(annE)%in%c("probeset_id","transcript_cluster_id","start","stop"))) {
    annE[,k]=as.integer(annE[,k])
}
annE$geneSym=sapply(annE$gene_assignment,function(x) {
    if (x=="---") {
        gene=NA
    } else {
        gene=strsplit(x," /// ")[[1]]
        gene=paste(unique(sapply(gene,function(x) {strsplit(x," // ")[[1]][2]},USE.NAMES=F)),collapse=" ")
    }
    gene
},USE.NAMES=F)
save(annE,file="annExprT.RData")
table(expr[,1]%in%annE$transcript_cluster_id)

if (F) {
    load("annExprP.RData")
    annEP=annE
    load("annExprT.RData")
    annET=annE
}

## -------------------------------------------
## -------------------------------------------
load("annExprT.RData")
annE$chr=sub("chr","",annE$seqname)
annE$chr[which(annE$chr=="X")]="23"
annE$chr[which(annE$chr=="Y")]="24"
annE$chr=as.integer(annE$chr)
i=match(expr$probe.ID,annE$transcript_cluster_id)
annE=annE[i,]
#expr=as.matrix(expr[,2:118])
expr=as.matrix(expr[,match(paste(phen$id,"_",phen$subtype,sep=""),colnames(expr))])
candGene=c("APOBEC1","APOBEC2","APOBEC3A","APOBEC3B","APOBEC3C","APOBEC3D","APOBEC3E","APOBEC3F","APOBEC3G","APOBEC3H","APOBEC4","AICDA","AID","ARP2","CDA2","HIGM2","HEL-S-284")
iC=grep(paste(candGene,collapse="|"),annE$geneSym)
sort(annE$geneSym[iC][! annE$geneSym[iC]%in%candGene])
sort(candGene)
iC=iC[which(!annE$geneSym[iC]%in%c("AIDA","ATRAID","FARP2","PARP2"))]
iCA=iC
candGene=c("IKZF1","IKZF2","IKZF3","IKZF4","IKZF5","ZPBP2","ORMDL3","GSDMB","GSDMA")
iC=grep(paste(candGene,collapse="|"),annE$geneSym)
annE[iC,c("seqname","start","stop","geneSym")]
iCI=iC
#phenE

#save(expr,phen,annE,iCA,iCI,file="tmp_expr2.RData")

## -------------------------------------------
## -------------------------------------------
load(file="tmp_expr2.RData")

load("annAll.RData")
#load("ann.RData")
ann=ann[match(rownames(methG),ann$IlmnID),]
if (F) {
    i1=c(); i21=c()
    for (i2 in which(annE$chr%in%unique(annE$chr[iC]))) {
        i=which(ann$CHR==annE$chr[i2])
        i1=c(i1,which.min(abs(ann$MAPINFO[i]-annE$start[i2])))
        i21=c(i21,i2)
    }
}
i1=c(); i21=c()
for (i2 in which(annE$chr%in%unique(annE$chr[iC]))) {
    i=which(ann$CHR==annE$chr[i2])
    i1=c(i1,i[which.min(abs(ann$MAPINFO[i]-annE$start[i2]))])
    i21=c(i21,i2)
}
tmp=rep(NA,length(i21))
tmpC=rep("",length(i21))
map=data.frame(probeId=annE$transcript_cluster_id[i21],cpgId=ann$IlmnID[i1],apobec=as.integer(i21%in%iC),stringsAsFactors=F)
methG=methG[i1,]
methL=methL[i1,]
#save(phenG,methG,phenL,methL,map,file="tmp_expr_apobec.RData")

i1=which(ann[,1]=="cg02051562"); i2=which(annE[,2]==7922804)
i1=which(ann[,1]=="cg06897860"); i2=which(annE[,2]==8119403)
i=1000
i1=which(ann[,1]==map$cpgId[i]); i2=which(annE[,2]==map$probeId[i])
table(ann$CHR[i1]==annE$chr[i2])
(ann$MAPINFO[i1]-annE$start[i2])

## -------------------------------------------
## -------------------------------------------
load("tmp_expr_apobec.RData")

library(qvalue)
pThres=0.05

transformFlag="_mVal"
transformFlag=""

iC=iCA
varName1=c("grp","meth")
parInfo=list(name="_apobec",mfcol=c(2,3),pv=T,outFormat="png")

iC=iCI; iC=iC[order(annE$geneSym[iC])]
varName1=c("grp")
parInfo=list(name="_ikzf",mfcol=c(1,1),pv=F,outFormat="pdf")

for (vId1 in 1:length(varName1)) {
    for (transformFlag in c("_mVal","")) {
        if (varName1[vId1]=="grp" & transformFlag=="_mVal") next
        load("tmp_expr_apobec.RData")
        j=match(phen$id,colnames(methL)); jE=which(!is.na(j)); jM=j[jE]
        
        if (transformFlag=="_mVal") {
            offset=10^-6
            meth=methG
            meth[which(meth==0)]=offset; meth[which(meth==1)]=1-offset
            #cat("Min meth:",min(meth),min(meth,na.rm=T),"\n")
            #cat("Max meth:",max(meth),max(meth,na.rm=T),"\n")
            meth=log2(meth/(1-meth))
            methG=meth
            meth=methL
            meth[which(meth==0)]=offset; meth[which(meth==1)]=1-offset
            #cat("Min meth:",min(meth),min(meth,na.rm=T),"\n")
            #cat("Max meth:",max(meth),max(meth,na.rm=T),"\n")
            meth=log2(meth/(1-meth))
            methL=meth
        }

        colId=2
        phen$grp=phen$subtype
        phen$grp[which(phen$grp=="t1221")]="telaml"
        phen$grp[which(!phen$grp%in%c("hyperdiploid","telaml","others","S1","S2","S3","S4"))]=NA
        phen$grp[which(phen$grp=="hyperdiploid")]="1:hyperdiploid"
        phen$grp[which(phen$grp=="telaml")]="2:telaml"
        phen$grp[which(phen$grp=="others")]="3:others"
        phen$grp[which(phen$grp=="S1")]="4:S1"
        phen$grp[which(phen$grp=="S2")]="5:S2"
        phen$grp[which(phen$grp=="S3")]="6:S3"
        phen$grp[which(phen$grp=="S4")]="7:S4"
        grpUniq=sort(unique(phen$grp))
        grpInfo=matrix("",nrow=length(grpUniq)*(length(grpUniq)-1)/2,ncol=2)
        k=1
        for (k1 in 1:(length(grpUniq)-1)) {
            for (k2 in (k1+1):length(grpUniq)) {
                grpInfo[k,1]=grpUniq[k1]
                grpInfo[k,2]=grpUniq[k2]
                k=k+1
            }
        }
        tmp=rep(NA,length(iC))
        tmpC=rep("",length(iC))
        tbl=matrix(nrow=length(iC),ncol=nrow(grpInfo))
        colnames(tbl)=paste("pv_",grpInfo[,2],"V",grpInfo[,1],sep="")
        tbl=data.frame(probeId=tmp,cpgId=tmp,geneSym=tmpC,model=tmpC,coef=tmp,tbl,stringsAsFactors=F)
        tbl2=matrix(nrow=length(iC),ncol=2*length(grpUniq))
        colnames(tbl2)=paste(rep(c("expr","meth"),each=length(sort(unique(phen$grp)))),"_",sort(unique(phen$grp)),sep="")
        for (subset1Flag in c("")) {
            cat("\n\n=================== ",varName1[vId1],", ",transformFlag," ============\n",sep="")
            if (parInfo$outFormat=="png") {
                png(paste("plot_expr",parInfo$name,ifelse(varName1[vId1]=="grp","",varName1[vId1]),subset1Flag,"_%02d.png",sep=""),width=3*240, height=2*240)
            } else {
                pdf(paste("plot_expr",parInfo$name,ifelse(varName1[vId1]=="grp","",varName1[vId1]),subset1Flag,".pdf",sep=""))
            }
            par(mfcol=parInfo$mfcol)
            xlim=NULL
            for (hId in 1:length(iC)) {
                #for (hId in 1) {
                #print(hId)
                clin=phen
                if (varName1[vId1]=="meth") {
                    clin$expr=expr[iC[hId],]
                    clin$meth=NA
                    k=which(map$probeId==annE$probeset_id[iC[hId]])
                    cpgId=which(rownames(methL)==map$cpgId[k])[1]
                    clin$meth[jE]=methL[cpgId,jM]
                    clin$subtype=clin$grp
                    tbl$model[hId]="expr ~ meth * subtype"
                    varNameThis="subtype"
                    grpUniq=sort(unique(clin[,"subtype"]))
                } else {
                    varNameThis="grp"
                    tbl$model[hId]="expr ~ subtype"
                    grpUniq=sort(unique(clin[,varName1[vId1]]))
                }
                j=1:nrow(clin)
                for (subsetFlag in c("")) {
                    switch(subsetFlag,
                    "case"={
                        j=which(clin$caco==1)
                    },
                    "ctrl"={
                        j=which(clin$caco==0)
                    }
                    )
                    #cat(hId,"\n")
                    #print(table(clin$subtype[j]))
                    tbl$probeId[hId]=annE$transcript_cluster_id[iC[hId]]
                    tbl$geneSym[hId]=strsplit(annE$geneSym[iC[hId]]," ")[[1]][1]
                    for (grpId1 in 1:(length(grpUniq)-1)) {
                        for (grpId2 in (grpId1+1):length(grpUniq)) {
                            jj=j[which(clin[j,varNameThis]%in%grpUniq[c(grpId1,grpId2)])]
                            if (varName1[vId1]%in%c("meth")) {
                                tbl$cpgId[hId]=rownames(methL)[cpgId]
                                res=summary(lm(as.formula(tbl$model[hId]),data=clin[jj,]))$coef
                                colId=nrow(res)
                            } else {
                                res=summary(lm(expr[iC[hId],jj]~clin[jj,varName1[vId1]]))$coef
                                colId=2
                            }
                            tbl$coef[hId]=res[colId,"Estimate"]
                            tbl[hId,paste("pv_",grpUniq[grpId2],"V",grpUniq[grpId1],sep="")]=res[colId,"Pr(>|t|)"]
                        }
                    }
                    if (varName1[vId1]%in%c("meth")) {
                        for (grpId1 in 1:length(grpUniq)) {
                            #cat(grpId1,grpUniq[grpId1],mean(clin$expr[j][which(clin$subtype[j]==grpUniq[grpId1])],na.rm=T),"\n")
                            tbl2[hId,paste("expr_",grpUniq[grpId1],sep="")]=mean(clin$expr[j][which(clin$subtype[j]==grpUniq[grpId1])],na.rm=T)
                            tbl2[hId,paste("meth_",grpUniq[grpId1],sep="")]=mean(clin$meth[j][which(clin$subtype[j]==grpUniq[grpId1])],na.rm=T)
                        }
                    } else {
                        k=hId
                        if (length(k)==1) {
                            if (is.numeric(clin[j,varName1[vId1]])) {
                                res=lm(as.formula(tbl$model[k]),data=clin[j,])
                                plot(clin[j,varName1[vId1]],clin$adult[j],xlim=xlim,main=paste(tbl$set[k],"\nmodel ",tbl$model[k],", pv ",tbl$pv[k],sep=""),xlab=varName1[vId1],ylab="Log2 expression")
                                abline(c(res$coef),lty="solid",col="red")
                            } else {
                                ttl=sort(unique(clin[j,varName1[vId1]])); ttl=substr(ttl,3,nchar(ttl))
                                ttl[which(ttl=="hyperdiploid")]="hyper\ndiploid"
                                if (varName1[vId1]=="hyperdipTelamlVsOther") ttl=c("other","hyperdiploid/tel-aml")
                                heading2=paste(tbl$probeId[k],", ",tbl$geneSym[k],sep="")
                                if (parInfo$name%in%c("apobec")) heading2=paste(heading2,"\npv ",signif(tbl$pv_12[k],2),"(HO) ",signif(tbl$pv_13[k],2),"(HT) ",signif(tbl$pv_23[k],2),"(OT)",sep="")
                                boxplot(expr[iC[hId],j]~clin[j,varName1[vId1]],names=ttl,notch=ifelse(parInfo$name%in%c("apobec"),T,F),ylim=xlim,main=heading2,xlab="",ylab="Log2 expression",las=0)
                                if (!parInfo$name%in%c("apobec")) {
                                    for (p in 1:length(grpUniq)) {
                                        jj=j[which(clin[j,varName1[vId1]]==grpUniq[p])]
                                        points(rep(p,length(jj)),expr[iC[hId],jj],col="red",pch=20)
                                    }
                                }
                            }
                        }
                    }
                }
            }
            dev.off()
        }
        tbl=cbind(tbl,tbl2)
        #names(tbl)=sub("_12","_hyperdipVsOther",names(tbl))
        #names(tbl)=sub("_13","_hyperdipVsTelaml",names(tbl))
        #names(tbl)=sub("_23","_telamlVsOther",names(tbl))
        for (k in grep("pv_",names(tbl))) {
            cat("\n\n========== ",names(tbl)[k]," ============\n",sep="")
            #print(tbl)
            kk=which(tbl[,k]<pThres)
            cat("PV<",pThres,": ",length(kk),"\n",sep="")
            if (length(kk)>0) print(tbl[kk,])
            #qv=qvalue(tbl[,k])$qvalues
            qv=p.adjust(tbl[,k],method="BH")
            cat("QV<",pThres,": ",sum(qv<pThres),"\n",sep="")
        }
        if (varName1[vId1]=="grp") {
            k=c(which(names(tbl)%in%c("probeId","geneSym","model","coef")),grep("pv_",names(tbl)))
            write.table(tbl[,k],file=paste("stat_expr",parInfo$name,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        } else {
            if (transformFlag=="_mVal") {
                tblM=tbl
                write.table(tbl,file=paste("stat_expr_methMval",parInfo$name,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            } else {
                tblB=tbl
                write.table(tbl,file=paste("stat_expr_methBval",parInfo$name,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            }
        }
    }
}

table(clin$subtype)
"
hyperdiploid       others       telaml
       31           25           24
"

#save(tblM,tblB,file="tmp_result_expr_apobec.RData")



## -------------------------------------------
## -------------------------------------------
load("tmp_expr_apobec.RData")

i21=which(map$probeId%in%c(8073062,8073056))
xlim=range(expr[which(annE$transcript_cluster_id%in%map$probeId[i21]),],na.rm=T)
ylim=c(0,1)
png(paste("plot_expr_",tolower("APOBEC3A_B"),".png",sep=""),width=3*240, height=2*240)
par(mfcol=c(2,1))
for (i2 in i21) {
    i=which(annE$transcript_cluster_id==map$probeId[i2])
    plot(density(expr[i,],na.rm=T),xlim=xlim,ylim=ylim,main=paste(annE$transcript_cluster_id[i],": ",gsub(" ",", ",annE$geneSym[i]),sep=""),xlab="log2(expression)")
}
dev.off()
grpUniq=sort(unique(phen$grp))
png(paste("plot_expr_",tolower("APOBEC3A_B"),"_%02d.png",sep=""),width=3*240, height=2*240)
par(mfrow=c(2,3))
for (i2 in i21) {
    i=which(annE$transcript_cluster_id==map$probeId[i2])
    for (grpId in 1:length(grpUniq)) {
        j=which(phen$grp==grpUniq[grpId])
        plot(density(expr[i,j],na.rm=T),xlim=xlim,ylim=ylim,main=paste(annE$transcript_cluster_id[i],": ",gsub(" ",", ",annE$geneSym[i]),sep=""),sub=grpUniq[grpId],xlab="log2(expression)")
    }
}
dev.off()

## -------------------------------------------
## -------------------------------------------
