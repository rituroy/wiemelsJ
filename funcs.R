## Create data objects

########################################################################
########################################################################
########################################################################

getClinData=function(setFlag,subsetFlag) {

    normFlag=""
    normFlag="_funNorm"
    normFlag="_bmiq"

    ## ---------------------------------

    if (computerFlag=="cluster") {
        datType2="allGuthSet2"
        datName2="ALL Guthrie, set2"
        dirClin2="data/set2/"
        dirMeth2=dirClin2
        #fileList2=c("mDat_funNorm_set2","beta_funNorm_set2","beta_bmiq_allGuthSet2","wbcAdjustedMDat_allGuthSet2_leukChip","wbcAdjustedData_allGuthSet2_leukChip")
        #fName22="clin_guthrieSet2_20140619"
        fileList2=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2","_allGuthSet2"),sep="")
        fName22="clin_allGuthSet2_20160928"
        
        datType1="allGuth"
        datName1="ALL Guthrie, set1"
        dirClin1="data/set1/"
        dirMeth1=dirClin1
        #fileList1=c("mDat_funNorm_set1","beta_funNorm_set1","beta_bmiq_allGuthSet1","wbcAdjustedMDat_allGuthSet1_leukChip","wbcAdjustedData_allGuthSet1_leukChip")
        #fName12="final"
        fileList1=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1","_allGuthSet1"),sep="")
        fName12="clin_allGuthSet1_20160928"
        
        datTypeL="allLeuk"
        datNameL="ALL Leukemia, set1"
        datadirL1="data/set1/"
        datadirL2=datadirL1
        fileListL=paste("beta",normFlag,"_leuk",sep="")
        fNameL2="i.LEU.v2"
        datadirL3=datadirL1
        fNameL3="0708011 Sample_Sheet (Fetal blood)"
        
        nRow=-1
    } else {
        dirSrc="/Users/royr/UCSF/"
        dirSrc2=dirSrc
        setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
        
        datType2="allGuthSet2"
        datName2="ALL Guthrie, set2"
        dirClin2="docs/all/set2/"
        dirMeth2=dirClin2
        #fileList=c("beta_funNorm_set2")
        #fName22="clin_guthrieSet2_20140619"
        fileList2=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2","_allGuthSet2"),sep="")
        fName22="clin_allGuthSet2_20160928"
        
        datType1="allGuth"
        datName1="ALL Guthrie, set1"
        dirClin1="docs/all/set1/"
        dirMeth1=dirClin1
        #fileList=c("beta_funNorm_set1")
        #fName12="final"
        fileList1=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1","_allGuthSet1"),sep="")
        fName12="clin_allGuthSet1_20160928"
        
        datTypeL="allLeuk"
        datNameL="ALL Leukemia, set1"
        datadirL1="docs/all/set1/"
        datadirL2="docs/all/set1/LEU.data/"
        fileListL=paste("beta",normFlag,"_leuk",sep="")
        fNameL2="i.LEU.v2"
        datadirL3="docs/all/set1/preBcell/"
        fNameL3="0708011 Sample_Sheet (Fetal blood)"
        
        nRow=1000
    }
    #betaFileId=3
    betaFileId=1

    ## -------------------------------------------
    ## All smoke variables from set1 & set2
    ## All smoke variables in set2 are in set1
    varSmoke=data.frame(varIn=c("passive_sm_home_post","smoke_mo_ever","smoke_mo_preg","M_SM_BF_DI","smoke_mo_3months","smoke_mo_3months_N","smoke_mo_preg_N","m_cignum_bf","smoke_mo_after","smoke_mo_after_N","smoke_fa_ever","smoke_fa_3months","smoke_fa_3months_N","smoke_mo_bf","smoke_mo_bf_N","smoke3","smoke2"),
    varOut=c("paSmHmPost","moEver","moPreg","moBf","mo3m","mo3mN","moPregN","moBfN","moAfter","moAfterN","faEver","fa3m","fa3mN","moBf","moBfN","smoke3","smoke2"),stringsAsFactors=F)

    ## Smoke variables in common between set1 & set2
    varSmoke=data.frame(varIn=c("smoke_mo_ever","smoke_mo_preg","smoke_mo_3months","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_after","smoke_mo_after_N","smoke_fa_ever","smoke_fa_3months","smoke_fa_3months_N","smoke_mo_bf","smoke_mo_bf_N","smoke3","smoke2"),
    varOut=c("moEver","moPreg","mo3m","mo3mN","moPregN","moAfter","moAfterN","faEver","fa3m","fa3mN","moBf","moBfN","smoke3","smoke2"),stringsAsFactors=F)

    ## -------------------------------------------
    clin1=read.table(paste(dirClin1,"final.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    clin14=read.table(paste(dirClin1,"matches.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    if (computerFlag=="cluster") {
        clin12=read.table(paste(dirClin1,"i.LEU.v2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        clin13=read.table(paste(dirClin1,"Ritu_birthweight_dataset - covariate data.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        bw=read.table(paste("data/pobw-2014-12-26.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    } else {
        clin12=read.table(paste(dirClin1,"LEU.data/i.LEU.v2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        clin13=read.table(paste(dirClin1,"birthWeight/Ritu_birthweight_dataset - covariate data.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        bw=read.table(paste("docs/birthWeight/pobw-2014-12-26.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    }
    names(clin12)[match("Subject_ID",names(clin12))]=c("subjectId")
    names(clin12)[match(c("Plate","Sentrix_ID","Sentrix_Position"),names(clin12))]=c("Batch","Beadchip","Position")
    #clin2=read.table(paste(dirClin2,"clinGuthrieReplJune2012.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    clin2=read.table(paste(dirClin2,"clin_guthrieSet2_20140619.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    names(bw)[match("SubjectID",names(bw))]=c("subjectId")

    x=c("caco","birthWt","sex","int_ch_ethnicity","ch_hispanic_bc","int_ch_race","gestage")
    varib=c("caco","birthWt","sex","ethnicity","hispanicStatus","race","gestage")
    tbl=data.frame(varib,names(clin1)[match(x,names(clin1))],names(clin2)[match(x,names(clin2))],stringsAsFactors=F)
    names(tbl)=c("variable","nameInSet1","nameInSet2")
    #rownames(tbl)=tbl$variable
    tbl[which(tbl$variable=="caco"),"nameInSet1"]="Leukemia"
    tbl[which(tbl$variable=="birthWt"),"nameInSet1"]="birth_weight"
    tbl[which(tbl$variable=="birthWt"),"nameInSet2"]="birth_wt"
    #write.table(tbl,file="birthWtRelatedVariables.txt", sep="\t", col.names=T, row.names=F, quote=F)

    names(clin1)[match(c("TargetID","Subject_ID","Leukemia","mo_race_bc","fa_race_bc","birth_weight"),names(clin1))]=c("guthrieId","subjectId","caco","mo_race","fa_race","birthWt")
    names(clin13)[match(c("birth_wt"),names(clin13))]=c("birthWt")
    names(clin14)[match(c("Array.ID","Old.match","New.match"),names(clin14))]=c("arrayId","guthrieIdO","guthrieIdN")
    clin14$Beadchip=substr(clin14$arrayId,1,10)
    clin14$Position=substr(clin14$arrayId,12,17)
    clin14$type=substr(clin14$guthrieIdN,5,5)
    clin14[!clin14$type%in%c("G","C"),]
    clin14$type[!clin14$type%in%c("G","C")]=NA

    j=grep("DFE",names(clin2))
    names(clin2)[j]=sub("_new","",names(clin2)[j])
    names(clin2)[match(c("subjectID","Plate","birth_wt"),names(clin2))]=c("subjectId","Batch","birthWt")
    clin1$int_ch_ethnicity=clin1$ch_ethnicity

    clin1$id=paste("X",clin1$guthrieId,sep="")
    clin2$id=paste("X",clin2$guthrieId,sep="")

    j=match(clin1$subjectId,clin12$subjectId)
    j1=which(!is.na(j)); j2=j[j1]
    #tmp=vector(mode="character",length=nrow(clin1))
    tmp=rep("",nrow(clin1))
    clin1$BeadchipL=clin1$PositionL=tmp
    clin1$BeadchipL[j1]=clin12$Beadchip[j2]
    clin1$BeadchipL[clin1$BeadchipL==""]=NA
    clin1$PositionL[j1]=clin12$Position[j2]
    clin1$PositionL[clin1$PositionL==""]=NA
    clin1$Subtype=rep(NA,nrow(clin1))
    clin1$Subtype[j1]=clin12$Subtype[j2]
    clin1$Subtype[which(clin1$caco==0)]="control"
    clin2$Subtype=rep("",nrow(clin2))
    clin2$Subtype[which(clin2$smhyper==1)]="hyperdiploid"
    clin2$Subtype[which(clin2$smtelaml==1)]="telaml"
    clin2$Subtype[which(clin2$smhyper==0 & clin2$smtelaml==0)]="nonHypTelaml"
    clin2$Subtype[which(clin2$caco==0)]="control"

    clin1$hispNoHispWt=clin1$int_ch_ethnicity
    clin1$hispNoHispWt[!clin1$int_ch_ethnicity%in%c(1,2)]=NA
    clin2$hispNoHispWt=clin2$int_ch_ethnicity
    clin2$hispNoHispWt[!clin2$int_ch_ethnicity%in%c(1,2)]=NA

    clin1$hyperdipCtrl=as.integer(clin1$Subtype=="hyperdiploid")
    clin1$hyperdipCtrl[!clin1$Subtype%in%c("hyperdiploid","control")]=NA
    clin2$hyperdipCtrl=as.integer(clin2$Subtype=="hyperdiploid")
    clin2$hyperdipCtrl[!clin2$Subtype%in%c("hyperdiploid","control")]=NA

    clin1$telamlCtrl=as.integer(clin1$Subtype=="t1221")
    clin1$telamlCtrl[!clin1$Subtype%in%c("t1221","control")]=NA
    clin2$telamlCtrl=as.integer(clin2$Subtype=="telaml")
    clin2$telamlCtrl[!clin2$Subtype%in%c("telaml","control")]=NA

    clin1$nonHypTelamlCtrl=as.integer(clin1$Subtype%in%c("mll","others","t119"))
    clin1$nonHypTelamlCtrl[!clin1$Subtype%in%c("mll","others","t119","control")]=NA
    clin2$nonHypTelamlCtrl=as.integer(clin2$Subtype=="nonHypTelaml")
    clin2$nonHypTelamlCtrl[!clin2$Subtype%in%c("nonHypTelaml","control")]=NA

    clin1$hyperdipTelaml=as.integer(clin1$Subtype=="hyperdiploid")
    clin1$hyperdipTelaml[!clin1$Subtype%in%c("hyperdiploid","t1221")]=NA
    clin2$hyperdipTelaml=as.integer(clin2$Subtype=="hyperdiploid")
    clin2$hyperdipTelaml[!clin2$Subtype%in%c("hyperdiploid","telaml")]=NA

    clin1$hyperdipNonHypTelaml=as.integer(clin1$Subtype=="hyperdiploid")
    clin1$hyperdipNonHypTelaml[!clin1$Subtype%in%c("hyperdiploid","mll","others","t119")]=NA
    clin2$hyperdipNonHypTelaml=as.integer(clin2$Subtype=="hyperdiploid")
    clin2$hyperdipNonHypTelaml[!clin2$Subtype%in%c("hyperdiploid","nonHypTelaml")]=NA

    clin1$telamlNonHypTelaml=as.integer(clin1$Subtype=="t1221")
    clin1$telamlNonHypTelaml[!clin1$Subtype%in%c("t1221","mll","others","t119")]=NA
    clin2$telamlNonHypTelaml=as.integer(clin2$Subtype=="telaml")
    clin2$telamlNonHypTelaml[!clin2$Subtype%in%c("telaml","nonHypTelaml")]=NA

    clin1$hyperdipTelamlOther=as.integer(clin1$Subtype=="hyperdiploid")
    clin1$hyperdipTelamlOther[!clin1$Subtype%in%c("hyperdiploid","t1221","others")]=NA
    clin1$hyperdipTelamlOther=as.character(clin1$hyperdipTelamlOther)
    clin2$hyperdipTelamlOther=as.integer(clin2$Subtype=="hyperdiploid")
    clin2$hyperdipTelamlOther[!clin2$Subtype%in%c("hyperdiploid","telaml","others")]=NA
    clin2$hyperdipTelamlOther=as.character(clin2$hyperdipTelamlOther)

    clin1$hyperdipTelamlVsOther=as.integer(clin1$Subtype%in%c("hyperdiploid","t1221"))
    clin1$hyperdipTelamlVsOther[!clin1$Subtype%in%c("hyperdiploid","t1221","others")]=NA
    clin1$hyperdipTelamlVsOther=as.character(clin1$hyperdipTelamlVsOther)
    clin2$hyperdipTelamlVsOther=as.integer(clin2$Subtype%in%c("hyperdiploid","t1221"))
    clin2$hyperdipTelamlVsOther[!clin2$Subtype%in%c("hyperdiploid","telaml","others")]=NA
    clin2$hyperdipTelamlVsOther=as.character(clin2$hyperdipTelamlVsOther)

    j=match(clin1$id,clin13$targetID)
    j1=which(!is.na(j)); j2=j[j1]
    clin1$birthWtPerc=rep(NA,nrow(clin1))
    clin1$birthWtPerc[j1]=clin13$Percent[j2]
    clin2$birthWt=as.numeric(clin2$birthWt)

    j=match(clin1$subjectId,bw$subjectId)
    j1=which(!is.na(j)); j2=j[j1]
    clin1$pobw=clin1$pred_btw=clin1$dbirwt=rep(NA,nrow(clin1))
    clin1$dbirwt[j1]=bw$dbirwt[j2]
    clin1$pred_btw[j1]=bw$pred_btw[j2]
    clin1$pobw[j1]=bw$pobw[j2]
    j=match(clin2$subjectId,bw$subjectId)
    j1=which(!is.na(j)); j2=j[j1]
    clin2$pobw=clin2$pred_btw=clin2$dbirwt=rep(NA,nrow(clin2))
    clin2$dbirwt[j1]=bw$dbirwt[j2]
    clin2$pred_btw[j1]=bw$pred_btw[j2]
    clin2$pobw[j1]=bw$pobw[j2]

    clin1$Beadchip=clin1$Position=NA
    j=match(clin14$guthrieIdN,clin1$guthrieId); j2=which(!is.na(j)); j1=j[j2]
    clin1$Beadchip[j1]=clin14$Beadchip[j2]
    clin1$Position[j1]=clin14$Position[j2]

    k=match(names(clin1),names(clin2)); k1=which(!is.na(k)); k2=k[k1]
    tbl=rbind(cbind(clin1[,k1],set=rep("set1",nrow(clin1))),cbind(clin2[,k2],set=rep("set2",nrow(clin2))))
    write.table(tbl,file="clin_guthrieSet1Set2_20140619.txt", sep="\t", col.names=T, row.names=F, quote=F)

    beta1=read.table(paste(dirMeth1,fileList1[1],".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=10)
    beta2=read.table(paste(dirMeth2,fileList2[1],".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=10)
    clin1=clin1[match(colnames(beta1)[-1],clin1$id),]
    clin2=clin2[match(colnames(beta2)[-1],clin2$id),]

    if (subsetFlag=="case") {
        clin1=clin1[which(clin1$caco==1),]
        clin2=clin2[which(clin2$caco==1),]
    } else if (subsetFlag=="ctrl") {
        clin1=clin1[which(clin1$caco==0),]
        clin2=clin2[which(clin2$caco==0),]
        
    }

    clin=clin1
    clin$DFE_foodCat=as.integer(clin$DFE_Food>=400)
    clin$DFE_natCat=as.integer(clin$DFE_nat>=200)
    clin$DFE_totCat=as.integer(clin$DFE_tot>=400)
    clin$DFE_supCat=as.integer(clin$DFE_sup>0)
    clin$smoke_mo_3months=as.integer(clin$smoke_mo_3months_N>0)
    clin$smoke_mo_after=as.integer(clin$smoke_mo_after_N>0)
    clin$smoke_mo_bf=as.integer(clin$smoke_mo_bf_N>0)
    #names(clin)[match(varSmoke$varIn,names(clin))]=varSmoke$varOut
    ## Smoking 3 & 2 categories
    clin$smoke3=NA
    clin$smoke3[which(clin$smoke_mo_preg==0 & clin$smoke_mo_3months==0 & clin$smoke_fa_3months==0)]=0
    clin$smoke3[which(clin$smoke_mo_preg==1)]=1
    clin$smoke3[which(clin$smoke_mo_3months>0 & clin$smoke_mo_preg==0)]=2
    clin$smoke3[which(clin$smoke_fa_3months==1 & clin$smoke_mo_3months==0 & clin$smoke_mo_preg==0)]=3
    clin$smoke3[which(clin$smoke3==2)]=1
    clin$smoke3=factor(clin$smoke3)
    clin$smoke2=factor(as.integer(as.integer(as.character(clin$smoke3))>0))
    clin$race3=as.integer(as.character(clin$int_ch_race))
    clin$race3[which(clin$race3%in%c(2,3))]=5
    clin$race3=as.factor(clin$race3)
    clin1=clin

    clin=clin2
    clin$DFE_foodCat=as.integer(clin$DFE_Food>=400)
    clin$DFE_natCat=as.integer(clin$DFE_nat>=200)
    clin$DFE_totCat=as.integer(clin$DFE_tot>=400)
    clin$DFE_supCat=as.integer(clin$DFE_sup>0)
    clin$smoke_mo_3months=as.integer(clin$smoke_mo_3months_N>0)
    clin$smoke_mo_after=as.integer(clin$smoke_mo_after_N>0)
    clin$smoke_mo_bf=as.integer(clin$smoke_mo_bf_N>0)
    #names(clin)[match(varSmoke$varIn,names(clin))]=varSmoke$varOut
    ## Smoking 3 & 2 categories
    clin$smoke3=NA
    clin$smoke3[which(clin$smoke_mo_preg==0 & clin$smoke_mo_3months==0 & clin$smoke_fa_3months==0)]=0
    clin$smoke3[which(clin$smoke_mo_preg==1)]=1
    clin$smoke3[which(clin$smoke_mo_3months>0 & clin$smoke_mo_preg==0)]=2
    clin$smoke3[which(clin$smoke_fa_3months==1 & clin$smoke_mo_3months==0 & clin$smoke_mo_preg==0)]=3
    clin$smoke3[which(clin$smoke3==2)]=1
    clin$smoke3=factor(clin$smoke3)
    clin$smoke2=factor(as.integer(as.integer(as.character(clin$smoke3))>0))
    clin$race3=as.integer(as.character(clin$int_ch_race))
    clin$race3[which(clin$race3%in%c(2,3))]=5
    clin$race3=as.factor(clin$race3)
    clin2=clin

    ## -------------------------------------------
    ## Make the ethnicity proportions in Set 1 the same as in Set2

    if (subsetFlag=="propOfCacoEthnAsInSet2") {
        table(clin2$caco,clin2$int_ch_ethnicity)
        clin2$caco_ethn=clin2$caco*3+clin2$int_ch_ethnicity
        #table(clin2$caco_ethn,clin2$int_ch_ethnicity,clin2$caco,exclude=NULL)
        
        clin1$caco_ethn=clin1$caco*3+clin1$int_ch_ethnicity
        x1=table(clin1$caco_ethn)
        x2=sum(!is.na(clin1$caco_ethn))*table(clin2$caco_ethn)/sum(!is.na(clin2$caco_ethn))
        x2-x1 ## If there are any values > 0, adjust by the group with highest value
        x3=round((x1[1]/x2[1])*sum(!is.na(clin1$caco_ethn))*table(clin2$caco_ethn)/sum(!is.na(clin2$caco_ethn)))
        x3-x1 ## ALL <= 0
        
        
        j1=order(clin1$caco_ethn,clin1$guthrieId)
        j2=1:nrow(clin1)
        set.seed(54534)
        j2=sample(1:length(j1),length(j1),replace=F)
        j3=which(is.na(clin1$caco_ethn))
        grp=clin1$caco_ethn
        grpUniq=sort(unique(grp[!is.na(grp)]))
        for (grpId in 1:length(grpUniq)) {
            j3=c(j3,j2[which(grp[j2]==grpUniq[grpId])][1:x3[grpId]])
        }
        ## Same proportions in the 2 sets
        round(table(grp[j3],exclude=NULL)/sum(!is.na(grp[j3])),2)
        round(table(clin2$caco_ethn,exclude=NULL)/sum(!is.na(clin2$caco_ethn)),2)
        
        clin1=clin1[j3,]
        
        tbl=data.frame(guthrieId=paste("X",clin1$guthrieId,sep=""))
        write.table(tbl,file=paste("guthrieId_",subsetFlag,"_set1.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
    }

    invisible(list(clin1=clin1,clin2=clin2,dirClin1=dirClin1,dirClin2=dirClin2,dirMeth1,dirMeth2))
}

##############################################
