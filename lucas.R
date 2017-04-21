##############################################

dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
dirSrc3=paste("code/",sep="")
setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))

##############################################

computerFlag="cluster"
computerFlag=""

##############################################

source(paste(dirSrc3,"funcs.R",sep=""))
res=getClinData()
clin1=res$clin1
clin2=res$clin2
rm(res)

##############################################

##############################################
## ----------------------------------
## Compare fetal/adult proportions from Lucas's signature

library(coin)


if (setFlag=="") {
    header1=c("Set2: Guthrie","Set1: Guthrie","Set1: Guthrie (good samples)","Set1: Leukemia")
} else {
    header1=c("Set2: Guthrie")
}

datadir="docs/lucas/allGuthSet1/"
tbl1=read.table(paste(datadir,"Fetal signature Guthrie scaled.csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
names(tbl1)[match(c("X..","X.Adult.","X.Fetal."),names(tbl1))]=c("id","adult","fetal")
for (k in 1:ncol(tbl1)) if (is.character(tbl1[,k])) tbl1[,k]=gsub("\"","",tbl1[,k])
j=match(paste(clin1$Beadchip,clin1$Position,sep="_"),tbl1$id); j1=which(!is.na(j)); j2=j[j1]
#tmp=vector(mode="numeric",length=nrow(clin1))
tmp=rep(NA,nrow(clin1))
tbl=data.frame(adult=tmp,fetal=tmp,stringsAsFactors=F)
for (k in c("adult","fetal")) {
    tbl[j1,k]=tbl1[j2,k]
}
lucas1=tbl

datadir="docs/lucas/allGuthSet1/"
tbl1=read.table(paste(datadir,"Fetal signature Guthrie scaled QN without bad samples.csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
names(tbl1)[match(c("X..","X.Adult.","X.Fetal."),names(tbl1))]=c("id","adult","fetal")
for (k in 1:ncol(tbl1)) if (is.character(tbl1[,k])) tbl1[,k]=gsub("\"","",tbl1[,k])
j=match(paste(clin1$Beadchip,clin1$Position,sep="_"),tbl1$id); j1=which(!is.na(j)); j2=j[j1]
#tmp=vector(mode="numeric",length=nrow(clin1))
tmp=rep(NA,nrow(clin1))
tbl=data.frame(adult=tmp,fetal=tmp,stringsAsFactors=F)
for (k in c("adult","fetal")) {
    tbl[j1,k]=tbl1[j2,k]
}
lucas1=cbind(lucas1,tbl)

datadir="docs/lucas/allLeukSet1/"
tbl1=read.table(paste(datadir,"Fetal signature.csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
names(tbl1)[match(c("X","Adult","Fetal"),names(tbl1))]=c("id","adult","fetal")
j=match(paste(clin1$BeadchipL,clin1$PositionL,sep="_"),tbl1$id); j1=which(!is.na(j)); j2=j[j1]
#tmp=vector(mode="numeric",length=nrow(clin1))
tmp=rep(NA,nrow(clin1))
tbl=data.frame(adult=tmp,fetal=tmp,stringsAsFactors=F)
for (k in c("adult","fetal")) {
    tbl[j1,k]=tbl1[j2,k]
}
lucas1=cbind(lucas1,tbl)

names(lucas1)=paste(c("adult","fetal"),rep(c("G","GG","L"),each=2),sep="")

datadir="docs/lucas/allGuthSet2/"
lucas2=read.table(paste(datadir,"Fetal prop.csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
names(lucas2)[match(c("X","Adult","Fetal"),names(lucas2))]=c("id","adult","fetal")
lucas2=lucas2[match(paste(clin2$Beadchip,clin2$Position,sep="_"),lucas2$id),]

clin1=cbind(clin1,lucas1[,which(!names(lucas1)%in%names(clin1))])
clin2=cbind(clin2,lucas2[,which(!names(lucas2)%in%names(clin2))])

clin=clin1; datadir=datadir11
tbl=read.table(paste(datadir,"epistructure_allGuthSet1.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
id=clin$id[!clin$id%in%tbl$id]
tmp=tbl[1:length(id),]
for (k in 1:ncol(tmp)) tmp[,k]=NA
tmp$id=id
tbl=rbind(tbl,tmp)
clin=cbind(clin,tbl[match(clin$id,tbl$id),grep("epistr",names(tbl))])
clin1=clin

clin=clin2; datadir=datadir21
tbl=read.table(paste(datadir,"epistructure_allGuthSet2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
id=clin$id[!clin$id%in%tbl$id]
tmp=tbl[1:length(id),]
for (k in 1:ncol(tmp)) tmp[,k]=NA
tmp$id=id
tbl=rbind(tbl,tmp)
clin=cbind(clin,tbl[match(clin$id,tbl$id),grep("epistr",names(tbl))])
clin2=clin

for (plotType in c("qqplot_","densityPlot_")) {
    for (varType in c("","logit_")) {
        for (varName in c("Adult")) {
            png(paste(plotType,varType,varName,"Lucas.png",sep=""))
            par(mfrow=c(2,2))
            for (hId in 1:length(header1)) {
                switch(header1[hId],
                "Set2: Guthrie"={
                    clin=clin2
                },
                "Set1: Guthrie"={
                    clin=clin1; names(clin)[match(c("adultG","fetalG"),names(clin))]=c("adult","fetal")
                },
                "Set1: Guthrie (good samples)"={
                    clin=clin1; names(clin)[match(c("adultGG","fetalGG"),names(clin))]=c("adult","fetal")
                },
                "Set1: Leukemia"={
                    clin=clin1; names(clin)[match(c("adultL","fetalL"),names(clin))]=c("adult","fetal")
                }
                )
                if (varType=="logit_") {
                    if (F) {
                        par(mfrow=c(2,1))
                        offset=10^-16/100 ## < than the min for adult/fetal in set1/set2
                        x=c(clin1[,tolower(varName)],clin2[,tolower(varName)])
                        offset=min(x[which(x!=0)])
                    }
                    if (varName=="Adult") offset=0.001
                    header2=paste("Lucas: logit(",varName," proportion)",sep="")
                    varThis=paste(tolower(varName),"Logit",sep="")
                    clin[,varThis]=clin[,tolower(varName)]/100
                    j=which(clin[,tolower(varName)]==0); if (length(j)!=0) clin[,varThis][j]=offset
                    j=which(clin[,tolower(varName)]==100); if (length(j)!=0) clin[j,varThis]=1-offset
                    clin[,varThis]=log2(clin[,varThis]/(1-clin[,varThis]))
                    switch(header1[hId],
                    "Set2: Guthrie"={
                        clin2[,varThis]=clin[,varThis]
                    },
                    "Set1: Guthrie"={
                        clin1[,paste(varThis,"G",sep="")]=clin[,varThis]
                    },
                    "Set1: Guthrie (good samples)"={
                        clin1[,paste(varThis,"GG",sep="")]=clin[,varThis]
                    },
                    "Set1: Leukemia"={
                        clin1[,paste(varThis,"L",sep="")]=clin[,varThis]
                    }
                    )
                } else {
                    header2=paste("Lucas: ",varName," proportion",sep="")
                    varThis=tolower(varName)
                }
                if (plotType=="qqplot_") {
                    qqnorm(clin[,varThis],main=paste(header1[hId],"\n",header2,sep="")); qqline(clin[,varThis])
                } else {
                    #lim=range(c(as.matrix(cbind(clin1[,grep(varName,names(clin1))],clin2[,tolower(varName)]))),na.rm=T)
                    xlim=c(0,100)
                    ylim=c(0,.11)
                    j=which(clin[,varThis]<=xlim[2])
                    plot(density(clin[j,varThis],na.rm=T),xlim=xlim,ylim=ylim,main=header1[hId],xlab=paste(header2," - median ",round(median(clin[j,varThis],na.rm=T),2),sep=""))
                    #axis(side=3)
                    print(summary(clin[,varThis]))
                }
            }
        }
        dev.off()
    }
}


#for (subsetFlag in c("","_noPreterm")) {
for (subsetFlag in c("")) {
    varInfo=as.data.frame(matrix(
    c("sex","sex","ctrl","set2",
    "hispNoHispWt","hispNoHispWt","ctrl","set2",
    "epistr1","epistr1","ctrl","set2",
    "epistr2","epistr2","ctrl","set2",
    "epistr3","epistr3","ctrl","set2",
    "epistr4","epistr4","ctrl","set2",
    "epistr5","epistr5","ctrl","set2",
    "pobw","pobw","ctrl","set2",
    "pobw","pobw + epistr1 + epistr2 + epistr3 + epistr4 + epistr5","ctrl","set2",
    "birthWt","birthWt","ctrl","set2",
    "birthWt","birthWt + gestage","ctrl","set2",
    "ch_ageref","ch_ageref","ctrl","set2",
    "smoke_mo_ever","smoke_mo_ever","ctrl","set2",
    "DFE_nat","DFE_nat","ctrl","set2",
    "DFE_Food","DFE_Food","ctrl","set2",
    "DFE_fort","DFE_fort","ctrl","set2",
    "DFE_sup","DFE_sup","ctrl","set2",
    "DFE_tot","DFE_tot","ctrl","set2",
    "gestage","gestage","ctrl","set2",
    
    "sex","sex","case","set2",
    "gestage","gestage","case","set2",
    "pobw","pobw","case","set2",
    "pobw","pobw + epistr1 + epistr2 + epistr3 + epistr4 + epistr5","case","set2",
    "birthWt","birthWt","case","set2",
    "birthWt","birthWt + gestage","case","set2",
    "ch_ageref","ch_ageref","case","set2",
    "hyperdipTelamlVsOther","hyperdipTelamlVsOther","case","set2",
    "hyperdipTelaml","hyperdipTelaml","case","set2",
    
    "sex","sex","case+ctrl","set2",
    "hispNoHispWt","hispNoHispWt","case+ctrl","set2",
    "epistr1","epistr1","case+ctrl","set2",
    "epistr2","epistr2","case+ctrl","set2",
    "epistr3","epistr3","case+ctrl","set2",
    "epistr4","epistr4","case+ctrl","set2",
    "epistr5","epistr5","case+ctrl","set2",
    "caco","caco + ch_ageref","case+ctrl","set2",
    "caco","caco + sex + epistr1 + epistr2 + epistr3 + epistr4 + epistr5","case+ctrl","set2",
    "caco","adult * pobw","case+ctrl","set2",
    
    "sex","sex","case","set1",
    "gestage","gestage","case","set1",
    "pobw","pobw","case","set1",
    "pobw","pobw + epistr1 + epistr2 + epistr3 + epistr4 + epistr5","case","set1",
    "birthWt","birthWt","case","set1",
    "birthWt","birthWt + gestage","case","set1",
    "ch_ageref","ch_ageref","case","set1",
    "hyperdipTelamlVsOther","hyperdipTelamlVsOther","case","set1",
    "hyperdipTelaml","hyperdipTelaml","case","set1"
    ),
    ncol=4,byrow=T),stringsAsFactors=F)
    names(varInfo)=c("variable","model","subset","set")
    k=which(varInfo$set=="set2" & !varInfo$variable%in%c("caco"))
    tbl=varInfo[k,]; tbl$set="Set1: Guthrie"; varInfo=rbind(varInfo,tbl)
    tbl=varInfo[k,]; tbl$set="Set1: Guthrie (good samples)"; varInfo=rbind(varInfo,tbl)
    varInfo$set[which(varInfo$set=="set1")]="Set1: Leukemia"
    varInfo$set[which(varInfo$set=="set2")]="Set2: Guthrie"
    varName1=varInfo$variable
    varName2=c("adultLogit")
    varName2=c("adult")
    tbl=NULL
    for (k2 in 1:length(varName2)) {
        #cat("================ Linear regression test \n")
        for (k1 in 1:length(varName1)) {
            if (length(grep(varName2[k2],varInfo$model[k1]))) {
                testFlag="Logistic regression"
                modelC=paste(varName1[k1],"~",varInfo$model[k1])
            } else {
                testFlag="Linear regression"
                modelC=paste(varName2[k2],"~",varInfo$model[k1])
            }
            #cat("\n================",modelC,"================\n")
            modelF=as.formula(modelC)
            for (hId in 1:length(header1)) {
                if (varName1[k1]%in%c("caco") & header1[hId]%in%c("Set1: Leukemia")) next
                switch(header1[hId],
                "Set2: Guthrie"={
                    clin=clin2
                },
                "Set1: Guthrie"={
                    clin=clin1; names(clin)[match(c("adultG","fetalG"),names(clin))]=c("adult","fetal")
                },
                "Set1: Guthrie (good samples)"={
                    clin=clin1; names(clin)[match(c("adultGG","fetalGG"),names(clin))]=c("adult","fetal")
                },
                "Set1: Leukemia"={
                    clin=clin1; names(clin)[match(c("adultL","fetalL"),names(clin))]=c("adult","fetal")
                }
                )
                if (subsetFlag=="_noPreterm") {
                    j=which(clin$gestage>=37)
                } else {
                    j=1:nrow(clin)
                }
                switch(varInfo$subset[k1],
                "case"={
                    j=j[which(clin$caco[j]==1)]
                },
                "ctrl"={
                    j=j[which(clin$caco[j]==0)]
                }
                )
                clin=clin[j,]
                x=sum(!duplicated(clin[!is.na(clin[,varName1[k1]]),varName1[k1]]))
                x2=any(!is.na(clin[,varName2[k2]]))
                if (x>1 & x2) {
                    if (testFlag=="Linear regression") {
                        res=lm(modelF,data=clin)
                        coef=NA
                        if (is.numeric(clin[,varName1[k1]]) | x==2) {
                            res=summary(res)$coef
                            coef=res[2,"Estimate"]
                            pv=res[2,"Pr(>|t|)"]
                        } else {
                            pv=anova(res)[1,"Pr(>F)"]
                        }
                    } else {
                        res=summary(glm(caco ~ adult * pobw,data=clin,family="binomial"))$coef
                        k=nrow(res)
                        coef=res[k,"Estimate"]
                        pv=res[k,"Pr(>|z|)"]
                    }
                    ttl="P-value "
                    suf=""
                    if (pv<0.05) {
                        suf=" ****"
                    } else if (pv<0.1) {
                        suf=" **"
                    }
                    #cat(ttl,signif(pv,2),suf,"\n")
                    tbl1=c(header1[hId],varName1[k1],varName2[k2],testFlag,coef,pv,"",paste(signif(pv,2),suf,sep=""),modelC,ifelse(header1[hId]%in%varInfo$set[k1],1,0),varInfo$subset[k1])
                    tbl=rbind(tbl,tbl1)
                }
            }
        }
    }
    rownames(tbl)=NULL
    colnames(tbl)=c("set","variable1","variable2","test","coef","pValue","corrCoef","pv","model","keep","subset")
    tbl=as.data.frame(tbl,stringsAsFactors=F)
    #tbl=cbind(varInfo,tbl)
    for (colId in which(names(tbl)%in%c("coef","pValue","corrCoef"))) {
        tbl[,colId]=as.numeric(tbl[,colId])
    }
    for (colId in which(names(tbl)%in%c("pValue","corrCoef"))) {
        tbl[,colId]=round(tbl[,colId],2)
    }
    for (colId in which(names(tbl)%in%c("coef"))) {
        tbl[,colId]=signif(tbl[,colId],2)
    }
    tbl$model=sub("adultLogit","adult",tbl$model)
    k=which(tbl$subset!="")
    if (length(k)!=0) {
        tbl$set[k]=paste(tbl$set[k],", ",tbl$subset[k],sep="")
    }
    tbl$variableDescription=""
    tbl=tbl[which(tbl$keep==1),c("set","variable1","variableDescription","model","coef","pv")]
    k=which(tbl$variable1=="hispNoHispWt"); tbl$variableDescription[k]="Hispanic vs. non-hispanic white"
    k=which(tbl$variable1=="hyperdipTelamlOther"); tbl$variableDescription[k]="Hyperdiploid vs. tel/aml vs. other"
    k=which(tbl$variable1=="hyperdipTelaml"); tbl$variableDescription[k]="Hyperdiploid vs. tel/aml"
    names(tbl)[match(c("variable1"),names(tbl))]=c("variableOfInterest")
    tbl=tbl[order(tbl$model,tbl$set),]
    
    write.table(as.matrix(tbl),file=paste("stat_adultLucas",subsetFlag,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F, append=F)
    
    tbl[tbl$model!="adult ~ caco + sex + epistr1 + epistr2 + epistr3 + epistr4 + epistr5",c("set","variableOfInterest","model","coef","pv")]
    
    if (subsetFlag=="_noPreterm") {
        tblNPT=tbl
    } else {
        tbl0=tbl
    }
}

colId=c("set","model","coef","pv")

tbl2=read.table(paste("results/lucas/Archive/stat_adultLucas.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl2$set=sub("Leukemia","Set1: Leukemia",tbl2$set)
unique(tbl$model)[!(unique(tbl$model)%in%unique(tbl2$model))]
unique(tbl2$model)[!(unique(tbl2$model)%in%unique(tbl$model))]
sort(unique(tbl$set)[!(unique(tbl$set)%in%unique(tbl2$set))])
sort(unique(tbl2$set)[!(unique(tbl2$set)%in%unique(tbl$set))])
paste(tbl2$set,tbl2$model)[!(paste(tbl2$set,tbl2$model)%in%paste(tbl$set,tbl$model))]
paste(tbl$set,tbl$model)[!(paste(tbl$set,tbl$model)%in%paste(tbl2$set,tbl2$model))]


tbl1=tbl

tbl1=tbl0; tbl2=tblNPT
tbl2=read.table(paste("results/lucas/wrong/stat_logitAdultLucas.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl2=read.table(paste("results/lucas/Archive/stat_adultLucas.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl2$set=sub("Leukemia","Set1: Leukemia",tbl2$set)
tbl1$pv=as.numeric(sub(" **","",sub(" ****","",tbl1$pv,fixed=T),fixed=T))
tbl2$pv=as.numeric(sub(" **","",sub(" ****","",tbl2$pv,fixed=T),fixed=T))
pThres=0.05
k=match(paste(tbl1$set,tbl1$model),paste(tbl2$set,tbl2$model)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
    #if ((tbl1$pv[k1[k]]<pThres & tbl2$pv[k2[k]]>pThres) | (tbl1$pv[k1[k]]>pThres & tbl2$pv[k2[k]]<pThres)) {
    if (tbl1$pv[k1[k]]!=tbl2$pv[k2[k]]) {
        cat("\n\n========== ",tbl1$set[k1[k]],": ",tbl1$model[k1[k]],"\n",sep="")
        print(cbind(tbl1[k1[k],c("set","model","coef","pv")],tbl2[k2[k],c("coef","pv")]))
    }
}
k=tbl$model!="adult ~ caco + sex + epistr1 + epistr2 + epistr3 + epistr4 + epistr5"
cbind(tbl1[k,c("set","model","coef","pv")],tbl2[k,c("coef","pv")])



varName1=c("birthWt","gestage","pobw","ch_ageref","hyperdipTelamlVsOther")
for (subset1Flag in c("")) {
    if (subset1Flag=="_noPreterm") {
        tbl=tblNPT
    } else {
        tbl=tbl0
    }
    for (vId1 in 1:length(varName1)) {
        png(paste("plot_lucas_",varName1[vId1],subset1Flag,".png",sep=""),width=2*480, height=480)
        par(mfcol=c(2,4))
        xlim=NULL
        for (hId in 1:length(header1)) {
            switch(header1[hId],
            "Set2: Guthrie"={
                clin=clin2
            },
            "Set1: Guthrie"={
                clin=clin1; names(clin)[match(c("adultG","fetalG"),names(clin))]=c("adult","fetal")
            },
            "Set1: Guthrie (good samples)"={
                clin=clin1; names(clin)[match(c("adultGG","fetalGG"),names(clin))]=c("adult","fetal")
            },
            "Set1: Leukemia"={
                clin=clin1; names(clin)[match(c("adultL","fetalL"),names(clin))]=c("adult","fetal")
            }
            )
            if (subsetFlag=="_noPreterm") {
                j=which(clin$gestage>=37)
            } else {
                j=1:nrow(clin)
            }
            clin=clin[j,]
            for (subsetFlag in c("case","ctrl")) {
                switch(subsetFlag,
                "case"={
                    j=which(clin$caco==1)
                },
                "ctrl"={
                    j=which(clin$caco==0)
                }
                )
                k=which(tbl$set==paste(header1[hId],", ",subsetFlag,sep="") & tbl$model==paste("adult ~ ",varName1[vId1],sep=""))
                if (length(k)==1) {
                    if (is.numeric(clin[j,varName1[vId1]])) {
                        res=lm(as.formula(tbl$model[k]),data=clin[j,])
                        plot(clin[j,varName1[vId1]],clin$adult[j],xlim=xlim,main=paste(tbl$set[k],"\nmodel ",tbl$model[k],", pv ",tbl$pv[k],sep=""),xlab=varName1[vId1],ylab="Lucas: Adult proportion")
                        abline(c(res$coef),lty="solid",col="red")
                    } else {
                        ttl=sort(unique(clin[j,varName1[vId1]]))
                        if (varName1[vId1]=="hyperdipTelamlVsOther") ttl=c("other","hyperdiploid/tel-aml")
                        boxplot(clin$adult[j]~clin[j,varName1[vId1]],names=ttl,notch=T,ylim=xlim,main=paste(tbl$set[k],"\nmodel ",tbl$model[k],", pv ",tbl$pv[k],sep=""),xlab="",ylab="Lucas: Adult proportion")
                    }
                }
            }
        }
        dev.off()
    }
}

varList=c("adultG","adultGG","adultL")
varName=paste(c("Set1 Guthrie","Set1 Guthrie (good samples)","Get1 Leukemia"),": Adult proportion",sep="")
png(paste("plot_lucas_adult.png",sep=""),width=480, height=480)
par(mfcol=c(2,2))
lim=c(0,100)
lim=NULL
for (subset1Flag in c("")) {
    if (subset1Flag=="_noPreterm") {
        tbl=tblNPT
    } else {
        tbl=tbl0
    }
    for (vId1 in 1:(length(varList)-1)) {
        xlim=NULL
        for (vId2 in (vId1+1):length(varList)) {
            clin=clin1
            if (subsetFlag=="_noPreterm") {
                j=which(clin$gestage>=37)
            } else {
                j=1:nrow(clin)
            }
            res=lm(as.formula(paste(varList[vId2],"~",varList[vId1])),data=clin[j,])
            plot(clin[j,varList[vId1]],clin[j,varList[vId2]],xlim=lim,ylim=lim,main=paste("pv ",signif(summary(res)$coef[2,"Pr(>|t|)"],4),sep=""),xlab=varName[vId1],ylab=varName[vId2])
            abline(c(res$coef),lty="solid",col="red")
            print(summary(res)$coef)
        }
    }
}
dev.off()


clin=clin1
summary(lm(adultL~ch_ageref,data=clin))
summary(lm(adultL~adultG*ch_ageref,data=clin))
summary(lm(adultL~adultGG*ch_ageref,data=clin))
summary(lm(adultL~adultG+ch_ageref,data=clin))
summary(lm(adultL~adultG,data=clin))
summary(lm(adultL~adultGG+ch_ageref,data=clin))
summary(lm(adultL~adultGG,data=clin))
clin=clin2
summary(lm(adult~caco+ch_ageref,data=clin))
summary(lm(adultLogit~caco+ch_ageref,data=clin))

library(MASS)
colId=c("adult","caco","sex","mo_age","fa_age","gestage","smoke_mo_ever","smoke_mo_preg","smoke_fa_ever","smoke_fa_3months","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_bf_N","smoke_mo_after_N","smoke_fa_3months_N","rs1801131_1","rs1801131_2","rs1801133_1","rs1801133_2","DFE_Food","DFE_fort","DFE_nat","DFE_sup","DFE_tot","t1221","mll","hyperdiploid_51to67","wbc","Batch","Position","fa_race","mo_race","int_ch_ethnicity","Beadchip","ch_ageref","income","ch_hispanic_bc","int_ch_race","birthWt","Subtype","hispNoHispWt","hyperdipCtrl","telamlCtrl","nonHypTelamlCtrl","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","hyperdipTelamlOther","hyperdipTelamlVsOther","dbirwt","pred_btw","pobw","DFE_foodCat","DFE_natCat","DFE_totCat","DFE_supCat","smoke_mo_3months","smoke_mo_after","smoke_mo_bf","smoke3","smoke2","race3","epistr1","epistr2","epistr3","epistr4","epistr5")
colId=c("adult","sex","mo_age","fa_age","gestage","smoke_mo_ever","smoke_mo_preg","smoke_fa_ever","smoke_fa_3months","smoke_mo_3months_N","smoke_mo_preg_N","smoke_mo_bf_N","smoke_mo_after_N","smoke_fa_3months_N","rs1801131_1","rs1801131_2","rs1801133_1","rs1801133_2","DFE_Food","DFE_fort","DFE_nat","DFE_sup","DFE_tot","t1221","mll","hyperdiploid_51to67","wbc","Batch","Position","fa_race","mo_race","int_ch_ethnicity","Beadchip","ch_ageref","income","ch_hispanic_bc","int_ch_race","birthWt","Subtype","hispNoHispWt","hyperdipCtrl","telamlCtrl","nonHypTelamlCtrl","hyperdipTelaml","hyperdipNonHypTelaml","telamlNonHypTelaml","hyperdipTelamlOther","hyperdipTelamlVsOther","dbirwt","pred_btw","pobw","DFE_foodCat","DFE_natCat","DFE_totCat","DFE_supCat","smoke_mo_3months","smoke_mo_after","smoke_mo_bf","smoke3","smoke2","race3","epistr1","epistr2","epistr3","epistr4","epistr5")
colId=c("adult","sex","mo_age","fa_age","gestage","smoke_mo_ever","smoke_fa_ever","rs1801131_1","rs1801131_2","rs1801133_1","rs1801133_2","DFE_Food","DFE_fort","DFE_nat","DFE_sup","DFE_tot","wbc","Batch","Position","fa_race","mo_race","int_ch_ethnicity","Beadchip","ch_ageref","income","int_ch_race","birthWt","Subtype","dbirwt","pred_btw","pobw","epistr1","epistr2","epistr3","epistr4","epistr5")
j=which(clin2$caco==0)
clin=clin2[,colId]
j=1:nrow(clin)
for (k in 1:ncol(clin)) {
    j=j[which(!is.na(clin[j,k]))]
}
clin=clin[j,]
fit <- lm(adult ~ ., data = clin)
fit2 <- stepAIC(fit, k=log(nrow(clin)), trace = FALSE)
fit2$anova
##############################################
