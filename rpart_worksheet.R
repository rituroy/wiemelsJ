
####################################################################
####################################################################

library(rpart)

varFlag="pobw"
varFlag="pobwBi"
varFlag="flagged"

datType="_allGuthSet2"
subsetFlag="ctrl"

subset2Flag="_pobwLow0.9High1.1"
subset2Flag=""

subsetFFlag="_qv0.01"
subsetFFlag="_set2qv0.05_set1pv0.05"
subsetFFlag="_qv0.05"
subsetFFlag="_pv0.001"

rpartCtrlFlag="_cp0.05"
rpartCtrlFlag="_minbct10"
rpartCtrlFlag=""

subsetName=ifelse(subsetFlag=="","",paste("_",subsetFlag,"Subset",sep=""))
subsetName2=subset2Flag
fNameOut=paste("_",sub("~","Resp_",gsub("*","X",gsub("+","_",modelFlag,fixed=T),fixed=T)),subsetName,covFlag,covPCFlag,datType,subsetName2,normFlag,transformFlag,subsetFFlag,sep="")

x=strsplit(subsetFFlag,"_")[[1]]
if (subsetFFlag=="") {
    i=1:nrow()
} else if (subsetFFlag=="_set2qv0.05_set1pv0.05") {
    i=which(stat2[,substr(x[2],5,6)]<as.numeric(substr(x[2],7,nchar(x[2]))) & stat1[,substr(x[3],5,6)]<as.numeric(substr(x[3],7,nchar(x[3]))))
} else {
    i=which(stat2[,substr(subsetFFlag,2,3)]<as.numeric(substr(subsetFFlag,4,nchar(subsetFFlag))))
}

prId=match(stat2$cpgId[i],rownames(meth2))
samId=1:nrow(phen2)
switch(subsetFlag,
"ctrl"={
    samId=samId[which(phen2$caco==0)]
}
)
switch(subset2Flag,
"_pobwLow0.9High1.1"={
    samId=samId[which(phen2$pobw[samId]<=0.9 | phen2$pobw[samId]>=1.1)]
}
)

iA2=match(rownames(meth2)[prId],ann$IlmnID)
if (varFlag=="pobw") {
    varId="pobw"
    grp=phen2[samId,varId]
} else if (varFlag=="pobwBi") {
    varId="pobwBi"
    grp=as.character(phen2[samId,varId])
    grp[which(grp=="0")]="low pobw"
    grp[which(grp=="1")]="high pobw"
}
grpUniq=sort(unique(grp))
nm=ann$geneSym[iA2]
i=which(nm=="")
nm[i]=ann$IlmnID[iA2][i]
id=unique(nm[duplicated(nm)])
for (k in 1:length(id)) {
    i=which(nm==id[k])
    nm[i]=paste(nm[i],"_",1:length(i),sep="")
}


dat=as.data.frame(t(meth2[prId,samId]))
colnames(dat)=nm
dat=cbind(class=grp,dat)
dat2=dat

if (length(grep("_minbct",rpartCtrlFlag))==1) {
    fit=rpart(class ~ ., control=rpart.control(minbucket=as.integer(sub("_minbct","",rpartCtrlFlag))), data = dat)
} else if (length(grep("_cp",rpartCtrlFlag))==1) {
    fit=rpart(class ~ ., control=rpart.control(cp=as.numeric(sub("_cp","",rpartCtrlFlag))), data = dat)
} else {
    fit=rpart(class ~ ., data = dat)
}
fit2=fit
png(paste("rpart",fNameOut,rpartCtrlFlag,".png",sep=""),width=2.5*480,height=480)
par(xpd = NA) # otherwise on some devices the text is clipped
plot(fit)
text(fit, use.n = TRUE,cex=.9)
dev.off()
#write.table(print(fit),file="tmp.txt", sep="\t", col.names=F, row.names=F, quote=F, append=F)
#fit=fit2; tbl=read.table(paste("rpart",fNameOut,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)

dat=as.data.frame(t(meth2[prId,samId]))
dat=cbind(class=grp,dat)
dat2=dat
if (length(grep("_minbct",rpartCtrlFlag))==1) {
    fit=rpart(class ~ ., control=rpart.control(minbucket=as.integer(sub("_minbct","",rpartCtrlFlag))), data = dat)
} else if (length(grep("_cp",rpartCtrlFlag))==1) {
    fit=rpart(class ~ ., control=rpart.control(cp=as.numeric(sub("_cp","",rpartCtrlFlag))), data = dat)
} else {
    fit=rpart(class ~ ., data = dat)
}
fit1=fit
png(paste("rpart",fNameOut,"_cpgId",rpartCtrlFlag,".png",sep=""),width=2.5*480,height=480)
par(xpd = NA) # otherwise on some devices the text is clipped
plot(fit)
text(fit, use.n = TRUE,cex=.9)
dev.off()

datadir="rpart/"
datadir=""
sink(paste(datadir,"rpart",fNameOut,"_cpgId",rpartCtrlFlag,".txt",sep=""),append=F)
fit1
sink()

## ----------------------------------
## Get the prediction criteria

datadir="rpart/"
datadir=""

fit=fit1; tbl=read.table(paste(datadir,"rpart",fNameOut,"_cpgId",rpartCtrlFlag,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
x=sapply(tbl[5:nrow(tbl),1],function(x) {
    offset=2
    z=strsplit(strsplit(x,")")[[1]][1],"")[[1]]
    y=gsub(" +"," ",sub(" +","",sub(" +$", "", x)))
    for (gId in 1:length(grpUniq)) {
        y=sub(grpUniq[gId],sub(" ","",grpUniq[gId]),y)
    }
    y=strsplit(y," ")[[1]]
},USE.NAMES=F)


xx=tbl[5:nrow(tbl),1]
pp=c()
ppLeaf=c()
grpLeaf=c()
kk=1:16
kk=1:length(xx)
#for (k in 1:length(xx)) {
nodeFlag=c()
nodePrevFlag=0
for (k in kk) {
    x=xx[k]
    #pId=pPrevId+1
    offset=2
    z=strsplit(strsplit(x,")")[[1]][1],"")[[1]]
    pId=length(grep(" ",z))-offset
    z=nchar(strsplit(x,")")[[1]][1])+1
    nodeThisFlag=z
    
    y=gsub(" +"," ",sub(" +","",sub(" +$", "", x)))
    y=sub("> ",">",sub("< ","<",y))
    for (gId in 1:length(grpUniq)) {
        y=sub(grpUniq[gId],sub(" ","",grpUniq[gId]),y)
    }
    y=strsplit(y," ")[[1]]
    k2=grep("-",y[2])
    if (length(k2)==1) {
        y[2]=paste(sub("-","(-",y[2]),")",sep="")
    }
    #ppThis=paste(y[2:3],collapse="")
    ppThis=paste("dat$",y[2],sep="")
    
    if (nodePrevFlag<nodeThisFlag) {
        nodeFlag=c(nodeFlag,z)
        pp=c(pp,ppThis)
    } else {
        for (k2 in length(nodeFlag):0) {
            if (k2==0) {
                pp=c()
                nodeFlag=c()
                break
            }
            if (nodeFlag[k2]<nodeThisFlag) {
                pp=pp[1:k2]
                nodeFlag=nodeFlag[1:k2]
                break
            }
        }
        nodeFlag=c(nodeFlag,z)
        pp=c(pp,ppThis)
    }
    if (y[length(y)]=="*") {
        ppLeaf=c(ppLeaf,paste(pp,collapse=" & "))
        grpLeaf=c(grpLeaf,y[5])
    }
    nodePrevFlag=nodeThisFlag
}

## ----------------------------------
## Check observed vs. predicted class

if (T) {
    cat("\n\n=================== ",paste("rpart",fNameOut,"_cpgId",rpartCtrlFlag,sep="")," ==================",sep="")
    cat("\n\n=================== Set2 (training set) ==================\n\n",sep="")
    dat=dat2
    classObs=sub(" ","",dat$class)
    classObs=dat$class
    if (varFlag=="pobw") {
        classPred=rep(NA,nrow(dat))
    } else {
        classPred=rep("",nrow(dat))
    }
    for (k in 1:length(ppLeaf)) {
        a=paste("classPred[which(",ppLeaf[k],")]=grpLeaf[k]",sep="")
        eval(parse(text=a))
    }
    if (varFlag=="pobw") {
        classPred=as.numeric(classPred)
        cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
    } else {
        classPred[which(classPred=="lowpobw")]="low pobw"
        classPred[which(classPred=="highpobw")]="high pobw"
        print(table(classObs,classPred))
        cat("\nMisclassification rate: ",round(mean(classObs!=classPred),2),sep="")
    }
    classPred2=classPred
    
    ## -----------------------
    samId=1:nrow(phen1)
    switch(subsetFlag,
    "ctrl"={
        samId=samId[which(phen1$caco==0)]
    }
    )
    switch(subset2Flag,
    "_pobwLow0.9High1.1"={
        samId=samId[which(phen1$pobw[samId]<=0.9 | phen1$pobw[samId]>=1.1)]
    }
    )
    
    if (varFlag=="pobw") {
        varId="pobw"
        grp=phen1[samId,varId]
    } else if (varFlag=="pobwBi") {
        varId="pobwBi"
        grp=as.character(phen1[samId,varId])
        grp[which(grp=="0")]="low pobw"
        grp[which(grp=="1")]="high pobw"
    }
    dat=as.data.frame(t(meth1[prId,samId]))
    dat=cbind(class=grp,dat)
    dat1=dat
    
    cat("\n\n=================== Set1 (validation set) ==================\n\n",sep="")
    dat=dat1
    classObs=sub(" ","",dat$class)
    classObs=dat$class
    if (varFlag=="pobw") {
        classPred=rep(NA,nrow(dat))
    } else {
        classPred=rep("",nrow(dat))
    }
    for (k in 1:length(ppLeaf)) {
        a=paste("classPred[which(",ppLeaf[k],")]=grpLeaf[k]",sep="")
        eval(parse(text=a))
    }
    if (varFlag=="pobw") {
        classPred=as.numeric(classPred)
        cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
    } else {
        classPred[which(classPred=="lowpobw")]="low pobw"
        classPred[which(classPred=="highpobw")]="high pobw"
        print(table(classObs,classPred))
        cat("\nMisclassification rate: ",round(mean(classObs!=classPred),2),sep="")
    }
    classPred1=classPred
    
    ## -----------------------
    samId=which(phenV$tissue=="normal placenta tissue")
    samId=which(phenV$tissue=="normal whole cord blood")
    grp=phenV[samId,"birthWt"]
    cutoff=median(grp,na.rm=T)
    cutoff=3600
    grp[which(grp<=cutoff)]=0
    grp[which(grp>cutoff)]=1
    grp=as.character(grp)
    grp[which(grp=="0")]="low pobw"
    grp[which(grp=="1")]="high pobw"
    dat=as.data.frame(t(methV[prId,samId]))
    dat=cbind(class=grp,dat)
    datV=dat
    
    if (F) {
        cat("\n\n=================== Mulligan 2014 (validation set) ==================\n\n",sep="")
        dat=datV
        classObs=sub(" ","",dat$class)
        classObs=dat$class
        varFlag2="_continuous"
        varFlag2="_categorical"
        if (varFlag2=="_continuous") {
            classPred=rep(NA,nrow(dat))
        } else {
            classPred=rep("",nrow(dat))
        }
        for (k in 1:length(ppLeaf)) {
            a=paste("classPred[which(",ppLeaf[k],")]=grpLeaf[k]",sep="")
            eval(parse(text=a))
        }
        if (varFlag2=="_continuous") {
            classPred=as.numeric(classPred)
            cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
        } else {
            classPred[which(classPred=="lowpobw")]="low pobw"
            classPred[which(classPred=="highpobw")]="high pobw"
            print(table(classObs,classPred))
            cat("\nMisclassification rate: ",round(mean(classObs!=classPred),2),sep="")
        }
        classPredV=classPred
    }
}


