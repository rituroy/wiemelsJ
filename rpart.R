
####################################################################
####################################################################

library(rpart)

datType="_allGuthSet2"

varResp="classObs"
varSamId="beadPos"
rpartCtrlFlag=""

nRep=2
nRep=10

dat=dat2
fNameOut=paste(datType,sep="")
fName2Out="_tmp"
modelThis=as.formula(paste(varResp,"~",paste(varList,collapse="+"),sep=""))
if (length(grep("_minbct",rpartCtrlFlag))==1) {
} else {
    ctrlThis=rpart.control(minsplit = 20, minbucket = round(20/3), cp = 0.01,maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval=nRep,surrogatestyle = 0, maxdepth = 30)
    ctrlThis=rpart.control(minsplit = 2, minbucket = 1, cp = 0.01,maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval=nRep,surrogatestyle = 0, maxdepth = 30)
    ctrlThis=rpart.control(minsplit = 10, minbucket = 2, cp = 0.01,maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval=nRep,surrogatestyle = 0, maxdepth = 30)
    ctrlThis=rpart.control(minsplit = 5, minbucket = 2, cp = 0.01,maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval=nRep,surrogatestyle = 0, maxdepth = 30)
}
fit=rpart(modelThis,data=dat,method = "class",control=ctrlThis)
fit2=fit
png(paste("rpart",fNameOut,rpartCtrlFlag,".png",sep=""),width=2.5*480,height=480)
par(xpd = NA) # otherwise on some devices the text is clipped
plot(fit)
text(fit, use.n = TRUE,cex=1.75)
dev.off()
#write.table(print(fit),file="tmp.txt", sep="\t", col.names=F, row.names=F, quote=F, append=F)
#fit=fit2; tbl=read.table(paste("rpart",fNameOut,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)

fit1=fit

dat=dat2
datadir="rpart/"
datadir=""
sink(paste(datadir,"rpart",fNameOut,fName2Out,rpartCtrlFlag,".txt",sep=""),append=F)
fit1
sink()

## ----------------------------------
## Get the prediction criteria

datadir="rpart/"
datadir=""

fit=fit1; tbl=read.table(paste(datadir,"rpart",fNameOut,fName2Out,rpartCtrlFlag,".txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
grpUniq=sort(unique(dat[,varResp]))
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
    cat("\n\n=================== ",paste("rpart",fNameOut,fName2Out,rpartCtrlFlag,sep="")," ==================",sep="")
    cat("\n\n=================== Set2 (training set) ==================\n\n",sep="")
    dat=dat22
    classObs=dat$classObs
    if (is.numeric(classObs)) {
        classPred=rep(NA,nrow(dat))
    } else {
        classPred=rep("",nrow(dat))
    }
    for (k in 1:length(ppLeaf)) {
        a=paste("classPred[which(",ppLeaf[k],")]=grpLeaf[k]",sep="")
        eval(parse(text=a))
    }
    cat("\nN: ",length(classObs),"\n",sep="")
    if (is.numeric(dat[,varResp])) {
        classPred=as.numeric(classPred)
        cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
    } else {
        print(table(classObs,classPred))
        cat("\nMisclassification rate: ",round(mean(classObs!=classPred,na.rm=T),2),sep="")
    }
    dat$classPred=classPred
    dat22=dat
    
    ## -----------------------
    cat("\n\n=================== Set1 (validation set) ==================\n\n",sep="")
    dat=dat12
    classObs=dat$classObs
    if (is.numeric(classObs)) {
        classPred=rep(NA,nrow(dat))
    } else {
        classPred=rep("",nrow(dat))
    }
    for (k in 1:length(ppLeaf)) {
        a=paste("classPred[which(",ppLeaf[k],")]=grpLeaf[k]",sep="")
        eval(parse(text=a))
    }
    cat("\nN: ",length(classObs),"\n",sep="")
    if (is.numeric(dat[,varResp])) {
        classPred=as.numeric(classPred)
        cat("\nCor: ",round(cor(classObs,classPred,use="complete.obs"),2),"\n")
    } else {
        print(table(classObs,classPred))
        cat("\nMisclassification rate: ",round(mean(classObs!=classPred,na.rm=T),2),sep="")
    }
    dat$classPred=classPred
    dat12=dat
}

"
=================== rpart_allGuthSet2_tmp ==================



=================== Set2 (training set) ==================

Without missing data:

classPred
classObs   0   1
0 409   2
1   8  16

Misclassification rate: 0.02

With missing data:
N: 451
classPred
classObs   0   X
0 423   3
X   9  16

Misclassification rate: 0.03

=================== Set1 (validation set) ==================


N: 479
classPred
classObs   0   X
0 452   8
X  11   8

Misclassification rate: 0.04

"
