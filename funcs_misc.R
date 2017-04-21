
############################################

inferWBCbyLme = function(Y, pheno, modelFix, modelRand, coefWBC, 
contrastWBC=NULL, detailLevel=1, silentErrors=TRUE){
	
	
	Y=targetDataHNSCC_Assay[prId,]; pheno=targetDataHNSCC_Covariates; modelFix=y~caco; modelRnd=~1|Beadchip; coefWBC=coefEsts1[CpGSelection,]; contrastWBC=Lwbc; detailLevel=1; silentErrors=TRUE
	Y=targetDataHNSCC_Assay[prId,]; pheno=targetDataHNSCC_Covariates; modelFix=y~gestage; modelRnd=~1|Beadchip; coefWBC=coefEsts1[CpGSelection,]; contrastWBC=Lwbc; detailLevel=1; silentErrors=TRUE

	M = dim(coefWBC)[1]
	n = dim(Y)[2]
	if(dim(Y)[1] != M) stop("Y does not match coefWBC in dimension\n")
	
	if(detailLevel == -1){
		lmeVcomp = list(
						mu = matrix(NA,M,n),
						e = matrix(NA,M,n),
						a = list()
						)
	}
	
	sigmaResid = sigmaIcept = nObserved = nClusters = rep(NA, M)
	coefEsts = list()
	for(j in 1:M){
		ii = !is.na(Y[j,])
		nObserved[j] = sum(ii)
		pheno$y = Y[j,]
		
		try({
			fit = try(lme(modelFix, random=modelRand, data=pheno[ii,]), 
					  silent=silentErrors)
			
			if (inherits(fit,"try-error")) {
			fit = lm(modelFix, data=pheno[ii,])
			fitCoef = fit$coef
			sigmaResid[j] = summary(fit)$sigma
			sigmaIcept[j] = 0
			nClusters[j] = 0
			if(detailLevel == -1){
			lmeVcomp$mu[j,ii] = predict(fit)
			lmeVcomp$e[j,ii] = residuals(fit)
			lmeVcomp$a[[j]] = list()
			}
			} else {
			fitCoef = fit$coef$fixed
			sigmaResid[j] = fit$sigma
			sigmaIcept[j] = sqrt(getVarCov(fit)[1])
			nClusters[j] = length(fit$coef$random[[1]])
			if(detailLevel == -1){
			lmeVcomp$mu[j,ii] = predict(fit,level=0)
			lmeVcomp$e[j,ii] = residuals(fit,level=1)
			lmeVcomp$a[[j]] = fit$coef$random
			}
			}
			coefEsts[[j]] = fitCoef
			})
	}
	coefEsts = t(matrix(unlist(coefEsts), length(fitCoef), M))
	rownames(coefEsts) = rownames(coefWBC)
	colnames(coefEsts) = names(fitCoef)
	
	if(is.null(contrastWBC)) {
		Z = cbind(1,coefWBC)
	} else {
		Z = cbind(1, coefWBC %*% t(contrastWBC))
	}
	colnames(Z)[1] = "<Intercept>"
	
	ZtZ = t(Z) %*% Z
	ZtZqr = qr(ZtZ)
	G = solve(ZtZqr, t(Z) %*% coefEsts)
	
	out = list(method="lme", GammaMatrix=G, Beta1Matrix=coefEsts, 
			   sigma=cbind(resid=sigmaResid, intercept=sigmaIcept),
			   N=cbind(obs=nObserved,clusters=nClusters))
	
#	if(detailLevel>=1){
		out$var.unscaled=solve(ZtZqr)
		
		d = dim(Z)[2]
		out$predicted = Z %*% G
		residual2 = coefEsts - out$predicted
		out$Sigma = (1/(M-d)) * t(residual2)%*%residual2
		
		X = model.matrix(modelFix, pheno)
		Xbar = apply(X,2,mean)
		Xctr = (X - matrix(1,n,1) %*% Xbar)
		residual2z = Xctr %*% t(residual2)
		out$ssu = apply(residual2z*residual2z, 2, sum)
		
		residual1 = Y - coefEsts %*% t(X) 
		out$sse = apply(residual1*residual1, 1, sum, na.rm=TRUE)
		
		residual0 = Y - apply(Y,1,mean,na.rm=TRUE) %*% matrix(1,1,n)
		out$ss0 = apply(residual0*residual0, 1, sum, na.rm=TRUE)
		
		if(detailLevel>=1){
			out$residuals = list(stage1=residual1, stage2=residual2, inner=t(residual2z))
		}
	}
	if(detailLevel==-1){
		trueLme = which(sapply(lmeVcomp$a, length)>0)
		ch = sort(unique(unlist(lapply(lmeVcomp$a[trueLme], function(u)rownames(u[[1]])))))
		nch = length(ch)
		a = matrix(0, M, nch)
		colnames(a) = ch
		for(j in trueLme){
			aa = lmeVcomp$a[[j]][[1]]
			a[j,rownames(aa)] = aa
		}
		lmeVcomp$ch = names(lmeVcomp$a[[trueLme[1]]])
		lmeVcomp$a = a
		out$lmeVcomp = lmeVcomp
	}
	
	class(out) = "inferWBC"
	out
}

