############################################################################################3
## GLLVM estimation variational approximation (VA). 
## VA distribution is normal with vameans as mean vector and vacov covariance matrix
## Allows Poisson, negative.binomial, such that Var = mu + mu^2/phi, Bernoulli using probit regression

## maxit is not the maximum number of total iterations, but rather the max number of sub-iterations for each Quasi-Newton method
##############################################################################################
library(mvtnorm)
library(mvabund)
library(MASS)
library(Matrix)

glvm.va <- function(y, # Response matrix n*m
					X = NULL, # Optional, n*p matrix => factor structure?
					family = "poisson", 
					num.lv = 2, # n of factor
					max.iter = 200, # mac iter of update iteration to perform
					eps = 1e-4, # convergence criterion
					row.eff = FALSE, # optional row effect
					covmat.struc = "unstructured", # or "diagonal"
					trace = TRUE, # Do you want to indicate estimation log?
					plot = FALSE, 
					sd.errors = FALSE, # se for model
					maxit = 100 # Max iter to parform for all updates involving Quasi-Newton
					) {
     n<-dim(y)[1]; p<-dim(y)[2]; # nは受験者数，pは観測変数（項目）数
	 num.X <- 0
	 # predictor が指定する場合には行列に変換してnum.Xにその数（列数）を入れる。そうでなければ0が指定される。
     if(!is.null(X)) { X <- as.matrix(X); num.X <- dim(X)[2] } ## no of predictors

	 # 使えるfamilyはpoisson, neg.binomial, binomial, ordinal （これは論文で提案されている手法に一致しているようだ）
     if(!(family %in% c("poisson","negative.binomial","binomial","ordinal"))) 
		stop("Inputed family not allowed...sorry =(")
	if(!(covmat.struc %in% c("unstructured","diagonal"))) 
		stop("A_i, i.e. covariance of vartiational distribution for latent variable, must be either unstructured or diagonal in form. Thanks.")
	if(num.lv == 1) covmat.struc <- "diagonal" ## Prevents it going to "unstructured" loops and causing chaos
	trial.size <- 1
	
	y <- as.matrix(y)
	if(!is.numeric(y)) stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
	if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
	if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")
	if(!is.null(X)) if(is.null(colnames(X))) colnames(X) <- paste("x",1:ncol(X),sep="")

	if(!is.null(X)) X <- as.matrix(X)# 不要？
	if(family == "ordinal") { 
		max.levels <- apply(y,2,max); #最大カテゴリ数を取得
		# どうやら現状だと複数のfamilyを同時にはサポートしていない。 => inefficientなだけで，2パラとGRMの混合推定はできそう？
		if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks") 
	}
	

	## Get initial starting values by fitting a bunch of GLMs
	res <- start.values.va(y = y, X = X, family = family, trial.size = trial.size, num.lv = num.lv)
	new.beta0 <- beta0 <- res$params[,1]
	new.beta <- beta <- NULL; if(!is.null(X)) new.beta <- beta <- res$params[,2:(num.X+1)];
	new.row.params <- row.params <- NULL; if(row.eff) new.row.params <- row.params <- rep(0,n)
	# means <- 変分分布の平均
	# lambda <- 識別力
	# vacov <- 
	new.vameans <- vameans <- new.lambda <- lambda <- new.vacov <- vacov <- NULL

	if(num.lv > 0) {
		new.vameans <- vameans <- res$index; 
		new.lambda <- lambda <- as.matrix(res$params[,(ncol(res$params)-num.lv+1):ncol(res$params)]);
		new.lambda[upper.tri(new.lambda)] <- lambda[upper.tri(lambda)] <- 0 # おそらく識別のための制約
		if(covmat.struc == "unstructured") { 
			new.vacov <- vacov <- array(NA,dim=c(n,num.lv,num.lv)); # 受験者ごとにcovmatを作る。
			for(i in 1:n) { 
				new.vacov[i,,] <- vacov[i,,] <- diag(rep(1,num.lv)) 
				}
			}
		if(covmat.struc == "diagonal") { 
			new.vacov <- vacov <- matrix(1,n,num.lv); 
			}
		zero.cons <- which(new.lambda == 0) 
		}

	new.zeta <- zeta <- NULL; if(family == "ordinal") { new.zeta <- zeta <- res$zeta }
	new.phi <- phi <- NULL; if(family == "negative.binomial") { new.phi <- phi <- res$phi }

	
	
	current.loglik <- -1e6; iter <- 1; err <- 10; diag.iter <- 5
	while((err > (1 + eps) || err < (1 - eps)) && iter <= max.iter) {
#	while(err > eps && iter <= max.iter) {
		if(trace) cat("Iteration:", iter, "\n")

		
		## Use a diagonal covariance matrix for the first diag.iter iterations, then switch to unstructured covariance matrix. 
		if(covmat.struc == "unstructured" & iter <= diag.iter) { 
			tmp.covmat.struc <- "diagonal"; 
			if(iter == 1) { new.vacov <- vacov <- matrix(1,n,num.lv) }
			if(iter == diag.iter) { tmp.covmat.struc <- "unstructured"; new.vacov <- vacov <- vacov.convert(vacov,type=2) } 
			}
		else { tmp.covmat.struc <- covmat.struc }
		
	
		## Calculate VA log-likelihood 
		ll0 <- function(x,v=NULL,vacov=NULL,phi,zeta=NULL) {
			x2 <- x
			new.phi <- phi
			new.vacov <- vacov
			new.vameans <- new.lambda <- NULL
			if(num.lv > 0) { 
				new.vameans <- matrix(v,n,num.lv); 
				new.lambda <- matrix(c(x2[1:(p*num.lv)]),p,num.lv); x2 <- x2[-(1:(p*num.lv))] }
			new.zeta <- zeta
			new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
			new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] }
			new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }

			mu.mat <- matrix(new.beta0,n,p,byrow=T) 
			if(!is.null(X)) mu.mat <- mu.mat + X%*%t(new.beta)
			if(row.eff) mu.mat <- mu.mat + matrix(new.row.params,n,p,byrow=F)
			if(num.lv > 0) mu.mat <- mu.mat  + new.vameans%*%t(new.lambda) 
			eta.mat <- mu.mat 
			if(num.lv > 0) eta.mat <- eta.mat + calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat

			
			if(family=="poisson") { out1 <- sum(y * mu.mat - lfactorial(y) - exp(eta.mat)) }

			if(family=="negative.binomial") { 
				phi.mat <- matrix(new.phi,n,p,byrow=T)
				if(num.lv > 0) eta.mat <- mu.mat - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat
				out1 <- sum(-phi.mat * mu.mat - lfactorial(y) + phi.mat * log(phi.mat) - lgamma(phi.mat) - (y+phi.mat) * log(1+phi.mat/exp(eta.mat)) + lgamma(y + phi.mat)) }
				
			if(family=="binomial") {        
				out1 <- dbinom(as.matrix(y),size=trial.size,prob=pnorm(mu.mat),log=T)
				out1 <- sum(out1[is.finite(out1)]) 
				if(num.lv > 0) out1 <- out1 - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat.sum }
	
			if(family == "ordinal") { 
				out1 <- matrix(NA,n,p); 
				for(j in 1:p) { 
					out1[y[,j]==1,j] <- pnorm(new.zeta[j,1]-mu.mat[y[,j] == 1,j],log=T)
					out1[y[,j]==max(y[,j]),j] <- log(1-pnorm(new.zeta[j,max(y[,j])-1]-mu.mat[y[,j]==max(y[,j]),j]))
					if(max(y[,j]) > 2) { 
						j.levels <- 2:(max(y[,j])-1)
						for(k in j.levels) { out1[y[,j] == k,j] <- log(pnorm(new.zeta[j,k]-mu.mat[y[,j]==k,j]) - pnorm(new.zeta[j,k-1]-mu.mat[y[,j]==k,j])) } }
					}
				out1 <- sum(out1[is.finite(out1)]) 					
				if(num.lv > 0) out1 <- out1 - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat.sum }

	
			if(num.lv > 0) {
				if(tmp.covmat.struc == "unstructured") { foo2 <- function(i) { 0.5*(log(det(new.vacov[i,,])) - sum(diag(new.vacov[i,,])) - sum(new.vameans[i,]^2)) } }
				
				if(tmp.covmat.struc == "diagonal") { foo2 <- function(i) { 0.5*(sum(log(t(new.vacov)[,i])) - sum(t(new.vacov)[,i]) - sum(new.vameans[i,]^2)) } }

				out1 <- out1 + sum(sapply(1:n,foo2)) 
				}
					
					
			ll <- out1
			return(ll) }
			

		if(num.lv > 0) {
			## Update vameans (a_i) by maximizing lower bound for loglik
			grad.var <- function(x,v=NULL,vacov=NULL,phi,zeta=NULL) {
				x2 <- x
				new.phi <- phi
				new.vacov <- vacov
				new.vameans <- new.lambda <- NULL
				if(num.lv > 0) { 
					new.vameans <- matrix(v,n,num.lv); 
					new.lambda <- matrix(c(x2[1:(p*num.lv)]),p,num.lv); x2 <- x2[-(1:(p*num.lv))] }
				new.zeta <- zeta
				new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
				new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] }
				new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }

				grad.vameans <- grad.vacov <- NULL
				eta.mat <- matrix(new.beta0,n,p,byrow=T) + new.vameans %*% t(new.lambda) + calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat 
				if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
				if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)

				
				if(family=="poisson") { 
					sum1 <- (y-exp(eta.mat))
					for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(sum1,2,new.lambda[,l],"*"))-new.vameans[,l]) } }

				if(family=="negative.binomial") {        
					phi.mat <- matrix(new.phi,n,p,byrow=T)
					eta.mat <- matrix(new.beta0,n,p,byrow=T) + new.vameans %*% t(new.lambda) - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat 
					if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
					if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
	
					sum1 <- -phi.mat - (y+phi.mat)*(-phi.mat/(exp(eta.mat)+phi.mat))
					for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(sum1,2,new.lambda[,l],"*"))-new.vameans[,l]) } } 
			
				if(family=="binomial") { 
					eta.mat <- eta.mat - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat 
					sum1 <- dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10)
					for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(sum1,2,new.lambda[,l],"*"),na.rm=T)-new.vameans[,l]) } }

				if(family=="ordinal") { 
					eta.mat <- eta.mat - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat 
					deriv.trunnorm <- matrix(NA,n,p)
					for(j in 1:p) { 
						deriv.trunnorm[y[,j]==1,j] <- -dnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])/pnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])
						deriv.trunnorm[y[,j]==max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j])/(1-pnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j]))
						if(max(y[,j])>2) {
							j.levels <- 2:(max(y[,j])-1)
							for(k in j.levels) { deriv.trunnorm[y[,j]==k,j] <- (-dnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])+dnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) / (pnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])-pnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) } }
						}
					deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0	

					for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,rowSums(sweep(deriv.trunnorm,2,new.lambda[,l],"*"))-new.vameans[,l]) } }
				
			return(c(grad.vameans))
			}

			q <- try(optim(c(vameans), method = "BFGS", fn = ll0, gr=grad.var, x=c(lambda,beta0,beta,row.params), vacov = vacov, phi=phi, zeta = zeta, control = list(trace=0, fnscale=-1, maxit = maxit, reltol=1e-2)),silent = T)
		
			if(!inherits(q, "try-error")) {
				if(trace) cat("Variational parameters updated", "\n")
				new.vameans <- q$par[1:(num.lv*n)]; new.vameans <- matrix(new.vameans,n,num.lv) }
			else { new.vameans <- vameans }
             
		
			## Update covariance matrix A_i via iterative fixed-point algorithm
			eta.mat <- cbind(rep(1,n),X)%*%t(cbind(beta0,beta)) + new.vameans%*%t(lambda) + calc.quad(vacov,lambda,tmp.covmat.struc)$mat
			if(!is.null(row.params)) eta.mat <- eta.mat + matrix(row.params,n,p)
			if(family == c("poisson")) { mu.mat <- exp(eta.mat) }
			if(family == c("negative.binomial")) { 
				phi.mat <- matrix(phi,n,p,byrow=T)
				eta.mat <- eta.mat - 2*calc.quad(vacov,lambda,tmp.covmat.struc)$mat
				mu.mat <- (y+phi.mat)*(phi.mat/(exp(eta.mat)+phi.mat)) }

			for(i in 1:n) {
				error <- 1; vacov.iter <- 0; 
				if(tmp.covmat.struc == "unstructured") new.vacov.mat <- vacov[i,,]
				if(tmp.covmat.struc == "diagonal") new.vacov.mat <- diag(x=vacov[i,],nrow=num.lv)

				while(error > 1e-2 && vacov.iter < 100) {	
					cw.vacov.mat <- new.vacov.mat
		
					if(tmp.covmat.struc == "unstructured") { 
						lambda2 <- sapply(1:p,function(j,lambda) lambda[j,]%*%t(lambda[j,]), lambda = new.lambda)
						lambda2 <- t(lambda2) 
					
						if(family %in% c("poisson","negative.binomial")) { new.vacov.mat <- solve(diag(rep(1,num.lv)) + matrix(apply(mu.mat[i,]*lambda2,2,sum),nrow=num.lv)) }
						if(family %in% c("binomial","ordinal")) { new.vacov.mat <- solve(diag(rep(1,num.lv)) + matrix(apply(lambda2,2,sum),nrow=num.lv)) } }
					
					if(tmp.covmat.struc == "diagonal") {
						lambda2 <- new.lambda^2
						if(family %in% c("poisson","negative.binomial")) { new.vacov.mat <- solve(diag(rep(1,num.lv)) + diag(apply(mu.mat[i,]*lambda2,2,sum),num.lv,num.lv)) }
						if(family %in% c("binomial","ordinal")) { new.vacov.mat <- solve(diag(rep(1,num.lv)) + diag(apply(lambda2,2,sum),num.lv,num.lv)) } }
					
					error <-  sum((new.vacov.mat-cw.vacov.mat)^2) 
					if(tmp.covmat.struc == "unstructured") vacov[i,,] <- new.vacov.mat 
					if(tmp.covmat.struc == "diagonal") vacov[i,] <- diag(new.vacov.mat) 

					eta.mat <- cbind(rep(1,n),X)%*%t(cbind(beta0,beta)) + new.vameans%*%t(lambda) + calc.quad(vacov,lambda,tmp.covmat.struc)$mat
					if(!is.null(row.params)) eta.mat <- eta.mat + matrix(row.params,n,p)
					if(family == c("poisson")) { mu.mat <- exp(eta.mat) }
					if(family == c("negative.binomial")) { 
						phi.mat <- matrix(phi,n,p,byrow=T)
						eta.mat <- eta.mat - 2*calc.quad(vacov,lambda,tmp.covmat.struc)$mat
						mu.mat <- (y+phi.mat)*(phi.mat/(exp(eta.mat)+phi.mat)) }
				
					vacov.iter <- vacov.iter + 1
					}

				if(tmp.covmat.struc == "unstructured") new.vacov[i,,] <- new.vacov.mat 
				if(tmp.covmat.struc == "diagonal") new.vacov[i,] <- diag(new.vacov.mat) 
				
				if((family %in% c("binomial","ordinal")) & i == 1) break; ## Since same A_i matrix for all i, then do one and then break 
				}
				
			if(family %in% c("binomial","ordinal")) { 
				if(tmp.covmat.struc == "diagonal") for(i2 in 2:n) new.vacov[i2,] <- new.vacov[1,]			
				if(tmp.covmat.struc == "unstructured") for(i2 in 2:n) new.vacov[i2,,] <- new.vacov[1,,] }
			
               } ## Ends loop for if num.lv > 0


		## Update model params (lambda_j, beta_j, alpha_i) by maximizing lower bound for loglik
		grad.mod <- function(x,v=NULL,vacov=NULL,phi,zeta=NULL) {
			x2 <- x
			new.phi <- phi
			new.vacov <- vacov
			new.vameans <- new.lambda <- NULL
			if(num.lv > 0) { 
				new.vameans <- matrix(v,n,num.lv); 
				new.lambda <- matrix(c(x2[1:(p*num.lv)]),p,num.lv); x2 <- x2[-(1:(p*num.lv))] }
			new.zeta <- zeta
			new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
			new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] }
			new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }

			grad.lambda <- grad.beta0 <- grad.beta <- grad.row.params <- NULL 

			if(tmp.covmat.struc == "unstructured" & num.lv > 0) { covmat.lambda <- array(NA,dim=c(p,n,num.lv)); for(i in 1:n) { covmat.lambda[,i,] <- new.lambda%*%new.vacov[i,,] } }
			
			eta.mat <- matrix(new.beta0,n,p,byrow=T) 
			if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
			if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
			if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.lambda) + calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat 
			
			
			if(family=="poisson") { 
				grad.beta0 <- colSums(y - exp(eta.mat))

				if(num.lv > 0) {
					for(l in 1:num.lv) {
						if(tmp.covmat.struc == "unstructured") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep(t(covmat.lambda[,,l]),1,new.vameans[,l],"+")*exp(eta.mat)
						if(tmp.covmat.struc == "diagonal") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep((new.vacov[,l]%*%t(new.lambda[,l])),1,new.vameans[,l],"+")*exp(eta.mat)
						grad.lambda <- c(grad.lambda,colSums(sum1)) } }
					
				if(!is.null(X)) { for (l in 1:num.X) { sum1 <- sweep(y-exp(eta.mat),1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1)) } } 
				if(row.eff) { grad.row.params <- rowSums(y-exp(eta.mat)) }  
				}

			if(family=="negative.binomial") { 
				phi.mat <- matrix(new.phi,n,p,byrow=T)
				eta.mat <- matrix(new.beta0,n,p,byrow=T) 
				if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
				if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
				if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.lambda) - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat 

				grad.beta0 <- colSums(-phi.mat - (y+phi.mat)*(-phi.mat/(exp(eta.mat)+phi.mat)))

				if(num.lv > 0) {
					for(l in 1:num.lv) {
						if(tmp.covmat.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(covmat.lambda[,,l]),1,-new.vameans[,l],"+")*(-(y+phi.mat)*(phi.mat/(exp(eta.mat)+phi.mat)))
						if(tmp.covmat.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.vacov[,l] %*% t(new.lambda[,l])),1,-new.vameans[,l],"+")*(-(y+phi.mat)*(phi.mat/(exp(eta.mat)+phi.mat)))
						grad.lambda <- c(grad.lambda,colSums(sum1)) } }
 
				if(!is.null(X)) {          
					for(l in 1:num.X) { sum1 <- sweep(-phi.mat-(y+phi.mat)*(-phi.mat/(exp(eta.mat)+phi.mat)),1,X[,l],"*")
						grad.beta <- c(grad.beta,colSums(sum1)) } }	 
				if(row.eff) { grad.row.params <- rowSums(-phi.mat-(y+phi.mat)*(-phi.mat/(exp(eta.mat)+phi.mat))) }  
				}

			if(family=="binomial") { 
				eta.mat <- matrix(new.beta0,n,p,byrow=T) 
				if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
				if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
				if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.lambda)
				
				grad.beta0 <- colSums(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),na.rm=T) 

				if(num.lv > 0) {
					for(l in 1:num.lv) {
						if(tmp.covmat.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,new.vameans[,l],"*") - t(covmat.lambda[,,l])
						if(tmp.covmat.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,new.vameans[,l],"*") - (new.vacov[,l] %*% t(new.lambda[,l]))
						grad.lambda <- c(grad.lambda,colSums(sum1,na.rm=T)) } }
					
				if(!is.null(X)) {          
					for (l in 1:num.X) { sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1,na.rm=T)) } } 
				if(row.eff) { 
					grad.row.params <- rowSums(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),na.rm=T) }  
				}
		
			if(family=="ordinal") { 
				eta.mat <- matrix(new.beta0,n,p,byrow=T) 
				if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
				if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
				if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.lambda)
				deriv.trunnorm <- matrix(0,n,p)
				
				for(j in 1:p) { 
					deriv.trunnorm[y[,j]==1,j] <- -dnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])/pnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])
					deriv.trunnorm[y[,j]==max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j])/(1-pnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j]))
					if(max(y[,j])>2) {
						j.levels <- 2:(max(y[,j])-1)
						for(k in j.levels) { deriv.trunnorm[y[,j]==k,j] <- (-dnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])+dnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) / (pnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])-pnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) } } }
				deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0	

				grad.beta0 <- colSums(deriv.trunnorm) 
						
				if(num.lv > 0) { 
					for(l in 1:num.lv) {
						if(tmp.covmat.struc == "unstructured") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - t(covmat.lambda[,,l])
						if(tmp.covmat.struc == "diagonal") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - (new.vacov[,l]%*%t(new.lambda[,l]))
						grad.lambda <- c(grad.lambda,colSums(sum1,na.rm=T)) }
					}

				if(!is.null(X)) { for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1,na.rm=T)) } } 
				if(row.eff) { grad.row.params <- rowSums(deriv.trunnorm,na.rm=T) }  
				}

				
			return(c(grad.lambda, grad.beta0, grad.beta, grad.row.params))
			}

		lower.vec <- rep(-30,length(c(lambda,beta0,beta,row.params))); upper.vec <- rep(30,length(c(lambda,beta0,beta,row.params)))
		lower.vec[which(upper.tri(lambda))] <- -1e-2; upper.vec[which(upper.tri(lambda))] <- 1e-2
		q <- try(optim(c(lambda,beta0,beta,row.params), v=c(new.vameans), vacov = new.vacov, phi=phi, zeta = zeta, method="L-BFGS-B", lower = lower.vec, upper = upper.vec, fn=ll0, gr=grad.mod, control=list(trace=0, fnscale=-1, maxit = maxit, factr=1e-2)), silent=T)

		if(!inherits(q, "try-error")) {
			if(trace) cat("Model parameters updated", "\n")
			x2 <- q$par
			if(num.lv > 0) { new.lambda <- matrix(x2[1:(p*num.lv)],p,num.lv); new.lambda[upper.tri(new.lambda)] <- 0; x2 <- x2[-(1:(p*num.lv))] }
			new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
			new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] } 
			new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
			}
		else { new.lambda <- lambda; new.beta0 <- beta0; new.beta <- beta; new.row.params <- row.params }
		
		
		## Update overdispersion parameter if necessary
		if(family=="negative.binomial") {
			grad.phi <- function(x,v=NULL,vacov=NULL,phi,zeta = NULL) {
				x2 <- x
				new.phi <- phi
				new.vacov <- vacov
				new.vameans <- new.lambda <- NULL
				if(num.lv > 0) { 
					new.vameans <- matrix(v,n,num.lv); 
					new.lambda <- matrix(c(x2[1:(p*num.lv)]),p,num.lv); x2 <- x2[-(1:(p*num.lv))] }
				new.zeta <- zeta
				new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
				new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] }
				new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
				phi.mat <- matrix(new.phi,n,p,byrow=T)
			
			
				eta.mat <- matrix(new.beta0,n,p,byrow=T) 
				if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
				if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
				if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.lambda) 
				mu.mat <- eta.mat
				if(num.lv > 0) mu.mat <- eta.mat - calc.quad(new.vacov,new.lambda,tmp.covmat.struc)$mat 

				out <- colSums(-eta.mat + 1 + log(phi.mat) - digamma(phi.mat) - log(1+phi.mat/exp(mu.mat)) - (y+phi.mat)/(phi.mat+exp(mu.mat)) + digamma(y+phi.mat))
				return(out)				
				}
				
			q <- try(optim(phi, x = c(new.lambda,new.beta0,new.beta,new.row.params), v=c(new.vameans), vacov = new.vacov, zeta = zeta, method="L-BFGS-B", lower = 1e-3, upper = 1e3, fn=ll0, gr=grad.phi, control=list(trace=0, fnscale=-1, factr = 1e-3)), silent=T)
			ll0( x = c(new.lambda,new.beta0,new.beta,new.row.params), v=c(vameans), vacov = vacov, phi = phi)
			
			if(!inherits(q, "try-error")) { if(trace) cat("Dispersion parameters updated", "\n"); new.phi <- q$par[1:p] } else { new.phi <- phi }
			}
					
		
		## Update cutoffs for ordinal data (zeta_j) if required
		if(family == "ordinal") {
			eta.mat <- cbind(rep(1,n),X)%*%t(cbind(new.beta0,new.beta)) 
			if(!is.null(row.params)) eta.mat <- eta.mat + matrix(new.row.params,n,p)
			if(num.lv > 0) eta.mat <- eta.mat + new.vameans%*%t(new.lambda) 
			
			func.zetaj <- function(cw.zeta, j) { ## Exclude the first column of zeta, which is fixed at 0
				zeta0 <- c(0,cw.zeta); out <- 0; 
				out <- out + sum(pnorm(zeta0[1]-eta.mat[which(y[,j] == 1),j],log=T))
				out <- out + sum(log(1-pnorm(zeta0[max(y[,j])-1]-eta.mat[which(y[,j] == max(y[,j])),j]))) 
				if(max(y[,j])>2) {
					j.levels <- 2:(max(y[,j])-1)
					for(k in j.levels) { out <- out + sum(log(pnorm(zeta0[k]-eta.mat[y[,j]==k,j])-pnorm(zeta0[k-1]-eta.mat[y[,j]==k,j]))) } }
				out }

			grad.zetaj <- function(cw.zeta, j) { ## Exclude the first column of zeta, which is fixed at 0
				zeta0 <- c(0,cw.zeta); deriv.trunnorm <- numeric(length(zeta0)); ## L-1 length
				deriv.trunnorm[length(deriv.trunnorm)] <- deriv.trunnorm[length(deriv.trunnorm)] - sum(dnorm(zeta0[max(y[,j])-1]-eta.mat[which(y[,j] == max(y[,j])),j])/(1-pnorm(zeta0[max(y[,j])-1]-eta.mat[which(y[,j] == max(y[,j])),j]))) 
				if(max(y[,j])>2) {
					j.levels <- 2:(max(y[,j])-1)
					for(k in j.levels) {
						deriv.trunnorm[k] <- deriv.trunnorm[k] + sum(dnorm(zeta0[k]-eta.mat[y[,j]==k,j])/(pnorm(zeta0[k]-eta.mat[y[,j]==k,j])-pnorm(zeta0[k-1]-eta.mat[y[,j]==k,j])),na.rm=T) 
						deriv.trunnorm[k-1] <- deriv.trunnorm[k-1] - sum(dnorm(zeta0[k-1]-eta.mat[y[,j]==k,j])/(pnorm(zeta0[k]-eta.mat[y[,j]==k,j])-pnorm(zeta0[k-1]-eta.mat[y[,j]==k,j])),na.rm=T) } }
				deriv.trunnorm[-1] }
			
			for(j in 1:p) {
				if(max(y[,j]) == 2) { new.zeta[j,] <- zeta[j,] }
				if(max(y[,j]) > 2) {
					constraint.mat <- matrix(0,max(y[,j])-2,max(y[,j])-2); constraint.mat[1,1] <- 1 ## Constrains zeta2 > 0
					if(nrow(constraint.mat) > 1) { 
						for(k in 2:nrow(constraint.mat)) constraint.mat[k,(k-1):k] <- c(-1,1) }

					update.zeta <- constrOptim(theta = zeta[j,2:(max(y[,j])-1)], f = func.zetaj, grad = grad.zetaj, ui = constraint.mat, ci = rep(0,max(y[,j])-2),j = j,outer.eps = 1e-3,control=list(trace=0, fnscale = -1))
					if(!inherits(update.zeta, "try-error")) new.zeta[j,2:(max(y[,j])-1)] <- update.zeta$par
					if(inherits(update.zeta, "try-error")) new.zeta[j,] <- zeta[j,] } 				
				}	
			cat("Cutoffs updated \n")
			}
		
			
		## Take values of loglik from optim -function to define stopping rule
		q <- list(value = ll0(c(new.lambda,new.beta0,new.beta,new.row.params), v=c(new.vameans), vacov = new.vacov, phi=new.phi, zeta = new.zeta))
		new.loglik <- q$value
 		err <- abs(new.loglik/current.loglik); 
  		if(trace) cat("New Loglik:", new.loglik,"Current Loglik:", current.loglik, "Ratio", err, "\n")

		
		## Plot old and new ordination points for spp's and sites
		if(trace == TRUE && num.lv <= 2 && num.lv > 0 && plot == TRUE) { 
			par(mfrow=c(1,2))
			plot(new.vameans[!duplicated(new.vameans),],type="n",xlab=ifelse(num.lv==1,"Row index (unique elements only)","LV1"),ylab="LV2",main="Ordination of rows")
			text(new.vameans[!duplicated(new.vameans),],labels=which(!duplicated(new.vameans)==1))

			plot(new.lambda,type="n",xlab=ifelse(num.lv==1,"Column index","Coefs for LV1"),ylab="Coefs for LV2",main="Ordination of columns")
			text(new.lambda,labels=seq(1,dim(y)[2])) }
	
	
		current.loglik <- new.loglik; 
		beta0 <- new.beta0; 
		beta <- new.beta; 
		lambda <- new.lambda; 
		phi <- new.phi; 
		vameans <- new.vameans; 
		vacov <- new.vacov; 
		row.params <- new.row.params
		zeta <- new.zeta
		
		iter = iter + 1
		} 

		
     ## Bling up the output 
	spp.coefs <- matrix(beta0,ncol=1); rownames(spp.coefs) <- colnames(y); colnames(spp.coefs) <- "Intercept"
	if(!is.null(X)) { spp.coefs <- cbind(beta0,beta); colnames(spp.coefs) <- c("Intercept",colnames(X)); rownames(spp.coefs) <- colnames(y) } 
	 
	out.list <-  list(y = y, X = X, num.lv = num.lv, row.eff = row.eff, beta = spp.coefs, iter=iter, logLik=new.loglik, family = family, covmat.struc = covmat.struc)
     if(num.lv > 0) { 
		rownames(lambda) <- colnames(y); colnames(lambda) <- 1:num.lv
		lambda[upper.tri(lambda)] <- 0
		rownames(vameans) <- rownames(y); colnames(vameans) <- 1:num.lv 
		out.list$A = vacov; out.list$lambda = lambda; out.list$lvs = vameans
		out.list$lambda[!is.finite(out.list$lambda)] <- 0
		out.list$lvs[!is.finite(out.list$lvs)] <- 0 }

	if(family == "negative.binomial") { out.list$phi <- 1/phi; names(out.list$phi) <- colnames(y); } ## Flip it back so that it is V = mu + phi*mu^2
	if(row.eff) { out.list$row.params <- row.params; names(row.params) <- rownames(y) }

	if(family == "ordinal") { 
		out.list$zeta <- zeta; rownames(out.list$zeta) <- colnames(y); 
		colnames(out.list$zeta) <- paste(1:(max(y)-1),"|",2:max(y),sep="") }

		
	if(sd.errors) {
		cat("Calculating standard errors for parameters...\n")
		get.sds <- calc.infomat(lambda, beta0, beta, row.params, vameans, vacov, phi, zeta, num.lv = num.lv, family = family, covmat.struc = covmat.struc, row.eff = row.eff, y = y, X = X)
		out.list$se <- get.sds
		}

		
     return(out.list)
     }

	
	
## Starting values by fitting marginal GLMs
start.values.va <- function(y, X = NULL, family, trial.size = 1, num.lv = 0) {
     N <- nrow(y); p <- ncol(y)
     y <- as.matrix(y)
     if(!is.numeric(y)) 
		stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
	if(!(family %in% c("poisson","negative.binomial","binomial","ordinal"))) 
		stop("inputed family not allowed...sorry =(")

     if(num.lv > 0) { 
		unique.ind <- which(!duplicated(y))
		rs <- as.vector(rowSums(y, na.rm = TRUE))
		len.uni <- length(unique(rs))
		rs <- factor(rs, labels = 1:len.uni)
		rs <- as.numeric(levels(rs))[as.integer(rs)]
		index <- matrix(seq(-3, 3, len = len.uni)[rs], ncol=1)
		if(num.lv > 1) { index <- cbind(index,rmvnorm(nrow(index),rep(0,num.lv-1))) }		
		unique.index <- as.matrix(index[unique.ind,]) ## Construct the starting latent variables in a ``clever" way
		}
	if(num.lv == 0) { index <- NULL }

	
     y <- as.matrix(y)
	if (family == "ordinal") {
	    max.levels <- apply(y, 2, max)
	    if (any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
	}
	if(is.null(rownames(y))) rownames(y) <- paste("row",1:n,sep="")
	if(is.null(colnames(y))) colnames(y) <- paste("col",1:p,sep="")

	
     options(warn = -1)

	if(family != "ordinal") { ## Using logistic instead of probit regession here for binomial
		if(!is.null(X) & num.lv > 0) fit.mva <- manyglm(y~X+index,family=family,K=trial.size)  
		if(is.null(X) & num.lv > 0) fit.mva <- manyglm(y~index,family=family,K=trial.size)  
		if(!is.null(X) & num.lv == 0) fit.mva <- manyglm(y~X,family=family,K=trial.size)  
		if(is.null(X) & num.lv == 0) fit.mva <- manyglm(y~1,family=family,K=trial.size)  
		params <- t(fit.mva$coef) 
		}
	if(family == "negative.binomial") { phi <- fit.mva$phi } else { phi <- NULL }
	
	if(family == "ordinal") {
		max.levels <- max(y)
		params <- matrix(NA,p,ncol(cbind(1,X,index))) # lambda?
		zeta <- matrix(NA,p,max.levels-1) ## max.levels # beta?
		zeta[,1] <- 0 ## polr parameterizes as no intercepts and all cutoffs vary freely. Change this to free intercept and first cutoff to zero 

		for(j in 1:p) {  
			# 素点を低い順に並べて，そこに-3から3の仮のパラメタを与え(index)，orderd logisticを実行して，それを項目の初期値計算に使う。
			# 戻り値のパラメトリゼーションは，`logit P(Y <= k | x) = zeta_k - eta`
			y.fac <- factor(y[,j]); 
			if(length(levels(y.fac)) > 2) { 
				if(!is.null(X) & num.lv > 0) cw.fit <- polr(y.fac ~ X+index, method = "probit")
				if(is.null(X) & num.lv > 0) cw.fit <- polr(y.fac ~ index, method = "probit") 
				if(!is.null(X) & num.lv == 0) cw.fit <- polr(y.fac ~ X, method = "probit")
				if(is.null(X) & num.lv == 0) cw.fit <- polr(y.fac ~ 1, method = "probit") 
				params[j,] <- c(cw.fit$zeta[1],-cw.fit$coefficients) 
				zeta[j,2:length(cw.fit$zeta)] <- cw.fit$zeta[-1]-cw.fit$zeta[1] 
				}
			if(length(levels(y.fac)) == 2) {				
				if(!is.null(X) & num.lv > 0) cw.fit <- glm(y.fac ~ X+index, family = binomial(link = "probit"))
				if(is.null(X) & num.lv > 0) cw.fit <- glm(y.fac ~ index, family = binomial(link = "probit")) 
				if(!is.null(X) & num.lv == 0) cw.fit <- glm(y.fac ~ X, family = binomial(link = "probit"))
				if(is.null(X) & num.lv == 0) cw.fit <- glm(y.fac ~ 1, family = binomial(link = "probit")) 
				params[j,] <- cw.fit$coef 
				}
			} 
		}

     out <- list(params=params,phi=phi) # phiはordinalではNULL(only overdispersed counts model)
	if(num.lv > 0) out$index <- index
     if(family == "ordinal") out$zeta <- zeta
     options(warn = 0)
     
     return(out) 
}

     
     
## Calculates the commonly used (1/2) lambda'_j A_i lambda_j 
## vacov is a n x num.lv x nmu.lv array, lambda is p x num.lv matrix
calc.quad <- function(vacov,lambda,covmat.struc) {
	if(covmat.struc == "diagonal") out <- 0.5*(vacov)%*%t(lambda^2) 

	if(covmat.struc == "unstructured") {
		if(class(vacov) == "array") { n <- dim(vacov)[1]; num.lv <- dim(vacov)[2] }
		if(class(vacov) == "matrix") { num.lv <- dim(vacov)[2]; n <- 1 }
		if(class(lambda) == "matrix") { p <- dim(lambda)[1] }
		if(class(lambda) == "numeric") { p <- 1; lambda <- matrix(lambda,1) }	
	
		out <- matrix(NA,n,p)
		if(n == 1) { for(j in 1:p) { out[1,j] <- 0.5*t(lambda[j,])%*%vacov%*%lambda[j,] } }
		if(n > 1) { 		
			if(n < p) out <- t(sapply(1:n, function(x) 0.5*rowSums(lambda*(lambda%*%vacov[x,,])))) 
			if(n > p) { 
				vacov.mat <- aperm(vacov,c(3,2,1)); dim(vacov.mat) <- c(num.lv,num.lv*n) 
				f <- function(x) 0.5*rowSums(matrix((t(vacov.mat)*lambda[x,])%*%lambda[x,],ncol=num.lv,byrow=T))
				out <- sapply(1:p,f) } }
		}		

	return(list(mat = out, mat.sum = sum(out))) 
	}
	
	
	
## Calculate info mat based on scores of VA log-likelihood
calc.infomat <- function(theta = NULL, beta0, beta = NULL, row.params = NULL, vameans = NULL, lambda=NULL, phi = NULL, zeta=NULL, num.lv, family, Lambda.struc, row.eff, y, X = NULL) { 
	n <- nrow(y); p <- ncol(y)
	if(!is.null(X)) num.X <- ncol(X)
	trial.size <- 1
	if(family == "ordinal") { 
		max.levels <- max(y); 
		if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks") 
		}
	
	
	## score of VA logL wrt model parameters
 	#A.mat <- -nd2(x0=c(theta,beta0,beta,row.params,phi,zeta2), f = grad.all.mod, vameans=vameans, lambda = lambda) 
 	#x=c(theta,beta0,beta,row.params,phi,t(zeta2))
	grad.all.mod <- function(x,vameans=NULL,lambda=NULL) {
		x2 <- x
		new.vameans <- new.lambda <- new.theta <- NULL
		if(num.lv > 0) { 
			new.lambda <- lambda
			new.vameans <- vameans
			new.theta <- matrix(c(x2[1:(p*num.lv)]),p,num.lv); x2 <- x2[-(1:(p*num.lv))] 
			}
		new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
		new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] }
		new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
		new.phi <- NULL; if(family == "negative.binomial") { new.phi <- x2[1:p]; x2 <- x2[-(1:p)] }
		new.zeta <- NULL; if(family == "ordinal") { 
			new.zeta <- matrix(NA,p,max.levels-2); 
			for(j in 1:p) { if(max(y[,j]) > 2) { new.zeta[j,1:(max(y[,j])-2)] <- x2[1:(max(y[,j])-2)]; x2 <- x2[-(1:(max(y[,j])-2))] } }
			new.zeta <- cbind(0,new.zeta); 
			}

		
		eta.mat <- matrix(new.beta0,n,p,byrow=T) 
		if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
		if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta) + calc.quad(new.lambda,new.theta,Lambda.struc)$mat 				
		if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
		
		grad.theta <- grad.beta0 <- grad.beta <- grad.phi <- grad.row.params <- grad.zeta <- NULL;

		if(Lambda.struc == "unstructured" & num.lv > 0) { 
			Lambda.theta <- array(NA,dim=c(p,n,num.lv)); 
			for(i2 in 1:n) { Lambda.theta[,i2,] <- new.theta%*%new.lambda[i2,,] } 
			}
			
		if(family=="poisson") { 
			grad.beta0 <- colSums(y-exp(eta.mat))

			if(num.lv > 0) {
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep(t(Lambda.theta[,,l]),1,new.vameans[,l],"+")*exp(eta.mat)
					if(Lambda.struc == "diagonal") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep((new.lambda[,l]%*%t(new.theta[,l])),1,new.vameans[,l],"+")*exp(eta.mat)
					grad.theta <- c(grad.theta,colSums(sum1)) } }

			if(!is.null(X)) { for (l in 1:num.X) { sum1 <- sweep(y-exp(eta.mat),1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1)) } } 
			if(row.eff) { grad.row.params <- rowSums(y-exp(eta.mat)) } }

			
		if(family=="negative.binomial") { 
			phi.mat <- matrix(new.phi,n,p,byrow=T)
			eta.mat <- matrix(new.beta0,n,p,byrow=T) 
			if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
			if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
			if(num.lv > 0) mu.mat <- eta.mat + new.vameans %*% t(new.theta) - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
			if(num.lv > 0) mu.mat.noquad <- eta.mat + new.vameans %*% t(new.theta) 
			if(num.lv == 0) mu.mat <- eta.mat  

			grad.beta0 <- colSums(-phi.mat - (y+phi.mat)*(-phi.mat/(exp(mu.mat)+phi.mat)))

			if(num.lv > 0) {
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(Lambda.theta[,,l]),1,-new.vameans[,l],"+")*(-(y+phi.mat)*(phi.mat/(exp(mu.mat)+phi.mat)))
					if(Lambda.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.lambda[,l] %*% t(new.theta[,l])),1,-new.vameans[,l],"+")*(-(y+phi.mat)*(phi.mat/(exp(mu.mat)+phi.mat)))
					grad.theta <- c(grad.theta,colSums(sum1)) } }
 
			if(!is.null(X)) {          
				for(l in 1:num.X) { sum1 <- sweep(-phi.mat-(y+phi.mat)*(-phi.mat/(exp(mu.mat)+phi.mat)),1,X[,l],"*")
					grad.beta <- c(grad.beta,colSums(sum1)) } }	 
			if(row.eff) { grad.row.params <- rowSums(-phi.mat-(y+phi.mat)*(-phi.mat/(exp(mu.mat)+phi.mat))) }  

			grad.phi <- colSums(-mu.mat.noquad + 1 + log(phi.mat) - digamma(phi.mat) - log(1+phi.mat/exp(mu.mat)) - (y+phi.mat)/(phi.mat+exp(mu.mat)) + digamma(y+phi.mat)) }

			
		if(family=="binomial") { 
			if(num.lv > 0) eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
				
			grad.beta0 <- colSums(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))),na.rm=T)

			if(num.lv > 0) {
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
					if(Lambda.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
					grad.theta <- c(grad.theta,colSums(sum1,na.rm=T)) } }
					
			if(!is.null(X)) { 
				for(l in 1:num.X) { 
					sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1,na.rm=T)) } } 
			if(row.eff) { grad.row.params <- rowSums((dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat)))),na.rm=T) } }
		
		
		if(family=="ordinal") { 
			grad.zeta <- matrix(0,p,ncol(new.zeta))
			if(num.lv > 0) eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 

			deriv.trunnorm <- matrix(0,n,p)
			for(j in 1:p) { 
				deriv.trunnorm[y[,j]==1,j] <- -dnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])/pnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])
				deriv.trunnorm[y[,j]==max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j])/(1-pnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j]))
				j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
				for(k in j.levels) {  
					deriv.trunnorm[y[,j]==k,j] <- (-dnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])+dnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) / (pnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])-pnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) } 
				}
			deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0	

			grad.beta0 <- colSums(deriv.trunnorm,na.rm=T)
						
			if(num.lv > 0) { 
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
					if(Lambda.struc == "diagonal") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - (new.lambda[,l]%*%t(new.theta[,l]))
					grad.theta <- c(grad.theta,colSums(sum1,na.rm=T)) }
				}

			if(!is.null(X)) { for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1,na.rm=T)) } } 
			if(row.eff) { grad.row.params <- rowSums(deriv.trunnorm,na.rm=T) }  
			
			for(j in 1:p) { 
				zeta0 <- new.zeta[j,]; 
				grad.zeta[j,max(y[,j])-1] <- grad.zeta[j,max(y[,j])-1] - sum(dnorm(zeta0[max(y[,j])-1]-eta.mat[which(y[,j] == max(y[,j])) ,j])/(1-pnorm(zeta0[max(y[,j])-1]-eta.mat[which(y[,j] == max(y[,j])),j]))) 
				j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
				for(k in j.levels) {
					grad.zeta[j,k] <- grad.zeta[j,k] + sum(dnorm(zeta0[k]-eta.mat[y[,j]==k,j])/(pnorm(zeta0[k]-eta.mat[y[,j]==k,j])-pnorm(zeta0[k-1]-eta.mat[y[,j]==k,j]))) 
					grad.zeta[j,k-1] <- grad.zeta[j,k-1] - sum(dnorm(zeta0[k-1]-eta.mat[y[,j]==k,j])/(pnorm(zeta0[k]-eta.mat[y[,j]==k,j])-pnorm(zeta0[k-1]-eta.mat[y[,j]==k,j]))) } }
				
			grad.zeta <- grad.zeta[,-1]		
			if(length(unique(apply(y,2,max))) > 1) { grad.zeta <- c(grad.zeta[-which(grad.zeta == 0)]) }
			}
			
		return(c(grad.theta, grad.beta0, grad.beta, grad.row.params, grad.phi, grad.zeta))
		}

		
		
	## score of VA logL wrt vameans_sel.i and lambda_sel.i
	grad.all.var <- function(x,vameans=NULL,lambda=NULL,mod.x,sel.i) {
		x2 <- mod.x
		if(num.lv > 0) { 
			new.vameans <- vameans; new.lambda <- lambda
			new.theta <- matrix(c(x2[1:(p*num.lv)]),p,num.lv); x2 <- x2[-(1:(p*num.lv))] }
		new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
		new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] }
		new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
		new.phi <- NULL; if(family == "negative.binomial") { new.phi <- x2[1:p]; x2 <- x2[-(1:p)] }
		new.zeta <- NULL; if(family == "ordinal") { 
			new.zeta <- matrix(NA,p,max.levels-2); 
			for(j in 1:p) { if(max(y[,j]) > 2) { new.zeta[j,1:(max(y[,j])-2)] <- x2[1:(max(y[,j])-2)]; x2 <- x2[-(1:(max(y[,j])-2))] } }
			new.zeta <- cbind(0,new.zeta); }

		## Replace vameans and lambda with elements in x
		if(num.lv > 0) { 
			new.vameans[sel.i,] <- x[1:num.lv]; va.x <- x[-(1:num.lv)]
			if(Lambda.struc == "diagonal") { new.lambda[sel.i,] <- va.x[1:(num.lv)]; va.x <- va.x[-(1:(num.lv))] }
			if(Lambda.struc == "unstructured") { ## Rebuilt lambda array from unique elements
				new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,],diag=T)] <- va.x[1:(num.lv+num.lv*(num.lv-1)/2)]; va.x <- va.x[-(1:(num.lv+num.lv*(num.lv-1)/2))] 
				new.lambda[sel.i,,][upper.tri(new.lambda[sel.i,,])] <- new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,])]
				} 
			}
			
		grad.vameans <- NULL
		grad.lambda <- matrix(0,num.lv,num.lv)
		eta.mat <- matrix(new.beta0,n,p,byrow=T) + new.vameans %*% t(new.theta) + calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
		if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
		if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
		
		if(family=="poisson") { 
			sum1 <- (y[sel.i,]-exp(eta.mat[sel.i,]))
			for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(sum1*new.theta[,l])-new.vameans[sel.i,l]) } }

		if(family=="negative.binomial") {        
			eta.mat <- matrix(new.beta0,n,p,byrow=T) + new.vameans %*% t(new.theta) - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
			if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
			if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
	
			sum1 <- -new.phi - (y[sel.i,]+new.phi)*(-new.phi/(exp(eta.mat[sel.i,])+new.phi))
			for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(sum1*new.theta[,l])-new.vameans[sel.i,l]) } }
			
		if(family=="binomial") { 
			eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
			sum1 <- dnorm(eta.mat[sel.i,])*(y[sel.i,]-trial.size*pnorm(eta.mat[sel.i,]))/(pnorm(eta.mat[sel.i,])*(1-pnorm(eta.mat[sel.i,]))+1e-10)
			for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(sum1*new.theta[,l])-new.vameans[sel.i,l]) } }

		if(family=="ordinal") { 
			eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
			deriv.trunnorm <- matrix(NA,n,p)
			for(j in 1:p) { 
				deriv.trunnorm[y[,j]==1,j] <- -dnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])/pnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])
				deriv.trunnorm[y[,j]==max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j])/(1-pnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j]))
				if(max(y[,j])>2) {
					j.levels <- 2:(max(y[,j])-1)
					for(k in j.levels) { deriv.trunnorm[y[,j]==k,j] <- (-dnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])+dnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) / (pnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])-pnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) } }
				}
				deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0	

			for(l in 1:num.lv) { grad.vameans <- c(grad.vameans,sum(deriv.trunnorm[sel.i,]*new.theta[,l])-new.vameans[sel.i,l]) } }

			
		if(family %in% c("poisson","negative.binomial")) {
			if(family == "poisson") { deriv1 <- exp(eta.mat[sel.i,]) }
			if(family == "negative.binomial") {
				eta.mat <- matrix(new.beta0,n,p,byrow=T) + new.vameans %*% t(new.theta) - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
				if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
				if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
				deriv1 <- (y[sel.i,]+new.phi)*(new.phi/(exp(eta.mat[sel.i,])+new.phi)) }
			
			if(Lambda.struc == "unstructured") {	
				theta2 <- sapply(1:p,function(j,theta) theta[j,]%*%t(theta[j,]), theta = new.theta)
				theta2 <- t(theta2)
				grad.lambda <- -0.5*matrix(apply(deriv1*theta2,2,sum),nrow=num.lv) + 0.5*solve(new.lambda[sel.i,,]) - 0.5*diag(nrow=num.lv) }
			if(Lambda.struc == "diagonal") { 
				theta2 <- new.theta^2
				grad.lambda <- -0.5*diag(apply(deriv1*theta2,2,sum),nrow=num.lv) + 0.5*solve(diag(x=new.lambda[sel.i,],nrow=num.lv)) - 0.5*diag(nrow=num.lv) }
			}
		
		if(family %in% c("binomial","ordinal")) { 
			if(Lambda.struc == "unstructured") {	
				theta2 <- sapply(1:p,function(j,theta) theta[j,]%*%t(theta[j,]), theta = new.theta)
				theta2 <- t(theta2)
				grad.lambda <- -0.5*matrix(apply(theta2,2,sum),nrow=num.lv) + 0.5*solve(new.lambda[sel.i,,]) - 0.5*diag(nrow=num.lv) }
			if(Lambda.struc == "diagonal") { 
				theta2 <- new.theta^2
				grad.lambda <- -0.5*diag(apply(theta2,2,sum),nrow=num.lv) + 0.5*solve(diag(x=new.lambda[sel.i,],nrow=num.lv)) - 0.5*diag(nrow=num.lv) }
			}
			
		grad.lambda2 <- NULL
		if(Lambda.struc == "diagonal") grad.lambda2 <- diag(x=as.matrix(grad.lambda)) ## Only extract diagonal elements
		if(Lambda.struc == "unstructured") grad.lambda2 <- grad.lambda[lower.tri(grad.lambda,diag=T)]
				
		return(c(grad.vameans,grad.lambda2)) }

		
	## cross derivatives of VA logL - taking VA parameters vameans_sel.i and lambda_sel.i and diffing wrt to model parameters
	grad.all.cross <- function(x,vameans=NULL,lambda=NULL,mod.x,sel.i) {
		x2 <- mod.x
		if(num.lv > 0) { 
			new.vameans <- vameans; new.lambda <- lambda
			new.theta <- matrix(c(x2[1:(p*num.lv)]),p,num.lv); x2 <- x2[-(1:(p*num.lv))] }
		new.beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
		new.beta <- NULL; if(!is.null(X)) { new.beta <- matrix(x2[1:(p*num.X)],p,num.X); x2 <- x2[-(1:(p*num.X))] }
		new.row.params <- NULL; if(row.eff) { new.row.params <- x2[1:n]; x2 <- x2[-(1:n)] }
		new.phi <- NULL; if(family == "negative.binomial") { new.phi <- x2[1:p]; x2 <- x2[-(1:p)] }
		new.zeta <- NULL; if(family == "ordinal") { 
			new.zeta <- matrix(NA,p,max.levels-2); 
			for(j in 1:p) { if(max(y[,j]) > 2) { new.zeta[j,1:(max(y[,j])-2)] <- x2[1:(max(y[,j])-2)]; x2 <- x2[-(1:(max(y[,j])-2))] } }
			new.zeta <- cbind(0,new.zeta); }

		## Replace vameans and lambda with elements in x
		if(num.lv > 0) { 
			new.vameans[sel.i,] <- x[1:num.lv]; va.x <- x[-(1:num.lv)]
			if(Lambda.struc == "diagonal") { new.lambda[sel.i,] <- va.x[1:(num.lv)]; va.x <- va.x[-(1:(num.lv))] }
			if(Lambda.struc == "unstructured") { ## Rebuilt lambda array from unique elements
				new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,],diag=T)] <- va.x[1:(num.lv+num.lv*(num.lv-1)/2)]; va.x <- va.x[-(1:(num.lv+num.lv*(num.lv-1)/2))] 
				new.lambda[sel.i,,][upper.tri(new.lambda[sel.i,,])] <- new.lambda[sel.i,,][lower.tri(new.lambda[sel.i,,])]
				} 
			}
			
		eta.mat <- matrix(new.beta0,n,p,byrow=T) 
		if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
		if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
		if(num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.theta) + calc.quad(new.lambda,new.theta,Lambda.struc)$mat 				
		
		grad.theta <- grad.beta0 <- grad.beta <- grad.phi <- grad.row.params <- grad.zeta <- NULL;

		if(row.eff) grad.row.params <- numeric(n)
		if(Lambda.struc == "unstructured" & num.lv > 0) { Lambda.theta <- array(NA,dim=c(p,n,num.lv)); for(i2 in 1:n) { Lambda.theta[,i2,] <- new.theta%*%new.lambda[i2,,] } }
			
		if(family=="poisson") { 
			grad.beta0 <- colSums(y-exp(eta.mat))

			if(num.lv > 0) {
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep(t(Lambda.theta[,,l]),1,new.vameans[,l],"+")*exp(eta.mat)
					if(Lambda.struc == "diagonal") sum1 <- sweep(y,1,new.vameans[,l],"*") - sweep((new.lambda[,l]%*%t(new.theta[,l])),1,new.vameans[,l],"+")*exp(eta.mat)
					grad.theta <- c(grad.theta,colSums(sum1)) } }

			if(!is.null(X)) { for (l in 1:num.X) { sum1 <- sweep(y-exp(eta.mat),1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1)) } } 
			if(row.eff) { grad.row.params <- rowSums(y-exp(eta.mat)) } }

		if(family=="negative.binomial") { 
			phi.mat <- matrix(new.phi,n,p,byrow=T)
			eta.mat <- matrix(new.beta0,n,p,byrow=T) 
			if(!is.null(X)) eta.mat <- eta.mat + X %*% t(new.beta)
			if(row.eff) eta.mat <- eta.mat + matrix(new.row.params,n,p)
			if(num.lv > 0) mu.mat <- eta.mat + new.vameans %*% t(new.theta) - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
			if(num.lv > 0) mu.mat.noquad <- eta.mat + new.vameans %*% t(new.theta) 
			if(num.lv == 0) mu.mat <- eta.mat  

			grad.beta0 <- colSums(-phi.mat - (y+phi.mat)*(-phi.mat/(exp(mu.mat)+phi.mat)))

			if(num.lv > 0) {
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep(t(Lambda.theta[,,l]),1,-new.vameans[,l],"+")*(-(y+phi.mat)*(phi.mat/(exp(mu.mat)+phi.mat)))
					if(Lambda.struc == "diagonal") sum1 <- sweep(-phi.mat,1,new.vameans[,l],"*") + sweep((new.lambda[,l] %*% t(new.theta[,l])),1,-new.vameans[,l],"+")*(-(y+phi.mat)*(phi.mat/(exp(mu.mat)+phi.mat)))
					grad.theta <- c(grad.theta,colSums(sum1)) } }
 
			if(!is.null(X)) {          
				for(l in 1:num.X) { sum1 <- sweep(-phi.mat-(y+phi.mat)*(-phi.mat/(exp(mu.mat)+phi.mat)),1,X[,l],"*")
					grad.beta <- c(grad.beta,colSums(sum1)) } }	 
			if(row.eff) { grad.row.params <- rowSums(-phi.mat-(y+phi.mat)*(-phi.mat/(exp(mu.mat)+phi.mat))) }  

			grad.phi <- colSums(-mu.mat.noquad + 1 + log(phi.mat) - digamma(phi.mat) - log(1+phi.mat/exp(mu.mat)) - (y+phi.mat)/(phi.mat+exp(mu.mat)) + digamma(y+phi.mat)) }

			
		if(family=="binomial") { 
			if(num.lv > 0) eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 
				
			grad.beta0 <- colSums(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),na.rm=T)

			if(num.lv > 0) {
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
					if(Lambda.struc == "diagonal") sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,new.vameans[,l],"*") - (new.lambda[,l] %*% t(new.theta[,l]))
					grad.theta <- c(grad.theta,colSums(sum1,na.rm=T)) } }
					
			if(!is.null(X)) { 
				for(l in 1:num.X) { 
					sum1 <- sweep(dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10),1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1,na.rm=T)) } } 
			if(row.eff) { grad.row.params <- rowSums((dnorm(eta.mat)*(y-trial.size*pnorm(eta.mat))/(pnorm(eta.mat)*(1-pnorm(eta.mat))+1e-10)),na.rm=T) } }
		
		if(family=="ordinal") { 
			grad.zeta <- matrix(0,p,ncol(new.zeta))
			if(num.lv > 0) eta.mat <- eta.mat - calc.quad(new.lambda,new.theta,Lambda.struc)$mat 

			deriv.trunnorm <- matrix(0,n,p)
			for(j in 1:p) { 
				deriv.trunnorm[y[,j]==1,j] <- -dnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])/pnorm(new.zeta[j,1]-eta.mat[y[,j]==1,j])
				deriv.trunnorm[y[,j]==max(y[,j]),j] <- dnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j])/(1-pnorm(new.zeta[j,max(y[,j])-1]-eta.mat[y[,j]==max(y[,j]),j]))
				j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
				for(k in j.levels) {  
					deriv.trunnorm[y[,j]==k,j] <- (-dnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])+dnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) / (pnorm(new.zeta[j,k]-eta.mat[y[,j]==k,j])-pnorm(new.zeta[j,k-1]-eta.mat[y[,j]==k,j])) } 
				}
			deriv.trunnorm[!is.finite(deriv.trunnorm)] <- 0	

			grad.beta0 <- colSums(deriv.trunnorm,na.rm=T)
						
			if(num.lv > 0) { 
				for(l in 1:num.lv) {
					if(Lambda.struc == "unstructured") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - t(Lambda.theta[,,l])
					if(Lambda.struc == "diagonal") sum1 <- sweep(deriv.trunnorm,1,new.vameans[,l],"*") - (new.lambda[,l]%*%t(new.theta[,l]))
					grad.theta <- c(grad.theta,colSums(sum1,na.rm=T)) }
				}

			if(!is.null(X)) { for(l in 1:num.X) { sum1 <- sweep(deriv.trunnorm,1,X[,l],"*"); grad.beta <- c(grad.beta,colSums(sum1,na.rm=T)) } } 
			if(row.eff) { grad.row.params <- rowSums(deriv.trunnorm,na.rm=T) }  
			
			for(j in 1:p) { 
				zeta0 <- new.zeta[j,]; 
				grad.zeta[j,max(y[,j])-1] <- grad.zeta[j,max(y[,j])-1] - sum(dnorm(zeta0[max(y[,j])-1]-eta.mat[which(y[,j] == max(y[,j])) ,j])/(1-pnorm(zeta0[max(y[,j])-1]-eta.mat[which(y[,j] == max(y[,j])),j]))) 
				j.levels <- (1:max(y[,j]))[-c(1,max(y[,j]))]
				for(k in j.levels) {
					grad.zeta[j,k] <- grad.zeta[j,k] + sum(dnorm(zeta0[k]-eta.mat[y[,j]==k,j])/(pnorm(zeta0[k]-eta.mat[y[,j]==k,j])-pnorm(zeta0[k-1]-eta.mat[y[,j]==k,j]))) 
					grad.zeta[j,k-1] <- grad.zeta[j,k-1] - sum(dnorm(zeta0[k-1]-eta.mat[y[,j]==k,j])/(pnorm(zeta0[k]-eta.mat[y[,j]==k,j])-pnorm(zeta0[k-1]-eta.mat[y[,j]==k,j]))) } 
				}
				
			grad.zeta <- grad.zeta[,-1]		
			if(length(unique(apply(y,2,max))) > 1) { grad.zeta <- c(grad.zeta[-which(grad.zeta == 0)]) }
			}
				
		return(c(grad.theta, grad.beta0, grad.beta, grad.row.params, grad.phi, grad.zeta))
		}

		
	zero.cons <- which(theta == 0)
	if (family == "ordinal") {
	    zeta2 <- as.vector(t(zeta[, -1]))
	    find.na.zeta <- which(is.na(zeta2))
	    if (length(find.na.zeta) > 0) zeta2 <- zeta2[-find.na.zeta]
	}
	if (family != "ordinal") {
	    zeta2 <- NULL
	}
	
	A.mat <- -nd2(x0 = c(theta, beta0, beta, row.params, phi, zeta2), f = grad.all.mod, vameans = vameans, lambda = lambda)
	if (length(zero.cons) > 0) {
	    A.mat <- A.mat[-zero.cons, ]
	    A.mat <- A.mat[, -zero.cons]
	}
	dim(A.mat)
 	
 	
 	if(num.lv > 0) {
        D.mat.fun <- function(i3) {
            if (Lambda.struc == "unstructured") {
                lambda.i <- lambda[i3, , ][lower.tri(lambda[i3, , ], diag = T)]
            }
            if (Lambda.struc == "diagonal") {
                lambda.i <- lambda[i3, ]
            }
        
            return(-nd2(f = grad.all.var, x0 = c(vameans[i3, ], lambda.i), vameans = vameans, lambda = lambda, mod.x = c(theta, beta0, beta, row.params, phi, zeta2), sel.i = i3))
        }
			
		unique.ind <- which(!duplicated(y))
		D.mat <- vector("list", n)
		for (k in 1:length(unique.ind)) {
		    subD <- D.mat.fun(i3 = unique.ind[k])
		    match.seq <- which(apply(y, 1, identical, y[unique.ind[k], ]) == 1) ## Unfortunately replace() only works if you want to replace elements with numerics
		    for (k2 in match.seq) {
		        D.mat[[k2]] <- subD
		    }
		}
		# D.mat <- lapply(1:n,D.mat.fun)
		D.mat <- as.matrix(bdiag(D.mat))
		rm(subD)

		B.mat <- vector("list", n)
		for (i3 in 1:length(unique.ind)) {
		    # cat("Onto row",i3,"\n")
		    if (Lambda.struc == "unstructured") {
		        lambda.i <- lambda[unique.ind[i3], , ][lower.tri(lambda[unique.ind[i3], , ], diag = T)]
		    }
		    if (Lambda.struc == "diagonal") {
		        lambda.i <- lambda[unique.ind[i3], ]
		    }
		
		    # grad.all.var(x=c(vameans[i3,],lambda.i), vameans = vameans, lambda = lambda, mod.x = c(theta,beta0,beta,row.params,phi,zeta2), sel.i = i3)
		    subB.mat <- -nd2(f = grad.all.cross, x0 = c(vameans[unique.ind[i3], ], lambda.i), vameans = vameans, lambda = lambda, mod.x = c(theta, beta0, beta, row.params, phi, zeta2), sel.i = unique.ind[i3])
		    if (length(zero.cons) > 0) {
		        subB.mat <- subB.mat[-zero.cons, ]
		    }
		    match.seq <- which(apply(y, 1, identical, y[unique.ind[i3], ]) == 1)
		    for (k2 in match.seq) {
		        B.mat[[k2]] <- subB.mat
		    }
		}
		B.mat <- do.call(cbind, B.mat)
		
		cov.mat.mod <- ginv(A.mat - B.mat %*% solve(D.mat) %*% t(B.mat))
	}
		
	if(num.lv == 0) cov.mat.mod <- ginv(A.mat)
		
		
	#cov.mat.mod[abs(cov.mat.mod) < 1e-5] <- 0
	sd.errors <- sqrt(diag(abs(cov.mat.mod))) 

	
	theta.errors <- NULL; 
	if(num.lv > 0) { 
		theta.vec <- sd.errors[1:(p*num.lv-length(zero.cons))]; theta.errors <- matrix(0,p,num.lv); 
		if(length(zero.cons) > 0) theta.errors[-zero.cons] <- theta.vec
		if(length(zero.cons) == 0) theta.errors <- matrix(theta.vec,p,num.lv)
		sd.errors <- sd.errors[-(1:(p*num.lv-length(zero.cons)))] 
		if(num.lv > 1) theta.errors[upper.tri(theta.errors)] <- 0 
		rownames(theta.errors) <- colnames(y); 
		colnames(theta.errors) <- 1:num.lv
		}
	beta.errors <- sd.errors[1:p]; sd.errors <- sd.errors[-(1:p)]
	names(beta.errors) <- colnames(y)

	if(!is.null(X)) { 
		beta.errors <- cbind(beta.errors,matrix(sd.errors[1:(p*ncol(X))],p,ncol(X))); 
		sd.errors <- sd.errors[-(1:(p*ncol(X)))] 
		}

	row.params.errors <- NULL; 
	if(row.eff) { row.params.errors <- sd.errors[1:n]; sd.errors <- sd.errors[-(1:n)]; names(row.params.errors) <- rownames(y) }
	
	phi.errors <- NULL; 
	if(family == "negative.binomial") { phi.errors <- sd.errors[1:p]; sd.errors <- sd.errors[-(1:p)]; names(phi.errors) <- colnames(y) }
	
	zeta.errors <- NULL;
	if(family == "ordinal") { 
		zeta.errors <- matrix(NA,p,max.levels-2); 
		for(j in 1:p) { if(max(y[,j]) > 2) { zeta.errors[j,1:(max(y[,j])-2)] <- sd.errors[1:(max(y[,j])-2)]; sd.errors <- sd.errors[-(1:(max(y[,j])-2))] } }
		zeta.errors <- cbind(0,zeta.errors) 
		rownames(zeta.errors) <- colnames(y); colnames(zeta.errors) <- paste(1:(max(y)-1),"|",2:max(y),sep="")
		}

		
	return(list(beta = beta.errors, row.params = row.params.errors, theta = theta.errors, phi = phi.errors, zeta = zeta.errors)) 
	}
	

	
## A function to compute highly accurate first-order derivatives
## From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)
nd2 <- function(x0, f, m = NULL, D.accur = 2, ...) {
    D.n <- length(x0)
    if (is.null(m)) {
        D.f0 <- f(x0, ...)
        m <- length(D.f0)
    }
    if (D.accur == 2) {
        D.w <- tcrossprod(rep(1, m), c(-1 / 2, 1 / 2))
        D.co <- c(-1, 1)
    } else {
        D.w <- tcrossprod(rep(1, m), c(1 / 12, -2 / 3, 2 / 3, -1 / 12))
        D.co <- c(-2, -1, 1, 2)
    }
    D.n.c <- length(D.co)
    macheps <- .Machine$double.eps
    D.h <- macheps^(1 / 3) * abs(x0)
    D.deriv <- matrix(NA, nrow = m, ncol = D.n)

    for (ii in 1:D.n) {
        D.temp.f <- matrix(0, m, D.n.c)
        for (jj in 1:D.n.c) {
            D.xd <- x0 + D.h[ii] * D.co[jj] * (1:D.n == ii)
            D.temp.f[, jj] <- f(D.xd, ...)
        }
        D.deriv[, ii] <- rowSums(D.w * D.temp.f) / D.h[ii]
    }
    return(D.deriv)
}


vacov.convert <- function(vacov, type = 1) {
    if (type == 1) { ## Convert unstructured to diagonal
        n <- dim(vacov)[1]
        vacov2 <- matrix(1, dim(vacov)[1], dim(vacov)[3])
        for (i in 1:n) {
            vacov2[i, ] <- diag(vacov[i, , ])
        }
        return(vacov2)
    }
    if (type == 2) { ## Convert diagonal to unstructured
        n <- nrow(vacov)
        vacov2 <- array(0, dim = c(nrow(vacov), ncol(vacov), ncol(vacov)))
        for (i in 1:n) {
            diag(vacov2[i, , ]) <- vacov[i, ]
        }
        return(vacov2)
    }
}