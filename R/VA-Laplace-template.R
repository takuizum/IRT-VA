###################################################
## Fitting GLLVMs using the Laplace approximation (Huber et al., 2004)
## Please note the information matrix is currently only produced for some the model parameters, i.e. no standard errors are produced for the predictions of the LVs, betas, and phis
## maxit is not the maximum number of total iterations, but rather the max number of sub-iterations for each Quasi-Newton method

## NOTE: For binary data the Laplace approximation can be a somewhat unstable and "bobble" without reaching convergence. We recommend you set maxit = 1 for binary data
#######################
library(mvabund)
library(mvtnorm)
library(graphics)
# library(numDeriv)
library(MASS)

glvm.laplace <- function(y, X = NULL, num.lv = 2, family = "poisson", row.eff = FALSE, max.iter = 100, eps = 1e-4, trace = FALSE, seed = NULL, plot = FALSE, maxit = 100, info = FALSE) {
	n <- dim(y)[1]; p <- dim(y)[2]	
     num.lv <- num.lv
     if(!is.null(colnames(y))) colnames(y) <- 1:p
	if(is.null(rownames(y))) rownames(y) <- paste("Row",1:n,sep="")
	if(is.null(colnames(y))) colnames(y) <- paste("Col",1:p,sep="")
	if(!is.null(X)) { X <- as.matrix(X); if(is.null(colnames(X))) colnames(X) <- paste("x",1:ncol(X),sep="") }

	
     ## Set initial values for model parameters (including dispersion prm) and latent variables
	if(!is.null(seed)) set.seed(seed)
	fit <- start.values(y, X, family,num.lv)
  
	beta0 <- fit$params[,1] ## column intercepts
	betas=NULL; if(!is.null(X))betas <- c(fit$X.params) ## covariates coefficients
	lambdas=NULL; if(num.lv>0) lambdas <- c(fit$params[,2:(num.lv+1)]) ## LV coefficients
	row.params <- NULL; if(row.eff) row.params <- rnorm(n) ## row parameters
	phis=NULL; if(family=="negative.binomial"){ phis <- runif(p) }  ## dispersion params

	new.params <- params <- c(row.params, beta0, lambdas)
	lvs=NULL; if(num.lv>0) lvs <- matrix(fit$index,ncol=num.lv); ## LVs

     
	current.loglik <- -1e6; iter <- 1; err <- 10; 
	while((err > (1 + eps) || err < (1 - eps)) && iter <= max.iter) {
		## LA-likelihood
		ll <- function(params, lvs, y, phis, betas = NULL, type = family) {
			if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
			if(!row.eff) { row.params <- 0; x2 <- params }
			beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
			if(num.lv>0) { 
				lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
				LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
				if(num.lv==1) LambdaLambdaT <- t(LambdaLambdaT) }

			eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params 
			if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat)
			if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas,nrow=p))
			I.numlv <- diag(x=num.lv)
			SGam <- 0
				
			switch(type,
				bernoulli = {
				if(num.lv>0) {
					V <- binomial()$variance(binomial()$linkinv(eta.mat))			
					G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
					if(num.lv==1) SGam <- SGam - sum(0.5*log(G));
					if(num.lv==2) SGam <- SGam - sum(0.5*log(G[,1]*G[,4]-G[,2]*G[,3]))
					if(num.lv>2) { for(i in 1:n) {
						Ga <- matrix(G[i,],nrow=num.lv, ncol=num.lv)
						SGam <- SGam - 0.5 * log(det(Ga)) } 
						}
					}
				SGam <- SGam + sum(dbinom(y, 1, prob=binomial()$linkinv(eta.mat), log = T))
				}, 
				poisson = {
				if(num.lv>0) {
					V <- poisson()$variance(poisson()$linkinv(eta.mat))		
					G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
					if(num.lv==1) SGam <- SGam - sum(0.5 * log(G));
					if(num.lv==2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
					if(num.lv>2){ for(i in 1:n) {
						Ga <- matrix(G[i,],nrow=num.lv, ncol=num.lv)
						SGam <- SGam - 0.5 * log(det(Ga)) } 
						}
					}
				SGam <- SGam + sum(dpois(y, lambda = exp(eta.mat), log = T))
				},
				negative.binomial = {
				if(num.lv>0) {
					V <- exp(eta.mat) * (1 + y * matrix(phis,n,p,byrow=TRUE)) / (1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))^2
					G<-t(t(V %*% LambdaLambdaT) + c(I.numlv))
					if(num.lv==1) SGam <- SGam - sum(0.5 * log(G));
					if(num.lv==2) SGam <- SGam - sum(0.5 * log(G[,1] * G[,4] - G[,2] * G[,3]))
					if(num.lv>2) { for(i in 1:n) {
						Ga <- matrix(G[i,],nrow=num.lv, ncol=num.lv)
						SGam <- SGam - 0.5 * log(det(Ga)) } 
						}
					}
				SGam <- SGam + sum(dnbinom(y, mu = exp(eta.mat), size = matrix(1/phis,n,p,byrow=TRUE), log = T))
				})
				
			if(num.lv>0) SGam <- SGam - 0.5 * sum(diag(lvs %*% t(lvs))) 
			return(SGam)
			}

		
		## Return n x num.lv^2 matrix with each row as inverse of Gamma(\theta,z_i) matrix 
		## as in eq. 9 in Huber et al. 
		Hub.Gamma <- function(y, eta, lambdas.mat, phis) {
			LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
			I.numlv <- diag(x=num.lv)
			if(num.lv==1) LambdaLambdaT <- t(LambdaLambdaT) 
			if(family == "bernoulli")  V <- binomial()$variance(binomial()$linkinv(eta)) 
			if(family == "poisson")  V <- poisson()$variance(poisson()$linkinv(eta)) 
			if(family == "negative.binomial")  V <- exp(eta) * (1 + y * matrix(phis,n,p,byrow=TRUE)) / (1 + matrix(phis,n,p,byrow=TRUE) * exp(eta))^2 

			s <- matrix(0,n,num.lv^2)
			for(i in 1:n) s[i,] <- c(solve(matrix(colSums(V[i,] * LambdaLambdaT), num.lv, num.lv) + I.numlv))
			return(s) 
			}	

			
		## Gradients of parameters row.params, beta0 and lambdas
		## Should simplify code by using switch
		param.gradient <- function(params, lvs, y, phis, betas = NULL, type = family) {
			if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
			if(!row.eff) { row.params <- 0; x2 <- params }
			beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
			if(num.lv>0) { 
				lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
				LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
				if(num.lv==1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) } }
			
			eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params 
			if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat) 
			if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas,nrow=p))

			grad.beta0 <- numeric(p); 
			grad.lambdas <- NULL; if(num.lv>0) grad.lambdas <- matrix(0, p, num.lv) 
			grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n) 
			
			if(num.lv>0) { 
				GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis) 
				Kron <- list()
				es<-diag(num.lv)
				for(s in 1:num.lv) {
					Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x))
					}
				if(num.lv==1) Kron[[1]] <- matrix(Kron[[1]],nrow=1)
				}
				
			switch(type,
				bernoulli = {
				V <- (exp(eta.mat) * (1 - exp(eta.mat))/(1 + exp(eta.mat))^3) ## third derivative
				V2 <- binomial()$variance(binomial()$linkinv(eta.mat))
        
				if(row.eff) grad.row.params <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,1,sum) + rowSums(y - binomial()$linkinv(eta.mat))
				grad.beta0 <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,2,sum) + colSums(y - binomial()$linkinv(eta.mat)) 
				if(num.lv>0){ 
					for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * apply((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V2),2,sum) + colSums((y - binomial()$linkinv(eta.mat)) * lvs[,s]) }
					}  
			         },				
				poisson = {
					V <- exp(eta.mat) ## third derivative and variance
					if(row.eff) grad.row.params <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,1,sum) + rowSums(y - poisson()$linkinv(eta.mat))
					grad.beta0 <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,2,sum) + colSums(y - poisson()$linkinv(eta.mat)) 
					if(num.lv>0) {
						for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * apply((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V),2,sum) + colSums((y - poisson()$linkinv(eta.mat)) * lvs[,s]) }
						} 
					},
				negative.binomial = {
					V <- exp(eta.mat) * (1 + y * matrix(phis,n,p,byrow=TRUE)) * (1 - matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))/(1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))^3 ## third derivative
					V2 <- exp(eta.mat) * (1 + y * matrix(phis,n,p,byrow=TRUE))/(1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))^2 ## variance
					yminuseta <- (y - exp(eta.mat)) / (1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))
		
					if(row.eff) grad.row.params <-  -0.5 * apply(GA %*% t(LambdaLambdaT) * V,1,sum) + rowSums(yminuseta)
					grad.beta0 <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,2,sum) + colSums(yminuseta)
					if(num.lv>0) { 
						for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * apply((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V2),2,sum) + colSums(yminuseta * lvs[,s]) }
						}
					})

			if(num.lv>0) grad.lambdas[upper.tri(grad.lambdas)] <- 0	
			score.out <- c(grad.row.params, grad.beta0, grad.lambdas)
			return(score.out)
			}

		update.params <- optim(params, fn = ll, gr = param.gradient, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit = maxit), lvs = lvs, y = y, phis = phis, betas = betas, type = family) 
		new.params <- update.params$par

     
		## Update covariate parameters (betas) if appropriate 
		if(is.null(X)) betas <- NULL
		if(!is.null(X)) {
			betas.gradient <- function(params, lvs, y, phis, betas, type = family) {
				 if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
				 if(!row.eff) { row.params <- 0; x2 <- params }
				 beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
				 if(num.lv>0) {
					lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
					LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
					if(num.lv==1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) } 
					}
				
				 eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params 
				 if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat) 
				 eta.mat <- eta.mat + X %*% t(matrix(betas,nrow=p))
				
				 grad.betas <- matrix(0, p, ncol(X))
				
				 if(num.lv>0) { GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis) }  
				
 				 switch(type,
					bernoulli = {
						V <- (exp(eta.mat) * (1 - exp(eta.mat)) / (1 + exp(eta.mat))^3) ## third derivative
						grad.betas <- -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X + t(y - binomial()$linkinv(eta.mat)) %*% X 
						}, 
						
					poisson = {
						V <- exp(eta.mat) ## third derivative and variance
						grad.betas <- -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X  + t(y - poisson()$linkinv(eta.mat)) %*% X 
						},

					negative.binomial = {
						V <- exp(eta.mat) * (1 + y * matrix(phis,n,p,byrow=TRUE)) * (1 - matrix(phis,n,p,byrow=TRUE) * exp(eta.mat)) / (1 + matrix(phis,n,p,byrow=T) * exp(eta.mat))^3 ## third derivative
						yminuseta <- (y - exp(eta.mat)) / (1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))
						grad.betas <- -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X + t(yminuseta) %*% X
						})

				score.out <- c(grad.betas)
				return(score.out)
				}

			update.betas <- optim(betas, fn = ll, gr = betas.gradient, method = "BFGS", control = list(trace = 0, fnscale = -1, maxit=maxit), params = new.params, lvs = lvs, y = y, phis = phis, type = family) 
			betas <- update.betas$par
			}
		
		
		## Update dispersion params for NB distribution
		if(family == "negative.binomial") {
			## LA-likelihood
			ll.nb <- function(x, params, lvs, y, betas = NULL) {
				if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
				if(!row.eff) { row.params <- 0; x2 <- params }
				beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
				if(num.lv>0) {
					lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
					LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x%*%t(x) ))
					if(num.lv==1) {LambdaLambdaT=t(LambdaLambdaT); lvs <- matrix(lvs) }
					}
				phis <- exp(x)

				eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params 
				if(num.lv>0) eta.mat <- eta.mat + lvs%*%t(lambdas.mat)
				if(!is.null(X)) eta.mat <- eta.mat + X%*%t(matrix(betas,nrow=p))
				I.numlv <- diag(x=num.lv)
				SGam <- 0
				
				if(num.lv>0){
					V <- exp(eta.mat)*(1 + y*matrix(phis,n,p,byrow=TRUE)) / (1 + matrix(phis,n,p,byrow=TRUE)*exp(eta.mat))^2
					G<-t(t(V%*%LambdaLambdaT) + c(I.numlv))
					if(num.lv==1) SGam <- SGam - sum(0.5*log(G));
					if(num.lv==2) SGam <- SGam - sum(0.5*log(G[,1]*G[,4]-G[,2]*G[,3]));
					if(num.lv>2){ for(i in 1:n) { 
						Ga <- matrix(G[i,],nrow=num.lv, ncol=num.lv); 
						SGam <- SGam - 0.5*log(det(Ga))} 
						}
					}
				SGam <- SGam + sum(dnbinom(y, mu= exp(eta.mat), size=matrix(1/phis,n,p,byrow=TRUE), log = T))
					
				if(num.lv>0) SGam <- SGam - 0.5*sum(diag(lvs%*%t(lvs))) ## Clever!
				return(SGam)
				}

 			Lphi.gradient <- function(x, params, lvs, y, betas) {
				if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
				if(!row.eff) { row.params <- 0; x2 <- params }
				beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
				if(num.lv>0) {
					lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
					LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x) ))
					if(num.lv==1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
					}
 			  phis <- exp(x)
 			  
				eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params 
				if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat) 
				if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas,nrow=p))
				# 
				V <- exp(eta.mat)
				yf <- y + matrix(1 / phis,n,p,byrow=TRUE)
				V2 <- V *matrix(phis,n,p,byrow=TRUE)* (y/(1+matrix(phis,n,p,byrow=TRUE)*V)^2 - 2*V*(1+y*matrix(phis,n,p,byrow=TRUE))/(1+matrix(phis,n,p,byrow=TRUE)*V)^3)
				V3 <- matrix(1 / phis,n,p,byrow=TRUE) * (log(1 + matrix(phis,n,p,byrow=TRUE) * V) - digamma(yf) + digamma(matrix(1 / phis,n,p,byrow=TRUE)))
				if(num.lv>0) { GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis) }   
				
				grad.phis <- colSums(-0.5*GA%*%t(LambdaLambdaT)*V2 + V3 - matrix(phis,n,p,byrow=TRUE)*V*yf/(1+matrix(phis,n,p,byrow=TRUE)*V) + y)
				return(grad.phis)
				}

 			update.phis <- optim(par = log(phis), fn = ll.nb, gr = Lphi.gradient, method = "BFGS", control = list(trace = 0, fnscale = -1), lvs = lvs, y = y, params = new.params, betas = betas)
			phis <- exp(update.phis$par)
			}
     
     
		if(row.eff) { row.params <- new.params[1:n]; x2 <- new.params[-(1:n)] }
		if(!row.eff) { row.params <- NULL; x2 <- new.params }
		beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
		lambdas.mat=NULL; if(num.lv>0) lambdas.mat <- matrix(x2,ncol=num.lv);
		
		if(num.lv>0) {
			## Update latent variables lvs
			ll.i <- function(lvs.i, params, y, phis, betas = NULL, i) {
				if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
				if(!row.eff) { row.params <- rep(0,n); x2 <- params }
				beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
				lambdas.mat <- matrix(x2, ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0

				eta.i <- beta0 + row.params[i] + lambdas.mat%*%lvs.i 
				if(!is.null(betas)) eta.i <- eta.i + matrix(betas,nrow=p)%*%X[i,]
				
				if(family == "binomial") { SGam <- sum(dbinom(y[i,], 1, prob=binomial()$linkinv(eta.i), log = T)) }
				if(family == "poisson") { SGam <- sum(dpois(y[i,], lambda = exp(eta.i), log = T)) } 
				if(family == "negative.binomial") { SGam <- sum(dnbinom(y[i,], mu= exp(eta.i), size=1/phis, log = T)) }
						
				SGam <- SGam - 0.5*t(lvs.i)%*%lvs.i
				return(SGam)
				}
					
			index <- matrix(0,n,num.lv); 
			for(i in 1:n) {
				update.lvs <- optim(lvs[i,], fn = ll.i, gr = NULL, method = "BFGS", control = list(trace = 0, fnscale = -1,maxit=maxit), y = y, params = new.params, phis = phis, betas = betas, i = i)
				index[i,] <- update.lvs$par
				}

				
			## Plot ordinations
			if(plot) {
				plot(rbind(index), type="n", main="Ordination points", xlab = "LV1", ylab = "LV2"); 
				text(index, labels = 1:n) ## new
				}

			lvs <- index
			}

		params <- new.params

		new.loglik <- ll(params = new.params, lvs = lvs, y = y, phis = phis, betas = betas)
		err <- abs(new.loglik/current.loglik); 
		if(trace) { cat("Iterations #", iter, "\n"); cat("New Loglik:", new.loglik, "Current Loglik:", current.loglik, "Ratio", err,"\n") }
		current.loglik <- new.loglik; 
		iter <- iter + 1
		}
		
 
	if(!is.null(X)) betas <- matrix(betas,ncol=ncol(X))
	out <- list(y = y, X = X, num.lv = num.lv, row.eff = row.eff, family = family, beta0 = beta0, beta = betas, iter = iter, logLik = current.loglik, lambda = lambdas.mat, lvs = lvs)
	names(out$beta0)  <- colnames(y); 
	if(!is.null(X)) { rownames(out$betas) <- colnames(y); colnames(out$betas) <- colnames(X); }
	if(family == "negative.binomial") { out$phi = phis; names(out$phi) <- colnames(y) }
	if(num.lv>0) {  
		colnames(out$lambda) <- colnames(out$lvs) <- paste("LV",1:num.lv,sep=""); 
		rownames(out$lambda) <- colnames(y);
		rownames(out$lvs) <- 1:n }
	if(row.eff) names(out$row.params) <- rownames(y)
     
     if(info) {
		cat("Calcuating information matrix...\n")
		
		## Gradients of all model parameters in order: row.params, beta0, lambdas, betas, phis
		full.param.gradient <- function(params, lvs, y, X = NULL, type = family, phis) {
			if(row.eff) { row.params <- params[1:n]; x2 <- params[-(1:n)] }
			if(!row.eff) { row.params <- 0; x2 <- params }
			beta0 <- x2[1:p]; x2 <- x2[-(1:p)]
			if(num.lv>0) { 
				lambdas.mat <- matrix(x2[1:(p*num.lv)], ncol=num.lv); lambdas.mat[upper.tri(lambdas.mat)] <- 0
				LambdaLambdaT <- t(apply(lambdas.mat, 1, function(x) x %*% t(x)))
				if(num.lv==1) { LambdaLambdaT <- t(LambdaLambdaT); lvs <- matrix(lvs) }
				x2 <- x2[-(1:(p*num.lv))] }

			if(!is.null(X)) { betas <- x2[1:(p*ncol(X))]; x2 <- x2[-(1:(p*ncol(X)))] }
# 			if(family == "negative.binomial") { phis <- x2[1:p]; x2 <- x2[-(1:p)] }
			
			eta.mat <- matrix(beta0, n, p, byrow=TRUE) + row.params 
			if(num.lv>0) eta.mat <- eta.mat + lvs %*% t(lambdas.mat) 
			if(!is.null(X)) eta.mat <- eta.mat + X %*% t(matrix(betas,nrow=p))

			grad.beta0 <- numeric(p); 
			grad.lambdas <- NULL; if(num.lv>0) grad.lambdas <- matrix(0, p, num.lv) 
			grad.row.params <- NULL; if(row.eff) grad.row.params <- numeric(n) 
			grad.betas <- NULL; if(!is.null(X)) grad.betas <- matrix(0, p, ncol(X)) 
			grad.phis <- NULL; 
			
			if(num.lv>0) { 
				GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis) 
				Kron <- list()
				es<-diag(num.lv)
				for(s in 1:num.lv) { Kron[[s]]<-apply(lambdas.mat, 1, function(x) kronecker(es[s,],t(x))) + apply(lambdas.mat, 1, function(x) kronecker(t(es[s,]),x)) }
				if(num.lv==1) Kron[[1]] <- matrix(Kron[[1]],nrow=1) }


			switch(type,
				bernoulli = {
					V <- (exp(eta.mat) * (1 - exp(eta.mat))/(1 + exp(eta.mat))^3) ## third derivative
					V2 <- binomial()$variance(binomial()$linkinv(eta.mat))
		
					if(row.eff) grad.row.params <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,1,sum) + rowSums(y - binomial()$linkinv(eta.mat))
					grad.beta0 <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,2,sum) + colSums(y - binomial()$linkinv(eta.mat)) 

					if(num.lv>0) { 
						for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * apply((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V2),2,sum) + colSums((y - binomial()$linkinv(eta.mat)) * lvs[,s]) }
						}
					if(!is.null(X)) { grad.betas <- -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X + t(y - binomial()$linkinv(eta.mat)) %*% X }
					},

				poisson = {
					V <- exp(eta.mat) ## third derivative and variance
					if(row.eff) grad.row.params <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,1,sum) + rowSums(y - poisson()$linkinv(eta.mat))
					grad.beta0 <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,2,sum) + colSums(y - poisson()$linkinv(eta.mat)) 

					if(num.lv>0) {
						for(s in 1:num.lv) { grad.lambdas[,s] <- -0.5 * apply((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V),2,sum) + colSums((y - poisson()$linkinv(eta.mat)) * lvs[,s]) }
						} 
					if(!is.null(X)) { grad.betas <- -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X  + t(y - poisson()$linkinv(eta.mat)) %*% X }
					},

				negative.binomial = {
					V <- exp(eta.mat) * (1 + y * matrix(phis,n,p,byrow=TRUE)) * (1 - matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))/(1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))^3 ## third derivative
					V2 <- exp(eta.mat) * (1 + y * matrix(phis,n,p,byrow=TRUE))/(1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))^2 ## variance
					yminuseta <- (y - exp(eta.mat)) / (1 + matrix(phis,n,p,byrow=TRUE) * exp(eta.mat))
		
					if(row.eff) grad.row.params <-  -0.5 * apply(GA %*% t(LambdaLambdaT) * V,1,sum) + rowSums(yminuseta)
					grad.beta0 <- -0.5 * apply(GA %*% t(LambdaLambdaT) * V,2,sum) + colSums(yminuseta)
					if(num.lv>0) { 
						for(s in 1:num.lv){ grad.lambdas[,s] <- -0.5 * apply((GA %*% t(LambdaLambdaT) * V * matrix(lvs[,s],n,p,byrow=FALSE) + GA %*% Kron[[s]] * V2),2,sum) + colSums(yminuseta * lvs[,s]) }
						}
					if(!is.null(X)) { grad.betas <- -0.5 * t(GA %*% t(LambdaLambdaT) * V) %*% X + t(yminuseta) %*% X }

					V <- exp(eta.mat)
					yf <- y + matrix(1 / phis,n,p,byrow=TRUE)
					V2 <- V * (y/(1+matrix(phis,n,p,byrow=TRUE)*V)^2 - 2*V*(1+y*matrix(phis,n,p,byrow=TRUE))/(1+matrix(phis,n,p,byrow=TRUE)*V)^3)
					V3 <- matrix(1 / phis^2,n,p,byrow=TRUE) * (log(1 + matrix(phis,n,p,byrow=TRUE) * V) - digamma(yf) + digamma(matrix(1 / phis,n,p,byrow=TRUE)))
					if(num.lv>0) { GA <- Hub.Gamma(y = y, eta = eta.mat, lambdas.mat, phis) }   
# 					grad.phis <- colSums(-0.5*GA%*%t(LambdaLambdaT)*V2 + V3 - V*yf/(1+matrix(phis,n,p,byrow=TRUE)*V) + y/matrix(phis,n,p,byrow=TRUE))
					}
				)

			if(num.lv>0) grad.lambdas[upper.tri(grad.lambdas)] <- 0	
			score.out <- c(grad.row.params, grad.beta0, grad.lambdas, grad.betas)#, grad.phis)
			return(score.out)
			}
			
			
		get.info <- try(-nd2.la(x0 = c(row.params, beta0, lambdas.mat, betas), func = full.param.gradient, lvs = lvs, y = y, X = X, type = family, phis = phis), silent = TRUE) 
		if(!inherits(get.info, "try-error")) {
  			find.const.lambdas <- which(is.finite(colSums(get.info)))
  			get.info <- get.info[find.const.lambdas,find.const.lambdas]
			out$info <- get.info
  			out$se <- try(sqrt(diag(ginv(get.info))),silent=T)
			}
		#out$param.vec <- c(row.params, beta0, lambdas.mat, betas, phis)
		}
          
	return(out)
	}

	
## Starting values by fitting marginal glm's
start.values <- function(y, X = NULL, family = "poisson", num.lv = 2) {
	N <- nrow(y); p <- ncol(y)
	y <- as.matrix(y)
	
	unique.ind <- which(!duplicated(y)) ## Get good lvs to start off
	rs <- as.vector(rowSums(y, na.rm = TRUE))
	len.uni <- length(unique(rs))
	rs <- factor(rs, labels = 1:len.uni)
	rs <- as.numeric(levels(rs))[as.integer(rs)]
	index <- matrix(seq(-3, 3, len = len.uni)[rs], ncol=1)
	if(num.lv > 1) { index <- cbind(index,rmvnorm(nrow(index),rep(0,num.lv-1))) }		
	unique.index <- as.matrix(index[unique.ind,]) 
	
	index <- cbind(index,X)
	if(num.lv==0) index <- cbind(index=rep(1,N),X)
	
	fit <- manyglm(y ~ index, family = family)  
	if(num.lv==0) { params <- cbind(fit$coef[1,], rep(0,p))}
	if(num.lv>0) { params <- cbind(t(fit$coef[1:(num.lv+1),]), rep(0,p)); }
	if(family %in% c("normal","negative.binomial")) params[,ncol(params)] <- fit$phi + 1e-5
	X.params <- NULL; if(!is.null(X)) X.params <- as.matrix(t(fit$coef[-(1:(num.lv+1)),]))
	lvs=NULL; 
	if(num.lv>0) { 
		lvs <- index[,1:num.lv];
		params[,2:(num.lv+1)][upper.tri(params[,2:(num.lv+1)])] <- 0 }
  
	return(list(params=params, index=lvs, X.params = X.params))
	}

     
####################
## Let's test with Spider data set

#data(spider)
#y <- as.matrix(spider$abund)
# fit <- GLVM.LA(y, X = as.matrix(spider$x), num.lv = 2, family = "negative.binomial", row.eff = TRUE, eps=1e-3)

