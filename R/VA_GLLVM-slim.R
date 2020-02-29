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

source("R/VA-GLLVM-subfunctions.R")

#' Variational approximation of generalized linear latent variavle models
#' 
#' @param 
#' @param X A covariate matrix,, i.e. Testlet effect construction, LLTM.
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @return 
#' @export
#' 
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
	n <- dim(y)[1]
	p <- dim(y)[2]
	# nは受験者数，pは観測変数（項目）数
	num.X <- 0
	# predictor を指定する場合には行列に変換してnum.Xにその数（列数）を入れる。そうでなければ0が指定される。
	if (!is.null(X)) {
		X <- as.matrix(X)
		num.X <- dim(X)[2]
	} ## no of predictors

	# 使えるfamilyはpoisson, neg.binomial, binomial, ordinal （これは論文で提案されている手法に一致しているようだ）
	if (!(family %in% c("poisson", "negative.binomial", "binomial", "ordinal"))) {
		stop("Inputed family not allowed...sorry =(")
	}
	if (!(covmat.struc %in% c("unstructured", "diagonal"))) {
		stop("A_i, i.e. covariance of vartiational distribution for latent variable, must be either unstructured or diagonal in form. Thanks.")
	}
	if (num.lv == 1) covmat.struc <- "diagonal" ## Prevents it going to "unstructured" loops and causing chaos
	trial.size <- 1
	
	y <- as.matrix(y)
	if (!is.numeric(y)) stop("y must a numeric. If ordinal data, please convert to numeric with lowest level equal to 1. Thanks")
	if (is.null(rownames(y))) rownames(y) <- paste("Row", 1:n, sep = "")
	if (is.null(colnames(y))) colnames(y) <- paste("Col", 1:p, sep = "")
	if (!is.null(X)) if (is.null(colnames(X))) colnames(X) <- paste("x", 1:ncol(X), sep = "")
	
	if (!is.null(X)) X <- as.matrix(X) # 不要？
	if (family == "ordinal") {
	    max.levels <- apply(y, 2, max) # 最大カテゴリ数を取得
	    # どうやら現状だと複数のfamilyを同時にはサポートしていない。 => inefficientなだけで，2パラとGRMの混合推定はできそう？
	    if (any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks")
	}
	

	## Get initial starting values by fitting a bunch of GLMs
	res <- start.values.va(y = y, X = X, family = family, trial.size = trial.size, num.lv = num.lv)
	new.beta0 <- beta0 <- res$params[, 1]
	new.beta <- beta <- NULL
	if (!is.null(X)) new.beta <- beta <- res$params[, 2:(num.X + 1)]
	new.row.params <- row.params <- NULL
	if (row.eff) new.row.params <- row.params <- rep(0, n)
	# means <- 変分分布の平均
	# lambda <- 識別力
	# vacov <- A
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

	new.zeta <- zeta <- NULL
	if (family == "ordinal") {
	    new.zeta <- zeta <- res$zeta
	}
	new.phi <- phi <- NULL
	if (family == "negative.binomial") {
	    new.phi <- phi <- res$phi
	}

	
	
	current.loglik <- -1e6; iter <- 1; err <- 10; diag.iter <- 5
    while((err > (1 + eps) || err < (1 - eps)) && iter <= max.iter) {
#	while(err > eps && iter <= max.iter) {
		if(trace) cat("Iteration:", iter, "\n")

		
		## Use a diagonal covariance matrix for the first diag.iter iterations, then switch to unstructured covariance matrix. 
		if (covmat.struc == "unstructured" & iter <= diag.iter) {
		    tmp.covmat.struc <- "diagonal"
		    if (iter == 1) {
		        new.vacov <- vacov <- matrix(1, n, num.lv)
		    }
		    if (iter == diag.iter) {
		        tmp.covmat.struc <- "unstructured"
		        new.vacov <- vacov <- vacov.convert(vacov, type = 2) # 3d array!
		    }
		} else { 
            tmp.covmat.struc <- covmat.struc 
        }			
		# Update a(vamean)
		if (num.lv > 0) {
			## Update vameans (a_i) by maximizing lower bound for loglik
		    q <- try(optim(c(vameans), method = "BFGS", fn = ll0, gr = grad.var, x = c(lambda, beta0, beta, row.params), vacov = vacov, phi = phi, zeta = zeta, control = list(trace = 0, fnscale = -1, maxit = maxit, reltol = 1e-2)), silent = T)

		    if (!inherits(q, "try-error")) {
		        if (trace) cat("Variational parameters updated", "\n")
		        new.vameans <- q$par[1:(num.lv * n)]
		        new.vameans <- matrix(new.vameans, n, num.lv)
		    } else {
		        new.vameans <- vameans
			}	
		
		
			## Update covariance matrix A_i via iterative fixed-point algorithm
			eta.mat <- cbind(rep(1, n), X) %*% t(cbind(beta0, beta)) + new.vameans %*% t(lambda) + calc.quad(vacov, lambda, tmp.covmat.struc)$mat
			if (!is.null(row.params)) eta.mat <- eta.mat + matrix(row.params, n, p)
			if (family == c("poisson")) {
				mu.mat <- exp(eta.mat)
			}
			if (family == c("negative.binomial")) {
				phi.mat <- matrix(phi, n, p, byrow = T)
				eta.mat <- eta.mat - 2 * calc.quad(vacov, lambda, tmp.covmat.struc)$mat
				mu.mat <- (y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat))
			}

			# Update a for i in 1:n by fitting a penalized probit GLM
			for (i in 1:n) {
				error <- 1
				vacov.iter <- 0
				if (tmp.covmat.struc == "unstructured") new.vacov.mat <- vacov[i, , ]
				if (tmp.covmat.struc == "diagonal") new.vacov.mat <- diag(x = vacov[i, ], nrow = num.lv)
				# Update vacov(A), specifically, A = (I + ∑λλᵀ)⁻¹
				while (error > 1e-2 && vacov.iter < 100) {
					cw.vacov.mat <- new.vacov.mat
					
					if (tmp.covmat.struc == "unstructured") {
						lambda2 <- sapply(1:p, function(j, lambda) lambda[j, ] %*% t(lambda[j, ]), lambda = new.lambda)
						lambda2 <- t(lambda2)
		
						if (family %in% c("poisson", "negative.binomial")) {
							new.vacov.mat <- solve(diag(rep(1, num.lv)) + matrix(apply(mu.mat[i, ] * lambda2, 2, sum), nrow = num.lv))
						}
						if (family %in% c("binomial", "ordinal")) {
							new.vacov.mat <- solve(diag(rep(1, num.lv)) + matrix(apply(lambda2, 2, sum), nrow = num.lv))
						}
					}
		
					if (tmp.covmat.struc == "diagonal") {
						lambda2 <- new.lambda^2
						if (family %in% c("poisson", "negative.binomial")) {
							new.vacov.mat <- solve(diag(rep(1, num.lv)) + diag(apply(mu.mat[i, ] * lambda2, 2, sum), num.lv, num.lv))
						}
						if (family %in% c("binomial", "ordinal")) {
							new.vacov.mat <- solve(diag(rep(1, num.lv)) + diag(apply(lambda2, 2, sum), num.lv, num.lv))
						}
					}
		
					error <- sum((new.vacov.mat - cw.vacov.mat)^2)
					if (tmp.covmat.struc == "unstructured") vacov[i, , ] <- new.vacov.mat
					if (tmp.covmat.struc == "diagonal") vacov[i, ] <- diag(new.vacov.mat)
		
					eta.mat <- cbind(rep(1, n), X) %*% t(cbind(beta0, beta)) + new.vameans %*% t(lambda) + calc.quad(vacov, lambda, tmp.covmat.struc)$mat
					if (!is.null(row.params)) eta.mat <- eta.mat + matrix(row.params, n, p)
					if (family == c("poisson")) {
						mu.mat <- exp(eta.mat)
					}
					if (family == c("negative.binomial")) {
						phi.mat <- matrix(phi, n, p, byrow = T)
						eta.mat <- eta.mat - 2 * calc.quad(vacov, lambda, tmp.covmat.struc)$mat
						mu.mat <- (y + phi.mat) * (phi.mat / (exp(eta.mat) + phi.mat))
					}
		
					vacov.iter <- vacov.iter + 1
				} ## End of updating A
		
				if (tmp.covmat.struc == "unstructured") new.vacov[i, , ] <- new.vacov.mat
				if (tmp.covmat.struc == "diagonal") new.vacov[i, ] <- diag(new.vacov.mat)
		
				if ((family %in% c("binomial", "ordinal")) & i == 1) break ## Since same A_i matrix for all i, then do one and then break
			}
			
			if (family %in% c("binomial", "ordinal")) {
				if (tmp.covmat.struc == "diagonal") for (i2 in 2:n) new.vacov[i2, ] <- new.vacov[1, ]
				if (tmp.covmat.struc == "unstructured") for (i2 in 2:n) new.vacov[i2, , ] <- new.vacov[1, , ]
			}
		} ## Ends loop for if num.lv > 0
		
		lower.vec <- rep(-30, length(c(lambda, beta0, beta, row.params)))
		upper.vec <- rep(30, length(c(lambda, beta0, beta, row.params)))
		lower.vec[which(upper.tri(lambda))] <- -1e-2; upper.vec[which(upper.tri(lambda))] <- 1e-2
		q <- try(optim(c(lambda,beta0,beta,row.params), v=c(new.vameans), vacov = new.vacov, phi=phi, zeta = zeta, method="L-BFGS-B", lower = lower.vec, upper = upper.vec, fn=ll0, gr=grad.mod, control=list(trace=0, fnscale=-1, maxit = maxit, factr=1e-2)), silent=T)

		if (!inherits(q, "try-error")) {
		    if (trace) cat("Model parameters updated", "\n")
		    x2 <- q$par
		    if (num.lv > 0) {
		        new.lambda <- matrix(x2[1:(p * num.lv)], p, num.lv)
		        new.lambda[upper.tri(new.lambda)] <- 0
		        x2 <- x2[-(1:(p * num.lv))]
		    }
		    new.beta0 <- x2[1:p]
		    x2 <- x2[-(1:p)]
		    new.beta <- NULL
		    if (!is.null(X)) {
		        new.beta <- matrix(x2[1:(p * num.X)], p, num.X)
		        x2 <- x2[-(1:(p * num.X))]
		    }
		    new.row.params <- NULL
		    if (row.eff) {
		        new.row.params <- x2[1:n]
		        x2 <- x2[-(1:n)]
		    }
		} else { 
            new.lambda <- lambda; new.beta0 <- beta0; new.beta <- beta
            new.row.params <- row.params 
        }
		
		
		## Update overdispersion parameter if necessary
		if (family == "negative.binomial") {
		    q <- try(optim(phi, x = c(new.lambda, new.beta0, new.beta, new.row.params), v = c(new.vameans), vacov = new.vacov, zeta = zeta, method = "L-BFGS-B", lower = 1e-3, upper = 1e3, fn = ll0, gr = grad.phi, control = list(trace = 0, fnscale = -1, factr = 1e-3)), silent = T)
		    ll0(x = c(new.lambda, new.beta0, new.beta, new.row.params), v = c(vameans), vacov = vacov, phi = phi)
		
		    if (!inherits(q, "try-error")) {
		        if (trace) cat("Dispersion parameters updated", "\n")
		        new.phi <- q$par[1:p]
		    } else {
		        new.phi <- phi
		    }
		}
					
		
		## Update cutoffs for ordinal data (zeta_j) if required
		if (family == "ordinal") {
		    eta.mat <- cbind(rep(1, n), X) %*% t(cbind(new.beta0, new.beta))
		    if (!is.null(row.params)) eta.mat <- eta.mat + matrix(new.row.params, n, p)
		    if (num.lv > 0) eta.mat <- eta.mat + new.vameans %*% t(new.lambda)
		
		    for (j in 1:p) {
		        if (max(y[, j]) == 2) {
		            new.zeta[j, ] <- zeta[j, ]
		        }
		        if (max(y[, j]) > 2) {
		            constraint.mat <- matrix(0, max(y[, j]) - 2, max(y[, j]) - 2)
		            constraint.mat[1, 1] <- 1 ## Constrains zeta2 > 0
		            if (nrow(constraint.mat) > 1) {
		                for (k in 2:nrow(constraint.mat)) constraint.mat[k, (k - 1):k] <- c(-1, 1)
		            }
		
		            update.zeta <- constrOptim(theta = zeta[j, 2:(max(y[, j]) - 1)], f = func.zetaj, grad = grad.zetaj, ui = constraint.mat, ci = rep(0, max(y[, j]) - 2), j = j, outer.eps = 1e-3, control = list(trace = 0, fnscale = -1))
		            if (!inherits(update.zeta, "try-error")) new.zeta[j, 2:(max(y[, j]) - 1)] <- update.zeta$par
		            if (inherits(update.zeta, "try-error")) new.zeta[j, ] <- zeta[j, ]
		        }
		    }
		    cat("Cutoffs updated \n")
		}
		
			
		## Take values of loglik from optim -function to define stopping rule
		q <- list(value = ll0(c(new.lambda,new.beta0,new.beta,new.row.params), v=c(new.vameans), vacov = new.vacov, phi=new.phi, zeta = new.zeta))
		new.loglik <- q$value
 		err <- abs(new.loglik/current.loglik); 
  		if(trace) cat("New Loglik:", new.loglik,"Current Loglik:", current.loglik, "Ratio", err, "\n")

		
		## Plot old and new ordination points for spp's and sites
		if (trace == TRUE && num.lv <= 2 && num.lv > 0 && plot == TRUE) {
		    par(mfrow = c(1, 2))
		    plot(new.vameans[!duplicated(new.vameans), ], type = "n", xlab = ifelse(num.lv == 1, "Row index (unique elements only)", "LV1"), ylab = "LV2", main = "Ordination of rows")
		    text(new.vameans[!duplicated(new.vameans), ], labels = which(!duplicated(new.vameans) == 1))
		
		    plot(new.lambda, type = "n", xlab = ifelse(num.lv == 1, "Column index", "Coefs for LV1"), ylab = "Coefs for LV2", main = "Ordination of columns")
		    text(new.lambda, labels = seq(1, dim(y)[2]))
		}
	
	
		current.loglik <- new.loglik
		beta0 <- new.beta0
		beta <- new.beta
		lambda <- new.lambda
		phi <- new.phi
		vameans <- new.vameans
		vacov <- new.vacov
		row.params <- new.row.params
		zeta <- new.zeta
		
		iter = iter + 1
	}

		
    ## Bling up the output
    spp.coefs <- matrix(beta0, ncol = 1)
    rownames(spp.coefs) <- colnames(y)
    colnames(spp.coefs) <- "Intercept"
    if (!is.null(X)) {
        spp.coefs <- cbind(beta0, beta)
        colnames(spp.coefs) <- c("Intercept", colnames(X))
        rownames(spp.coefs) <- colnames(y)
    }
	 
	out.list <-  list(y = y, X = X, num.lv = num.lv, row.eff = row.eff, beta = spp.coefs, iter=iter, logLik=new.loglik, family = family, covmat.struc = covmat.struc)
    if (num.lv > 0) {
        rownames(lambda) <- colnames(y)
        colnames(lambda) <- 1:num.lv
        lambda[upper.tri(lambda)] <- 0
        rownames(vameans) <- rownames(y)
        colnames(vameans) <- 1:num.lv
        out.list$A <- vacov
        out.list$lambda <- lambda
        out.list$lvs <- vameans
        out.list$lambda[!is.finite(out.list$lambda)] <- 0
        out.list$lvs[!is.finite(out.list$lvs)] <- 0
    }

	if (family == "negative.binomial") {
	    out.list$phi <- 1 / phi
	    names(out.list$phi) <- colnames(y)
	} ## Flip it back so that it is V = mu + phi*mu^2
	if (row.eff) {
	    out.list$row.params <- row.params
	    names(row.params) <- rownames(y)
	}

	if (family == "ordinal") {
	    out.list$zeta <- zeta
	    rownames(out.list$zeta) <- colnames(y)
	    colnames(out.list$zeta) <- paste(1:(max(y) - 1), "|", 2:max(y), sep = "")
	}
	
	
	if (sd.errors) {
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
	if(family == "ordinal") { 
		max.levels <- apply(y,2,max); 
		if(any(max.levels == 1) || all(max.levels == 2)) stop("Ordinal data requires all columns to have at least has two levels. If all columns only have two levels, please use family == binomial instead. Thanks") }
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

     
