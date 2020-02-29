################
## Template for simulation Setting 1 (binary data) for assessing the performance of GLLVM using VA
################
rm(list = ls())
library(ltm)
library(mirt)
library(vegan)
source("VA-GLLVM-template.R")
	
## Dataset generation
create.binary.life <- function(true.fit, link = "probit") {
	n <- nrow(true.fit$lv); p <- nrow(true.fit$lv.coefs)
	sim.y <- matrix(0,n,p)
	eta <- cbind(1,as.matrix(true.fit$lv))%*%t(as.matrix(true.fit$lv.coefs)) ## Assumes lv.coefs has beta_0j as first column
	if(!is.null(true.fit$X)) eta <- eta + true.fit$X%*%t(true.fit$X.coefs)
	
	for(j in 1:p) {  z.j <- rnorm(n,eta[,j],1); sim.y[which(z.j > 0),j] <- 1 }
	while(is.list(apply(sim.y,2,table))) { ## Ensures all items have at least one success
		sim.y <- matrix(0,n,p)
		for(j in 1:p) {  z.j <- rnorm(n,eta[,j],1); sim.y[which(z.j > 0),j] <- 1 } }
	return(list(y = sim.y, true.ft = true.fit)) 
	}


simpa.noX <- vector("list",3)
n.multiplier <- c(1,2,4)
for(t in 1:3) { 
	n.items <- 10 ## This can be adjusted 
	true.fit <- list(lv = rbind(rmvnorm(25*n.multiplier[t],c(-2,2)),rmvnorm(15*n.multiplier[t],c(0,-1)),rmvnorm(10*n.multiplier[t],c(1,1))), z = rep(1:3,c(25,15,10)*n.multiplier[t]))
	true.fit$lv.coefs <- cbind(runif(n.items,-1,1),seq(-2,2,length=n.items),seq(1,-1,length=n.items))
	simpa.noX[[t]] <- replicate(1000,create.binary.life(true.fit=true.fit)) }
save(simpa.noX, file = "simpa10itemlv2noX.RData")


##############
## Example scatter plot of scores
rm(list = ls())
load("simpa10itemlv2noX.RData")
par(las = 1, cex = 1.5)
plot(simpa.noX[[2]][2,2]$true.ft$lv, pch = simpa.noX[[2]][2,2]$true.ft$z, xlab = "LV1", ylab = "LV2", main = "Scatterplot of true latent variables", cex = 2)


##################
## Testing sims
load("simpa10itemlv2noX.RData")

all.times <- matrix(NA,1000,ncol=6)
all.lv.errors <- all.lambda.errors <- matrix(NA,1000,6) 
colnames(all.times) <- colnames(all.lv.errors) <- colnames(all.lambda.errors) <- c("LTM-hybrid","MIRT-EM","MIRT-MHRM","VA-Diag","VA-unstructured","Laplace")

for(t in 1:1000) {
	cat("Up to iteration",t,"\n")
	cw.dat <- simpa.noX[[comp]][,t]

	
	tic <- proc.time()
	fit <- try(ltm(cw.dat$y ~ z1 + z2, IRT.param = F),silent=T) ## control = list(GHk = xxx) to control quadrature points
	if(!inherits(fit,"try-error")) { 
		fit.scores <- factor.scores(fit, resp.patterns = cw.dat$y, method = "EB"); 
		all.lv.errors[t,1] <- procrustes(cw.dat$true.ft$lv, fit.scores$score.dat[,c("z1","z2")]+rnorm(nrow(cw.dat$y),0,1e-4), symmetric=T)$ss ## Add an little amount to the scores so that the procrustes error remains calculate-able
		all.lambda.errors[t,1] <- procrustes(cw.dat$true.ft$lv.coefs[,2:3], fit$coef[,2:3], symmetric = T)$ss 		
		}
	toc <- proc.time()
	all.times[t,1] <- (toc-tic)[3]

	
	tic <- proc.time()
	fit <- try(mirt(cw.dat$y,2,SE=T,SE.type="Louis",verbose=F),silent=T) ## quadpts=xxx to control quadrature points
	if(!inherits(fit,"try-error")) { 
		fit.scores <- as.data.frame(fscores(fit, response.pattern = cw.dat$y, method = "MAP")); 
		all.lv.errors[t,2] <- procrustes(cw.dat$true.ft$lv, fit.scores[,c("F1","F2")]+rnorm(nrow(cw.dat$y),0,1e-4), symmetric=T)$ss  
		probit.coefs <- matrix(NA,nrow(cw.dat$true.ft$lv.coefs),3); for(l in 1:nrow(probit.coefs)) { probit.coefs[l,] <- coef(fit)[[l]][1,1:3] }
		all.lambda.errors[t,2] <- procrustes(cw.dat$true.ft$lv.coefs[,2:3], probit.coefs[,1:2], symmetric = T)$ss 						
		}
	toc <- proc.time()
	all.times[t,2] <- (toc-tic)[3]

	
	tic <- proc.time()
	fit <- try(mirt(cw.dat$y,2,method="MHRM",SE=T,SE.type="Louis",verbose=F),silent=T)
	if(!inherits(fit,"try-error")) { 
		fit.scores <- as.data.frame(fscores(fit, response.pattern = cw.dat$y, method = "MAP")); 
		all.lv.errors[t,3] <- procrustes(cw.dat$true.ft$lv, fit.scores[,c("F1","F2")]+rnorm(nrow(cw.dat$y),0,1e-4), symmetric=T)$ss 
		probit.coefs <- matrix(NA,nrow(cw.dat$true.ft$lv.coefs),3); for(l in 1:nrow(probit.coefs)) { probit.coefs[l,] <- coef(fit)[[l]][1,1:3] }
		all.lambda.errors[t,3] <- procrustes(cw.dat$true.ft$lv.coefs[,2:3], probit.coefs[,1:2], symmetric = T)$ss 						
		}
	toc <- proc.time()
	all.times[t,3] <- (toc-tic)[3]
	
	
	tic <- proc.time()
	fit <- try(glvm.va(y=cw.dat$y, family = "binomial", num.lv = 2, row.eff = F, covmat.struc = "diagonal", trace = F, plot = F, sd = T),silent=T)
	if(!inherits(fit,"try-error")) { 
		all.lv.errors[t,4] <- procrustes(cw.dat$true.ft$lv, fit$lvs+rnorm(nrow(cw.dat$y),0,1e-4), symmetric=T)$ss
		all.lambda.errors[t,4] <- procrustes(cw.dat$true.ft$lv.coefs[,2:3], fit$theta + rnorm(ncol(cw.dat$y),0,1e-4), symmetric = T)$ss 
		} 
	toc <- proc.time()
	all.times[t,4] <- (toc-tic)[3]
	
	
	tic <- proc.time()
	fit <- try(glvm.va(y=cw.dat$y, family = "binomial", num.lv = 2, row.eff = F, covmat.struc = "unstructured", trace = F, plot = F, sd = T),silent=T)
	if(!inherits(fit,"try-error")) { 
		all.lv.errors[t,5] <- procrustes(cw.dat$true.ft$lv, fit$lvs+rnorm(nrow(cw.dat$y),0,1e-4), symmetric=T)$ss
		all.lambda.errors[t,5] <- procrustes(cw.dat$true.ft$lv.coefs[,2:3], fit$theta + rnorm(ncol(cw.dat$y),0,1e-4), symmetric = T)$ss 
		} 
	toc <- proc.time()
	all.times[t,5] <- (toc-tic)[3]

	
	tic <- proc.time()
	fit <- try(glvm.laplace(y=cw.dat$y, family = "bernoulli", num.lv = 2, row.eff = FALSE, eps = 1e-4, maxit = 1, max.iter = 50, plot = FALSE, info = TRUE),silent=T) 
	if(!inherits(fit,"try-error")) { 
		all.pro.errors[t,6] <- procrustes(cw.dat$true.ft$lv, fit$lvs + rnorm(nrow(cw.dat$y),0,1e-4), symmetric=T)$ss
		all.mse[t,6] <- procrustes(cw.dat$true.ft$lv.coefs[,2:3], fit$lambdas + rnorm(ncol(cw.dat$y),0,1e-4), symmetric = T)$ss 	
		} 
	toc <- proc.time()
	all.times[t,6] <- (toc-tic)[3]
	}
	
	
	
	