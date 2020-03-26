################
## Application of GLLVM to ordinal data to Youth Gratitude dataset
################
rm(list = ls())
library(psychotools)
#library(mirt)
source("R/VA-GLLVM-template.R")
	
data(YouthGratitude)
head(YouthGratitude)
grat <- YouthGratitude[,-(1:3)]
sel.sub <- which(apply(grat,1,sum) == round(apply(grat,1,sum))) ## Remove some weird non-integer ratings
grat <- grat[sel.sub,]
dim(grat)


fit.va <- glvm.va(y = grat, family = "ordinal", num.lv = 3, row.eff = FALSE, eps = 0.01, covmat.struc = "unstructured", plot = FALSE, maxit = 50) ## A larger eps is acceptable here given the size of the dataset. LVs don't change that much after 20 iterations anyway.
fit.va$beta
fit.va$lambda
fit.va$zeta
fit.va$A[1,,]

plot(fit.va$lvs, col = as.numeric(YouthGratitude$agegroup), xlab = "LV1", ylab = "LV2", main = "A: Unconstrained ordination of youths")
legend("topleft", col = unique(as.numeric(YouthGratitude$agegroup)), pch = 1, legend = levels(YouthGratitude$agegroup))

write.csv(fit.va$lvs, "data/va_theta_original.csv")

agegroup2 <- YouthGratitude[sel.sub,"agegroup"]; levels(agegroup2) <- c(0,0,1,1,1,1)
fit.va2 <- glvm.va(y = grat, X = as.matrix(as.numeric(agegroup2)-1), family = "ordinal", num.lv = 2, row.eff = FALSE, eps = 5, covmat.struc = "unstructured", plot = FALSE, maxit = 10) ## A larger eps is acceptable here given the size of the dataset. LVs don't change that much after 20 iterations anyway.

plot(fit.va2$lvs, col = as.numeric(YouthGratitude$agegroup), xlab = "LV1", ylab = "LV2", main = "B: Residual ordination of youths")
legend("topleft", col = unique(as.numeric(YouthGratitude$agegroup)), pch = 1, legend = levels(YouthGratitude$agegroup))

# Test

sample(c(0, 1), nrow(grat), replace = TRUE)

fit.va <- glvm.va(y = grat, X = matrix(sample(c(0, 1), nrow(grat), replace = TRUE), nrow = nrow(grat), 1), family = "ordinal", num.lv = 3, row.eff = FALSE, eps = .1, covmat.struc = "diagonal", plot = FALSE, maxit = 10)
str(fit.va)

fit.va$beta
fit.va$lambda
fit.va$zeta
fit.va$lvs # thetaのスケールを固定できない...

library(mirt)
fit <- mirt(grat, 3, method = "MHRM")
item <- coef(fit, simplify = TRUE)$item
theta <- fscores(fit)
head(theta)
write.csv(item, "data/mirt_mhrm_para.csv")
write.csv(theta, "data/mirt_mhrm_theta.csv")