library(psychotools)
#library(mirt)
source("R/VA-GLLVM-template.R")
	
data(YouthGratitude)
head(YouthGratitude)
grat <- YouthGratitude[,-(1:3)]
sel.sub <- which(apply(grat,1,sum) == round(apply(grat,1,sum))) ## Remove some weird non-integer ratings
grat <- grat[sel.sub,]
dim(grat)

fit.va <- glvm.va(y = grat, family = "ordinal", num.lv = 3, row.eff = FALSE, eps = 0.01, covmat.struc = "unstructured", plot = FALSE, maxit = 10) ## A larger eps is acceptable here given the size of the dataset. LVs don't change that much after 20 iterations anyway.
str(fit.va)

y <- grat
X <- NULL
family <- "ordinal"
num.lv <- 3
row.eff <- FALSE
max.iter <- 100
eps <- 5 # It is for test, too norestrictive. => default is 1e-4
covmat.struc = "unstructured"
plot <- FALSE
maxit <- 10
trace <- TRUE
seed <- 1234
info <- TRUE

n <- dim(y)[1]
p <- dim(y)[2]
# nは受験者数，pは観測変数（項目）数
num.X <- 0
# predictor が指定する場合には行列に変換してnum.Xにその数（列数）を入れる。そうでなければ0が指定される。
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
res$index # a(theta)
res$params # lambda(discrimination)
res$zeta # beta (cutoff parameter)
res$phi # NULL

new.beta0 <- beta0 <- res$params[, 1]
new.row.params <- row.params <- NULL
# means <- 変分分布の平均
# lambda <- 識別力
# vacov <-
if (row.eff) new.row.params <- row.params <- rep(0, n) # roweffect = 0 * n_person
new.vameans <- vameans <- new.lambda <- lambda <- new.vacov <- vacov <- NULL

if (num.lv > 0) {
    new.vameans <- vameans <- res$index
    new.lambda <- lambda <- as.matrix(res$params[, (ncol(res$params) - num.lv + 1):ncol(res$params)])
    new.lambda[upper.tri(new.lambda)] <- lambda[upper.tri(lambda)] <- 0 # おそらく識別のための制約
    if (covmat.struc == "unstructured") { # 受験者ごとにcovmatを作る。
        new.vacov <- vacov <- array(NA, dim = c(n, num.lv, num.lv))
        for (i in 1:n) {
            new.vacov[i, , ] <- vacov[i, , ] <- diag(rep(1, num.lv))
        }
    }
    if (covmat.struc == "diagonal") {# こっちは全受験者で共通？
        new.vacov <- vacov <- matrix(1, n, num.lv)
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

# init
current.loglik <- -1e6
iter <- 1
err <- 10
diag.iter <- 5

# Start iteration (in the while loop)
# Convergence criterion
(err > (1 + eps) || err < (1 - eps)) && iter <= max.iter)
# indicator
if (trace) cat("Iteration:", iter, "\n")
## Use a diagonal covariance matrix for the first diag.iter iterations, then switch to unstructured covariance matrix.
# 最初の不安定なパラメタでcovarianceをunstructureに推定することが難しいから，最初の5回はdiagで推定する。
if(covmat.struc == "unstructured" & iter <= diag.iter) { 
    tmp.covmat.struc <- "diagonal"; 
    if(iter == 1) { new.vacov <- vacov <- matrix(1,n,num.lv) }
    if(iter == diag.iter) { tmp.covmat.struc <- "unstructured"; new.vacov <- vacov <- vacov.convert(vacov,type=2) } 
} else { 
    tmp.covmat.struc <- covmat.struc 
}

