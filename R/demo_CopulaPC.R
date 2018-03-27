#############################################################################################
# Goal: This is a demo to show how the Copula PC algorithm works for incomplete data.
#############################################################################################

#### 0. Dependent Packages and Parameters ####

## dependencies
require(pcalg)
require(sbgcop)
require(infotheo)
source('R/gaussCItestLocal.R')
source('R/inferCopulaModel.R')
source('R/addMAR.R')

## parameters
# No. of variables
p <- 10
# sample size
n <- 1000
# significance level
alpha <- 0.05


#### 1. Simulate Data ####

## generate a random DAG
g <- randomDAG(p, 2/(p-1))
var.names <- g@nodes
## simulate Gaussian data given the DAG
Z <- rmvDAG(n, g)
## discretize some variables into ordinal
# index of discrete variables
iy <- sample(1:p, round(p/2))
# YY, mixed complete data
YY <- Z
YY[,iy] <- matrix(unlist(apply(Z[,iy], 2, discretize, nbins = 2)), nrow = n, byrow = F)
## add missing values
# the expected percentage of missing values
beta <- 0.25
# please choose MCAR or MAR
# Y, mixed incomplete data
# MCAR
Y <- apply(YY, 2, function(x){x[sample(n, round(n*runif(1,0,2*beta)))] = NA; x})
# MAR
ZZ <- addMAR(Z, beta)
Y <- YY
Y[is.na(ZZ)] <- NA

#### 2. Causal Discovery ####

## estimate underlying correlation matrix and (effective) sample size
# copula object
cop.obj <- inferCopulaModel(Y, nsamp = 1000, S0 = diag(p)/n, verb = T)
# correlation matrix samples
C_samples <- cop.obj$C.psamp[,, 501:1000]
# average correlation matrix
corr.cop <- apply(C_samples, c(1,2), mean)
# local effective sample size
less.cop <- ((1-corr.cop^2)^2)/apply(C_samples,c(1,2), var)
# global effective sample size
gess.cop <- mean(less.cop[upper.tri(less.cop)])

## call the PC algorithm for causal discovery
# CoPC + SS
graph.cpc.ss <- pc(suffStat = list(C = corr.cop, n = n), 
                   indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
# CoPC + GESS
graph.cpc.gess <- pc(suffStat = list(C = corr.cop, n = gess.cop), 
                     indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
# CoPC + LESS
graph.cpc.less <- pc(suffStat = list(C = corr.cop, n = n, ESS.Mat = less.cop), 
                     indepTest = gaussCItestLocal, labels = var.names, alpha = alpha, conservative = T)
