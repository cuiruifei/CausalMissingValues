#############################################################################################
# Goal: This is a demo to show how the Rank PC algorithm for incomplete data
#############################################################################################

#### 0. Dependent Packages and Parameters ####

## dependencies
require(pcalg)
source('R/gaussCItestLocal.R')
source('R/inferCopulaModel.R')

## parameters
# No. of variables
p <- 15
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
## nonparanormal transformation
# here we assume the observed distribution of each margin is chi-squared with df;
# one could change it to the desired distribution.
YY <- apply(Z, 2, function(x) qchisq(pnorm(x), df = sample(2:10, 1)))
## add missing values
# the expected percentage of missing values
beta <- 0.25
# please choose MCAR or MAR
# MCAR
Y <- apply(YY, 2, function(x){x[sample(n, round(n*runif(1,0,2*beta)))] = NA; x})
# MAR
Y <- addMAR(YY, beta)


#### 2. Causal Discovery ####

## estimate underlying correlation matrix and (effective) sample size
# No. of missing values in each margin
count.NA <- apply(Y, 2, function(x){ sum(is.na(x))})
# rank correlations (here, transformed Kendall's tau. one could use others.) 
corr.rank <- sin(pi/2 * cor(Y, use = 'pairwise.complete.obs', method = 'kendall'))
# local effective sample size (LESS)
less.rank <- (1-count.NA/n) %*% t(1-count.NA/n) * n
# global effective sample size (GESS)
gess.rank <- mean(less.rank[upper.tri(less.rank)])

## call the standard PC algorithm for causal discovery
# RPC (rank PC) + SS (original sample size)
graph.rpc.ss <- pc(suffStat = list(C = corr.rank, n = n), 
                   indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
# RPC + GESS
graph.rpc.gess <- pc(suffStat = list(C = corr.rank, n = gess.rank), 
                     indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
# RPC + LESS
graph.rpc.less <- pc(suffStat = list(C = corr.rank, n = n, ESS.Mat = less.rank), 
                     indepTest = gaussCItestLocal, labels = var.names, alpha = alpha, conservative = T)

