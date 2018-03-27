#############################################################################################
# Goal: This script aims to illustrate the 'Copula PC' algorithm on CFS data.
# Notations: CFS1, all available data; CFS1_0, the resulting data after adding some missing 
#            values manually.
#############################################################################################

#### 0. Dependent Packages and Parameters ####

## dependencies
require(pcalg)
require(sbgcop)
require(dplyr)
require(polycor)
source('R/gaussCItestLocal.R')
source('R/inferCopulaModel.R')

## parameters
#significance level in conditional independence tests
alpha <- 0.05


#### 1. Load Data ####

## read data
CFS1 <- read.csv('data/CFS1.csv', colClasses = c('NULL', rep(NA, 6)))
## basic information
n <- nrow(CFS1)
p <- ncol(CFS1)
var.names <- colnames(CFS1)


#### 2. Causal Discovery ####

#### get 'true' structure (here it means the resulting graph based on all available data)
graph.true <- pc(suffStat = list(C = cor(CFS1, use = 'pairwise.complete.obs'), n = n), 
                 indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
plot(graph.true, main = "Pseudo True Graph")


#### fill in some missing values manually
CFS1_0 <- CFS1
beta <- 0.2
# phy1 missing given obj1
CFS1_0[rank(CFS1_0$phy1) < n*beta,'obj1'] <- NA
# func1 missing given fat1
CFS1_0[rank(CFS1_0$func1) < n*beta,'fat1'] <- NA
# sens1 missing given foc1
CFS1_0[rank(CFS1_0$foc1) < n*beta,'sens1'] <- NA

#### recover the structure from CFS1_0

### Method 1: PC + LD (list-wise deletion)
graph.ld <- pc(suffStat = list(C = sin(pi/2 * cor(CFS1_0, use = 'complete.obs', method = 'kendall')), n = n),
   indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)

### Method 2: Pearson PC

## estimate underlying correlation matrix and (effective) sample size
count.NA <- apply(CFS1_0, 2, function(x){ sum(is.na(x))})
# pearson correlations
corr.rank <- cor(CFS1_0, use = 'pairwise.complete.obs')
# local effective sample size
less.rank <- (1-count.NA/n) %*% t(1-count.NA/n) * n
# global effective sample size
gess.rank <- mean(less.rank[upper.tri(less.rank)])

## call the PC algorithm for causal discovery
# RPC + SS
graph.rpc.ss <- pc(suffStat = list(C = corr.rank, n = n), 
                  indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
# RPC + GESS
graph.rpc.gess <- pc(suffStat = list(C = corr.rank, n = gess.rank), 
                    indepTest = gaussCItest, labels = var.names, alpha = alpha, conservative = T)
# RPC + LESS
graph.rpc.less <- pc(suffStat = list(C = corr.rank, n = n, ESS.Mat = less.rank), 
                    indepTest = gaussCItestLocal, labels = var.names, alpha = alpha, conservative = T)

### Method 3: Copula PC

## estimate underlying correlation matrix and (effective) sample size
# copula object
cop.obj <- inferCopulaModel(CFS1_0, nsamp = 1000, S0 = diag(p)/n, verb = F)
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

