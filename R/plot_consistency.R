###############################################################################
# Goal: Plot the behaviour of various correlation estimators for bivariate
#       normal data with missing values.
#
# Methods: 1, mean substitution (MS);
#          2, hot deck (HD);
#          3, Kendall's tau (Tau);
#          2, copula estimator (Cop);
#
# Setting: beta = 0.25, n = {100, 500, 1000}
#
###############################################################################


#### 0. Denpendencies ####

## packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(psych)
library(latex2exp)

## help function
data2stat <- function(data){
  
  describe(data) %>%
    #mutate(ci = qt(0.975,n-1) * se) %>%
    mutate(Methods = factor(c(rep('MS',3), rep('HD',3), rep('Tau',3), rep('Cop',3)), 
                            levels = c('MS', 'HD', 'Tau', 'Cop'), labels = c('MS', 'HD', 'Tau', 'Cop'))) %>%
    mutate(sample_size = factor(rep(c('n = 100', 'n = 500', 'n = 1000'), 4),levels = c('n = 100', 'n = 500', 'n = 1000'), labels = c('n = 100', 'n = 500', 'n = 1000'))) %>%
    select(Methods, sample_size, mean, sd)
}



#### 1. Read Data ####

## for MCAR
data.rho0 = read.table("results/consistency_rho0_MCAR")
data.rho0.3 = read.table("results/consistency_rho0.3_MCAR")
data.rho0.6 = read.table("results/consistency_rho0.6_MCAR")
data.rho0.9 = read.table("results/consistency_rho0.9_MCAR")
# ## for MAR
# data.rho0 = read.table("results/consistency_rho0_MAR")
# data.rho0.3 = read.table("results/consistency_rho0.3_MAR")
# data.rho0.6 = read.table("results/consistency_rho0.6_MAR")
# data.rho0.9 = read.table("results/consistency_rho0.9_MAR")


#### 2. Plot ####

#
pd <- position_dodge(0.4) # move them to the left and right

## plot 1: rho = 0
# statistical values
stat = data2stat(data.rho0) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
rho0 <- ggplot(stat, aes(x=Methods, y=mean, group=sample_size, color=sample_size)) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_point(aes(shape=sample_size), position=pd) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=pd) + xlab('') + ylab('') + ggtitle(TeX('$\\rho = 0$')) +
  theme_bw() +
  theme(legend.position='none') +
  #theme(legend.justification=c(0.01,0.01),legend.position=c(0.01,0.01), legend.text=element_text(size=6.5), legend.title = element_blank(), legend.key.size = unit(0.27, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 2: rho = 0.3
# statistical values
stat = data2stat(data.rho0.3) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
rho0.3 <- ggplot(stat, aes(x=Methods, y=mean, group=sample_size, color=sample_size)) +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(yintercept=0.3, linetype="dotted") +
  geom_point(aes(shape=sample_size), position=pd) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=pd) + xlab('') + ylab('') + ggtitle(TeX('$\\rho = 0.3$')) +
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 3: rho = 0.6
# statistical values
stat = data2stat(data.rho0.6) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
rho0.6 <- ggplot(stat, aes(x=Methods, y=mean, group=sample_size, color=sample_size)) +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(yintercept=0.6, linetype="dotted") +
  geom_point(aes(shape=sample_size), position=pd) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=pd) + xlab('') + ylab('') + ggtitle(TeX('$\\rho = 0.6$')) +
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 4: rho = 0.9
# statistical values
stat = data2stat(data.rho0.9) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
rho0.9 <- ggplot(stat, aes(x=Methods, y=mean, group=sample_size, color=sample_size)) +
  #coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_hline(yintercept=0.9, linetype="dotted") +
  geom_point(aes(shape=sample_size), position=pd) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=pd) + xlab('') + ylab('') + ggtitle(TeX('$\\rho = 0.9$')) +
  theme_bw() +
  #theme(legend.position='none') +
  theme(legend.justification=c(0.99,0.01),legend.position=c(0.99,0.01), legend.text=element_text(size=6.5), legend.title = element_blank(), legend.key.size = unit(0.27, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))


grid.arrange(rho0, rho0.3, rho0.6, rho0.9, nrow = 1, ncol = 4)

## save files for MCAR
# pdf(file = 'results/consistency_rho_MCAR.pdf', width = 8, height = 2)
# grid.arrange(rho0, rho0.3, rho0.6, rho0.9, nrow = 1, ncol = 4)
# dev.off()
# ## save files for MAR
# pdf(file = 'results/consistency_rho_MAR.pdf', width = 8, height = 2)
# grid.arrange(rho0, rho0.3, rho0.6, rho0.9, nrow = 1, ncol = 4)
# dev.off()
