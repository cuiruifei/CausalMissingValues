#####################################################################################
# Goal: Plot the property of Rank/Copula PC to for nonparanormal data missing values.
#
# Metrics: TPR, FPR, SHD against percentage of missing values.
#
# Setting: beta = {0.1, 0.2, 0.3}, n = {100, 500, 1000}
#
#####################################################################################


#### 0. Load Packages ####

## packages
library(ggplot2)
library(gridExtra)
library(dplyr)
library(psych)
library(ggpubr)

## help function
data2stat <- function(data){
  
  describe(data) %>%
    mutate(ci = qt(0.975,n-1) * se) %>%
    mutate(Methods = factor(c(rep('PC+LD',3), rep('PC+MS', 3), rep('RPC+SS',3), rep('RPC+GESS',3), rep('RPC+LESS',3), rep('CoPC+SS',3), rep('CoPC+GESS',3), rep('CoPC+LESS',3)), 
                            levels = c('PC+LD', 'PC+MS', 'RPC+SS', 'RPC+GESS','RPC+LESS', 'CoPC+SS','CoPC+GESS','CoPC+LESS'), 
                            labels = c('PC+LD', 'PC+MS', 'RPC+SS', 'RPC+GESS','RPC+LESS', 'CoPC+SS','CoPC+GESS','CoPC+LESS'))) %>%
    mutate(Beta.range = factor(rep(c(0.1,0.2,0.3),8))) %>%
    select(Methods, Beta.range, mean, ci)
}


#### 1. Read Data ####

## MCAR
data.n100 = read.table("results/TFPRSHD_ASIS_MCAR_n100")
data.n500 = read.table("results/TFPRSHD_ASIS_MCAR_n500")
data.n1000 = read.table("results/TFPRSHD_ASIS_MCAR_n1000")

# ## MAR
# data.n100 = read.table("results/TFPRSHD_ASIS_MAR_n100")
# data.n500 = read.table("results/TFPRSHD_ASIS_MAR_n500")
# data.n1000 = read.table("results/TFPRSHD_ASIS_MAR_n1000")


t(matrix(apply(data.n500, 2, mean), 3))
p = ncol(data.n100)

## data: 100 * 63; 
# column names
# 1:(p/3) for TPR (3*7)
# (p/3 + 1):(2*p/3) for FPR
# (2*p/3 + 1):p for SHD 


#### 2. Plot ####

#
pd <- position_dodge(0.7) # move them to the left and right

## plot 1: TPR (n=100)
# statistical values
stat = data2stat(data.n100[,1:(p/3)]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
TPR_n100 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 1)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) +
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('TPR') + ggtitle('TPR (n=100)') +
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 2: FPR (n=100)
# statistical values
stat = data2stat(data.n100[,(p/3 + 1):(2*p/3)]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
FPR_n100 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 0.11)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) +
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('FPR') + ggtitle('FPR (n=100)') +
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 3: SHD (n=100)
# statistical values
stat = data2stat(data.n100[,(2*p/3 + 1):p]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
SHD_n100 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color = Methods)) + coord_cartesian(ylim = c(1, 8)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) +
  #geom_bar(position=pd, stat="identity") +
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) +
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('SHD') + ggtitle('SHD (n=100)') +
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 4: TPR (n=500)
# statistical values
stat = data2stat(data.n500[,1:(p/3)]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
TPR_n500 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 1)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) + 
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('TPR') + ggtitle('TPR (n=500)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 5: FPR (n=500)
# statistical values
stat = data2stat(data.n500[,(p/3 + 1):(2*p/3)]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
FPR_n500 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 0.11)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) + 
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('FPR') + ggtitle('FPR (n=500)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 6: SHD (n=500)
# statistical values
stat = data2stat(data.n500[,(2*p/3 + 1):p]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
SHD_n500 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(1, 8)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) + 
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('SHD') + ggtitle('SHD (n=500)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 7: TPR (n=1000)
# statistical values
stat = data2stat(data.n1000[,1:(p/3)]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
TPR_n1000 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 1)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) + 
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('TPR') + ggtitle('TPR (n=1000)') + 
  theme_bw() +
  theme(legend.justification=c(0.01,0.01),legend.position=c(0.01,0.01), legend.text=element_text(size=6.5), legend.title = element_blank(), legend.key.size = unit(0.27, "cm")) +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 8: FPR (n=1000)
# statistical values
stat = data2stat(data.n1000[,(p/3 + 1):(2*p/3)]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
FPR_n1000 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(0, 0.11)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) + 
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('FPR') + ggtitle('FPR (n=1000)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))

## plot 9: SHD (n=1000)
# statistical values
stat = data2stat(data.n1000[,(2*p/3 + 1):p]) #coord_cartesian(ylim = c(0.7, 0.9))
# ggplot
SHD_n1000 <- ggplot(stat, aes(x=Beta.range, y=mean, group=Methods, color=Methods)) + coord_cartesian(ylim = c(1, 8)) +
  geom_line(position=pd) + scale_color_manual(values = c("purple", "#F8766D", "#C49A00", "#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) + 
  geom_point(aes(shape=Methods), position=pd) +scale_shape_manual(values=c(9, seq(0,7))) + 
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=pd) + xlab('beta') + ylab('SHD') + ggtitle('SHD (n=1000)') + 
  theme_bw() +
  theme(legend.position='none') +
  theme(text = element_text(size=8), plot.title = element_text(hjust = 0.5,  size = 8))


## save Figure
plot.Asia = ggarrange(TPR_n100, FPR_n100, SHD_n100, TPR_n500, FPR_n500, SHD_n500, TPR_n1000, FPR_n1000, SHD_n1000, nrow = 3, ncol = 3, common.legend = TRUE, legend="right")
plot.Asia
# # MCAR
# pdf(file = 'results/Asia_MCAR_nonparanormal.pdf', width = 7.5, height = 5)
# plot.Asia
# dev.off()
# # MAR
# pdf(file = 'results/Asia_MAR_nonparanormal.pdf', width = 7.5, height = 5)
# plot.Asia
# dev.off()


