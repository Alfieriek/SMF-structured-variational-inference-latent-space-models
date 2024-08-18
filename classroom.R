#application of SMF on the classroom dataset, visualizing the latent positions given the dynamic network

rm(list=ls())
set.seed(2023)
require('networkDynamicData')
library('rmatio')
library(MASS)
library(mvtnorm)
library(igraph)
library(DiagrammeR)
library(pROC)
library("ggplot2")
library(reshape2)
library("RColorBrewer")
library(gridExtra)
library(grid)
library(ggpubr)
library(cluster)
library(ndtv)
library(animation)


source('helper.R')
source('mix_DN_adaptive.R')
source('mean_field_DN_adaptive.R')


data(McFarland_cls33_10_16_96)


rawVerts<-read.table(paste(path.package('networkDynamic'),
                           "/extdata/cls33_10_16_96_vertices.tsv",sep=''),header=TRUE,
                     stringsAsFactors = FALSE)


rawEdges<-read.table(paste(path.package('networkDynamic'),
                           "/extdata/cls33_10_16_96_edges.tsv",sep=''),header=TRUE,
                     stringsAsFactors = FALSE)
n = 20
T = 8
Y = vector("list", T)
for ( t in 1:T){
  Y[[t]] = matrix(rep(0,n^2),nrow=n)
}
for (k in 1:length(rawEdges$end_minute)){
  t = ceiling(rawEdges$end_minute[k] /(44/T))
  if (rawEdges$interaction_type[k] == 'task')
  {Y[[t]][rawEdges$from_vertex_id[k],rawEdges$to_vertex_id[k]] = 1}
}
for ( t in 1:T){
  Y[[t]] = t(Y[[t]])+Y[[t]]
  Y[[t]][Y[[t]]>1] = 1
}
rawVerts$role[rawVerts$role == 'grade_12'] = 'grade_11'

rawVerts$role[rawVerts$role == 'grade_11'] = 'Student'


rawVerts$role[rawVerts$role == 'instructor'] = 'Instructor'

d = 2


########### predictive part #########

# simulation_iter = 25
# auc_mf =matrix(rep(0,simulation_iter*6),ncol=6)
# auc_smf =matrix(rep(0,simulation_iter*6),ncol=6)
# for (iter in 1:simulation_iter){
# 
# r0=1
# 
# for(t0 in 3:T){
# 
# Mix_list_pred = mix_DN_adaptive(Y = Y[1:(t0-1)], mean_beta_prior = 0, sigma_beta_prior = sqrt(10),gap=1e-2,global_prior = 'Gamma')
# 
# MF_list_pred = mean_field_DN_adaptive(Y = Y[1:(t0-1)], mean_beta_prior = 0, sigma_beta_prior = sqrt(10),gap=0.01)
# 
# 
# 
# auc_test_mf = rep((n-1)*n/2,0)
# auc_test_mix = rep((n-1)*n/2,0)
# auc_res = rep((n-1)*n/2,0)
# r=1
#   for (i in 1:n){
#     for (j in 1:n){
#       if (j<i){
#         auc_test_mf[r] = 1/(1+exp(-MF_list_pred$mean_beta-t(MF_list_pred$Mean_X[[t0-1]][i,])%*% MF_list_pred$Mean_X[[t0-1]][j,]))
#         auc_test_mix[r] = 1/(1+exp(-Mix_list_pred$mean_beta-t(Mix_list_pred$Mean_X[[t0-1]][i,])%*% Mix_list_pred$Mean_X[[t0-1]][j,]))
#         auc_res[r] = Y[[t0]][i,j]
#         r=r+1
#       }
#     }
#   }
# 
# auc_mf[iter,r0] =  auc(roc(auc_res,auc_test_mf))
# auc_smf[iter,r0] = auc(roc(auc_res,auc_test_mix))
# r0 = r0+1
# }
# }
# 
# cat('test AUC for MF',auc_mf,'\n')
# cat('test AUC for SMF',auc_smf,'\n')
# 
# 
# df = NULL
# 
# df$AUC = c(as.vector(auc_mf),as.vector(auc_smf))
# 
# df$iter = rep(rep(1:simulation_iter,6),2)
# 
# df$time_point = rep(as.vector(outer(rep(1,simulation_iter),3:8)),2)
# 
# df$time_point = factor(df$time_point)
# 
# df$method = c(rep('MF',simulation_iter*6),rep('SMF',simulation_iter*6))
# 
# 
# df = as.data.frame(df)
# 
# ggplot(df, aes(x = time_point, y = AUC, fill = method)) +
#   geom_boxplot(outlier.shape = NA) +
#   labs(x = "Time Point", y = "AUC", fill = "Method",
#        title = "Boxplot of AUC by Method") +
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5))


# ########## visualization part #######
Mix_list = mix_DN_adaptive(Y = Y, mean_beta_prior = 0, sigma_beta_prior = sqrt(10),
                  gap = 0.01,trans_sd_init =0.5,global_prior='Gamma')

MF_list = mean_field_DN_adaptive(Y = Y, mean_beta_prior = 0, sigma_beta_prior = sqrt(10),
                           gap = 0.01)


X_pos =  vector("list", T)

X_pos_MF = vector("list", T)

X_pos[[1]] = Mix_list$Mean_X[[1]]

X_pos_MF [[1]] = MF_list$Mean_X[[1]]

X_pos_MF [[1]] = t(procrustes_r(t(X_pos[[1]]),t(X_pos_MF[[1]]))$B.transformed)

for (t in 2:T){
  X_pos[[t]]= t(procrustes_r(t(X_pos[[t-1]]),t(Mix_list$Mean_X[[t]]))$B.transformed)
  X_pos_MF[[t]]= t(procrustes_r(t(X_pos_MF[[t-1]]),t(MF_list$Mean_X[[t]]))$B.transformed)
}


### Estmated Latent Positions at Different Time Points via SMF
m <- matrix(c(1,2,3,4,5,6,7,8,9,9,9,9),nrow = 3,ncol = 4,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.2))
par(mar=rep(1.5,4))
for (t in 1:T){
   plot_clus(X_pos[[t]],Y[[t]],as.factor(rawVerts$role),t)
}
pal <- brewer.pal(4,"Accent")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title(main='Estmated Latent Positions at Different Time Points via SMF',cex.main=1.5, line = -1, outer = TRUE,cex=1.8)
legend('top',inset = 0,legend=levels(as.factor(rawVerts$role)),bty = "n",cex=1.6,col=pal[c(1,2)],pch=c(1,0), horiz = F,pt.cex=c(2,1))



### Estmated Latent Positions at Different Time Points via MF
dev.new()
m <- matrix(c(1,2,3,4,5,6,7,8,9,9,9,9),nrow = 3,ncol = 4,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.2))
par(mar=rep(1.5,4))
for (t in 1:T){
  plot_clus(X_pos_MF[[t]],Y[[t]],as.factor(rawVerts$role),t)
}
pal <- brewer.pal(4,"Accent")
plot(1, type = "n",axes=FALSE, xlab="", ylab="")
title(main='Estmated Latent Positions at Different Time Points via MF',cex.main=1.5, line = -1, outer = TRUE,cex=1.8)
legend('top',inset = 0,legend=levels(as.factor(rawVerts$role)),pch=c(16,16),bty = "n",cex=1.6,col=pal[c(1,2)], horiz = F,pt.cex=c(2,1))






# ### Visualization via ndtv
# dev.new()
# net_list = vector("list", T)
# for (t in 1:T){
#   net_list[[t]] = as.network.matrix(Y[[t]],matrix.type = 'adjacency',directed = F)
# }
# net_dynamic  = networkDynamic(network.list = net_list)
# m <- matrix(c(1,2,3,4,5,6,7,8,9,9,9,9),nrow = 3,ncol = 4,byrow = TRUE)
# layout(mat = m,heights = c(0.4,0.4,0.2))
# par(mar=rep(2,4))
# slice.par = list(start=0,end=7,interval=1,aggregate.dur=1, rule="latest")
# filmstrip(net_dynamic,slice.par=slice.par,mfrow=c(2,4),displaylabels = T)
# title(main='Dynamic Networks Visualization via ndtv',cex.main=1.5, line = -2, outer = TRUE,cex=1.8)


# ### Animated visualization via SMF
# dev.new()
# ani.record(reset = TRUE)
# for (t in 1:T){
#   plot_clus(X_pos[[t]],Y[[t]],as.factor(rawVerts$role),t)
#   ani.record()
# }
# oopts = ani.options(interval = 1)
# saveLatex(ani.replay(), img.name = "record_plot")







