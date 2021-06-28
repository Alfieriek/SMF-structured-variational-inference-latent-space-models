#the plot function is incompatible with dplyr package
#application of PMP on the classroom dataset, visualizing the latent positions given the dynamic network 

rm(list=ls())
set.seed(1234)
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

source('mix_DN.R')
source('helper.R')
source('mix_GaussianDN.R')

data(McFarland_cls33_10_16_96)


rawVerts<-read.table(paste(path.package('networkDynamic'),
                           "/extdata/cls33_10_16_96_vertices.tsv",sep=''),header=TRUE,
                     stringsAsFactors = FALSE)


rawEdges<-read.table(paste(path.package('networkDynamic'),
                           "/extdata/cls33_10_16_96_edges.tsv",sep=''),header=TRUE,
                     stringsAsFactors = FALSE)
n = 20 
T = 4
Y = vector("list", T)
for ( t in 1:T){
  Y[[t]] = matrix(rep(0,n^2),nrow=n)
}
for (k in 1:length(rawEdges$end_minute)){
  t = ceiling(rawEdges$end_minute[k] /11)
  if (rawEdges$interaction_type[k] == 'task')
  {Y[[t]][rawEdges$from_vertex_id[k],rawEdges$to_vertex_id[k]] = 1}
}
for ( t in 1:T){
  Y[[t]] = t(Y[[t]])+Y[[t]]
  Y[[t]][Y[[t]]>1] = 1
}
rawVerts$role[rawVerts$role == 'grade_12'] = 'grade_11'

rawVerts$role[rawVerts$role == 'grade_11'] = 'student'


d = 2
Mix_list = mix_DN(Y = Y, mean_beta_prior = 0, sigma_beta_prior = 1, sigma_X1 =0.5, trans_sd = 0.5,
                  gap = 1, rho =1)

X_pos =  vector("list", T)

X_pos[[1]] = Mix_list$Mean_X[[1]]


for (t in 2:T){
  X_pos[[t]]= t(procrustes_r(t(X_pos[[t-1]]),t(Mix_list$Mean_X[[t]]))$B.transformed)
}
par(mfrow=c(2,2))
par(cex=0.8, mai=c(0.3,0.3,0.3,0.3))
for (t in 1:T){
  plot_clus(X_pos[[t]],Y[[t]],as.factor(rawVerts$role),t)
}