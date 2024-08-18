#application of SMF on the enron email dataset, compare with inverse gamma prior

rm(list=ls())
set.seed(2021)
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
library(networkDynamicData)
library(lubridate)
library(ggraph)
library(animation)
library(Bessel)
source('helper.R')
source('mix_DN_adaptive_miss.R')
source('mix_DN_adaptive_invgamma_miss.R')

data(enronEmails)

data = as.data.frame(enronEmails)

data_date = as.Date(as.POSIXct(data$onset, origin="1970-01-01"))

day(data_date) =1 

data_date = data_date -min(data_date)+1

 data_month = data_date

 
T = length(unique(data_month))
 
 for (t in 1:T){
   data_month[which(data_date==(unique(data_date))[t])] = rank(unique(data_date))[t]
 }



n = length(unique(data$head))


Y = vector("list", T)
for ( t in 1:T){
  Y[[t]] = matrix(rep(0,n^2),nrow=n)
}
Ind = Y



for (k in 1:length(data$edge.id)){
  Y[[data_month[k]]][data$tail[k],data$head[k]] = 1
}

detect_edges_gamma = NULL

detect_edges_inv = NULL

auc_gamma = NULL
auc_inv_gamma =NULL

length_out = 10

for(missing_prob in  seq(0.01,0.1,length.out =length_out)){

for (t in 1:T){
  Ind[[t]] = matrix(rep(0,n^2),nrow = n)
  for (i in 1:n){
    for (j in 1:n){
      if (j>i){
        Ind[[t]][i,j] = sample(x=c(0,1),size=1, prob =c(missing_prob,1-missing_prob))
      }
    }
  }
}


for ( t in 1:T){
  Y[[t]] = t(Y[[t]])+Y[[t]]
  Y[[t]][Y[[t]]>1] = 1
  Ind[[t]] = t(Ind[[t]])+Ind[[t]]
  Ind[[t]][Ind[[t]]>1]=1
}

d = 4

MF_list = mix_DN_adaptive (Y,Ind=Ind)

MF_invg_list = mix_DN_adaptive_invgamma(Y,Ind=Ind)



pred_gamma = NULL
res = NULL
pred_inv_gamma =NULL

for (t in 1:T){
  for (i in 1:n){
    for (j in 1:n){
      if (Ind[[t]][i,j]==0 && j<i){
        pred_gamma =c(pred_gamma, 1/(1+exp(-MF_list$mean_beta-t( MF_list$Mean_X[[t]][i,])%*% MF_list$Mean_X[[t]][j,])))
        pred_inv_gamma =c(pred_inv_gamma, 1/(1+exp(-MF_invg_list$mean_beta-t( MF_invg_list$Mean_X[[t]][i,])%*% MF_invg_list$Mean_X[[t]][j,])))
          res = c(res,Y[[t]][i,j])
      }
    }
  }
}

auc_gamma = c(auc_gamma ,as.numeric(auc(roc(res, pred_gamma,quiet=TRUE))))

auc_inv_gamma = c(auc_inv_gamma,as.numeric(auc(roc(res, pred_inv_gamma,quiet=TRUE))))


detect_edges_gamma=  c(detect_edges_gamma,sum(pred_gamma[res==1]>0.5)/sum(res==1))

detect_edges_inv = c(detect_edges_inv,sum(pred_inv_gamma[res==1]>0.5)/sum(res==1))

}


df_plot = data.frame(list(missing_prob=rep(seq(0.01,0.05,length.out =length_out),2), AUC = c(auc_gamma,auc_inv_gamma), 
   detection_ratio=c(detect_edges_gamma,detect_edges_inv), varaince_prior = c(rep('Gamma',length_out),rep('Inv_Gamma',length_out))))


p1 = ggplot(data=df_plot, aes(x=missing_prob, y=AUC, group=varaince_prior)) + theme(legend.position="none")+
  geom_line(size=1.5,aes(color = varaince_prior, linetype = varaince_prior)) + ylab("Test AUC")+
  ggtitle("AUC Comparison") + theme(plot.title = element_text(hjust = 0.5))
p2 =  ggplot(data=df_plot, aes(x=missing_prob, y=detection_ratio, group=varaince_prior)) + theme(legend.position="right")+
  geom_line(size=1.5,aes(color = varaince_prior, linetype = varaince_prior)) + ylab("Sucessfully detection ratio")+
  ggtitle("Detection Comparison") + theme(plot.title = element_text(hjust = 0.5))
grid.arrange(p1, p2, nrow = 1,widths = c(0.42,0.58))
