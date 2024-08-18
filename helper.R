# Some helper functions

library(pROC)
library(Bessel)
library(BayesLogit)

procrustes_r <- function(A, B, normalize=F ){
  # center and normalize A 
  # A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A, type = "F") 
  A.normalized <- A
  # 
  # # center and normalize B
  # B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B, type = "F")
  if (normalize == T){
    B.normalized <- B /B.size
  } else {B.normalized <- B }
  
  # Rotation matrix T 
  svd.results <- svd(B.normalized %*% t(A.normalized))
  U <- svd.results$u
  V <- svd.results$v
  T <- V %*% t(U)
  
  # B transformed
  B.transformed <- T %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.normalized,  type = "F")
  
  # Return
  return(list(rotation.mtx = T, B.transformed = B.transformed, RSS = RSS))
}


distance_squared_inner_prod = function(mu1,mu2,Sigma1,Sigma2){
  value = t(mu1)%*% mu2 %*% t(mu2) %*%mu1 +sum(diag(Sigma1 %*% Sigma2)) +  t(mu1)%*% Sigma2 %*%mu1 + t(mu2)%*% Sigma1 %*%mu2
  
  return(value)
}

positions_to_edges_Gaussian = function(Xt,beta,sigma){
  n = length(Xt[,1])
  edge_mat = matrix(rep(0,n*n),nrow=n)
  t=1
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {edge_mat[j,i] = rnorm(n=1,mean=beta+t(Xt[i,])%*%Xt[j,],sd=sigma)
      t= t+1}
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}

positions_to_edges = function(Xt,beta){
  n = length(Xt[,1])  
  edge_mat = matrix(rep(0,n*n),nrow=n)
  t=1
  for (i in 1:n){
    for(j in 1:n){
      if(j > i)
      {
        binary_mean = 1/(1+exp(-beta-t(Xt[i,])%*%Xt[j,]))
        edge_mat[j,i] = rbinom(n=1,size=1, prob = binary_mean)
        t= t+1}
    }
  }
  edge_mat = edge_mat+ t(edge_mat)
  return(edge_mat)
}


ggplot_dis= function(Xm,Ym){
  n0= length(Xm)
  data_fr = as.data.frame(cbind(c(1:n0),Xm,Ym))
  names(data_fr) = c("index","neg_inner_prod","connected")
  data_fr$connected = as.factor(data_fr$connected)
  levels(data_fr$connected) = c("unconnected","connected")
  gp = ggplot(data_fr,aes(x=index,y=neg_inner_prod))+geom_point(size=4,aes(shape=connected,color=connected))+ 
    theme(legend.position = "none")+theme(axis.text=element_text(size=16),
                                          axis.title=element_text(size=18,face="bold"))+theme(legend.title = element_blank())+theme(legend.text = element_text( size=16))
  return(gp)
}



plot_clus = function(Xm1,Y1,clus,t){
  n = length(Xm1[,1])
  diag(Y1) = 0
  Y1[lower.tri(Y1, diag = FALSE)] = 0
  NodeList <- data.frame((1:n), x=Xm1[,2] ,y=Xm1[,1])
  #  Edge_list = data.frame(matrix(ncol = 2, nrow = 0))
  Edge_list = igraph::as_data_frame(graph_from_adjacency_matrix(Y1,mode="undirected"),'edge')
  pal <- brewer.pal(4,"Accent")
  vertex.col = pal[factor(clus)]
  a<- graph_from_data_frame(vertices = NodeList, d= Edge_list, directed = FALSE)
  plot.igraph(a,vertex.color=vertex.col,vertex.shape=setdiff(shapes(), "")[clus], vertex.label=1:n,vertex.size=12,
              edge.width=0.1,axes=TRUE)
#  title(main=paste("Latent positions for time point",t,sep = ' '),cex.main=1.8)
#  legend('topright',legend=levels(factor(clus)),pch=c(1,0),col=vertex.col,bty = "n",cex=1.6)
  return(NULL)
}

