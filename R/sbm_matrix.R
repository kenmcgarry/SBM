# sbm_matrix.R
# For each data set, get the complex network statistics:
# 1. create r x r diagonal matrices for eigen(hubness), closeness, betweeness & degree
# 2. integrate the three matrices using linear algebra combinations 

# https://towardsdatascience.com/prototyping-a-recommender-system-step-by-step-part-2-alternating-least-square-als-matrix-4a76c58714a1

### Do it ###
## zachary complex net ##
gstats <- graphstats(zachary)
gstats[,5] <- row.names(gstats)
head(gstats,10)
by_hub <- gstats %>% arrange(desc(hubness))
xtable(by_hub)  # make table for the paper

## Fblog complex net ##
gstats <- graphstats(fblog)
gstats[,5] <- row.names(gstats)
head(gstats,10)
by_hub <- gstats %>% arrange(desc(hubness))
xtable(by_hub) # make table for the paper

## Lazega net
gstats <- graphstats(lazega)
gstats[,5] <- row.names(gstats)
head(gstats,10)
by_hub <- gstats %>% arrange(desc(hubness))
xtable(by_hub) # make table for the paper



n <- length(V(zachary)) # how many nodes?
A <- B <- C <- D <- E <- matrix(0, n, n) # create the matrices
A <- get.adjacency(zachary)
B <- gstats$between
C <- gstats$close
D <- gstats$degree
E <- gstats$hubness
E <- diag(E)

by_hub <- gstats %>% group_by(hubness)
by_hub <- gstats %>% arrange(desc(hubness))
xtable(by_hub)  # make table for the paper

## zachary SBM ##
# internals of sbm
zachary.sbm.output$q                # No of classes (n)
zachary.sbm.output$criterion        # ICL criterion used for model selection
zachary.sbm.output$Taus             # matrix of posterior probabilities (n x actors)
zachary.sbm.output$alphas           # vector of proportion for each class n
zachary.sbm.output$Pis              # class connectivity matrix (n x n)

# mergeclass <- function(Pis,alphas,q=NULL){
  q <- length(zachary.sbm.output$alphas)
  if (q==1) {D <- list(vector=data.frame(1,1)); a <- b <- 1} else {
    if (q==2) {a<-b<-1} else {a<-2; b<-3}
    D <- colSums(zachary.sbm.output$Pis)
    L <- diag(rep(1,q)) -  diag(D^(-1/2)) %*% zachary.sbm.output$Pis %*% diag(D^(-1/2))
    D <- eigen(L)  # D$vectors
    D$vectors  # plot(D$vector[,a],D$vector[,b])
  }

# class merging if necessary
q <- length(zachary.sbm.output$alphas)
q <- length(lazega.sbm.output$alphas)
q <- length(fblog.sbm.output$alphas)

U <- zachary.sbm.output$alphas
Z <- zachary.sbm.output$Pis  
X <- zachary.sbm.output$Taus

Z <- lazega.sbm.output$Pis  
X <- lazega.sbm.output$Taus
Z <- fblog.sbm.output$Pis
X <- fblog.sbm.output$Taus
  
X[X < .5] <- 0
X[X >=.5] <- 1 
rowSums(X)

for(i in 1:q){
  #cat("\n",which(Z[,i] > 0.5 )," ")
  merge <- which(Z[,i] > 0.7 )
  if(length(merge)>1){
    cat("\n merge classes ",merge)} else{cat("\n Nothing to merge")}
}

res <- NMF::nmf(X,3)   # version 1
res <- NMF::nmf(X)   # version 2
W <- basis(res)
H <- coef(res)
W <- res$W
H <- res$H
Y <- W%*%H   # should restore X (more or less)
D <- eigen(Y)

Z <- A %*% E

ev <- eigen(A)
K <- ev$values
V <- ev$vectors
J <- V %*% diag(K) %*% t(V)  # Matrix factorization,matrix A can be represented as the product.

Z <- U %*% diag(K) %*% t(U)  # test changes here
diag(Z)

### igraph csardi version ###
library(NMF)
basis <- matrix( c(1,2,3,4, 4,3,2,1),4,2)
A <- matrix(0,4,3)
A[,1] <- .9*basis[,1]+.1*basis[,2]
A[,2] <- .5*basis[,1]+.5*basis[,2]
A[,3] <- .1*basis[,1]+.9*basis[,2]
csNMF <- NMF::nmf(A, rank=2,method = 'brunet')
sum(abs(basis(csNMF) %*% coef(csNMF) - A ))              

# remember (Zitnik, 2013) paper for data integration
                                                                                                                                                

##############   NMF package  ##########################
# load the Golub data
data(esGolub)
X <- matrix(1:12,3,4)
res <- NMF::nmf(X,3)
W <- basis(res)
H <- coef(res)

# compute NMF for each method
# retrieve all the methods that have a secondary R version
meth <- NMF::nmfAlgorithm(version='R')
meth <- c(names(meth), meth)
meth
res <- NMF::nmf(X, 3, meth, seed=123456)

# nice plot of r estimation
estim.r <- NMF::nmf(esGolub, 2:6, nrun=10, seed=123456)
plot(estim.r)
estim.a <- NMF::nmf(A, 1:3,nrun=10, seed=123456)
plot(estim.a)
V.random <- randomize(A)
estim.a.random <- NMF::nmf(V.random, 1:5, nrun=10, seed=123456, method = 'lee')
plot(estim.a.random)
plot(estim.a,estim.a.random)

### better plot below ###
### plots only ####
# random graph example
g <- erdos.renyi.game(10, p=1/2) + erdos.renyi.game(10, p=1/2)
# Plot the adjacency matrix
A <- get.adjacency(g)
image(A)
# In this example, sorting the nodes by degree is not a good idea
i <- order( degree(g) )
image(A[i,i])
plot(g)

# Zachary pretty plots
z <- graph.famous("Zachary")
deg <- igraph::degree(z)
lay <- layout.fruchterman.reingold(z)
fine <- 500 # this will adjust the resolving power.
palette <- colorRampPalette(c('lightblue','pink'))
degCol <- palette(fine)[as.numeric(cut(deg,breaks = fine))]
plot(z, layout=lay, vertex.color=degCol, vertex.size=deg*1.5, vertex.label.cex=0.6, main="Degree centrality")

clos <- igraph::closeness(z)
# Plot the graph:
closCol = palette(fine)[as.numeric(cut(clos,breaks = fine))]
plot(z,layout = lay, vertex.color=closCol, vertex.size=clos*1500, vertex.label.cex=0.6, main="Closeness centrality")

betw <- igraph::betweenness(z)
#Plot the graph
betwCol = palette(fine)[as.numeric(cut(betw,breaks = fine))]
plot(z,layout = lay, vertex.color=betwCol, vertex.size=betw*0.2, vertex.label.cex=0.6, main="Betweeness centrality")

ev <- igraph::evcent(z)
ev <- igraph::evcent(z)$vector
# Produce the plot:
evCol = palette(fine)[as.numeric(cut(ev,breaks = fine))]
plot(z,layout = lay, vertex.size=ev*40, vertex.color=evCol, vertex.label.cex=0.6, main="Eigenvector centrality")

par(mfrow=c(2,2), mar = c(5,0,4,0), mai=c(0,0.1,0.1,0.5), plt=c(0.,0.6,0.2,0.9))

#op <- par(mfrow = c(2,2), omi=c(0.7,0,0,0), plt=c(0.1,0.6,0.1,0.9))
plot(z,layout = lay, vertex.color=degCol, vertex.size=deg*1.5, vertex.label.cex=0.8);
title(main="Degree centrality",cex.main=2);
plot(z,layout = lay, vertex.color=closCol, vertex.size=clos*1500, vertex.label.cex=0.8);
title(main="Closeness centrality",cex.main=2);
plot(z,layout = lay, vertex.color=betwCol, vertex.size=betw*0.2, vertex.label.cex=0.8);
title(main="Betweeness centrality",cex.main=2);
plot(z,layout = lay, vertex.size=ev*40, vertex.color=evCol, vertex.label.cex=0.8);
title(main="Eigenvector centrality",cex.main=2)

# next time simply plot four separate graphs (as above) and use a multiplot to join them!












