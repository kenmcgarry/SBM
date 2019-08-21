# sbm_matrix.R
# For each data set, get the complex network statistics:
# 1. create r x r diagonal matrices for eigen(hubness), closeness, betweeness & degree
# 2. integrate the three matrices using linear algebra combinations 
#
# For each data set, get the SBM statistics:
# 1.
# 2.
#
# For each data set, get the link communities
# 1.
# 2. 

# https://towardsdatascience.com/prototyping-a-recommender-system-step-by-step-part-2-alternating-least-square-als-matrix-4a76c58714a1

### Do it ###
## zachary complex net ##
gstats <- graphstats(zachary)
gstats[,5] <- row.names(gstats)
head(gstats)

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
res <- NNLM::nnmf(X,3)   # version 2
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
J <- V %*% diag(K) %*% t(V)  # Matrix factorization, the matrix A can 
                             # be represented as the product.

Z <- U %*% diag(K) %*% t(U)  # test changes here
diag(Z)

##################### MATRIX FACTORIZATION ########################
# integrate the complex network, the SBM and the link clustering network
# https://cran.r-project.org/web/packages/blocksdesign/vignettes/design_Vignette.pdf
#
library(blockmatrix)

A <- array(rnorm(9,mean=1),c(3,3))
B <- array(rnorm(9,mean=2),c(3,3))
C <- 0
D <- array(rnorm(9,mean=4),c(3,3))
F <- array(rnorm(9,mean=10),c(3,3))
M <- blockmatrix(names=c("A","0","D","0"),A=A,D=D,dim=c(2,2))
E <- blockmatrix(names=c("0","F","D","0"),F=F,D=D,dim=c(2,2))
R <- M+E
S <- solve(R)
P <- blockmatmult(R,E)
l <- list(A=A,B=B,C=C,D=D,F=F)
mv <- array(c("A","B","C","D","F","F"),c(3,2))
BB <- blockmatrix(value=mv,list=l)


A <- array(rnorm(9,mean=1),c(3,3))
B <- 0 #array(rnorm(9,mean=2),c(3,3))
C <- 0
D <- array(rnorm(9,mean=4),c(3,3))
F <- array(rnorm(9,mean=10),c(3,3))
M <- blockmatrix(names=c("A","0","D","0"),A=A,D=D,dim=c(2,2))
E <- blockmatrix(names=c("0","F","D","0"),F=F,D=D,dim=c(2,2))
E[,1] <- M[,1]

##### blockmodeling package #########
library(blockmodeling)
n <- 20
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(5, 15))
tclu <- table(clu)
net[clu == 1, clu == 1] <- rnorm(n = tclu[1] * tclu[1], mean = 0, sd = 1)
net[clu == 1, clu == 2] <- rnorm(n = tclu[1] * tclu[2], mean = 4, sd = 1)
net[clu == 2, clu == 1] <- rnorm(n = tclu[2] * tclu[1], mean = 0, sd = 1)
net[clu == 2, clu == 2] <- rnorm(n = tclu[2] * tclu[2], mean = 0, sd = 1)

# Ploting the network
plotMat(M = net, clu = clu, print.digits.cells = 3)
class(net) <- "mat"
plot(net, clu = clu)

# We select a random partition  and then optimize  it
all.par <- nkpartitions(n = n, k = length(tclu))# Forming the partitions
all.par <- lapply(apply(all.par, 1, list), function(x) x[[1]])# Optimizing one partition
res <- optParC(M = net,clu = all.par[[sample(1:length(all.par), size = 1)]],
               approaches = "hom", homFun = "ss" , blocks = "com")
plot(res) # Hopefully we get the original partition

# Optimizing 10 random partitions with optRandomParC
res <- optRandomParC(M = net, k = 2, rep = 10,approaches = "hom", homFun = "ss", 
                     blocks = "com")
plot(res) # Hopefully we get the original partition# Using indirect approach - 
# structural equivalence
D <- sedist(M = net)
plot.mat(net, clu = cutree(hclust(d = D, method = "ward.D"), k = 2))

############ NMFN package #############
#library(NMFN)
library(NNLM)

X <- matrix(1:16,4,4)
z.mm   <- NNLM::nnmf(X,3)             # 3 factors via multiplicative update
z.als  <- NNLM::nnmf(X,2,'nnmf_als')  # 3 factors via alternating least square
z.prob <- NNLM::nnmf(X,3,'nnmf_prob') # 3 factors via multinomial

k <- 15;
init <-list(W =matrix(runif(nrow(nsclc)*k), ncol = k),
            H =matrix(runif(ncol(nsclc)*k), nrow = k))


##############   NMF package  ##########################
library(NMF)
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


############# another package #############################



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












