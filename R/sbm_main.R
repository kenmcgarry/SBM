# sbm_main.R
# Commenced 27/07/19, including LaTex file.
# Using the Kolacyzk and Csardi book.  # Minimum Kullback-Leibler (MKL) positions
# https://cran.r-project.org/web/packages/heuristica/vignettes/README.html
# https://www.smartcat.io/blog/2017/fast-matrix-factorization-in-r/
# https://blog.acolyer.org/2019/02/18/the-why-and-how-of-nonnegative-matrix-factorization/
# https://blog.acolyer.org/2019/02/13/beyond-news-contents-the-role-of-social-context-for-fake-news-detection/

library(mixer)
library(dplyr)
library(sand)
library(network)
library(ggplot2)
#library(sna)
#library(ergm)
#library(latentnet)
library(igraph)
library(Matrix)
library(UserNetR)
#library(linkcomm)
library(data.table)
library(poweRlaw)
library(xtable)
library(gplots)
library(NMF)
library(blockmatrix)
library(blockmodeling)

setwd("C:/common_laptop/R-files/sbm")
source("sbm_functions.R")
source("sbm_buildmodels.R")
source("sbm_matrix.R")
source("sbm_others.R")  # compare/contrast with other methods


##################### MATRIX FACTORIZATION ########################
# integrate the complex network and the SBM
# https://cran.r-project.org/web/packages/blocksdesign/vignettes/design_Vignette.pdf


##### blockmatrix package examples ####

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
n <- 15
net <- matrix(NA, ncol = n, nrow = n)
clu <- rep(1:2, times = c(5, 10))
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

# generate a synthetic dataset with known classes: 20 features, 23 samples (10+5+8)
n <- 20; counts <- c(10, 5, 8);
p <- sum(counts)
x <- NMF::syntheticNMF(n, counts)
dim(x)

# build the true cluster membership
groups <- unlist(mapply(rep, seq(counts), counts))
# run on a data.frame
res <- NMF::nmf(data.frame(x), 3)
# missing method: use algorithm suitable for seed
res <- NMF::nmf(x, 2, seed=rnmf(2, x))
NMF::algorithm(res)

# compare some NMF algorithms (tracking the approximation error)
res <- NMF::nmf(x, 2, list('brunet', 'lee', 'nsNMF'), .options='t')
plot(res)

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
plot(res)


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












