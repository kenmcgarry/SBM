# sbm_functions.R

# Calculate some statistics for complex network, pass an igraph object 
# and return a dataframe of results.
graphstats <- function(gt) {
  #modularity <-modularity(gt, membership(cluster_walktrap(gt)))
  #nedges <- ecount(gt)
  #nverts <- vcount(gt)
  #transit <- transitivity(gt)
  degree <- igraph::degree(gt)
  close <- igraph::closeness(gt)
  between <- igraph::betweenness(gt,directed=FALSE)
  hubness <- igraph::hub_score(gt)$vector
  gs <- data.frame(degree, close, between, hubness)
  return(gs)
}

# calculate entropy
my.ent <- function(x){-sum(x*log(x,2))}

# produce a power law graph of degree
# needs powerLaw libraried in.
plot_power <- function(gs){
  m <- displ$new(gs$degree)
  ##Estimate the cut-off
  estimate_xmin(m)
  m$setXmin(105); m$setPars(2.644)
  
  temp <- plot(m,xlab="Degree",ylab="Percentiles") # only use this to grab coordinates for ggplot2
  
  ggplot(data=temp,aes(x=x,y=y*100))+
    geom_point()+
    geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)))+
    labs(title="",x="Degree distribution",y="Percentiles")+
    theme(axis.text.x=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.title.y = element_text(color="black", size=14, face="bold"))+
    theme(axis.title.x = element_text(color="black", size=14, face="bold"))+
    geom_vline(xintercept = 5, linetype="dashed", color = "red", size=1)
  
  cat("\nPercentiles",quantile(gs$degree, c(.70, .80, .90)) )
  
}

# ggplot2 version of power law plot
plot_power2 <- function(gs){
  
  G.degrees <- gs$degree
  G.degree.histogram <- as.data.frame(table(G.degrees))
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
    geom_point(alpha = 0.5, color = "blue",size = 3)+
    #labs(title="",x="Degree distribution",y="Percentiles")+
    theme(axis.text.x=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.title.y = element_text(color="black", size=14, face="bold"))+
    theme(axis.title.x = element_text(color="black", size=14, face="bold"))+
    
    scale_x_continuous("Degree",expand = c(0,0),
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16),
                       trans = "log10") +
    scale_y_continuous("hub Frequency",expand=c(0,0),
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,20),
                       trans = "log10") +
    ggtitle("Degree Distribution (log-log)") +
    geom_segment(aes(x = 1, y = 6, xend = 3, yend = 6),color="red",linetype="dashed",size=1)  + # Horiz
    geom_segment(aes(x = 3, y = 1, xend = 3, yend = 6),color="red",linetype="dashed",size=1) # Vert
}


############################################################
# Plot the Pis matrix and alphas vector using spectral decomposition
plotparam <- function(Pis,alphas,q=NULL){
  q <- length(alphas)
  if (q==1) {D <- list(vector=data.frame(1,1)); a <- b <- 1} else {
    if (q==2) {a<-b<-1} else {a<-2; b<-3}
    D <- colSums(Pis)
    L <- diag(rep(1,q)) -  diag(D^(-1/2)) %*% Pis %*% diag(D^(-1/2))
    D <- eigen(L)
  }
  
  plot(D$vector[,a],D$vector[,b],
       cex=1/min(alphas^(1/2))*alphas^(1/2)*3,
       axes=FALSE,xlab="",ylab="",
       main="Spectral view of the connection matrix",pch=19,col="red")
  points(D$vector[,a],D$vector[,b],cex=1/min(alphas^(1/2))*alphas^(1/2)*3)
  
  text(D$vector[,a],D$vector[,b],label=1:q)  
  # plot arrows
  gplot((Pis>median(Pis))*Pis,vertex.cex=1/min(alphas^(1/2))*alphas^(1/2)*3,edge.lwd=(Pis>median(Pis))*Pis*1/min(median(Pis)),label=1:length(alphas),label.pos=6)
}



############################################################
# Plot the icl criterion

ploticl <- function(x,q,...)
{
  if (x$method == "bayesian" ){
    title = "Bayesian criterion vs class number"
    y.lab = "Bayesian criterion"
  } else {
    title = "Integrated Classification Likelihood"
    y.lab = "ICL"
  }
  Q<-unlist(lapply(x$output,ICL<-function(x) length(x$alphas)))
  ICL<-unlist(lapply(x$output,ICL<-function(x) x$criterion))
  plot(Q,ICL,xlab="Number of classes",ylab=y.lab,main=title)
  lines(Q,ICL)
  abline(v=q,col="red",lty=2)
}


############################################################
# Plot the reorganized adjacency matrix

mixture <- function(x,alphas,lambdaq){
  fx<-0; for (q in 1:length(alphas)) {
    fx<-fx+alphas[q]*dpois(x,lambda=lambdaq[q])
  }
  return(fx)
}

plotmixture<-function(degrees,Pis,alphas,n, directed=FALSE){
  if( directed )
    colSums(Pis*alphas)*(2*n-2)->lambdaq
  else
    colSums(Pis*alphas)*(n-1)->lambdaq
  
  # Remove unconnected nodes
  degrees <- degrees[ which( degrees != 0) ]
  min(degrees):max(degrees)->x
  mixture(x,alphas,lambdaq)->y
  histo<-hist(degrees,plot=FALSE)
  plot(histo,ylim=c(0,max(histo$density,y)),freq=FALSE,col=7,main="Degree distribution",)
  lines(x,y,lwd=2,col="blue")
  points(x,y)
}

##############################################################
#  Spectral Clustering using normalized Laplacian
##############################################################
spectralkmeans <- function(x,q=2){
  #INPUT:
  #    x is an adjacency matrix
  #OUTPUT:
  #    An object of class "kmeans" which is a list with components:
  n<-dim(x)[1]
  D<-colSums(x)
  L<-diag(rep(1,n)) -  diag(D^(-1/2))%*% x %*% diag(D^(-1/2))
  eigen(L)->D
  kmeans(as.matrix(D$vectors[,max(1,(n-q)): (n-1)]),q) -> res         
}

##############################################################
#  Compute the rand index between two partition
##############################################################
randError<-function(x, y) {
  # function to calculate the adjusted rand statistic
  # x and y are vectors containing the two partitions to be compared
  # first, get crosstabs
  ctab <- table(x,y);
  
  # now calculate 4 intermediary sums
  cellsum <- sum(ctab*(ctab-1)/2)
  totsum <- sum(ctab)*(sum(ctab)-1)/2
  
  # use matrix multiplication to get row and column marginal sums
  rows <- ctab %*% rep(1,ncol(ctab))
  rowsum <- sum(rows*(rows-1)/2)
  cols <- rep(1,nrow(ctab)) %*% ctab
  colsum <- sum(cols*(cols-1)/2)
  # now put them together
  adj.rand <- (cellsum - (rowsum*colsum/totsum))/(.5*(rowsum +colsum)-(rowsum*colsum/totsum))
  return (adj.rand);
}

##########################################################
graph.affiliation<-function( n=100,
                             alphaVect=c(1/2,1/2), lambda=0.7, epsilon=0.05,
                             directed=FALSE) {
  # INPUT  n: number of vertex
  #           alphaVect : vecteur of class proportion
  #           lambda: proba of edge given  same classe
  #           epsilon: proba of edge given two different classes
  # OUTPUT x: adjacency matrix
  #              cluster: class vector
  #           
  
  x<-matrix(0,n,n);
  Q<-length(alphaVect);
  NodeToClass <- vector(length=n) 
  rmultinom(1, size=n, prob = alphaVect)->nq;
  Z<-class.ind(rep(1:Q,nq));
  Z<-Z[sample(1:n,n),];
  for (i in 1:n) {
    NodeToClass[i] <- which.max( Z[i,] )
  }
  for (i in 1:n) {
    if ( i != n) {
      for (j in (i+1):n) {
        # if i and j in same class
        if ( NodeToClass[i] ==  NodeToClass[j]) p<-lambda else  p<-epsilon
        if ( (rbinom(1,1,p) )) { x[i,j] <- 1 }
      }
      if ( directed ) {
        if ( i != 1) {
          for (j in 1:(i-1)) {
            if ( NodeToClass[i] ==  NodeToClass[j]) p<-lambda else  p<-epsilon
            if ( (rbinom(1,1,p) )) { x[i,j] <- 1 }
          }
        }
      }
    }
  }
  if ( ! directed ) {
    x <- x + t(x)
  }
  return(list(x=x,cluster=apply(Z,1,which.max)) )   
}
################################
class.ind<-function (cl)
{ 
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)))
  x[(1:n) + n * (unclass(cl) - 1)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  x
}
#################################

subgraph_density <-function(graph, gvert){
  graph %>% 
  igraph::induced_subgraph(gvert) %>%
  density()
}

##########################################################################
#' Calculate AICc (Akaike Information Criterion). Lower values of AICc indicate some 
#' combination of better fit to the data and more parsimony in the model (fewer free parameters).
#' @param LnL The log-likelihood (typically negative, but may not be for continuous data).
#' @param numparams The number of parameters for each model.
#' @param samplesize The number of data on which the model conferred likelihood.
#' @return \code{AICcval} A vector of AICc results.
#' @export
#' @seealso \code{\link{calc_AICc_column}}, \code{\link{calc_AICc_column}}
getAICc <- function(LnL, numparams, samplesize)
{
 if (numparams >= samplesize)
  {
    stop("ERROR!  You cannot have more parameters than samples in AICc, you get bizarre results.")
  }
  
  # Calculate AIC
  AICval = 2*numparams - 2*LnL
  
  # Correction for finite sample size
  correction_val = (2*numparams*(numparams+1)) / (samplesize - numparams - 1)
  AICc_val = AICval + correction_val
  
  return(AICc_val)
}



#####################################################
# The Girvan-Newman benchmarks can already be generated through igraph_preference_game, 
# or more generally a block-model. Hence, you can also generate clusters of varying sizes. 
# The only difference with LFR is that the degree distribution will be something like 
# a mixture of Poisson distributions and cannot be separately specified [Vincent Traag].
# However, here we try different approach.

girvan_benchmark <- function(g){
# let's see if we have communities here using the Girvan-Newman algorithm
# 1st we calculate the edge betweenness, merges, etc...
  ebc <- igraph::edge.betweenness.community(g, directed=F)

# Now we have the merges/splits and we need to calculate the modularity
# for each merge for this we'll use a function that for each edge
# removed will create a second graph, check for its membership and use
# that membership to calculate the modularity
  mods <- sapply(0:ecount(g), function(i){
  g2 <- igraph::delete.edges(g, ebc$removed.edges[seq(length=i)])
  cl <- igraph::clusters(g2)$membership
  igraph::modularity(g,cl)
})

# we can now plot all modularities
  plot(mods, pch=20, type="line", 
       col="blue",panel.first = grid(),ylab="Adjusted Rand Index",xlab="External degree")
}



