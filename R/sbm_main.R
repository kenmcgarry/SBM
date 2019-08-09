# sbm_main.R
# commenced 27/07/19
# using the Kolacyzk and Csardi book.  # Minimum Kullback-Leibler (MKL) positions
# https://cran.r-project.org/web/packages/heuristica/vignettes/README.html

library(mixer)
library(dplyr)
library(sand)
library(network)
library(sna)
library(ergm)
library(latentnet)
library(heuristica)
library(igraph)
library(Matrix)
library(UserNetR)

help(package="UserNetR")  # useful list of data sets from Doug Luke.

####################################################################
# French political bog dataset
fblog <- upgrade_graph(fblog)
fblog.sbm <- mixer(as.matrix(get.adjacency(fblog)),qmin=2, qmax=15)
fblog.sbm.output <- getModel(fblog.sbm)
names(fblog.sbm.output)

fblog.sbm.output$q
fblog.sbm.output$alphas
fblog.sbm.output$Taus

my.ent <- function(x){-sum(x*log(x,2))}
apply(fblog.sbm.output$Taus[,1:3], 2, my.ent)
log(fblog.sbm.output$q, 2)
summary(apply(fblog.sbm.output$Taus, 2, my.ent))
plot(fblog.sbm, classes=as.factor(V(fblog)$PolParty))

###################################################################
# Lazega lawyer data set
lazega <- upgrade_graph(lazega)
summary(lazega)
lazega.sbm <- mixer(as.matrix(get.adjacency(lazega)),qmin=2, qmax=15)
lazega.sbm.output <- getModel(lazega.sbm)
names(lazega.sbm.output)

lazega.sbm.output$q
lazega.sbm.output$alphas
lazega.sbm.output$Taus

my.ent <- function(x){-sum(x*log(x,2))}
apply(lazega.sbm.output$Taus[,1:3], 2, my.ent)
log(lazega.sbm.output$q, 2)
summary(apply(lazega.sbm.output$Taus, 2, my.ent))
plot(lazega.sbm, classes=as.factor(V(lazega)$Practice))

###################################################################
# Zachary data set
data(karate)
zachary <- upgrade_graph(karate)
summary(zachary)
zachary.sbm <- mixer(as.matrix(get.adjacency(zachary)),qmin=2, qmax=15)
zachary.sbm.output <- getModel(zachary.sbm)
names(zachary.sbm.output)

zachary.sbm.output$q
zachary.sbm.output$RR
zachary.sbm.output$Taus

my.ent <- function(x){-sum(x*log(x,2))}
apply(zachary.sbm.output$Taus[,1:3], 2, my.ent)
log(zachary.sbm.output$q, 2)
summary(apply(zachary.sbm.output$Taus, 2, my.ent))
plot(zachary.sbm, classes=as.factor(V(zachary)$Faction))

######################################################################
# 


######################################################################
# 


######################################################################
# 


######################################################################
# 




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












