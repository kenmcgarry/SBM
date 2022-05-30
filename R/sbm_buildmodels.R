# sbm_buildmodels.R

library(mixer)
library(dplyr)
library(sand)
library(network)
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

help(package="UserNetR")  # useful list of data sets from Doug Luke.

###############################################################################################
# random-small
randsmall <- sample_gnm(34,78,directed=FALSE, loops = FALSE) # random version of Zachary
# add groundtruth i.e. a class label
V(randsmall)$groundtruth <- sample(c(0, 1), vcount(randsmall), replace = TRUE, prob = c(0.5, 0.5))

randsmall.sbm <- mixer(as.matrix(get.adjacency(randsmall)),qmin=2, qmax=15)
randsmall.sbm.output <- getModel(randsmall.sbm)
names(randsmall.sbm.output)

randsmall.sbm.output$q
randsmall.sbm.output$alphas
randsmall.sbm.output$Taus

apply(randsmall.sbm.output$Taus[,1:3], 2, my.ent)
log(randsmall.sbm.output$q, 2)
summary(apply(randsmall.sbm.output$Taus, 2, my.ent))
plot(randsmall.sbm, classes=as.factor(V(randsmall)$groundtruth))

###############################################################################################
# random-large
randlarge <- sample_gnm(192,1431,directed=FALSE, loops = FALSE) # random version of french blog
# add groundtruth i.e. a class label  # unique(V(fblog)$PolParty)
V(randlarge)$groundtruth <- sample(unique(V(fblog)$PolParty),vcount(randlarge),
                                   replace=TRUE,prob = rep(0.11,9))

randlarge.sbm <- mixer(as.matrix(get.adjacency(randlarge)),qmin=2, qmax=15)
randlarge.sbm.output <- getModel(randlarge.sbm)
names(randlarge.sbm.output)

randlarge.sbm.output$q
randlarge.sbm.output$alphas
randlarge.sbm.output$Taus

apply(randlarge.sbm.output$Taus[,1:3], 2, my.ent)
log(randlarge.sbm.output$q, 2)
summary(apply(randlarge.sbm.output$Taus, 2, my.ent))
plot(randlarge.sbm, classes=as.factor(V(randlarge)$groundtruth))

###############################################################################################
# French political blog dataset
fblog <- upgrade_graph(fblog)  # igraph object
fblog.sbm <- mixer(as.matrix(get.adjacency(fblog)),qmin=2, qmax=15)
fblog.sbm.output <- getModel(fblog.sbm)
names(fblog.sbm.output)

fblog.sbm.output$q
fblog.sbm.output$alphas
fblog.sbm.output$Taus

apply(fblog.sbm.output$Taus[,1:3], 2, my.ent)
log(fblog.sbm.output$q, 2)
summary(apply(fblog.sbm.output$Taus, 2, my.ent))
plot(fblog.sbm, classes=as.factor(V(fblog)$PolParty))



###################################################################
# Lazega lawyer data set
lazega <- upgrade_graph(lazega)  # igraph object
summary(lazega)
lazega.sbm <- mixer(as.matrix(get.adjacency(lazega)),qmin=2, qmax=15)
lazega.sbm.output <- getModel(lazega.sbm)
names(lazega.sbm.output)

lazega.sbm.output$q
lazega.sbm.output$alphas
lazega.sbm.output$Taus

apply(lazega.sbm.output$Taus[,1:3], 2, my.ent)
log(lazega.sbm.output$q, 2)
summary(apply(lazega.sbm.output$Taus, 2, my.ent))
plot(lazega.sbm, classes=as.factor(V(lazega)$Practice))

###################################################################
# Zachary data set
data(karate)
zachary <- upgrade_graph(karate)  # igraph object
rm(karate)
summary(zachary)

# build SBM network
zachary.sbm <- mixer(as.matrix(get.adjacency(zachary)),qmin=2, qmax=15)
zachary.sbm.output <- getModel(zachary.sbm)
names(zachary.sbm.output)

apply(zachary.sbm.output$Taus[,1:3], 2, my.ent)
log(zachary.sbm.output$q, 2)
summary(apply(zachary.sbm.output$Taus, 2, my.ent))
plot(zachary.sbm, classes=as.factor(V(zachary)$Faction))

# build linked network
#edgelist <- as_edgelist(zachary, names = TRUE)  # get edgelist for linkcomm process
#zachary.link <-getLinkCommunities(edgelist)
#plot(zachary.link, type = "graph",layout = "spencer.circle")
#plot(zachary.link, type = "members")
#getAllNestedComm(zachary.link)
#plot(zachary.link, type = "graph", clusterids = c(5,12))
#zachr <- getClusterRelatedness(zachary.link, hcmethod = "ward.D")

#zachcc <- getCommunityCentrality(zachary.link)
#head(sort(zachcc, decreasing = TRUE))
#zachdense <- LinkDensities(zachary.link)
#zachcm <- getCommunityConnectedness(zachary.link)
#zachmod <- getCommunityConnectedness(zachary.link, conn = "mod")
#plot(zachary.link, type = "commsumm", summary = "modularity")

#zach_df <- data.frame(density=zachdense, connectivity=zachcm,modularity=zachmod)
#xtable(zach_df)

# the membership matrix needs to be told there are 34 nodes, otherwise it just defaults to 20
#zachmatrix <- getCommunityMatrix(zachary.link,nodes=head(names(zachary.link$numclusters),34))
#zachoverlap <- get.community.overlaps(zachary.link)
#plot(zachary.link, type = "members")
#summary(zachary.link)

######################################################################
# Currently too large for immdeiate SBM results
#data("yeast")
#yeast.sbm <- mixer(as.matrix(get.adjacency(yeast)),qmin=2, qmax=15)
#yeast.sbm.output <- getModel(yeast.sbm)
#names(yeast.sbm.output)

#yeast.sbm.output$q
#fblog.sbm.output$alphas
#fblog.sbm.output$Taus

#apply(fblog.sbm.output$Taus[,1:3], 2, my.ent)
#log(fblog.sbm.output$q, 2)
#summary(apply(fblog.sbm.output$Taus, 2, my.ent))
#plot(fblog.sbm, classes=as.factor(V(fblog)$PolParty))

######################################################################
# 


######################################################################
# 


######################################################################
# 








