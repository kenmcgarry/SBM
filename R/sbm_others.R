# sbm_others.R
# detect communities on data using existing methods, walktrap, louvain etc etc.
# Compare with SBM+ method.

library(igraph)

######## Zachary data ##################
karate <- make_graph("Zachary")
nz <- graphstats(karate)

wcz <- igraph::cluster_walktrap(karate)
lcz <- igraph::cluster_louvain(karate)
scz <- igraph::cluster_spinglass(karate, spins=2)
trz <- igraph::transitivity(karate)
modz <- igraph::modularity(karate,membership(wcz))
sbmz <- sbm_nmf(karate)
modularity(wcz)

######## French political blog (Fblog) data ###################
fblog <- upgrade_graph(fblog)  
nf <- graphstats(fblog)
wcf <- igraph::cluster_walktrap(fblog)
lcf <- igraph::cluster_louvain(fblog)
scf <- igraph::cluster_spinglass(fblog, spins=2)
trf <- igraph::transitivity(fblog)
modf <- igraph::modularity(fblog,membership(wcf))
sbmf <- sbm_nmf(fblog)
  
######## Lageza data #################
lazega <- upgrade_graph(lazega)  # igraph object
Isolated <- which(degree(lazega)==0)# isolates present,these screw up, so remove 
lazega <- igraph::delete.vertices(lazega, Isolated)

lf  <- graphstats(lazega)
wcl <- igraph::cluster_walktrap(lazega)
lcl <- igraph::cluster_louvain(lazega)
scl <- igraph::cluster_spinglass(lazega, spins=2)
trl <- igraph::transitivity(lazega)
modl <- igraph::modularity(lazega,membership(wcl))
sbml <- sbm_nmf(lazega)




########################################################################

op <- par(mfrow = c(3,2),
          tcl=-0.5, 
          mai=c(0.3,0.3,0.3,0.3),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)

#par(mfrow=c(3,2))
par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(2,2,2,2))

plot(wc, karate, main="walktrap")
plot(lc, karate, main="louvain")
plot(sc, karate, main="spinglass")
#plot(fc, karate, main="fast-greedy")
#plot(cl, karate, main="cluster label")
plot(sc, karate, main="SBM+")

par(op)

#membership(fc)
sizes(fc)
sizes(sc)

sizes(lc)
membership(lc)

modularity(lc)
is_hierarchical(lc)
igraph::compare(wc, lc, method="nmi")





