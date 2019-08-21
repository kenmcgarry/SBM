# sbm_others.R
# detect communities on data using existing methods, walktrap, louvain etc etc.
# Compare with SBM+ method.

library(igraph)

######## Zachary data ##################
karate <- make_graph("Zachary")
wc <- cluster_walktrap(karate)
modularity(wc)
membership(wc)
plot(wc, karate)


wc <- cluster_walktrap(karate)
modularity(wc)
membership(wc)
plot(wc, karate)

lc <- cluster_louvain(karate)
sc <- cluster_spinglass(karate, spins=1)
fc <- cluster_fast_greedy(karate)
cl <- cluster_label_prop(karate)

# plot communities
par(mfrow=c(2,2))
plot(wc, karate, main="walktrap")
plot(lc, karate, main="louvain")
plot(sc, karate, main="spinglass")
plot(fc, karate, main="fast-greedy")

#membership(fc)
sizes(fc)
sizes(lc)
sizes(sc)



