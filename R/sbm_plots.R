# sbm_plots.R

# produces some of the various plots and tables used in the paper
dendrogram <- cluster_edge_betweenness(zachary)
V(zachary)[Faction == 1]$shape <- "circle"
V(zachary)[Faction == 2]$shape <- "square"
plot(dendrogram, zachary)

##### multiple plots in one diagram #######################################
op <- par(mfrow = c(3,2),par(mar=c(5,4,4,2)))
x <- seq(-pi,pi,0.1)
op

#### 1 #################################################################
plot(x, sin(x),
     main="Girvan-Newman plots for Zachary data",
     ylab="normalized mutual information", xlab="External degree Kout", type="l", col="blue")
lines(x,cos(x), col="red")
lines(x+.2,sin(x), col="black")
legend("topleft",
       c("SBM+","SpinGlass","Louvain"),
       fill=c("blue","red","black"))
#### 2
plot(x, sin(x),
     main="Girvan-Newman plots for Fblog data",
     ylab="normalized mutual information", xlab="External degree Kout",  type="l", col="blue")
lines(x,cos(x), col="red")
lines(x+.2,sin(x), col="black")
legend("topleft",
       c("SBM+","SpinGlass","Louvain"),
       fill=c("blue","red","black"))
#### 3
plot(x, sin(x),
     main="Girvan-Newman plots for Lageza data",
     ylab="normalized mutual information", xlab="External degree Kout", type="l", col="blue")
lines(x,cos(x), col="red")
lines(x+.2,sin(x), col="black")
legend("topleft",
       c("SBM+","SpinGlass","Louvain"),
       fill=c("blue","red","black"))
#### 4
plot(x, sin(x),
     main="Girvan-Newman plots for Lazega data",
     ylab="normalized mutual information", xlab="External degree Kout", type="l", col="blue")
lines(x,cos(x), col="red")
lines(x+.2,sin(x), col="black")
legend("topleft",
       c("SBM+","SpinGlass","Louvain"),
       fill=c("blue","red","black"))
#### 5
plot(x, sin(x),
     main="Girvan-Newman plots for French blog data",
     ylab="normalized mutual information", xlab="External degree Kout", type="l", col="blue")
lines(x,cos(x), col="red")
lines(x+.2,sin(x), col="black")
legend("topleft",
       c("SBM+","SpinGlass","Louvain"),
       fill=c("blue","red","black"))
#### 6
plot(x, sin(x),
     main="Girvan-Newman plots for Yeast data",
     ylab="normalized mutual information", xlab="External degree Kout", type="l", col="blue")
lines(x,cos(x), col="red")
lines(x+.2,sin(x), col="black")
legend("topleft",
       c("SBM+","SpinGlass","Louvain"),
       fill=c("blue","red","black"))

#################################################


