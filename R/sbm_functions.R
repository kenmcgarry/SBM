# sbm_functions.R

# Calculate some statistics for complex network, pass an igraph object 
# and return a dataframe of results.
graphstats <- function(gt) {
  #modularity <-modularity(gt, membership(cluster_walktrap(gt)))
  #nedges <- ecount(gt)
  #nverts <- vcount(gt)
  #transit <- transitivity(gt)
  degree <- igraph::degree(gt)
  close <- closeness(gt)
  between <- betweenness(gt,directed=FALSE)
  hubness <-hub_score(gt)$vector
  gstats <- data.frame(degree, close, between, hubness)
  return(gstats)
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





