#Create segregated and non-segregated neighborhoods, using binary and cluster data-generating processes

#Preamble
library(ggplot2)
library(tidyr)
library(plyr);library(dplyr, warn.conflicts = F)
library(tictoc)
rm(list=ls())
s = function(x){summary(factor(x))}
####################################################################################################

####################################################################################################
#Set global parameters
xlim = 100 #neighborhood size x. call this meters. or km if this is city scale??
ylim = 100 #neighborhood size y. call this meters. or km if this is city scale??
ngrid = 10 #number of grid spaces
####################################################################################################

####################################################################################################
#Functions to create neighborhoods (Binary DGP)

#Nonsegregated neighborhoods: just a uniform distribution
create_Nonseg = function(n, xlim, ylim, pminority){
  fr = data.frame(long = runif(n,0,xlim),
                  lat = runif(n,0,xlim),
                  minor = rbinom(n,1,pminority) == 1)
  return(fr)
}

#Seg neighborhoods: Assume that majority is uniform across whole space
create_Seg = function(n, xlim, ylim, pminority, conc){ #concentration parameter
  fr = data.frame(minor = rbinom(n,1,pminority) == 1,
                  long = NA,
                  lat = NA,
                  ghe = rbinom(n,1,conc)) #concentration parm: proportion of each group on the "designated" side
  #lat longs for majority (not concentrated), ie ghe = F
  fr$long[fr$minor == F & fr$ghe == F] = runif(sum(fr$minor == F & fr$ghe == F), 0, xlim) #x: runif from 0 to xlim
  fr$lat[fr$minor == F & fr$ghe == F] = runif(sum(fr$minor == F & fr$ghe == F), 0, ylim) #y: runif from 0 to ylim
  #lat longs for majority (concentrated)
  fr$long[fr$minor == F & fr$ghe == T] = runif(sum(fr$minor == F & fr$ghe == T), min = 0, #x: runif from 0 to limit
                                               max = (1 - pminority)*xlim) #concentration limit: proportional to size of majority
  fr$lat[fr$minor == F & fr$ghe == T] = runif(sum(fr$minor == F & fr$ghe == T), 0, ylim) #y: runif from 0 to xlim
  #lat longs for minority (not concentrated)
  fr$long[fr$minor == T & fr$ghe == F] = runif(sum(fr$minor == T & fr$ghe == F), 0, xlim) #x: runif from 0 to xlim
  fr$lat[fr$minor == T & fr$ghe == F] = runif(sum(fr$minor == T & fr$ghe == F), 0, ylim) #y: runif from 0 to ylim
  #lat longs for minority (concentrated)
  fr$long[fr$minor == T & fr$ghe == T] = runif(sum(fr$minor == T & fr$ghe == T),  #x: runif from edge to xlim
                                               min = (1-pminority)*xlim, max = xlim) #concentration limit: prop. to size of minority
  fr$lat[fr$minor == T & fr$ghe == T] = runif(sum(fr$minor == T & fr$ghe == T), 0, ylim) #y: runif from 0 to ylim
  return(fr)
}

#make simulated neighborhood and plot: Segregated
d = create_Seg(n = 200, xlim = 100, ylim = 100, pminority = 0.4, conc = 0.9) %>% rename(Ethnicity = minor) %>%
  mutate(Ethnicity = mapvalues(Ethnicity, from = c('TRUE', 'FALSE'), to = c('Minority','Majority')))
ggplot(data = d) + geom_point(aes(x = long, y = lat, shape = Ethnicity), size = 2) + theme_bw() +
  ggtitle('Binary segregation example') + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(1, 2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
#path0 = paste0('insert_path_here',Sys.Date(),'/')
#dir.create(path0)
#ggsave(filename = paste0(path0, 'binary_seg_example.png'), height = 150, width = 150, units = 'mm'); rm(path0)

#make simulated neighborhood and plot: Non-Segregated
d = create_Nonseg(n = 200, xlim = 100, ylim = 100, pminority = 0.4) %>% rename(Ethnicity = minor) %>%
  mutate(Ethnicity = mapvalues(Ethnicity, from = c('TRUE', 'FALSE'), to = c('Minority','Majority')))
ggplot(data = d) + geom_point(aes(x = long, y = lat, shape = Ethnicity), size = 2) + theme_bw() +
  ggtitle('Binary non-segregated example') + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(1, 2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
#path0 = paste0('insert_path_here',Sys.Date(),'/')
#dir.create(path0)
#ggsave(filename = paste0(path0, 'binary_seg_example.png'), height = 150, width = 150, units = 'mm'); rm(path0)
####################################################################################################

####################################################################################################
#Functions to create neighborhoods (Clustered DGP)
library(spatstat)

#Function to generate segregated clusters
# multitype Neyman-Scott process (each cluster is a multitype process)
#[This is same as Matern cluster process, only with marks]
#see for more info: https://rdrr.io/cran/spatstat/man/rPoissonCluster.html
nclust2 <- function(x0, y0, radius, n, types=c("yes", "no"), p) {
  X <- runifdisc(n, radius, centre=c(x0, y0))
  M <- sample(types, 1, replace=TRUE, prob = c(p, 1-p)) %>% rep(n)#this gives everybody in the cluster the same mark
  marks(X) <- M
  return(X)
}

#Function to generate segregated clusters
#Nonsegregated version: Type is NOT determined cluster, ie each cluster can have minority and nonminority members
nclust1 <- function(x0, y0, radius, n) {
  X <- runifdisc(n, radius, centre=c(x0, y0))
  return(X)
}

ClusterSeg = function(n_samp = 100, p = 0.3, nclus = 13, ka = 0.05, r = 3){ #segregated neighborhood, cluster DGP
  kids = rPoissonCluster(kappa = ka, #density of clusters; determines expected number of clusters
                         expand = 10, #allows cluster `parent' to be outside main window
                         rcluster = nclust2, #choice of function to generate clusters
                         radius = r, #radius of clusters. smaller clusters -> higher segregation
                         n = nclus, #number of people per cluster
                         p = p, #minority proportion
                         win = owin(c(0,100),c(0,100)))
  
  kids_df = data.frame(long = kids$x, #data frame containing locations and group identity of simulated people
                       lat = kids$y,
                       minor = kids$marks) 
    ns = min(nrow(kids_df), n_samp)
    return(kids_df[sample(nrow(kids_df), ns), ]) #return sample of n rows
}

ClusterNonSeg = function(n_samp = 100, p = 0.3, nclus = 13, ka = 0.05, r = 3){ #non-segregated neighborhood, cluster DGP
  kids = rPoissonCluster(kappa = ka, #density of clusters; determines expected number of clusters
                         expand = 10, #allows cluster `parent' to be outside main window
                         rcluster = nclust1, #choice of function to generate clusters
                         radius = r, #radius of clusters. smaller clusters -> higher segregation
                         n = nclus, #number of people per cluster
                         win = owin(c(0,100),c(0,100)))
  
  marks(kids) = sample( c('yes', 'no'), kids$n, replace=TRUE, prob = c(p, 1-p)) #p is minority proportion
  
  kids_df = data.frame(long = kids$x, #data frame containing locations and group identity of simulated people
                       lat = kids$y,
                       minor = kids$marks) 
  ns = min(nrow(kids_df), n_samp)
  return(kids_df[sample(nrow(kids_df), ns), ]) #return sample of n rows
}

#Create a non-segregated cluster neighborhood and plot
d = ClusterNonSeg(ka = 0.002, nclus = 50, r = 20, p = 0.3, n_samp = 1000) %>% rename(Ethnicity = minor) %>%
  mutate(Ethnicity = mapvalues(Ethnicity, from = c('yes', 'no'), to = c('Minority','Majority')))

ggplot(data = d) + geom_point(aes(x = long, y = lat, shape = Ethnicity), size = 2) + theme_bw() +
  ggtitle('Cluster non-segregated\nexample') + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(1, 2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
#path0 = paste0('insert path',Sys.Date(),'/')
#dir.create(path0)
#ggsave(filename = paste0(path0, 'cluster_nonseg_example.png'), height = 150, width = 150, units = 'mm'); rm(path0)

#Create a segregated cluster neighborhood and plot
d = ClusterSeg(ka = 0.002, nclus = 50, r = 20, p = 0.3, n_samp = 1000) %>% rename(Ethnicity = minor) %>%
  mutate(Ethnicity = mapvalues(Ethnicity, from = c('yes', 'no'), to = c('Minority','Majority')))

ggplot(data = d) + geom_point(aes(x = long, y = lat, shape = Ethnicity), size = 2) + theme_bw() +
  ggtitle('Cluster segregation example') + theme(text = element_text(size=16), plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(1, 2)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
#path0 = paste0('insert path',Sys.Date(),'/')
#dir.create(path0)
#ggsave(filename = paste0(path0, 'cluster_seg_example.png'), height = 150, width = 150, units = 'mm'); rm(path0)
########################################################################################################################
