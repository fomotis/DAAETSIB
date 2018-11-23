# small R script for debugging

library(flowCore) #loads the package flowCore
library(ggcyto)
library(flowQ)
library(flowViz)

source("C:/Users/Jurg/Dropbox/Doktorat/Projects/P3_Costa Experiment/3_Programs/data_import/cluster_algorithm.R")
setwd("C:/Users/Jurg/Dropbox/Doktorat/Projects/P3_Costa Experiment/2_data")

experiment <- "invasion" # set to "growth", "invasion" or "exudates"

if (experiment =="growth"){
  exp.id <- "Exp2"
  n.samples <- 48
  n <- 20 # number of times measured
  org.files <- dir("C:/Users/Jurg/Desktop/Experimental Data/Exp2Growth/Originals", full.names = T)[1:n]
  samples <- rep(list(1:48),n)
  samples[[17]][c(19:21,34:36,40:45)] <- 73:84
  sample[[18]][c(25:27)] <- c(1,5,9)
  sample[[19]] <- c(46+1:24+0:23%/%3,12+1:24+0:23%/%3)
  sample[[20]] <- c(rep(c(6:1,2,6),each = 3),rep(c(12:7,8,12),each =3))
  
}
if (experiment =="invasion"){
  exp.id <- "Exp3"
  n.samples <- 48
  n <- 12 # number of times measured
  org.files <- dir("C:/Users/Jurg/Desktop/Experimental Data/Exp3invasion/Originals", full.names = T)[1:n]
  samples <- rep(list(1:48),n)
}
if (experiment =="exudates"){
  exp.id <- "Exp5"
  n.samples <- 36
  n <- 15 # number of times measured
  vec.times <- c(1,1:14)
  org.files <- dir("C:/Users/Jurg/Desktop/Experimental Data/Exp5exudates_redo/Originals", full.names = T)[1:n]
  org.files <- org.files[vec.times]
  samples <- list(1:36,c(37:72,85:96),1:36)
}
if (experiment =="variance"){
  exp.id <- paste("C:/Users/Jurg/Dropbox/Doktorat/Projects/Z2_smallerprojects",
                  "Variance of Flowcytometer, Stan/2_data/Exp_var", sep = "/")
  n.samples <- 42
  n <- 2 # number of times measured
  org.files <- dir("C:/Users/Jurg/Desktop/Experimental Data/Exp_variance_Stan/Originals"
                   , full.names = T)[1:n]
  samples <- rep(list(1:42),n)
}
mus <- data.matrix(read.csv(paste(exp.id,"mus.csv", sep = ","), sep = ","))
sigmas <- data.matrix(read.csv(paste(exp.id,"sigmas.csv", sep = ","), sep = ","))
percentages <- data.matrix(read.csv(paste(exp.id,
                                  "percentages.csv", sep = ","), sep = ","))
n.cells <- data.matrix(read.csv(paste(exp.id,"n_cells.csv", sep = ","), sep = ","))
cluster.index <- data.matrix(read.csv(paste(exp.id,
                                  "cluster.index.csv", sep = ","), sep = ","))

color_data1 <- function(data, tau, mu, sig, pars,clusters.pres= 1:5){
  # find probabilities

  lambda <- matrix(NA, nrow = nrow(data), ncol = 5)
  for (p in clusters.pres){
    lambda[,p] <- tau[p]*mvnorm(data,mu[[p]], sig[[p]])
  }
  
  lambda[,-clusters.pres] <- 0
  lambda <- lambda/rowSums(lambda[,clusters.pres],na.rm = T)
  plots <- list()
  counter <- 1
  for (p in pars){
    par1 <- cols[p[1]]
    par2 <- cols[p[2]]
    plots[[counter]] <- ggplot()+
      geom_point(data = as.data.frame(data),aes_string(par1, par2),colour = "grey")+
      geom_point(data = as.data.frame(data[which(lambda[,1]>0.7),]),aes_string(par1, par2),colour = "yellow")+
      geom_point(data = as.data.frame(data[which(lambda[,2]>0.7),]),aes_string(par1, par2),colour = "blue")+
      geom_point(data = as.data.frame(data[which(lambda[,3]>0.7),]),aes_string(par1, par2),colour = "orange")
      #geom_point(data = as.data.frame(data[which(lambda[,4]>0.7),]),aes_string(par1, par2),colour = "green")
      #geom_point(data = as.data.frame(data[which(lambda[,5]>0.7),]),aes_string(par1, par2),colour = "red")
    counter <- counter + 1
  }
  return(plots)
}

cyano.ranges <- matrix(c(1.8,4.0,0.5,2.0,1.5,2.5,3.0,4.0,1.5,3.5,
                         1.8,4.0,2.0,4.0,2.0,3.5,2.5,4.0,1.5,3.0),2,10)
# create noneg filter
cols.all <- c("FSC.HLin", "SSC.HLin", "GRN.B.HLin","YEL.B.HLin", "RED.B.HLin","NIR.B.HLin",
          "RED.R.HLin", "NIR.R.HLin","SSC.ALin", "SSC.W")
all.par = length(cols.all)
# remove boundary effects, all data are within 0 and 5 
mat<-matrix(rep(c(0.1,99990),all.par),ncol = all.par, dimnames = list(c("min", "max"),cols.all))
finite <- rectangleGate(filterId="positive", .gate=mat)

cols <- cols.all[c(1,4,5,7,8)]

# log transform all the data
logtrans <- transformList(cols, logTransform()) #log transform all parameters

fil <- 5
dataset <- 10
if (fil>15){ cyano.ranges <- matrix(c(1.8,2.8,0.0,1.5,1.5,2.5,1.5,2.5,1.2,2.5,
                                             1.8,3.0,1.5,3.0,1.5,3.0,1.5,2.5,1.0,2.5),2,10)}
data.set <- read.FCS(org.files[fil], emptyValue = F, dataset = samples[[fil]][dataset],
                     alter.names = T)
data.set <- Subset(data.set, finite)
data.set <- data.set[,c(1,4,5,7,8)]

data.set <- suppressWarnings(transform(data.set, logtrans))
data <- exprs(data.set)
ind <- (n.samples*(fil-1)+dataset-1)*5+1

ind2 <- n.samples*(fil-1)+dataset
tau.data <- c(percentages[ind2,3:7])
mu.data <- data.frame(t(data.matrix(mus[ind:(ind+4),])[,4:8]))
sig.data <- list()
for (i in 1:5){
  sig.data[[i]] <- matrix(sigmas[ind+i-1,4:28],5,5)
}
color_data1(data, tau.data, mu.data, sig.data,list(c(1,3)))

if (T){
  print(mus[ind:(ind+4),])
  print(cyano.ranges)
  print(n.cells[ind2,])
  print(cluster.index[ind:(ind+4),])
  }
