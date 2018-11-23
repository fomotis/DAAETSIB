rm(list = ls())

# working directory of data
str.wd <- "C:/Users/Jurg/Desktop/Experimental Data/Exp3invasion"
setwd("C:/Users/Jurg/Dropbox/Doktorat/Projects/P3_Costa Experiment/3_Programs/data_import")
source("cluster_algorithm.R")
source("data_import.R")
source("computing_density.R")
source("parameter_evaluation.R")
#############################################
# find which ones do actually have to be calculated
m.files <- dir(paste(str.wd,"Meta data",sep = "/"), full.names = T) # contain metadata
files.order <- 1:12
m.files <- m.files[files.order]

sample.list <- rep(list(1:48),length(files.order))
colClasses <- rep("NULL",67)
colClasses[c(4,7,9)]<- NA # Dilution factor, cells/ul

prelim.list <- create.conc.all(m.files,sample.list,colClasses)
conc.all <- prelim.list[[1]]
events.all <- prelim.list[[2]]

conc.all[,1:2] <- NA
################################################
# cluster the datasets

# Load the actual data
org.files <- dir(paste(str.wd,"Originals",sep = "/"), full.names = T)
org.files <- org.files[files.order]

cluster.data <- data.import(org.files,sample.list,conc.all, events.all)

save.datas(conc.all,cluster.data,"Exp3")
conc.list <- compute.density("../../2_Data/Exp3")

means.BS4 <- parameter.values(org.files, sample.list, conc.all, 
                              "../../2_Data/Exp3", "BS4.TRUE")
means.BS5 <- parameter.values(org.files, sample.list, conc.all, 
                              "../../2_Data/Exp3", "BS5.TRUE")

conc.BS4 <- conc.list[[1]]
conc.BS5 <- conc.list[[2]]

missing.data.handler<- function(conc.clus, conc.ref){
  no.cyano <- which(is.na(conc.clus/conc.clus) & !is.na(conc.ref))
  id.file <- (no.cyano-1)%/%nrow(conc.ref)+1
  id.sample <- no.cyano%%nrow(conc.ref)
  id.sample[which(id.sample==0)] = nrow(conc.ref)
  id.no.cyano <- data.frame(id.file,id.sample)
  print(id.no.cyano)
}

# check whether in all samples cyanos have been found

print("Samples, where no BS4 have been found:")
missing.data.handler(conc.BS4,conc.all)

print("Samples, where no BS5 have been found")
missing.data.handler(conc.BS5,conc.all)
# Explanation for missing data
# BS5, 1-3 in all files. Were not strong enough to invade
# File 7, 40: Sampling error, did not contain anything
# Sample 13: Got corupted during experiment
# File 8,17 and 9,17 very low densities of BS5
