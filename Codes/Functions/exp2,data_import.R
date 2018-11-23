rm(list = ls())

# working directory of data
str.wd <- "C:/Users/Jurg/Desktop/Experimental Data/Exp2growth"
setwd("C:/Users/Jurg/Dropbox/Doktorat/Projects/P3_Costa Experiment/3_Programs/data_import")
source("cluster_algorithm.R")
source("data_import.R")
source("computing_density.R")
source("parameter_evaluation.R")
#############################################
# find which ones do actually have to be calculated
m.files <- dir(paste(str.wd,"Meta data",sep = "/"), full.names = T) # contain metadata
files.order <- c(1:19,21)
m.files <- m.files[files.order]

sample.list <- rep(list(1:48),length(files.order))

# sample 20 contains the dates of starting the experiment
sample.list[[20]] <- c(rep(c(6:1,2,6),each = 3),rep(c(12:7,8,12),each =3))
colClasses <- rep("NULL",67)
colClasses[c(4,7,9)]<- NA # Dilution factor, cells/ul

prelim.list <- create.conc.all(m.files,sample.list,colClasses)
conc.all <- prelim.list[[1]]
events.all <- prelim.list[[2]]

# adapt wrong settings in FSC files:
conc.all[,1] <- conc.all[,1]/1*5
conc.all[,6] <- conc.all[,6]/100*2000
conc.all[41:42,15] <- NA # wrong pipetting
conc.all[,c(2:4,12)] <- NA # can't trust measurments/not needed
# Allow only one file per measurment
conc.all[c(1:6,25:30),6] <- NA
conc.all[c(7:24,31:48),5] <- NA
conc.all[c(1:24,28:48),c(7,9,11,14)] <- NA
conc.all[25:27,c(8,10,13,15)] <- NA

# The files 17:19 were done with different worklist
# they analysed samples from the experiment after freezing them
sample.list[[17]][c(19:21,34:36,40:45)] <- 73:84
sample.list[[18]][c(25:27)] <- c(1,5,9)
sample.list[[19]] <- c(46+1:24+0:23%/%3,12+1:24+0:23%/%3)

colClasses <- rep("NULL",66)
colClasses[c(3,6,8)]<- NA # Dilution factor, cells/ul
prelim.list <- create.conc.all(m.files, sample.list, colClasses)
conc.all[,17:19] <- prelim.list[[1]][,17:19]
events.all[,17:19] <- prelim.list[[2]][,17:19]

# not all measurments were done in file 17 and 19
conc.all[-c(19:21,34:36,40:45),17] <- NA
conc.all[!is.na(conc.all[,17]),16] <- NA
conc.all[-(25:27),18] <- NA
################################################
# cluster the datasets

# Load the actual data
org.files <- dir(paste(str.wd,"Originals",sep = "/"), full.names = T)
org.files <- org.files[files.order]

cluster.data <- data.import(org.files,sample.list,conc.all,events.all)

save.datas(conc.all,cluster.data,"Exp2")
conc.first <- compute.density("../../2_Data/Exp2")
means.BS4 <- parameter.values(org.files[1:15], sample.list, conc.all, 
                              "../../2_Data/Exp2", "BS4.TRUE")
means.BS5 <- parameter.values(org.files[1:15], sample.list, conc.all, 
                              "../../2_Data/Exp2", "BS5.TRUE")

cyano.ranges.redo <- matrix(c(1.8,2.8,0.0,1.5,1.5,2.5,1.5,2.5,1.2,2.5,
                              1.8,3.0,1.5,3.0,1.5,3.0,1.5,2.5,0.9,2.5),2,10)
conc.redo <- compute.density("../../2_Data/Exp2",cyano.ranges.redo)

conc.first[[1]][,16:19] <- conc.redo[[1]][,16:19]
conc.first[[2]][,16:19] <- conc.redo[[2]][,16:19]

conc.BS4 <- conc.first[[1]]
conc.BS5 <- conc.first[[2]]

missing.data.handler<- function(conc.clus, conc.ref,add.sample = 0){
  no.cyano <- which(is.na(conc.clus/conc.clus) & !is.na(conc.ref))
  id.file <- (no.cyano-1)%/%nrow(conc.ref)+1
  id.sample <- no.cyano%%nrow(conc.ref)
  id.sample[which(id.sample==0)] = nrow(conc.ref)
  id.sample <- id.sample +add.sample
  id.no.cyano <- data.frame(id.file,id.sample)
  print(id.no.cyano)
}

# check whether in all samples cyanos have been found

print("Samples, where no BS4 have been found:")
missing.data.handler(conc.BS4[1:24,],conc.all[1:24,])

print("Samples, where no BS5 have been found")
missing.data.handler(conc.BS5[25:48,],conc.all[25:48,],add.sample = 24)

write.csv(conc.first[[1]],"../../2_Data/Exp2,densitiesBS4,clustering.csv")
write.csv(conc.first[[2]],"../../2_Data/Exp2,densitiesBS5,clustering.csv")