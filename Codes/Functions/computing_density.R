
compute.density <- function(folder = NULL, cyano.ranges = NULL,
                            treshhold = 100){
  mus <- read.csv(paste(folder,"mus.csv", sep = ","))
  percentages <- read.csv(paste(folder,"percentages.csv", sep = ","))
  n.cells <- read.csv(paste(folder, "n_cells.csv",sep = ","))
  conc.all <- read.csv(paste(folder,"conc_all.csv",sep = ","),row.names = 1)
  n.samples <- nrow(conc.all)
  
  if (is.null(cyano.ranges)){
    cyano.ranges <- matrix(c(1.8,4.0,0.5,2.0,1.5,2.5,3.0,4.0,1.5,3.5,
                           1.8,4.0,2.0,4.0,2.0,3.5,2.5,4.0,1.5,3.0),2,10)}
  
  # which cluster respects which constraints of the clusters
  cluster.index <- data.frame(matrix(nrow = nrow(mus),ncol = 10))
  colnames(cluster.index) <- c("BS4_1", "BS4_2", "BS4_3", "BS4_4","BS4_5",
                               "BS5_1", "BS5_2", "BS5_3", "BS5_4","BS5_5")
  index <- c(4:8,4:8)
  names <- c("BS4_1", "BS4_2", "BS4_3", "BS4_4","BS4_5",
             "BS5_1", "BS5_2", "BS5_3", "BS5_4","BS5_5")
  for (i in 1:10){
    cluster.index[names[i]] <- ((mus[index[i]]>cyano.ranges[1,i]) * (mus[index[i]]<cyano.ranges[2,i]))[,1]
  }
  
  # which clusters fullfill all requirements
  cluster.index["BS4.TRUE"] <- apply(cluster.index[,1:5],1,prod)
  cluster.index["BS5.TRUE"] <- apply(cluster.index[,6:10],1,prod)
  cluster.index["cluster"] <- mus$cluster
  
  for (clu in 3:7){
    # to small cluster sizes, not trustworthy
    percentages[which(n.cells[,clu]<treshhold),clu] = NA
  }
  
  # add size of BS4 clusters (proportion of all cells)
  BS4.clusters <- data.frame(matrix(nrow = nrow(percentages)))
  BS4.clusters["CLU1"] <- percentages$CLU1*cluster.index[which(cluster.index$cluster==1),"BS4.TRUE"]
  BS4.clusters["CLU2"] <- percentages$CLU2*cluster.index[which(cluster.index$cluster==2),"BS4.TRUE"]
  BS4.clusters["CLU3"] <- percentages$CLU3*cluster.index[which(cluster.index$cluster==3),"BS4.TRUE"]
  BS4.clusters["CLU4"] <- percentages$CLU4*cluster.index[which(cluster.index$cluster==4),"BS4.TRUE"]
  BS4.clusters["CLU5"] <- percentages$CLU5*cluster.index[which(cluster.index$cluster==5),"BS4.TRUE"]
  conc.BS4 <- matrix(apply(BS4.clusters[,2:6],1,sum,na.rm = TRUE),nrow = n.samples)
  
  # add size of BS5 clusters (proportion of all cells)
  BS5.clusters <- data.frame(matrix(nrow = nrow(percentages)))
  BS5.clusters["CLU1"] <- percentages$CLU1*cluster.index[which(cluster.index$cluster==1),"BS5.TRUE"]
  BS5.clusters["CLU2"] <- percentages$CLU2*cluster.index[which(cluster.index$cluster==2),"BS5.TRUE"]
  BS5.clusters["CLU3"] <- percentages$CLU3*cluster.index[which(cluster.index$cluster==3),"BS5.TRUE"]
  BS5.clusters["CLU4"] <- percentages$CLU4*cluster.index[which(cluster.index$cluster==4),"BS5.TRUE"]
  BS5.clusters["CLU5"] <- percentages$CLU5*cluster.index[which(cluster.index$cluster==5),"BS5.TRUE"]
  conc.BS5 <- matrix(apply(BS5.clusters[,2:6],1,sum,na.rm = TRUE),nrow = n.samples)
  
  
  conc.BS4 <- conc.all[1:n.samples,]*conc.BS4
  conc.BS5 <- conc.all[1:n.samples,]*conc.BS5
  
  # save results
  write.csv(cluster.index[,11:13],paste(folder, "cluster.index.csv", sep = ","))
  write.csv(conc.BS4,paste(folder,"densitiesBS4,clustering.csv",sep = ","))
  write.csv(conc.BS5,paste(folder,"densitiesBS5,clustering.csv",sep = ","))
  return(list(conc.BS4, conc.BS5,cluster.index))
}