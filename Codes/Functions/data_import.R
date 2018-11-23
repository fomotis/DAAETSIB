library(flowCore) #loads the package flowCore
library(ggcyto)
library(flowQ)
library(flowViz)

# create noneg filter
cols.all <- c("FSC.HLin", "SSC.HLin", "GRN.B.HLin","YEL.B.HLin", "RED.B.HLin","NIR.B.HLin",
              "RED.R.HLin", "NIR.R.HLin","SSC.ALin", "SSC.W")
all.par = length(cols.all)
# remove boundary effects, all data are within 0 and 5 

cols <- cols.all[c(1,4,5,7,8)]
n.par <- length(cols) # number of parameters to use for clustering
n.clus <- 5 # number of clusters to use for clustering

# remove boundary effects, all data are within 0 and 5 
mat<-matrix(rep(c(0.1,99990),all.par),ncol = all.par, 
                  dimnames = list(c("min", "max"),cols.all))
finite <- rectangleGate(filterId="positive", .gate=mat)

# log transform all the data
logtrans <- transformList(cols, logTransform()) #log transform all parameters

create.conc.all <- function(m.files,sample.list,colClasses){
  n.files <- length(m.files)
  n.samples <- length(sample.list[[1]])
  m.datas <- list()
  col.names <- c("Events","Dilution factor","Cells/uL")
  for (i in 1:n.files){
    temp <- read.csv(m.files[i], skip = 8, header = F, 
                     colClasses = colClasses) #load meta data
    temp[3][temp[3]>1000]<- NA # to high densities, not accurate measurments
    m.datas[[i]] <- temp[sample.list[[i]],]
    colnames(m.datas[[i]]) <- col.names
  }
  # contains the concentrations in the original samples (dil*conc)
  concentrations <- list()
  events <- list()
  for(i in 1:n.files){
    conc <- as.double(m.datas[[i]]$'Dilution factor')*
      as.double(m.datas[[i]]$'Cells/uL')
    concentrations[[i]] <- conc
    events[[i]] <- as.numeric(as.character(m.datas[[i]]$"Events"))
  }
  conc.all <- do.call(cbind, concentrations)
  events.all <- do.call(cbind, events)
  return(list(conc.all, events.all))
}

data.import <- function(org.files, sample.list,conc.all,events.all){
  n.files <- length(org.files) # number of files to cluster
  n.samples <- length(sample.list[[1]]) # number of samples in each file
  
  ###################################################################################
  # prepare dataframes to save the data
  
  # name of the clusters
  clusters <- c("CLU1","CLU2","CLU3", "CLU4", "CLU5")
  
  # size of each cluster
  percentages <- data.frame(matrix(nrow = n.files*n.samples, ncol = (2+n.clus)))
  colnames(percentages) <- c("file", "dataset", clusters)
  percentages[,1] <- rep(1:n.files,each = n.samples)
  percentages[,2] <- rep(1:n.samples, times = n.files)
  
  # total number of cells in a cluster, to avoid statistical errors
  n.cells <- data.frame(matrix(nrow = n.files*n.samples, ncol = (2+n.clus)))
  colnames(n.cells) <- c("file", "dataset", clusters)
  n.cells[,1] <- rep(1:n.files,each = n.samples)
  n.cells[,2] <- rep(1:n.samples, times = n.files)
  
  # how many cells are "hard clustered"
  assigned <- data.frame(matrix(nrow = n.files*n.samples,ncol = 5))
  colnames(assigned) <- c("file", "dataset", "50%", "75%", "95%")
  assigned[,1] <- rep(1:n.files,each = n.samples)
  assigned[,2] <- rep(1:n.samples, times = n.files)
  
  
  # average of each cluster
  mus <- data.frame(matrix(nrow = n.files*n.samples*n.clus, ncol = 3+n.par))
  colnames(mus) <- c("file", "dataset","cluster", cols)
  mus[,1] <- rep(1:n.files,each = n.samples*n.clus)
  mus[,2] <- rep(1:n.samples, each = n.clus, times = n.files)
  mus[,3] <- rep(1:n.clus, times = n.files*n.samples)
  
  # standard deviation of each cluster
  sigmas <- data.frame(matrix(nrow = n.files*n.samples*n.clus, ncol = 3+n.par^2))
  sigmas[,1] <- rep(1:n.files,each = n.samples*n.clus)
  sigmas[,2] <- rep(1:n.samples, each = n.clus, times = n.files)
  sigmas[,3] <- rep(1:n.clus, times = n.files*n.samples)
  
  # the actual clustering
  for (id.file in 1:n.files){
    print(id.file) # progress report
    for (id.sample in 1:n.samples){
      # only cluster if density is low enough
      if(!is.na(conc.all[id.sample,id.file])){
        # load the datafile
        data.set <- read.FCS(org.files[id.file], emptyValue = F, alter.names = T,
                             dataset = sample.list[[id.file]][id.sample])
        # Check whether metadata and data.set are the same
        if (!(nrow(data.set)==events.all[id.sample,id.file])){
          return(data.frame(c(id.sample = id.sample, id.file = id.file, file = org.files[id.file],
                              events.org = nrow(data.set), events.meta = events.all[id.sample,id.file])))
        }
        
        # transform data
        data.set <- Subset(data.set, finite)
        data.set <- data.set[,c(1,4,5,7,8)] # only use important parameters
        data.set <- suppressWarnings(transform(data.set, logtrans)) # normal distribution in log space

        data <- exprs(data.set) # transform to dataframe
        
        # cluster the results
        save.vals <- clustering(data, 1:n.clus)
        
        # index for saving
        id1 <- (id.file-1)*n.samples + id.sample
        id2 <- (id1 - 1)*n.clus
        # multiply percentages with proportion that has been clustered
        percents <- save.vals$percentages*nrow(data)/events.all[id.sample, id.file]
        percentages[id1,] <- c(id.file,id.sample,percents)
        # approximate the number of cells in each cluster
        abs.cells <- round(save.vals$percentages*nrow(data))
        n.cells[id1,] <- c(id.file,id.sample,abs.cells)
        assigned[id1,] <- c(id.file,id.sample,save.vals$assigned)
        
        # save datas from clusters
        for (id.clus in 1:n.clus){
          mus[id2+id.clus,] <- unlist(c(id.file,id.sample,id.clus,save.vals$mus[id.clus]))
          sigmas[id2+id.clus,] <- unlist(c(id.file, id.sample, id.clus, save.vals$sigmas[id.clus]))
        }
      }
    }
  }
  return(list(percentages = percentages,n.cells = n.cells,
                   assigned = assigned,mus = mus, sigmas = sigmas))
}

save.datas <- function(conc.all,cluster.data,prefix,folder = NULL){
  if (is.null(folder)) {folder <- "C:/Users/Jurg/Dropbox/Doktorat/Projects/P4_Different_light/2_Data"}
  folder <- paste(folder,prefix,sep ="/")
  write.csv(cluster.data$mus,paste(folder,"mus.csv",sep = ",")
            ,row.names = F)
  write.csv(cluster.data$sigmas,paste(folder,"sigmas.csv",sep = ",")
            ,row.names = F)
  write.csv(cluster.data$n.cells,paste(folder,"n_cells.csv",sep = ",")
            ,row.names = F)
  write.csv(cluster.data$percentages,paste(folder,"percentages.csv",sep = ",")
            ,row.names = F)
  write.csv(cluster.data$assigned,paste(folder,"assigned.csv",sep = ",")
            ,row.names = F)
  write.csv(conc.all,paste(folder,"conc_all.csv",sep = ","))
}