library(magrittr)
#vsc_home
#vsc_home <- Sys.getenv("VSC_HOME")
#UNamur CECI HOME
ceci_home <- Sys.getenv("HOME")
ncore <- parallel::detectCores(as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE")))
ncl <- parallel::makeCluster(ncore)

###distances
##VSC
#source(paste0(vsc_home,"/Functions/distances.R"))

##My system
#source(paste0("Codes","/Functions/distances.R"))

##UNamur CECI
source(paste0(ceci_home,"/Functions/distances.R"))
snow::clusterExport(ncl, list("eucldis", "L1_inf","rescale"))

###TOP
##VSC
#source(paste0(vsc_home,"/Functions/TOP.R"))
#source(paste0(vsc_home,"/Functions/TOP_NEW.R"))

##My system
#source(paste0("Codes","/Functions/TOP.R"))
#source(paste0("Codes","/Functions/TOP_NEW.R"))

##UNamur CECI
source(paste0(ceci_home,"/Functions/TOP.R"))
source(paste0(ceci_home,"/Functions/TOP_NEW.R"))
snow::clusterExport(ncl, list("TOP", "TOP_new"))

###TED
##VSC
#source(paste0(vsc_home,"/Functions/TED.R"))
#source(paste0(vsc_home,"/Functions/TED_NEW.R"))

##My system
#source(paste0("Codes","/Functions/basetraitM.R"))
#source(paste0("Codes","/Functions/TED.R"))
#source(paste0("Codes","/Functions/TED_NEW.R"))

##UNamur CECI
source(paste0(ceci_home,"/Functions/basetraitM.R"))
source(paste0(ceci_home,"/Functions/TED.R"))
source(paste0(ceci_home,"/Functions/TED_NEW.R"))
snow::clusterExport(ncl, list("base_traitmatrix", "TED", "TED_new"))


#############Different Sample Size
n_dim <- 3
n_row <- seq(100, 10000, by = 100)
b_matrix <- matrix(NA, nrow = 100, ncol = n_dim)
b_matrix <- apply(b_matrix, 2, function(x)runif(100, 0, 1))
b_matrix2 <- base_traitmatrix(basetype = "cube", ncols = n_dim, dsize = 100)
snow::clusterExport(ncl, c("n_dim", "n_row", "b_matrix", "b_matrix2"))

TOPS <- TEDS <- TEDS_Classic <- vector("list", length = length(n_row))

for(i in 1:3){
    Tsdata1_2and3D <- matrix(runif(n_row[i] * n_dim, 0, 1), ncol = n_dim)
    middle <- mrfDepth::hdepthmedian(Tsdata1_2and3D, maxdir = 3000)$median
    dist_mean <- apply(Tsdata1_2and3D - middle,1,function(x){
      L1_inf(x)
    })
    mid_points <- which(dist_mean < quantile(dist_mean,0.20))
  
    #T1
    T1 <- Tsdata1_2and3D[-mid_points,]
    T1b <- Tsdata1_2and3D
    bb <- mid_points #sample(mid_points,length(mid_points)/4)
  
    #T1b
    T1b[bb, 1]  <- T1b[bb, 1] - min(T1b[bb, 1]) #min(T1b[,1])
    T1b[bb, 2]  <- T1b[bb,2] + 1.5
  
    #Scenario T_1c
    T1c <- T1b
    T1c[bb, 2] <- T1c[bb, 2] + 3.5
  
    #T2a
    pp1 <- base::which(Tsdata1_2and3D[,1] <= 0.3)
    T2a <- Tsdata1_2and3D[-pp1,]
    pp1b <- base::which(T2a[,1]<=0.35)
    pp1c <- base::which(T2a[pp1b,2]==max(T2a[pp1b,2]) )
    T2a[pp1b[pp1c], 1] <- 0.2; T2a[pp1b[pp1c], 2] <- max(T2a[, 2]); T2a[pp1b[pp1c], 3] <- max(T2a[, 3])
  
    #T2b
    T2b <- T2a
    pp2 <- base::sample(base::which(T2b[,1] > 0.3),12)
    T2b[pp2,1] <- 0.0
    print(i)
    parallel::clusterExport(ncl, 
                            c("i", "Tsdata1_2and3D","T1",
                              "T1b", "T1c", "T2a", "T2b"))
    
    sfc <- snow::clusterApply(cl = ncl, x = 1:10, fun = function(j){
      TOP_T1aT2 <- TOP(Tsdata1_2and3D); TOP_T1aT2_new <- TOP_new(Tsdata1_2and3D)
      TOP_T1 <- TOP(T1); TOP_T1_new <- TOP_new(T1)
      TOP_T1b <- TOP(T1b); TOP_new_T1b <- TOP_new(T1b)
      TOP_T1c <- TOP(T1c); TOP_new_T1c <- TOP_new(T1c)
      TOP_T2a <- TOP(T2a); TOP_new_T2a <- TOP_new(T2a)
      TOP_T2b <- TOP(T2b); TOP_new_T2b <- TOP_new(T2b)
    
      TOPP <- data.frame(TOP_T1aT2$TA, TOP_T1aT2_new$TA, TOP_T1aT2_new$N,
                         TOP_T1$TA, sum(TOP_T1_new$Areas)/TOP_T1aT2_new$N, TOP_T1aT2_new$N,
                         TOP_T1b$TA, TOP_new_T1b$TA, TOP_new_T1b$N,
                         TOP_T1c$TA, TOP_new_T1c$TA, TOP_new_T1c$N,
                         TOP_T2a$TA, TOP_new_T2a$TA, TOP_new_T2a$N,
                         TOP_T2b$TA, TOP_new_T2b$TA, TOP_new_T2b$N,
                         TOP_T1aT2$Areas[1], TOP_T1aT2_new$TA/TOP_T1aT2_new$N,
                         nrow(T1),nrow(T1b), nrow(T1c), nrow(T2a), nrow(T2b))
      names(TOPP) <- c("T1aT2","T1aT2_new","N_T1aT2","T1","T1_new","N_T1","T1b","T1b_new","N_T1b",
                       "T1c","T1c_new","N_T1c","T2a","T2a_new","N_T2a","T2b","T2b_new","N_T2b",
                       "T2c","T2c_new","n_T1","n_T1b","n_T1c","n_T2a","n_T2b")
    
      TEDD <- data.frame(TED(Tsdata1_2and3D, base_matrix = b_matrix), TED_new(Tsdata1_2and3D),
                         TED(T1, base_matrix = b_matrix), TED_new(T1),
                         TED(T1b, base_matrix = b_matrix), TED_new(T1b),
                         TED(T1c, base_matrix = b_matrix), TED_new(T1c),
                         TED(T2a, base_matrix = b_matrix), TED_new(T2a),
                         TED(T2b, base_matrix = b_matrix), TED_new(T2b))
      TEDD_classic <- data.frame(TED(Tsdata1_2and3D, base_matrix = b_matrix2),
                                 TED(T1, base_matrix = b_matrix2),
                                 TED(T1b, base_matrix = b_matrix2),
                                 TED(T1c, base_matrix = b_matrix2),
                                 TED(T2a, base_matrix = b_matrix2), 
                                 TED(T2b, base_matrix = b_matrix2))
      names(TEDD) <- c("T1aT2","T1aT2_new","T1","T1_new","T1b","T1b_new","T1c",
                       "T1c_new","T2a","T2a_new","T2b","T2b_new")
      names(TEDD_classic) <- c("T1aT2", "T1", "T1b", "T1c", "T2a", "T2b")
      
      return(list(TOPS = TOPP, TEDS = TEDD, TEDS_Classic = TEDD_classic))
  })
  
  #TOPS[[i]] <-  TOPP
  TOPS[[i]] <- lapply(sfc, "[[", 1) %>% do.call(rbind.data.frame, .)
  #TEDS[[i]] <-  TEDD
  TEDS[[i]] <- lapply(sfc, "[[", 2) %>% do.call(rbind.data.frame, .)
  TEDS_Classic[[i]] <- lapply(sfc, "[[", 3) %>% do.call(rbind.data.frame, .)
}

save.image("ThreeDimensionSimulation.RData")
