# -------------------------------------------------------
# Evolution of FEve-SR and FRic - SR relationships
# ---------------------------------------------------------

# STEP 2: RUN SIMULATIONS AND EMPIRICAL DATA ANALYSIS

# load functions & packages
rm(list=ls())
source("R/functions.R")
source("R/packages.R")

# load R data
load(here ( "Processed_data_sensitivity_2","image_rodents.RData"))
load(here ( "Processed_data","Spatial_covs_data.RData"))

dir.create("Output_sensitivity_2")
# -----------------------------------------
# Empirical traits

# run
empirical_FD <- lapply(match_comm_data, function (i) 
  
  dbFD(x=std_traits_subset,
       a=i$comm,
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       stand.FRic = T,
       m = "max",
       corr = "lingoes",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)
  )

# save
save (empirical_FD,
      file= here("Output_sensitivity_2", "empirical_FD_rodents.RData"))


# -------------------------------------------------------------------
## Simulate trait evolution according to a bivariate "BMM" model
# Number of traits
ntraits<-ncol(std_traits_subset)

# Number of simulated (pairs of) traits
nsim<-50
nc<-5

# simulate parameters
simul_param_BM <- lapply (match_data, function (i) 
  
  
  fitContinuous(phy=i$phy,  
                dat = (i$data), 
                model="BM", 
                #SE=NA,
                ncores = nc)
  
)

# ancestral states for each trait
theta<-rep(0,ntraits)

# Simulate

simul<-lapply (seq(1,length(tree_list)), function (i) 
  
        mvSIM(tree_list[[i]],
             nsim=nsim, 
             model="BM1",
             param=list(sigma=diag (c(simul_param_BM[[i]]$Body.mass.g$opt$sigsq,
                                      simul_param_BM[[i]]$LitSz$opt$sigsq,
                                      simul_param_BM[[i]]$MaxLifepsan$opt$sigsq
                                      
                                      )),
                        theta=theta,
                        ntraits=ntraits,
                        names_traits=c("Trait 1",
                                       "Trait 2",
                                       "Trait 3"))))


# reduce (per phylogeny) to have the average of multivariate traits
mean_simul <- lapply (simul, function (i)
  
            Reduce("+",i)/length(i))


# simulated FD
# run across simulations
simulated_FD <- lapply (seq (1,length (match_comm_data)), function (i)
  
                dbFD(x=mean_simul[[i]][which(rownames(mean_simul[[i]]) %in% colnames(match_comm_data[[i]]$comm)),],
                     a=data.matrix(match_comm_data[[i]]$comm),
                     w.abun=T,
                     stand.x=F,
                     calc.FRic = T,
                     stand.FRic = T,
                     m = "max",
                     corr = "lingoes",
                     calc.CWM = F,
                     calc.FDiv=F,
                     print.pco = T)
  
)

# save
save (simul_param_BM,
      simulated_FD,
      file= here("Output_sensitivity_2", "simulated_FD_BM_rodents.RData"))

# ----------------------------------------------
# niche filling (early burst)

# sigmas
# simulate parameters
simul_param_EB <- lapply (match_data, function (i) 
  
  
  fitContinuous(phy=i$phy,  
                dat = (i$data), 
                model="EB", 
               # SE=NA,
                ncores = nc)
  
  )

# Simulate traits based on the Early Burst model
# parallel
cl <- makeCluster(nc) ## number of cores

# export packages
clusterEvalQ(cl, library(mvMORPH))

# export your data and function
clusterExport(cl, c("simul_param_EB", 
                    "theta",
                    "ntraits",
                    "tree_list",
                    'nsim'))

simul_EB<-parLapply (cl, seq(1,length(tree_list)), function (i) 
  
  mvSIM(tree_list[[i]],
        nsim=nsim, 
        model="EB",
        param=list(sigma=diag (c(simul_param_EB[[i]]$Body.mass.g$opt$sigsq,
                                 simul_param_EB[[i]]$LitSz$opt$sigsq,
                                 simul_param_EB[[i]]$MaxLifepsan$opt$sigsq
                                 
                                 )), 
                   beta=diag (c(simul_param_EB[[i]]$Body.mass.g$opt$a,
                                simul_param_EB[[i]]$LitSz$opt$a,
                                simul_param_EB[[i]]$MaxLifepsan$opt$a)),
                   theta=theta,
                   ntraits=ntraits,
                   names_traits=c("Trait 1",
                                  "Trait 2",
                                  "Trait 3"))))

stopCluster(cl)

# reduce (per phylogeny) to have the average of multivariate traits
mean_simul_EB <- lapply (simul_EB, function (i)
  
  Reduce("+",i)/length(i))

# simulated FD
# run across simulations
simulated_FD_EB <- lapply (seq (1,length (match_comm_data)), function (i)
  
  dbFD(x=mean_simul_EB[[i]][which(rownames(mean_simul_EB[[i]]) %in% colnames(match_comm_data[[i]]$comm)),],
       a=data.matrix(match_comm_data[[i]]$comm),
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       stand.FRic = T,
       m = "max",
       corr = "lingoes",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)
  
)

# save
save (simul_param_EB,
      simulated_FD_EB,
      file= here("Output_sensitivity_2", "simulated_FD_EB_rodents.RData"))

#---------------------------------------------------------------
# OU
# simulate parameters
simul_param_OU <- lapply (match_data, function (i) 
  
  
  fitContinuous(phy=i$phy,  
                dat = (i$data),
                model="OU", 
               # SE=NA,
                ncores = nc)
  )


# parallel
cl <- makeCluster(nc) ## number of cores

# export packages
clusterEvalQ(cl, library(mvMORPH))

# export your data and function
clusterExport(cl, c("simul_param_OU", 
                    "theta",
                    "ntraits",
                    "tree_list",
                    'nsim'))


simul_OU<-parLapply (cl, seq(1,length(tree_list)), function (i) 
  
  tryCatch(
    mvSIM(tree_list[[i]],
          nsim=nsim, 
          model="OU1",
          param=list(sigma=diag (c(simul_param_OU[[i]]$Body.mass.g$opt$sigsq,
                                   simul_param_OU[[i]]$LitSz$opt$sigsq,
                                   simul_param_OU[[i]]$MaxLifepsan$opt$sigsq
                                   )), 
                     
                     alpha = diag (c(simul_param_OU[[i]]$Body.mass.g$opt$alpha,
                                     simul_param_OU[[i]]$LitSz$opt$alpha,
                                     simul_param_OU[[i]]$MaxLifepsan$opt$alpha
                                     )),
                     theta=theta,
                     ntraits=ntraits,
                     names_traits=c("Trait 1",
                                    "Trait 2",
                                    "Trait 3"))),
    error = function(e) return ("NULL"))
)

stopCluster (cl)

# reduce (per phylogeny) to have the average of multivariate traits
mean_simul_OU <- lapply (simul_OU, function (i)
  
  Reduce("+",i)/length(i))

# simulated FD
# run across simulations
simulated_FD_OU <- lapply (seq (1,length (match_comm_data)), function (i)
  
  dbFD(x=mean_simul_OU[[i]][which(rownames(mean_simul_OU[[i]]) %in% colnames(match_comm_data[[i]]$comm)),],
       a=data.matrix(match_comm_data[[i]]$comm),
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       m = "max",
       stand.FRic = T,
       corr = "lingoes",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)
  
)
# save
save (simul_param_OU,
      simulated_FD_OU,
      file= here("Output_sensitivity_2", "simulated_FD_OU_rodents.RData"))

# comparison of support of models
comp_mod <- lapply (seq(1,length(simul_param_BM)), function (i)
  
  data.frame(BM=AIC(simul_param_BM[[i]]$LitPerYear),
        EB=AIC(simul_param_EB[[i]]$LitPerYear),
        OU=AIC(simul_param_OU[[i]]$LitPerYear))

)
apply(do.call(rbind,comp_mod),2,mean)

# end