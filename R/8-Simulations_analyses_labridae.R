
# -------------------------------------------------------

# Evolution of FEve-SR and FRic - SR relationships - run models - Labridae family


# ---------------------------------------------------------
rm(list=ls())
# load functions & packages
source("R/functions.R")
source("R/packages.R")

# ----------------------------------------------------------------------
# load data
# load phylogenies
#fishtree_complete_phylogeny()
tree<- read.tree (here ("data_fish","TACT","Reef_fish_all_combined.trees"))#fishtree_complete_phylogeny()

# dataframe with taxa name
df_taxa <-data.frame (sp = (unique(tree[[1]]$tip.label[order(tree[[1]]$tip.label)])))
df_taxa$sp <- gsub ("_", " ", df_taxa$sp)

# load taxonomic validation of analysis for all spp
load (file = here ("output", "tax_validation_fish.RData"))

# set names
names (worms_record_fish) <- df_taxa$sp
test_taxa <- lapply (worms_record_fish,data.frame)
test_taxa <- test_taxa[unlist(lapply(test_taxa,nrow))>=1]
test_taxa <- test_taxa[unlist(lapply(test_taxa,ncol))>1]
test_taxa<-do.call(rbind,test_taxa)

# match with the table
df_taxa <- cbind (df_taxa,
                           (test_taxa[match (df_taxa$sp,rownames(test_taxa)),]))

# if missing, kepp the previous name
df_taxa$scientificname <-  ifelse (is.na(df_taxa$scientificname),
        df_taxa$sp,
        df_taxa$scientificname)

# filter labridae
df_taxa <- df_taxa[which(df_taxa$family %in% c("Labridae", "Scaridae")),]
df_taxa[grep("Cryptotomus", df_taxa$sp),]

# alter taxon names in the phlogeny
table(gsub (" ","_",df_taxa$sp [(match (tree[[2]]$tip.label,
                                  gsub (" ","_",df_taxa$sp)))]) == tree[[2]]$tip.label)


# select Labridae spp in the phylogeny
test_tree<- lapply (tree, function (i)
  
  drop.tip(i, 
           tip = i$tip.label[which(i$tip.label %in% gsub (" ","_",df_taxa$scientificname) == F)] )

)

# change tipnames by valid names
test_tree <- lapply (test_tree, function (i){

  i$tip.label <- df_taxa$scientificname [(match (i$tip.label,
                                            gsub (" ","_",df_taxa$sp)))]
  i

})
# table(gsub ("_"," ",tree[[3]]$tip.label) == test_tree[[3]]$tip.label)

# -------------------------------------------------

# load covariates
load (file = here ("Processed_data","Spatial_covs_data_reefs.RData"))

# load fish data
load(here ( "Processed_data","image_labridae.RData"))

# Empirical relatioship
# run
empirical_FD <- lapply (seq (1,length (subset_comm_data)), function (i)
                             
                             dbFD(x=subset_trait_data[[i]],
                                  a=subset_comm_data[[i]][which(rowSums(subset_comm_data[[i]]) >0),],
                                   w.abun=T,
                                   stand.x=F,
                                   calc.FRic = T,
                                   stand.FRic = T,
                                  m="max",
                                   corr = "lingoes",
                                  calc.CWM = F,
                                  calc.FDiv=F,
                                  print.pco = T)
                        )
# save
save (empirical_FD,
      file= here("Output", 
                 "empirical_FD_labridae.RData"))


# -------------------------------------

# Random traits

r_traits_sim <- lapply (seq(1,length(match_data)), function (i) {
  
  # matrix traits
  r_traits <- match_data[[i]]$data
  matrix_traits  <- apply (r_traits,2, function (x) {
  
            rnorm(nrow(r_traits),0,1)
  
  })
  dimnames(matrix_traits) <- dimnames(r_traits)
  matrix_traits
  
  }
)

# simulated FD
# run across simulations
simulated_FD_random <- lapply (seq (1,length (match_comm_data)), function (i)
  
                dbFD(x=r_traits_sim[[i]][which(rownames(r_traits_sim[[i]]) %in% colnames(subset_comm_data[[i]])),],
                     a=subset_comm_data[[i]][which(rowSums(subset_comm_data[[i]]) >0),],
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
save (simulated_FD_random,
      file= here("Output", "simulated_FD_random_labridae.RData"))

# -----------------------------------------------------------------
# trait simulation

## Simulate trait evolution according to a bivariate "BM" model
# Number of traits
ntraits<-ncol(subset_trait_data[[1]])
# Number of simulated (pairs of) traits
nsim<-50
# ncores
nc <- 5
# sigmas
#sigma<- (rbind(c(1,0.25,0.25),
#              c(0.25,1,0.25),
#              c(0.25,0.25,1)))
#
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

simul<-lapply (seq(1,length(match_data)), function (i) 
  
  mvSIM(match_data[[i]]$phy,
        nsim=nsim, 
        model="BM1",
        param=list(sigma=diag (c(simul_param_BM[[i]]$Body_size$opt$sigsq,
                                 simul_param_BM[[i]]$Trophic_level$opt$sigsq,
                                 simul_param_BM[[i]]$Aspect_ratio$opt$sigsq,
                                 simul_param_BM[[i]]$Depth_max$opt$sigsq,
                                 simul_param_BM[[i]]$TemPref_mean$opt$sigsq)),
                   theta=theta,
                   ntraits=ntraits,
                   names_traits=c("Trait 1",
                                  "Trait 2",
                                  "Trait 3",
                                  "Trait 4",
                                  "Trait 5"))))


# reduce (per phylogeny) to have the average of multivariate traits
mean_simul <- lapply (simul, function (i)
  
  Reduce("+",i)/length(i))


# simulated FD
# run across simulations
simulated_FD <- lapply (seq (1,length (subset_comm_data)), function (i)
  
  dbFD(x=mean_simul[[i]][which(rownames(mean_simul[[i]]) %in% colnames(subset_comm_data[[i]])),],
       a=subset_comm_data[[i]][which(rowSums(subset_comm_data[[i]]) >0),],
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       m="max",
       stand.FRic = T,
       corr = "lingoes",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)
  
)

# save
save (simul_param_BM,
      simulated_FD,
      file= here("Output", 
                 "simulated_FD_BM_labridae.RData"))

# ----------------------------------------------
# niche filling (early burst)
# estimating parameters
# simulate parameters
simul_param_EB <- lapply (match_data, function (i) 
  
  fitContinuous(phy=i$phy,  
                dat = (i$data), 
                model="EB", 
                #SE=NA,
                ncores = nc)
  
  )

# ancestral states for each trait
theta<-rep(0,ntraits)

#run trait simulation
simul_EB<-lapply (seq(1,length(match_data)), function (i) 
  
  tryCatch(
    mvSIM(match_data[[i]]$phy,
          nsim=nsim, 
          model="EB",
          param=list(sigma=diag (c(simul_param_EB[[i]]$Body_size$opt$sigsq,
                                   simul_param_EB[[i]]$Trophic_level$opt$sigsq,
                                   simul_param_EB[[i]]$Aspect_ratio$opt$sigsq,
                                   simul_param_EB[[i]]$Depth_max$opt$sigsq,
                                   simul_param_EB[[i]]$TemPref_mean$opt$sigsq)), 
                     beta=diag (c(simul_param_EB[[i]]$Body_size$opt$a,
                                  simul_param_EB[[i]]$Trophic_level$opt$a,
                                  simul_param_EB[[i]]$Aspect_ratio$opt$a,
                                  simul_param_EB[[i]]$Depth_max$opt$a,
                                  simul_param_EB[[i]]$TemPref_mean$opt$a)),
                     theta=theta,
                     ntraits=ntraits)),
  error = function(e) return ("NULL"))
  
  )

# rm error
#correct<-which(unlist(lapply (simul_EB,length)) == 50) # all successful simulations
#simul_EB <- (simul_EB[correct]) # remove

# reduce (per phylogeny) to have the average of multivariate traits
mean_simul_EB <- lapply (simul_EB, function (i)
  
  Reduce("+",i)/length(i))


# simulated FD
# run across simulations
#match_comm_data_sub <- match_comm_data [correct] # rm errors

simulated_FD_EB <- lapply (seq (1,length (subset_comm_data)), function (i)
  
  dbFD(x=mean_simul_EB[[i]][which(rownames(mean_simul_EB[[i]]) %in% colnames(subset_comm_data[[i]])),],
       a=subset_comm_data[[i]][which(rowSums(subset_comm_data[[i]]) >0),],
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       stand.FRic = T,
       corr = "lingoes",
       m="max",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)
  
)
# save
save (simul_param_EB,
      simulated_FD_EB,
      file= here("Output", 
                 "simulated_FD_EB_labridae.RData"))

# ------------------------------
# OU
# estimating parameters
simul_param_OU <- lapply (match_data, function (i) 
  
  fitContinuous(phy=i$phy,  
                dat = (i$data), 
                model="OU", 
                SE=NA)
  
  )
# ancestral states for each trait
theta<-rep(0,ntraits)

#run trait simulation
simul_OU<-lapply (seq(1,length(match_data)), function (i) 
  
  tryCatch(
    mvSIM(match_data[[i]]$phy,
          nsim=nsim, 
          model="OU1",
          param=list(sigma=diag (c(simul_param_OU[[i]]$Body_size$opt$sigsq,
                                   simul_param_OU[[i]]$Trophic_level$opt$sigsq,
                                   simul_param_OU[[i]]$Aspect_ratio$opt$sigsq,
                                   simul_param_OU[[i]]$Depth_max$opt$sigsq,
                                   simul_param_OU[[i]]$TemPref_mean$opt$sigsq)), 
                     alpha = diag (c(simul_param_OU[[i]]$Body_size$opt$alpha,
                                     simul_param_OU[[i]]$Trophic_level$opt$alpha,
                                     simul_param_OU[[i]]$Aspect_ratio$opt$alpha,
                                     simul_param_OU[[i]]$Depth_max$opt$alpha,
                                     simul_param_OU[[i]]$TemPref_mean$opt$alpha)),
                     theta=theta,
                     ntraits=ntraits,
                     names_traits=c("Trait 1",
                                    "Trait 2",
                                    "Trait 3",
                                    "Trait 4",
                                    "Trait 5"))),
    error = function(e) return ("NULL"))
  
)

# rm error
#correctOU<-which(unlist(lapply (simul_OU,length)) == 50) # all successful simulations
#simul_OU <- (simul_OU[correctOU]) # remove

# reduce (per phylogeny) to have the average of multivariate traits
mean_simul_OU <- lapply (simul_OU, function (i)
  
  Reduce("+",i)/length(i))

# simulated FD
# run across simulations
simulated_FD_OU <- lapply (seq (1,length (subset_comm_data)), function (i)
  
  tryCatch(
  dbFD(x=mean_simul_OU[[i]][which(rownames(mean_simul_OU[[i]]) %in% colnames(subset_comm_data[[i]])),],
       a=subset_comm_data[[i]][which(rowSums(subset_comm_data[[i]]) >0),],
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       stand.FRic = T,
       corr = "lingoes",
       m="max",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T),
  error = function(e) return ("NULL"))
  
)

# save
save (simul_param_OU,
      simulated_FD_OU,
      file= here("Output", "simulated_FD_OU_labridae.RData"))


# -----------------------------------
# simulate multiple optimum OU to compare with a single optimum OU

## make analysis input data.frame
regime <- match_data[[1]]$data$Body_size
regime <- cut(regime, breaks = c(-1,-0.5, 0.5,6))


data<-data.frame(Genus_species=rownames(match_data[[1]]$data),
                 Reg=as.factor (regime),
                 Body_size = match_data[[1]]$data$Body_size)



require("OUwie")
fitOU<-OUwie(match_data[[1]]$phy,
             data,
             model="OUM",
             simmap.tree = F,
             algorithm="invert")



data(tworegime)

#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree)
nodelabels(pch=21, bg=select.reg)



## Not run: 
#To see the first 5 lines of the data matrix to see what how to
#structure the data:
trait[1:5,]

#Now fit an OU model that allows different sigma^2:
OUwie(tree,trait,model=c("OUMV"))


# ------------------------------
# OU
# estimating parameters
simul_param_MOU <- lapply (match_data, function (i) 
  
  fitContinuous(phy=i$phy,  
                dat = (i$data), 
                model="OUM", 
                SE=NA)
  
  )
# ancestral states for each trait
theta<-rep(0,ntraits)

rm(list=ls())
# -----------------------------------------------------
# end