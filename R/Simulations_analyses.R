
# -------------------------------------------------------
# Evolution of FEve-SR and FRic - SR relationships
# ---------------------------------------------------------

# load functions & packages
source("R/functions.R")
source("R/packages.R")

# ----------------------------------------------------------------------
# load data
# load phylogenies
#fishtree_complete_phylogeny()
tree<- read.tree (here ("data","TACT","Reef_fish_all_combined.trees"))#fishtree_complete_phylogeny()

# load community data
# UVC fish data
peixes <- read.csv(here("data","UpdatedData_RMorais_et_al_2017.csv"))

## modify eventID to rm year
eventID_MOD  <- paste (peixes$Region, peixes$Locality,peixes$Site, peixes$eventDepth,
                       sep = "_")
peixes$eventID_MOD <- eventID_MOD
# number of belt transects
length(unique(peixes$Transect_id))

# number of localities
length(unique(peixes$Locality))

# obtain table 
tab_sp_site<-cast(formula = eventID_MOD ~ ScientificName,
     value="IndCounting",
     data=peixes,
     fun.aggregate=sum)
# transforming into DF
tab_sp_site<-(data.frame(tab_sp_site))
# site names
sites <- tab_sp_site$eventID_MOD

# effort
list_sites <- unique(peixes$eventID_MOD)
ntrans <- lapply (list_sites, function (i) {
  # number of transects per site
  length(unique(peixes[which(peixes$eventID_MOD %in% i),"Transect_id"]))
})
# df
effort_site <- data.frame(sites = list_sites,
                          effort = unlist(ntrans))
# match with table
effort_site <- effort_site[match(sites,effort_site$sites),]

# load trait data
traits_peixes <- read.csv(here("data","Atributos_especies_Atlantico_&_Pacifico_Oriental_2020_04_28.csv"),
                          h=T,sep=";")

# SW atlatic spp
traits_peixes <- traits_peixes[which(  traits_peixes$Province_13 == 1 | 
                                     traits_peixes$Province_14 == 1 | 
                                     traits_peixes$Province_47 ==1),]

# adjust names to match community, trait, and phylogeny
# trait
traits_peixes$Name <- tolower(gsub(" ",".",traits_peixes$Name)) 
traits_peixes<- traits_peixes[duplicated(traits_peixes$Name) !=T,]
rownames(traits_peixes) <- traits_peixes$Name

# calculate functional metrics
## subset of fish traits and communities
tab_sp_site<- tab_sp_site[,which(colnames(tab_sp_site) %in% traits_peixes$Name)] # spp in the trait dataset
#subset_traits_peixes <- traits_peixes[which(traits_peixes$Name %in% colnames(tab_sp_site)),] # traits in the community
## interesting traits
interesting_traits <- c("Body_size", "Trophic_level", "Aspect_ratio","Depth_max","TemPref_mean")
# subset
subset_traits_peixes <- traits_peixes[,interesting_traits]
## replacing comma by dot, and transforming into number
subset_traits_peixes$Body_size <- as.numeric(gsub (",",".",subset_traits_peixes$Body_size))
subset_traits_peixes$Trophic_level <- as.numeric(gsub (",",".",subset_traits_peixes$Trophic_level))
subset_traits_peixes$Aspect_ratio <- as.numeric(gsub (",",".",subset_traits_peixes$Aspect_ratio))
subset_traits_peixes$TemPref_mean <- as.numeric(gsub (",",".",subset_traits_peixes$TemPref_mean))
subset_traits_peixes$Depth_max <- as.numeric(gsub (",",".",subset_traits_peixes$Depth_max))

# standardize traits
std_traits <- apply (subset_traits_peixes, 2, scale) # scale trait values
std_traits<-data.frame(std_traits)# dataframe (to dbFD function)
rownames(std_traits)<- rownames(subset_traits_peixes) #lose names

# imputation without phylogeny
require(missForest)
std_traits <- missForest (std_traits, maxiter = 50,
                          ntree= 100,variablewise = T)
std_traits<-std_traits$ximp
  
# match spp names in trait and community dataset
#std_traits <- std_traits [match(colnames(tab_sp_site),rownames(std_traits)),]
#rownames(std_traits) == colnames(tab_sp_site)

## capture per trapping effort (number of transects per site)
tab_sp_site <- (tab_sp_site / effort_site$effort)

# adjust tip labels
# phylogeny
tree<-lapply (tree, function (i) {
  
  i$tip.label<-tolower(gsub("_",".", i$tip.label))
  ;
  i
  })


# table(colnames(tab_sp_site) %in% tree$tip.label)
# finally, match phylogeny, traits, and community
# match phylogenetic and trait data

match_data <- lapply (tree, function (i) 
  
                  match.phylo.data(i, 
                               std_traits)
)

# match phylgeny and community
match_comm_data<-lapply(match_data, function (i)
  
          match.phylo.comm(i$phy,
                           tab_sp_site)
          )

# subsetting tarit data
subset_trait_data <-lapply (seq(1,length(match_data)), function (i) 
  
  
  match_data[[i]]$data[which(rownames(match_data[[i]]$data) %in% colnames(match_comm_data[[i]]$comm)),]

)
# subset comm data  
subset_comm_data <- lapply (seq(1,length(match_comm_data)), function (i) 
  
  match_comm_data[[i]]$comm[,which(colnames(match_comm_data[[i]]$comm) %in% rownames(subset_trait_data[[i]]))]
  
)

# phylogenetic signal

psignal <- lapply (seq(1,ncol (match_data[[1]]$data)), function (i)
  
  phylosig(match_data[[1]]$phy, 
           match_data[[1]]$data[,i], 
           method="K", test=TRUE, nsim=999)
)

# df with res
psignal <- do.call (rbind, 
         
         lapply (psignal, function (i)
           
           data.frame (K=i$K,
                       pval=i$P)
         )
)

save.image(here ( "Output","image_fish.RData"))

# =======================================================
# run
empirical_FD <- lapply (seq (1,length (subset_comm_data)), function (i)
                             
                             dbFD(x=subset_trait_data[[i]],
                                  a=subset_comm_data[[i]],
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
      file= here("output", 
                 "empirical_FD_fish.RData"))

# -----------------------------------------------------------------
# trait simulation

## Simulate trait evolution according to a bivariate "BM" model
# Number of traits
ntraits<-ncol(std_traits)
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

# ancestral states for each traits
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
       a=subset_comm_data[[i]],
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
      file= here("output", 
                 "simulated_FD_BM_fish.RData"))

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
       a=subset_comm_data[[i]],
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
      file= here("output", 
                 "simulated_FD_EB_fish.RData"))

# ------------------------------
# OU
# estimating parameters
simul_param_OU <- lapply (match_data, function (i) 
  
  fitContinuous(phy=i$phy,  
                dat = (i$data), 
                model="OU", 
                SE=NA)
  
  )

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
       a=subset_comm_data[[i]],
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
      file= here("output", "simulated_FD_OU_fish.RData"))


# -----------------------------------------------------
# end