
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
peixes$eventID_MOD  <- substr(peixes$eventID, 1,nchar(as.character(peixes$eventID))-5) 

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
interesting_traits <- c("Body_size", "Trophic_level", "Aspect_ratio")#,"Depth_max","TemPref_mean")
# subset
subset_traits_peixes <- traits_peixes[,interesting_traits]
## replacing comma by dot, and transforming into number
subset_traits_peixes$Body_size <- as.numeric(gsub (",",".",subset_traits_peixes$Body_size))
subset_traits_peixes$Trophic_level <- as.numeric(gsub (",",".",subset_traits_peixes$Trophic_level))
subset_traits_peixes$Aspect_ratio <- as.numeric(gsub (",",".",subset_traits_peixes$Aspect_ratio))
#subset_traits_peixes$TemPref_mean <- as.numeric(gsub (",",".",subset_traits_peixes$TemPref_mean))
#subset_traits_peixes$Depth_max <- as.numeric(gsub (",",".",subset_traits_peixes$Depth_max))

# standardize traits
std_traits <- apply (subset_traits_peixes, 2, scale) # scale trait values
std_traits<-data.frame(std_traits)# dataframe (to dbFD function)
rownames(std_traits)<- rownames(subset_traits_peixes) #lose names

# remove NA
std_traits<- std_traits[is.na(std_traits$Aspect_ratio) !=T,]
std_traits<- std_traits[is.na(std_traits$Trophic_level) !=T,]

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
                           tab_sp_site )
          )

# subsetting tarit data
subset_trait_data <-lapply (seq(1,length(match_data)), function (i) 
  
  
  match_data[[i]]$data[which(rownames(match_data[[i]]$data) %in% colnames(match_comm_data[[i]]$comm)),]

)
# subset comm data  
subset_comm_data <- lapply (seq(1,length(match_comm_data)), function (i) 
  
  match_comm_data[[i]]$comm[,which(colnames(match_comm_data[[i]]$comm) %in% rownames(subset_trait_data[[i]]))]
  
)

# run
empirical_FD <- lapply (seq (1,length (match_data)), function (i)
                             
                             dbFD(x=subset_trait_data[[i]],
                                  a=subset_comm_data[[i]],
                                   w.abun=T,
                                   stand.x=F,
                                   calc.FRic = T,
                                   stand.FRic = T,
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
# sigmas
#sigma<- (rbind(c(1,0.25,0.25),
#              c(0.25,1,0.25),
#              c(0.25,0.25,1)))
#
# simulate parameters
simul_param_BM <- lapply (match_data, function (i) 
  
  mvgls(scale(i$data)~1, 
        tree=i$phy,  
        model="BM", penalty="RidgeArch",
        REML=T,
        method = "Mahalanobis")
)

# ancestral states for each traits
theta<-rep(0,ntraits)

# Simulate

simul<-lapply (seq(1,length(match_data)), function (i) 
  
  mvSIM(match_data[[i]]$phy,
        nsim=nsim, 
        model="BM1",
        param=list(sigma=simul_param_BM[[i]]$sigma$S, 
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
simulated_FD <- lapply (seq (1,length (subset_comm_data)), function (i)
  
  dbFD(x=mean_simul[[i]][which(rownames(mean_simul[[i]]) %in% colnames(subset_comm_data[[i]])),],
       a=subset_comm_data[[i]],
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
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
  
  mvgls(scale(i$data) ~ 1,  
        tree=i$phy, 
        model = "EB", 
        penalty="RidgeArch",
        method = "Mahalanobis",
        REML=T))

#run trait simulation
simul_EB<-lapply (seq(1,length(match_data)), function (i) 
  
  tryCatch(
    mvSIM(match_data[[i]]$phy,
          nsim=nsim, 
          model="EB",
          param=list(sigma=simul_param_EB[[i]]$sigma$S, 
                     beta = simul_param_EB[[i]]$coefficients,
                     theta=theta,
                     ntraits=ntraits,
                     names_traits=c("Trait 1",
                                    "Trait 2",
                                    "Trait 3"))),
  error = function(e) return ("NULL"))
  
  )

# rm error
correct<-which(unlist(lapply (simul_EB,length)) == 50) # all successful simulations
simul_EB <- (simul_EB[correct]) # remove

# reduce (per phylogeny) to have the average of multivariate traits
mean_simul_EB <- lapply (simul_EB, function (i)
  
  Reduce("+",i)/length(i))



# simulated FD
# run across simulations
match_comm_data_sub <- match_comm_data [correct] # rm errors

simulated_FD_EB <- lapply (seq (1,length (match_comm_data_sub)), function (i)
  
  dbFD(x=mean_simul_EB[[i]][which(rownames(mean_simul_EB[[i]]) %in% colnames(match_comm_data_sub[[i]]$comm)),],
       a=data.matrix(match_comm_data_sub[[i]]$comm),
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       stand.FRic = T,
       corr = "lingoes",
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
  
  mvgls(scale(i$data) ~ 1,  
        tree=i$phy, 
        model = "OU1", 
        penalty="RidgeArch",
        method = "Mahalanobis",
        REML=T))

# Simulate OU
alpha <- (rbind(c(1,1,1),
                c(1,1,1),
                c(1,1,1)))


alpha <- lapply (simul_param_OU, function (i) {
  
  # coefficients (alpha) in the diagonal
  diag(alpha) <- i$coefficients;
  alpha
  
})

#run trait simulation
simul_OU<-lapply (seq(1,length(match_data)), function (i) 
  
  tryCatch(
    mvSIM(match_data[[i]]$phy,
          nsim=nsim, 
          model="OU1",
          param=list(sigma=simul_param_OU[[i]]$sigma$S, 
                     alpha = alpha[[i]],
                     theta=theta,
                     ntraits=ntraits,
                     names_traits=c("Trait 1",
                                    "Trait 2",
                                    "Trait 3"))),
    error = function(e) return ("NULL"))
  
)
nsim=2

teste_garbage<-mvSIM(match_data[[i]]$phy,
                     nsim=nsim, 
                     model="OU1",
                     param=list(sigma=simul_param_OU[[i]]$sigma$S, 
                                alpha =(alpha[[i]]),#ifelse (alpha[[i]] == 0.5,0.2,alpha[[i]]),
                                theta=theta,
                                ntraits=ntraits,
                                names_traits=c("Trait 1",
                                               "Trait 2",
                                               "Trait 3")))

# rm error
correctOU<-which(unlist(lapply (simul_OU,length)) == 50) # all successful simulations
simul_OU <- (simul_OU[correctOU]) # remove

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
       stand.FRic = T,
       corr = "lingoes",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)
  
)

# save
save (simul_param_OU,
      simulated_FD_OU,
      file= here("output", "simulated_FD_OU_fish.RData"))


# -----------------------------------------------------
# empirical results
empirical_results <- data.frame (SR= apply(sapply(empirical_FD,"[[","nbsp"),1,mean),
                                 FRic= apply(sapply(empirical_FD,"[[","FRic"),1,mean),
                                 FEve=apply(sapply(empirical_FD,"[[","FEve"),1,mean),
                                 Dataset= "Empirical")

# average of simulated values (brownian motion)
simulated_results_BM <- data.frame (SR= apply(sapply(simulated_FD,"[[","nbsp"),1,mean),
                                 FRic= apply(sapply(simulated_FD,"[[","FRic"),1,mean),
                                 FEve=apply(sapply(simulated_FD,"[[","FEve"),1,mean),
                                 Dataset= "SimulatedBM")

# average of simulated values by EB
simulated_results_EB <- data.frame (SR= apply(sapply(simulated_FD_EB,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD_EB,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD_EB,"[[","FEve"),1,mean),
                                    Dataset = "SimulatedEB")
# average of simulated values by OU
simulated_results_OU <- data.frame (SR= apply(sapply(simulated_FD_OU,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD_OU,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD_OU,"[[","FEve"),1,mean),
                                    Dataset = "SimulatedOU")

# bind them
df_analyzes <- rbind(empirical_results,
                     simulated_results_BM,
                     simulated_results_EB,
                     simulated_results_OU)

##---------------------------------------------------------
# analyses
# MCMC settings
nc<-3
ni<-10000
nb<-5000
nt<-10

# run model (ancova)
model.ancova.FRic <- brm (FRic ~ poly(SR,2)*Dataset,
                     data=df_analyzes,
                     family = gaussian (link="identity"),
                     chains=nc,
                     iter = ni,
                     warmup = nb,
                     thin=nt)

# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)

# plotting
p1<-plot(conditional_effects(model.ancova.FRic,
                             method="fitted",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +

	 theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "top") ,
         points=T) [[1]] + 
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 


  xlab("Species richness gradient") + 
  
  ylab ("Functional Richness (FRic)")


# compare slopes

m.lst.FRic <- emtrends (model.ancova.FRic, "Dataset", var="SR")
# m.lst_tab.FRic <- summary(m.lst.FRic,point.est = mean)

# run model (ancova)
model.ancova.FEve <- brm (FEve ~ poly(SR,2)*Dataset,
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)

# plotting
p2<-plot(conditional_effects(model.ancova.FEve,
                             method="fitted",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "top") ,
         points=T) [[1]] + 
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 

  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (FEve)")

# organize plots
pdf(here("Output","Fig3_fish.pdf"), width=9,height=5)
grid.arrange(p1,p2,nrow=1)
dev.off()

# compare slopes
m.lst.FEve <- emtrends (model.ancova.FEve, "Dataset", var="SR")
# m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

# save results
save (model.ancova.FRic,
	m.lst.FRic,
	m.lst_tab.FRic,
	model.ancova.FEve,
	m.lst.FEve,
	m.lst_tab.FEve,
      file=here("Output", "GLM_test.RData"))

