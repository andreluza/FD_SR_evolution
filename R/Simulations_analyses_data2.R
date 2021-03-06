
# -------------------------------------------------------
# Evolution of FEve-SR and FRic - SR relationships
# ---------------------------------------------------------

# load functions & packages
source("R/functions.R")
source("R/packages.R")

# ----------------------------------------------------------------------
# load data

# load community data (luza et al.)
mammals_luza <- read.csv(here("data2",
                         "AppendixS1- Small_mammal_data.csv"),
                    h=T, sep=";")

# neotropical communities
NT_mammals_luza<- mammals_luza[which(mammals_luza$WWF_REALM2 == "Neotropic"),]
# rodents & marsupials
NT_mammals_luza<- NT_mammals_luza[which(NT_mammals_luza$ORDER %in% c("Rodentia")),]

# solve spp name (oxymycterus sp is now O quaestor (Pecanha et al. 2020))
NT_mammals_luza [grep("Luza",NT_mammals_luza$REFERENCE),"SPECIES"][which(NT_mammals_luza [grep("Luza",NT_mammals_luza$REFERENCE),"SPECIES"] == "Oxymycterus_sp.")] <- "Oxymycterus_quaestor"

# assemblage a dataset with abundance and effort 
sel_cols_luza <- c("REFERENCE", "SITE", "LAT", "LONG", "SPECIES","NUMBER_OF_RECORDS","EFFORT_PER_HABITAT")
NT_mammals_luza <- NT_mammals_luza [,which(colnames(NT_mammals_luza) %in% sel_cols_luza)]

# aggregate data per study and site
## effort (average per habitat (the scale of record in the database))
NT_mammals_luza$SAMPLEID <- paste(NT_mammals_luza$REFERENCE,NT_mammals_luza$SITE)
df_av_effort <- aggregate (NT_mammals_luza, by = list (NT_mammals_luza$SAMPLEID),
           FUN = mean,na.rm=T)
# select interesting cols & remove NA effort
df_av_effort <- df_av_effort[is.na(df_av_effort$EFFORT_PER_HABITAT) != T,
             c("Group.1", "EFFORT_PER_HABITAT","LAT","LONG")]

# community matrix
tab_sp_site<-cast(formula = SAMPLEID  ~ SPECIES,
                  value="NUMBER_OF_RECORDS",
                  data=NT_mammals_luza,
                  fun.aggregate=sum,
                  na.rm=T)

# rm study lacking effort data
tab_sp_site <- tab_sp_site[which(tab_sp_site$SAMPLEID %in% df_av_effort$Group.1),]

# check
# tab_sp_site$SAMPLEID == df_av_effort$Group.1
# captures per trapping unit
tab_sp_site <- tab_sp_site[,-1] / df_av_effort$EFFORT_PER_HABITAT

# -----------------------------------
# load community data (figueiredo et al.)
mammals_figueiredo <- read.csv(here("data2",
                              "Mammal_Communities.csv"),
                         h=T, sep=";")
# rm NAs
mammals_figueiredo<- mammals_figueiredo[is.na(mammals_figueiredo$SampleID)!=T,]

# set "_" in spp names
mammals_figueiredo$Valid_Species<- gsub (" ","_",mammals_figueiredo$Valid_Species)

# Deltamys_sp._nov. is Deltamys araucaria (Quintela et al. 2017)
mammals_figueiredo$Valid_Species[which(mammals_figueiredo$Valid_Species == "Deltamys_sp._nov.")] <- "Deltamys_araucaria"

# coordinates are here
localities_f <-  read.csv(here("data2",
                               "localities.csv"),
                          h=T, sep=";")
# rm NAs
localities_f<- localities_f[is.na(localities_f$SampleID)!=T,]

# interesting cols
sel_cols_figueiredo <- c("SampleID", "Latitude", "Longitude","Sampling_effort")
localities_f <- localities_f[which(localities_f$SampleID %in% unique(mammals_figueiredo$SampleID)),
                             sel_cols_figueiredo]

# matching
localities_f_match <- localities_f [match (mammals_figueiredo$SampleID, localities_f$SampleID),]

# bind
mammals_figueiredo <- cbind(mammals_figueiredo,
                            localities_f_match) 
# clean columns of garbage
mammals_figueiredo <- mammals_figueiredo[,-grep("X",colnames(mammals_figueiredo))]

# community table
tab_sp_site_f <- cast (SampleID ~ Valid_Species, 
                        data=mammals_figueiredo,
                        value = "Abundance",
                        fun.aggregate = sum,
                        na.rm=T)

# effort per effort unit
# check
# tab_sp_site_f$SampleID == localities_f$SampleID
tab_sp_site_f <- tab_sp_site_f[,-1] / localities_f$Sampling_effort

# ---------------------------------------------------
# bind things
colnames(df_av_effort)<- c("SampleID", "Sampling_effort","Latitude","Longitude")
# organize cols
df_av_effort[,match(colnames(localities_f),colnames(df_av_effort))]

# bind 
spatial_effort_data_LF <- rbind(df_av_effort,
                                localities_f)

# bind community matrix
# spp of figuiredo not in luza
not_in_luza<-colnames(tab_sp_site_f)[which(colnames(tab_sp_site_f) %in% colnames(tab_sp_site) == F)]
# table
not_in_luza_tab <- matrix(0, ncol = length(not_in_luza),
                             nrow=nrow(tab_sp_site),
                          dimnames=list(rownames(tab_sp_site),
                                        not_in_luza))
# not in figueiredo
not_in_f<-colnames(tab_sp_site)[which(colnames(tab_sp_site) %in% colnames(tab_sp_site_f) == F)]
# table
not_in_f_tab <- matrix(0, ncol = length(not_in_f),
                          nrow=nrow(tab_sp_site_f),
                          dimnames=list(rownames(tab_sp_site_f),
                                        not_in_f))

# cbind to each dataset
tab_sp_site<-cbind(tab_sp_site,not_in_luza_tab)
tab_sp_site_f<-cbind(tab_sp_site_f,not_in_f_tab)
# set in order
tab_sp_site<-tab_sp_site[,order(colnames(tab_sp_site))]
tab_sp_site_f<-tab_sp_site_f[,order(colnames(tab_sp_site_f))]

# collate datasets
combined_data <- rbind(tab_sp_site,
                       tab_sp_site_f)

# ------------------------------------------- #
# load phylogenies
# Load the fully resolved phylogenies
tree_list <- tree <- read.nexus(file=here("data2",
                                          "Sigmodontinae_413species100Trees.trees"))
# Adjusting the names

tree_list <- lapply (tree_list, function (i)
  
      {i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);
      
      i}
      
      )


# traits
traits<- read.csv (here("data2", "Penone_et_al_2016_mammal_trait_data_imputed.csv"),
                   sep=",")
# standradize traits
std_traits <- data.frame (Body.mass.g = scale(traits$Body.mass.g),
                          LitSz = scale(traits$LitSz),
                          LitPerYear = scale(traits$LitPerYear),
                          MaxLifepsan = scale(traits$MaxLifepsan.m),
                          PopDen.n.km2 = scale(traits$PopDen.n.km2))
rownames(std_traits)<- traits$IUCN.binomial                          

# rm spp not in trait dataset
combined_data <- combined_data[,which(colnames(combined_data) %in%
                                        rownames(std_traits))]
# rm spp not in phylogeny
combined_data <- combined_data[,which(colnames(combined_data) %in% 
                                        tree_list$UNTITLED$tip.label)]

# removing communities lacking spp
comm_enough_spp <- rowSums(combined_data>0)
# rm
combined_data <- combined_data [which(comm_enough_spp>=3),]

# do the same for spatial geo data
spatial_effort_data_LF<- spatial_effort_data_LF[which(comm_enough_spp>=3),]

# removing zero spp
spp_enough <- colSums(combined_data>0)
# rm
combined_data <- combined_data [,which(spp_enough>0)]

# matching datasets
# match phylogenetic and trait data

match_data <- lapply (tree_list, function (i) 
  
  match.phylo.data(i, 
                   std_traits)
)

# match phylgeny and community
match_comm_data<-lapply(match_data, function (i)
  
  match.phylo.comm(i$phy,
                   combined_data ) 
)

## all phylogenies have the same order of tiplabels
# thus we used the first community data
# check this out
# table(match_comm_data[[1]]$comm == match_comm_data[[2]]$comm)
# or 
# table(match_comm_data[[10]]$comm == match_comm_data[[100]]$comm)

# match trait and community
std_traits_subset <- std_traits[which(rownames(std_traits) %in% match_comm_data[[1]]$phy$tip.label),]

# ordering 
std_traits_subset<- std_traits_subset[order(rownames(std_traits_subset)),]

# run
empirical_FD <- dbFD(x=std_traits_subset,
       a=match_comm_data[[1]]$comm,
       w.abun=T,
       stand.x=F,
       calc.FRic = T,
       stand.FRic = T,
       corr = "lingoes",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)


## Simulate trait evolution according to a bivariate "BMM" model
# Number of traits
ntraits<-ncol(std_traits_subset)
# Number of simulated (pairs of) traits
nsim<-10
# sigmas
sigma<-(rbind(c(1,0.1,0.1,0.1,0.1),
              c(0.1,1,0.1,0.1,0.1),
              c(0.1,0.1,1,0.1,0.1),
              c(0.1,0.1,0.1,1,0.1),
              c(0.1,0.1,0.1,0.1,1)))

# ancestral states for each traits
theta<-rep(0,ntraits)

# Simulate

simul<-lapply (tree_list, function (i) 
  
        mvSIM(i,
             nsim=nsim, 
             model="BM1",
             param=list(sigma=sigma, 
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
simulated_FD <- lapply (seq (1,length (match_comm_data)), function (i)
  
                dbFD(x=mean_simul[[i]][which(rownames(mean_simul[[i]]) %in% colnames(match_comm_data[[i]]$comm)),],
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

# ----------------------------------------------
# niche filling (early burst)

# sigmas
betas<-(rbind(c(-0.5,0.1,0.1,0.1,0.1),
              c(0.1,-0.5,0.1,0.1,0.1),
              c(0.1,0.1,-0.5,0.1,0.1),
              c(0.1,0.1,0.-0.5,1,0.1),
              c(0.1,0.1,0.1,0.1,-0.5)))

# Simulate Eb

simul_EB<-lapply (tree_list, function (i) 
  
  mvSIM(i,
        nsim=nsim, 
        model="EB",
        param=list(sigma=sigma, 
                   beta=betas,
                   theta=theta,
                   ntraits=ntraits,
                   names_traits=c("Trait 1",
                                  "Trait 2",
                                  "Trait 3",
                                  "Trait 4",
                                  "Trait 5"))))


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
       corr = "lingoes",
       calc.CWM = F,
       calc.FDiv=F,
       print.pco = T)
  
)

# OU
# Simulate Eb
alpha <- (rbind(c(1,0.1,0.1,0.1,0.1),
                c(0.1,1,0.1,0.1,0.1),
                c(0.1,0.1,1,0.1,0.1),
                c(0.1,0.1,0.1,1,0.1),
                c(0.1,0.1,0.1,0.1,1)))

simul_OU<-lapply (tree_list, function (i) 
  
  mvSIM(i,
        nsim=nsim, 
        model="OU1",
        param=list(sigma=sigma, 
                   beta=betas,
                   theta=theta,
                   alpha=alpha,
                   ntraits=ntraits,
                   names_traits=c("Trait 1",
                                  "Trait 2",
                                  "Trait 3",
                                  "Trait 4",
                                  "Trait 5"))))


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

save.image("image_rodents.RData")

# -----------------------------------------------------
# empirical results
empirical_results <- data.frame (SR= empirical_FD$nbsp,
                                 FRic= empirical_FD$FRic,
                                 FEve=empirical_FD$FEve,
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
pdf(here("Output","Fig3_SM.pdf"), width=9,height=5)
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
      file=here("Output", "GLM_test_SM.RData"))

