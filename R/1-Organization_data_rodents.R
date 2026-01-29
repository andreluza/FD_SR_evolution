#test <- read.csv("https://raw.githubusercontent.com/andreluza/FD_SR_evolution/main/data/Atributos_especies_Atlantico_%26_Pacifico_Oriental_2020_04_28.csv",
#                 sep=";")

# -------------------------------------------------------

# Evolution of FEve-SR and FRic - SR relationships - organize rodent data

# ---------------------------------------------------------

# STEP 1: ORGANIZE RODENT DATA

# load functions & packages
rm(list=ls())
source("R/functions.R")
source("R/packages.R")

# Atlantic Rainforest Sites
Atl_For <- vect(here ("data_rodents", "limites_integradores_wgs84_v1_2_0", "ma_limite_integrador_muylaert_et_al_2018_wgs84.shp"))


# ----------------------------------------------------------------------
# Load community data (Luza et al. 2019)
mammals_luza <- read.csv(here("data_rodents",
                         "AppendixS1- Small_mammal_data.csv"),
                    h=T, sep=";",fileEncoding="latin1")

# Neotropical communities
NT_mammals_luza<- mammals_luza[which(mammals_luza$WWF_REALM2 == "Neotropic"),]
# Rodents only
NT_mammals_luza<- NT_mammals_luza[which(NT_mammals_luza$ORDER %in% c("Rodentia")),]

# solve spp name (oxymycterus sp is now O quaestor (Pecanha et al. 2020))
NT_mammals_luza [grep("Luza",NT_mammals_luza$REFERENCE),"SPECIES"][which(NT_mammals_luza [grep("Luza",NT_mammals_luza$REFERENCE),"SPECIES"] == "Oxymycterus_sp.")] <- "Oxymycterus_quaestor"

# assemblage a dataset with abundance and effort 
sel_cols_luza <- c("REFERENCE", "SITE", "LAT", "LONG", "SPECIES","NUMBER_OF_RECORDS","EFFORT_PER_HABITAT")
NT_mammals_luza <- NT_mammals_luza [,which(colnames(NT_mammals_luza) %in% sel_cols_luza)]

# Keep only sites from the Atlantic Forest
# with raster data, xy=TRUE works
# transform points into spatialVect
pts <- terra::vect(NT_mammals_luza,
                   geom=c("LONG", "LAT"),
                   crs= "+proj=longlat +datum=WGS84")
idx <- relate(pts, Atl_For, "intersects")[,1]  # TRUE se está dentro
Luza_coords_inside <- NT_mammals_luza[idx, ]
# check
plot(Atl_For)
points(Luza_coords_inside$LONG,Luza_coords_inside$LAT,pch=19,col="red")

# index for reference and site (many sites sampled per reference) 
Luza_coords_inside$SAMPLEID <- paste(Luza_coords_inside$REFERENCE,
                                     Luza_coords_inside$SITE)

# -----------------------------------
# load community data (Figueiredo et al.)
mammals_figueiredo <- read.csv(here("data_rodents",
                              "Mammal_Communities.csv"),
                         h=T, sep=";")
# rm NAs
mammals_figueiredo<- mammals_figueiredo[is.na(mammals_figueiredo$SampleID)!=T,]

# set "_" in spp names
mammals_figueiredo$Valid_Species<- gsub (" ","_",mammals_figueiredo$Valid_Species)

# Deltamys_sp._nov. is Deltamys araucaria (Quintela et al. 2017)
mammals_figueiredo$Valid_Species[which(mammals_figueiredo$Valid_Species == "Deltamys_sp._nov.")] <- "Deltamys_araucaria"

# coordinates are here
localities_f <-  read.csv(here("data_rodents",
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
# add points to the map
points(mammals_figueiredo$Longitude,mammals_figueiredo$Latitude,col="green",pch=1,cex=1.2)


# remove References in both data sets (remove from Luza) -----------------

Luza_coords_inside<-Luza_coords_inside [-grep("De_La_Sancha_2014",Luza_coords_inside$SAMPLEID),] # remove De La Sancha et al. 2014
Luza_coords_inside<-Luza_coords_inside [-grep("Passamani_&_Fernandez_2011",Luza_coords_inside$SAMPLEID),] # remove Passamani_&_Fernandez_2011
Luza_coords_inside<-Luza_coords_inside [-grep("Pardini_2004",Luza_coords_inside$SAMPLEID),] # remove Pardini_2004
Luza_coords_inside<-Luza_coords_inside [-grep("Stevens_&_Husband_1998",Luza_coords_inside$SAMPLEID),] # remove Stevens_&_Husband_1998
Luza_coords_inside<-Luza_coords_inside [-grep("De_Carvalho_Braga_et_al_2015 road",Luza_coords_inside$SAMPLEID),] # remove De_Carvalho_Braga_et_al_2015 road
unique(Luza_coords_inside$SAMPLEID)
points(Luza_coords_inside$LONG,Luza_coords_inside$LAT,col="blue",pch=1,cex=1.5)

# Obtain effort and community matrix for each data set -- site x species matrix ------------------

# aggregate data per study and site
## effort (average per habitat (the scale of record in the database))
df_av_effort <- aggregate (Luza_coords_inside, 
                           by = list (Luza_coords_inside$SAMPLEID),
                        FUN = mean,na.rm=T)

# select interesting cols & remove NA effort
df_av_effort <- df_av_effort[is.na(df_av_effort$EFFORT_PER_HABITAT) != T,
             c("Group.1", "EFFORT_PER_HABITAT","LAT","LONG")]


# community matrix from Luza et al.
tab_sp_site<-cast(formula = SAMPLEID  ~ SPECIES,
                  value="NUMBER_OF_RECORDS",
                  data=Luza_coords_inside,
                  fun.aggregate=sum,
                  na.rm=T)

# rm study lacking effort data
tab_sp_site <- tab_sp_site[which(tab_sp_site$SAMPLEID %in% df_av_effort$Group.1),]

# check
# tab_sp_site$SAMPLEID == df_av_effort$Group.1
# captures per trapping unit
tab_sp_site <- tab_sp_site[,-1] / df_av_effort$EFFORT_PER_HABITAT


# community table from Figueiredo et al. ------------
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


# plot Atlantic forest sites
plot(Atl_For)
points(Luza_coords_inside$LONG,Luza_coords_inside$LAT,col="red",cex=1,pch=19)
points(localities_f$Longitude,localities_f$Latitude)

# ------------------------------------------- #
# load phylogenies
# Load the fully resolved phylogenies
tree_list <- tree <- read.nexus(file=here("data_rodents",
                                          "Sigmodontinae_413species100Trees.trees"))
# Adjusting the names

tree_list <- lapply (tree_list, function (i)
  
      {i$tip.label <- gsub ("__CRICETIDAE__RODENTIA", "",i$tip.label);
      
      i}
      
      )

# traits
traits<- read.csv (here("data_rodents", "Penone_et_al_2016_mammal_trait_data_imputed.csv"),
                   sep=",")
cor (traits$MaxLifepsan,traits$AdultHeadBodyLen.mm,use="complete.obs")

# standradize traits
std_traits <- data.frame (Body.mass.g = scale(traits$Body.mass.g),
                          LitSz = scale(traits$LitSz),
                          LitPerYear = scale(traits$LitPerYear),
                          MaxLifepsan = scale(traits$MaxLifepsan.m),
                          AdultHeadBodyLen.mm = scale(traits$AdultHeadBodyLen.mm)
                          )
rownames(std_traits)<- traits$IUCN.binomial                          

# percentage of missing entries
table(is.na(std_traits))[2]/sum(table(is.na(std_traits)))

# imputation without phylogeny
#std_traits <- missForest (std_traits, maxiter = 50,
#                          ntree= 100,variablewise = T)
#std_traits<-std_traits$ximp

# rm spp missing in the trait dataset (could be non identified species, higher taxa, etc)
combined_data <- combined_data[,which(colnames(combined_data) %in%
                                        rownames(std_traits))]
# rm spp missing in the phylogeny (could be non identified species, higher taxa, etc)
combined_data <- combined_data[,which(colnames(combined_data) %in% 
                                        tree_list$UNTITLED$tip.label)]

# check colnames (only Sigmodontinae there)
colnames(combined_data)

# removing communities lacking spp
comm_enough_spp <- rowSums(combined_data>0)
# remove communities with less than 3 spp
combined_data <- combined_data [which(comm_enough_spp>=3),]

# do the same for spatial geo data
spatial_effort_data_LF<- spatial_effort_data_LF[which(comm_enough_spp>=3),]

# removing spp that were in the sites just removed # Not the case
(spp_enough <- colSums(combined_data>0))
# remove
combined_data <- combined_data [,which(spp_enough>0)]

# matching datasets -------------------------------------------------
# match phylogenetic and trait data

match_data <- lapply (tree_list, function (i) 
  
  match.phylo.data(i, 
                   std_traits)
)

# match phylogeny and community
match_comm_data<-lapply(match_data, function (i)
  
  match.phylo.comm(i$phy,
                   combined_data ) 
)

## all phylogenies have the same order of tiplabels
# thus we used the first community data
# to evaluate this run:
# table(match_comm_data[[1]]$comm == match_comm_data[[2]]$comm)
# and
# table(match_comm_data[[10]]$comm == match_comm_data[[100]]$comm)

# match trait and community
std_traits_subset <- std_traits[which(rownames(std_traits) %in% match_comm_data[[1]]$phy$tip.label),]

# ordering of species names
std_traits_subset<- std_traits_subset[order(rownames(std_traits_subset)),]
rownames(match_data$UNTITLED$data) == match_data$UNTITLED$phy$tip.label

# test of phylogenetic signal 
psignal <- lapply (seq(1,length(match_data)), function (k)
  
  lapply (seq(1,ncol (match_data[[k]]$data)), function (i)
    
    phylosig(match_data[[k]]$phy, 
             match_data[[k]]$data[,i], 
             method="K", test=TRUE, nsim=999)
  ))

# df with res
psignal <- lapply (psignal, function (k) do.call (rbind, 
                                                  
                                                  lapply (k, function (i)
                                                    
                                                    data.frame (K=i$K,
                                                                pval=i$P)
                                                  )
))

# signal K
apply(sapply (psignal, "[[", "K"),1,mean)
apply(sapply (psignal, "[[", "K"),1,sd)

# rownames(std_traits_subset) == match_comm_data$UNTITLED$phy$tip.label

# create dir for processed data
dir.create("Processed_data")

# save objects needed for the next step (calculation of Functional Indices, simulations)
save(match_comm_data,
     std_traits_subset,
     match_data,
     tree_list,
     spatial_effort_data_LF,
     combined_data,
     file = here ( "Processed_data","image_rodents.RData"))

# save.image(here ( "Output","image_rodents.RData"))
