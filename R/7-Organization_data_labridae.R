
# -------------------------------------------------------

# Evolution of FEve-SR and FRic - SR relationships - organizing data from Labridae family

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

# load community data
# UVC fish data
peixes <- read.csv(here("data_fish","UpdatedData_RMorais_et_al_2017.csv"))
# 
tapply (X = peixes$ScientificName,
              INDEX = peixes$Family,
              FUN=function (x) length(unique(x))) # haemulidae  et labridae
tapply (X = peixes$IndCounting,
              INDEX = peixes$Family,
              FUN=function (x) length(unique(x))) # haemulidae  et labridae (similar N records)

# subset Labridae
peixes <- peixes [which(peixes$Family %in% "labridae"),]

# dataframe with taxa name
df_taxa_community <-data.frame (sp = unique(peixes$ScientificName)[order(unique(peixes$ScientificName))])
df_taxa_community$sp <- gsub ("\\.", " ", df_taxa_community$sp)

# use validation from the analysis with the complete set of spp
# naming
#names (worms_record_fish_community) <- df_taxa_community$sp
test_taxa_comm <- lapply (worms_record_fish_community,data.frame)
test_taxa_comm <- test_taxa_comm[unlist(lapply(test_taxa_comm,nrow))>=1]
test_taxa_comm <- test_taxa_comm[unlist(lapply(test_taxa_comm,ncol))>1]
test_taxa_comm<-do.call(rbind,test_taxa_comm)


# match with the table
df_taxa_community$sp_worms <- (test_taxa_comm[match (df_taxa_community$sp,
                                                     tolower(test_taxa_comm$scientificname)),
                                              "scientificname"])

# if missing, keep the previous
df_taxa_community$names_sp <-  ifelse (is.na(df_taxa_community$sp_worms),
                                       df_taxa_community$sp,
                             df_taxa_community$sp_worms)


# match with dataset

peixes$validScientificName <- df_taxa_community$names_sp [match (gsub ("\\."," ", peixes$ScientificName),
                                  df_taxa_community$sp)]





# ======================================================


## modify eventID to rm year
eventID_MOD  <- paste (peixes$Region, peixes$Locality,peixes$Site, 
                       sep = "_")

peixes$eventID_MOD <- eventID_MOD


# number of belt transects
length(unique(peixes$Transect_id))

# number of localities
length(unique(peixes$Locality))

# obtain table 
tab_sp_site<-cast(formula = eventID_MOD ~ validScientificName,
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
traits_peixes <- read.csv(here("data_fish","Atributos_especies_Atlantico_&_Pacifico_Oriental_2020_04_28.csv"),
                          h=T,sep=";")


# select Labridae
traits_peixes<- traits_peixes[which(traits_peixes$Family %in% c("scaridae","labridae")),]

# adjust names to match community, trait, and phylogeny
# trait
#traits_peixes$Name <- tolower(gsub(" ",".",traits_peixes$Name)) 
#traits_peixes<- traits_peixes[duplicated(traits_peixes$Name) !=T,]
rownames(traits_peixes) <- traits_peixes$Name

# check all this!!!! 14 spp missing is too much

# calculate functional metrics
## subset of fish traits and communities
traits_peixes [grep ("Nicholsina", traits_peixes$Name),"Name"] <- "Nicholsina usta collettei" # adjust this species to match
colnames(tab_sp_site)[(which(gsub ("\\."," ", colnames(tab_sp_site)) %in% traits_peixes$Name == F))] # only higher taxa

# table
tab_sp_site <- tab_sp_site[,which(gsub ("\\."," ", colnames(tab_sp_site)) %in% traits_peixes$Name)] # spp in the trait dataset
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
# proportion of missing data
table(is.na(std_traits))[2]/sum(table(is.na(std_traits)))

# imput
require(missForest)
std_traits <- missForest (std_traits, maxiter = 50,
                          ntree= 100,variablewise = T)
std_traits<-std_traits$ximp
  
# match spp names in trait and community dataset
#std_traits <- std_traits [match(colnames(tab_sp_site),rownames(std_traits)),]
#rownames(std_traits) == colnames(tab_sp_site)

## capture per trapping effort (number of transects per site)
tab_sp_site <- (tab_sp_site / effort_site$effort)

# adjust trait names
rownames(std_traits) <- firstup (gsub ("\\.", " ",rownames(std_traits)))
colnames(tab_sp_site)<-firstup (gsub ("\\.", " ",colnames(tab_sp_site)))

# table(colnames(tab_sp_site) %in% tree$tip.label)
# finally, match phylogeny, traits, and community
# match phylogenetic and trait data

match_data <- lapply (test_tree, function (i) 
  
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
  
  
  match_data[[i]]$data[which(rownames(match_data[[i]]$data) %in% 
                               colnames(match_comm_data[[i]]$comm)),]

)

# subset comm data  
subset_comm_data <- lapply (seq(1,length(match_comm_data)), function (i) 
  
  match_comm_data[[i]]$comm[,which(colnames(match_comm_data[[i]]$comm) %in% rownames(subset_trait_data[[i]]))]
  
)

# phylogenetic signal (on all species from SW Atlantic)

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
(apply(sapply (psignal, "[[", "K"),1,mean))
apply(sapply (psignal, "[[", "K"),1,sd)

# save
save(match_comm_data,
     subset_trait_data,
     match_data,
     test_tree,
     subset_comm_data,
     effort_site,
     file = here ( "Processed_data","image_labridae.RData"))

# end
rm(list=ls())

