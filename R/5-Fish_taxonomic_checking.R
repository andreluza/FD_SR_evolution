
# -------------------------------------------------------
# Evolution of FEve-SR and FRic - SR relationships

# First to a taxonomic checking to correspond taxonomic nomenclature between phylogeny and fish data

# ---------------------------------------------------------

# load functions & packages
source("R/functions.R")
source("R/packages.R")

# ----------------------------------------------------------------------
# load data
# load phylogenies
#fishtree_complete_phylogeny()
tree<- read.tree (here ("data","TACT","Reef_fish_all_combined.trees"))#fishtree_complete_phylogeny()

# dataframe with taxa name
df_taxa <-data.frame (sp = (unique(tree[[1]]$tip.label[order(tree[[1]]$tip.label)])))
df_taxa$sp <- gsub ("_", " ", df_taxa$sp)

# worms's validation
worms_record_fish <- lapply (df_taxa$sp, function (i) 
  
  tryCatch (
    
    wm_records_taxamatch(i, fuzzy = TRUE, marine_only = TRUE)[[1]],
    
    error = function (e) print(NA)
    
    
  )
  
)

names (worms_record_fish) <- df_taxa$sp
test_taxa <- lapply (worms_record_fish,data.frame)
test_taxa <- test_taxa[unlist(lapply(test_taxa,nrow))>=1]
test_taxa <- test_taxa[unlist(lapply(test_taxa,ncol))>1]
test_taxa<-do.call(rbind,test_taxa)

# match with the table
df_taxa$sp_worms <- (test_taxa[match (df_taxa$sp,rownames(test_taxa)),"scientificname"])

# if missing, keep the previous name
df_taxa$names_sp <-  ifelse (is.na(df_taxa$sp_worms),
        df_taxa$sp,
        df_taxa$sp_worms)


# filter labridae
df_taxa <- df_taxa[which(df_taxa$family == "Labridae"),]

# alter taxon names in the phlogeny
table(gsub (" ","_",df_taxa$sp [(match (tree[[2]]$tip.label,
                                  gsub (" ","_",df_taxa$sp)))]) == tree[[2]]$tip.label)


# change tipnames by valid names
test_tree <- lapply (tree, function (i){

  i$tip.label <- df_taxa$names_sp [(match (i$tip.label,
                                        gsub (" ","_",df_taxa$sp)))]
  i

})
# table(gsub ("_"," ",tree[[3]]$tip.label) == test_tree[[3]]$tip.label)

# load community data
# UVC fish data
peixes <- read.csv(here("data","UpdatedData_RMorais_et_al_2017.csv"))

# dataframe with taxa name
df_taxa_community <-data.frame (sp = unique(peixes$ScientificName)[order(unique(peixes$ScientificName))])
df_taxa_community$sp <- gsub ("\\.", " ", df_taxa_community$sp)

# worms's validation
worms_record_fish_community <- lapply (df_taxa_community$sp, function (i) 
  
  tryCatch (
    
    wm_records_taxamatch(i, fuzzy = TRUE, marine_only = TRUE)[[1]],
    
    error = function (e) print(NA)
    
    
  )
  
)
# naming
names (worms_record_fish_community) <- df_taxa_community$sp
test_taxa_comm <- lapply (worms_record_fish_community,data.frame)
test_taxa_comm <- test_taxa_comm[unlist(lapply(test_taxa_comm,nrow))>=1]
test_taxa_comm <- test_taxa_comm[unlist(lapply(test_taxa_comm,ncol))>1]
test_taxa_comm<-do.call(rbind,test_taxa_comm)


# match with the table
df_taxa_community$sp_worms <- (test_taxa_comm[match (df_taxa_community$sp,rownames(test_taxa_comm)),
                                              "scientificname"])

table (df_taxa_community$sp_worms %in% df_taxa$sp_worms)

table(df_taxa_community$sp_worms %in% traits_peixes$Name )
table(df_taxa$sp_worms %in% traits_peixes$Name )


# if missing, kepp the previous
df_taxa_community$names_sp <-  ifelse (is.na(df_taxa_community$sp_worms),
                                       df_taxa_community$sp,
                             df_taxa_community$sp_worms)


# match with dataset

peixes$validScientificName <- df_taxa_community$names_sp [match (gsub ("\\."," ", peixes$ScientificName),
                                  df_taxa_community$sp)]

# save worms data
save (worms_record_fish, worms_record_fish_community , 
      file = here ("output", "tax_validation_fish.RData"))


