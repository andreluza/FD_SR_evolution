
# -------------------------------------------------------
# Evolution of FEve-SR and FRic - SR relationships

# Obtain and extract covariates from Bioracle and reef area map (Magris et al. 2020)

# ---------------------------------------------------------

# load functions & packages
source("R/functions.R")
source("R/packages.R")

# ----------------------------------------------------------------------

# load community data
# UVC fish data
peixes <- read.csv(here("data_fish","UpdatedData_RMorais_et_al_2017.csv"))

# save worms data
load (file = here ("output", "tax_validation_fish.RData"))

# -----------------------

#  BioOracle 
# code to download and save in the folder "environment"
# BiO Oracle - extracting covariate data
# Explore datasets in the package
devtools::install_github("lifewatch/sdmpredictors")
require(sdmpredictors)
layers <- list_layers()
# View (layers [grep ("Bio-ORACLE",layers$dataset_code),])
# Download specific layers to the current directory
# set preferred folder (to download data)
options(sdmpredictors_datadir=here ("Geo_data"))
# chlorophil has different extent - loading and extracting in two steps         
layers_oracle <- load_layers(c("BO2_tempmean_ss",
                               "BO2_ppmean_ss", 
                               "BO2_salinitymean_ss", 
                               "BO_damean"
                              ))

# list and process
require(raster)
biooracle_data <-  (list.files (here ("Geo_data"),pattern = ".tif"))
biooracle_data <- lapply (biooracle_data, function (i) 
  
    raster (here ("Geo_data",i))
)
biooracle_data<- stack (biooracle_data)#stack

## geo coordinates
require(dplyr)
coordinates_sites <- peixes  %>% 
  
  group_by(Region ,Locality, Site) %>%

    summarise(decimalLongitude = mean(Lon),
              decimalLatitude = mean(Lat))


# coord to sppoints df to extract
spdf <- SpatialPointsDataFrame(coords = coordinates_sites[,4:5], data = coordinates_sites,
                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

## extracting data
extracted_sea_data <- raster::extract (biooracle_data, 
                                       spdf,method='bilinear', 
                               fun=mean)

# Covariate data
reef_covariates <- cbind (coordinates_sites, 
                          extracted_sea_data)
# check
# temperature
plot(reef_covariates$decimalLatitude,
     reef_covariates$Present.Surface.Temperature.Mean)
# salinity
plot(reef_covariates$decimalLongitude,
     reef_covariates$Present.Surface.Salinity.Mean)
# productivity
plot(reef_covariates$decimalLatitude,
     reef_covariates$Present.Surface.Primary.productivity.Mean)
# turbidity
plot(reef_covariates$decimalLongitude,
     reef_covariates$BO_damean_lonlat)

# ------------------------------
# Reef area

shapefiles<-list.files(here ("Geo_data","magris_reef_map"),pattern=".shp")

# rm cumulative impacts
shapefiles<-shapefiles[-grep("pu", tolower (shapefiles))]

# load all at once
require(sf)
shapes <- lapply (shapefiles, function (shp) 
  
  read_sf (here ("Geo_data","magris_reef_map"), layer = gsub(".shp","",shp)) 
  
)

# subset of habitats
list_habitats <- list ("amazon" = "AO11",# amazon
                       "eastern"= c("EC11",# eastern
                                    "EC12",
                                    "EC13",
                                    "EC14",
                                    "EC15",
                                    "EC16"),
                       "noronha" = c("FS12", "FS13"),#noronha
                       "northeastern" = c("NC11",# northeastern
                                          "NC12",
                                          "NC13"),
                       "riogrande"="RC11", # rio grande
                       "southeastern"="SC11", #southeastern
                       "trindade" = "TS12"
)

# reef location
BR_reefs <- lapply (seq (1,length(shapes)), function (shp)
  # extract codes  
  shapes[[shp]][which(shapes[[shp]]$habitat %in% list_habitats[[shp]]),]
  
)

# penedos
spsp<-read_sf (here ("Geo_data","magris_reef_map","SPSP"),"PS11_2") 

# bind 
BR_reefs <- rbind(BR_reefs[[1]],
                 BR_reefs[[2]],
                 BR_reefs[[3]],
                 BR_reefs[[4]],
                 BR_reefs[[5]],
                 BR_reefs[[6]],
                 BR_reefs[[7]])

plot(spsp[1])
plot(BR_reefs[1])

#crs(BR) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
#BR <- spTransform(BR, CRS("+init=epsg:4326"))
BR_reefs <- st_transform(BR_reefs, crs = "+proj=longlat +datum=WGS84 +no_defs") # assign crs
spsp <- st_transform(spsp, crs = "+proj=longlat +datum=WGS84 +no_defs") # assign crs

# use dist2Line from geosphere - only works for WGS84 
sp_data <- st_as_sf(spdf, coords = c("decimalLongitude", "decimalLatitude"))
sf_use_s2(F) 
sp_data <- (st_buffer(sp_data,dist = 0.1))
plot(sp_data[1])

# extract
extract_area <- st_intersection (BR_reefs,sp_data) %>%
  mutate (area = st_area(.),
          prop = st_area(.)/mean(st_area(sp_data),na.rm=T)) %>%
  group_by(Region, Locality, Site) %>%
  reframe (prop = mean(prop,na.rm=T))

# check
extract_area [which(extract_area$prop == max(extract_area$prop)),]
dim(extract_area)
range(extract_area$prop)

# penedos which is different
extract_area_penedos <- st_intersection (spsp,sp_data)%>%
  mutate (area = st_area(.),
          prop = st_area(.)/mean(st_area(sp_data))) %>%
  group_by(Region, Locality, Site) %>%
  reframe (prop = mean(prop))

# reef area
reef_area <- rbind(extract_area,
                   extract_area_penedos)
                    # missing - likely zeroes
reef_area <- rbind (reef_area, 
                    st_drop_geometry (cbind (sp_data[which(paste (sp_data$Region, sp_data$Locality, sp_data$Site) %in% 
                                                            paste (reef_area$Region, reef_area$Locality, reef_area$Site) == F),1:3],
                                             prop = 0)))
# match order
reef_area <- reef_area [match(paste (sp_data$Region, sp_data$Locality, sp_data$Site),paste (reef_area$Region, reef_area$Locality, reef_area$Site)),]

# check order
reef_covariates$Site == reef_area$Site

# bind reef area 
reef_covariates$Reef.Area <- reef_area$prop

save (reef_covariates,
      file = here ("Processed_data","Spatial_covs_data_reefs.RData"))

# end

