
# Download Spatial Data Atlantic Forest
# from https://mauriciovancine.github.io/atlanticr/

# install.packages("remotes")
# remotes::install_github("mauriciovancine/atlanticr")
rm(list=ls())
require("atlanticr")
require(here)

# create dir to store rasters
dir.create ("Geo_data")
dir.create (here ("Geo_data", "spatial_atlantic"))

# read table with metrics
atlantic_spatial$metric_type
atlantic_spatial$metric_description[2]
atlantic_spatial$value[2]
atlantic_spatial$value_description[2]
atlantic_spatial$value_description[3]

# import data
atlanticr::atlantic_spatial_download(id = 458, path = here ("Geo_data", "spatial_atlantic")) # slope
atlanticr::atlantic_spatial_download(id = 457, path = here ("Geo_data", "spatial_atlantic")) # elevation
atlanticr::atlantic_spatial_download(id = 2, path = here ("Geo_data", "spatial_atlantic")) # land use

# load rasters with spatial data
require(raster)
slope<-raster("Geo_data/spatial_atlantic/458_atlantic_spatial_slope.tif")
plot(slope)
elevation<-raster("Geo_data/spatial_atlantic/457_atlantic_spatial_elevation.tif")
plot(elevation)
lulc<-raster("Geo_data/spatial_atlantic/002_atlantic_spatial_grouped_classes.tif")
plot(lulc==1)

# load community data to extract covariates
load( here ( "Processed_data","image_rodents.RData"))

# need to transform coordinates to match projection of AF data
spatial_data<-st_as_sf (spatial_effort_data_LF, coords = c("Longitude" ,"Latitude"), crs = 4326)
# reproject
spatial_data <- st_transform(spatial_data, crs = crs(elevation))
plot(elevation)
points(st_coordinates(spatial_data)[,1],st_coordinates(spatial_data)[,2])

# extract data for each coordinate
# slope
slope_ext <- terra::extract(slope, 
                            st_coordinates(spatial_data),
                            method = "simple",
                            cells=T,
                            buffer = 1000,
                            fun = mean ,
                            na.rm=T,
                            df=F)

# elevation
elevation_ext <- terra::extract(elevation, 
                         st_coordinates(spatial_data), 
                         method = "simple",
                            cells=T,
                            buffer = 1000,
                            fun = mean ,
                            na.rm=T,
                            df=F)

points(st_coordinates(spatial_data)[is.na(elevation_ext),1],st_coordinates(spatial_data)[is.na(elevation_ext),2],col="black",pch=19)

# spatial_effort_data_LF[is.na(elevation_ext),]
# We can remove these two sites : are inside biome Cerrado: Serra da Canastra and Serra do Cipé
spatial_effort_data_LF[95,]
# function to get the mode of land use in the buffer
get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}
lulc_ext <- terra::extract(lulc, 
                    st_coordinates(spatial_data),
                    method = "simple",
                    cells=T,
                    buffer = 1000,
                    #fun = get_mode ,
                    na.rm=T,
                    df=T)

# apply to all sites
sites <- unique(lulc_ext[,1])
lulc_extract <- lapply (sites, function (i) {
    # organize data
    res<-list (site = i ,
                forest = ifelse (is.na(table(lulc_ext [lulc_ext[,"ID"]==i,2])["1"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2]))),
                                0,table(lulc_ext [lulc_ext[,"ID"]==i,2])["1"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2]))),
                temp_crop = ifelse (is.na(table(lulc_ext [lulc_ext[,"ID"]==i,2])["5"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2]))),
                                0,table(lulc_ext [lulc_ext[,"ID"]==i,2])["5"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2]))),
                pastures = ifelse (is.na(table(lulc_ext [lulc_ext[,"ID"]==i,2])["4"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2]))),
                                0,table(lulc_ext [lulc_ext[,"ID"]==i,2])["4"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2]))),
                urban = ifelse (is.na(table(lulc_ext [lulc_ext[,"ID"]==i,2])["7"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2]))),
                                0,table(lulc_ext [lulc_ext[,"ID"]==i,2])["7"]/sum(table(lulc_ext [lulc_ext[,"ID"]==i,2])))
    )
  
    print(i)
    res
  }
                
)


# bind covariates to spatial data
spatial_data <- cbind (spatial_data,
                       forest = (sapply (lulc_extract, "[[", "forest")),
                       urban = (sapply (lulc_extract, "[[", "urban")),
                       pasture = (sapply (lulc_extract, "[[", "pastures")),
                       crops = (sapply (lulc_extract, "[[", "temp_crop")),
                       elevation = elevation_ext,
                       slope = slope_ext)
                       
# save spatial data
save (spatial_data, file = here ("Processed_data", "Spatial_covs_data.RData"))
rm(list=ls())



