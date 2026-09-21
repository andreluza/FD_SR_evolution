
# data organization
require(openxlsx)
require(here)
require(reshape)
# phylogeny

library(phytools)
library(geiger)
require(picante)

# trait simulation
require(mvMORPH)
require(parallel)

# trait imputation
require(missForest)

# functional diversity
require (FD)

#  fishtree 
require("fishtree")

# plot
require(ggplot2)
require(gridExtra)
require(scales)
require("tidybayes")

# analyses
require(brms)
require(emmeans)
require(cmdstanr) # install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# table
library("sjPlot")
library("sjmisc")
library("sjlabelled")

# maps & space
require(rnaturalearth)
require(rnaturalearthdata)
require(ggplot2)
require(gridExtra)
require(ggrepel)
require (scatterpie)
require(sf)
require(terra)

# taxonomic validation (Fish)
require(worrms)
