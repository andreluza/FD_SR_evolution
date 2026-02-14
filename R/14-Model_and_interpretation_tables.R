
# ---------------------------------------------

# Run Models : Relationship between functional richness and evenness, considering also traits produced by the three models of evolution
# Haemulidae

# ----------------------------------------------
rm(list=ls())
source("R/packages.R")

# -----------------------------------------------------
require(here)
##---------------------------------------------------------

# Load processed results
# environments
env_haemulidae <- new.env()
env_labridae <- new.env()
env_rodents <- new.env()

# load
load (file=here("Output", "GLM_test_haemulidae.RData"),envir = env_haemulidae)
load (file=here("Output", "GLM_test_labridae.RData"),envir = env_labridae)
load (file=here("Output", "GLM_test_rodents.RData"),envir = env_rodents)

# compare models
# Haemulidae
# FRic
require(kableExtra)
kable(
  rbind (
    loo_compare(env_haemulidae$model.ancova.FRic, 
                env_haemulidae$model.ancova.FRic_sq),
    # FEve
    loo_compare(env_haemulidae$model.ancova.FEve, 
                env_haemulidae$model.ancova.FEve_sq),
    # Labridae
    # FRic
    loo_compare(env_labridae$model.ancova.FRic, 
                env_labridae$model.ancova.FRic_sq),
    # FEve
    loo_compare(env_labridae$model.ancova.FEve, 
                env_labridae$model.ancova.FEve_sq),
    
    # Rodents
    # FRic
    loo_compare(env_rodents$model.ancova.FRic, 
                env_rodents$model.ancova.FRic_sq),
    # FEve
    loo_compare(env_rodents$model.ancova.FEve, 
                env_rodents$model.ancova.FEve_sq)
  ),format = "pipe"
  )

# FRIC ----------------------------------------------

# summary of results
# FRic
bayes_R2(env_rodents$model.ancova.FRic_sq,summary = T)
bayes_R2(env_labridae$model.ancova.FRic_sq,summary = T)
bayes_R2(env_haemulidae$model.ancova.FRic_sq,summary = T)

# FEve
bayes_R2(env_rodents$model.ancova.FEve_sq,summary = T)
bayes_R2(env_labridae$model.ancova.FEve,summary = T)
bayes_R2(env_haemulidae$model.ancova.FEve,summary = T)

# plots of coefficients
require(dplyr)
require(ggridges)
require(tidyr)

# Rodents
dat_rodents_FRic <- fixef(env_rodents$model.ancova.FRic_sq,summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "sigma_Intercept",
                                                        "logSR",
                                                        "IlogSRE2",
                                                        "DatasetSimulatedBM",
                                                        "DatasetSimulatedEB",
                                                        "DatasetSimulatedOU",
                                                        "DatasetSimulatedRandom",
                                                        "logEffort",
                                                        "logSR.DatasetSimulatedBM",
                                                        "logSR.DatasetSimulatedEB",
                                                        "logSR.DatasetSimulatedOU",
                                                        "logSR.DatasetSimulatedRandom",
                                                         "IlogSRE2.DatasetSimulatedBM",
                                                        "IlogSRE2.DatasetSimulatedEB",
                                                        "IlogSRE2.DatasetSimulatedOU",
                                                        "IlogSRE2.DatasetSimulatedRandom",
                                                        "sigma_logSR"
                                                        ))) %>%
  mutate ("Metric" = "FRic",
          "Taxon" = "Rodents")
  
# Wrasses

# Lbridae
dat_Labridae_FRic <- fixef(env_labridae$model.ancova.FRic_sq,summary = F) %>%
   data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "sigma_Intercept",
                                                        "logSR",
                                                        "IlogSRE2",
                                                        "DatasetSimulatedBM",
                                                        "DatasetSimulatedEB",
                                                        "DatasetSimulatedOU",
                                                        "DatasetSimulatedRandom",
                                                        "logEffort",
                                                        "logSR.DatasetSimulatedBM",
                                                        "logSR.DatasetSimulatedEB",
                                                        "logSR.DatasetSimulatedOU",
                                                        "logSR.DatasetSimulatedRandom",
                                                        "IlogSRE2.DatasetSimulatedBM",
                                                        "IlogSRE2.DatasetSimulatedEB",
                                                        "IlogSRE2.DatasetSimulatedOU",
                                                        "IlogSRE2.DatasetSimulatedRandom",
                                                        "sigma_logSR"
                                                        ))) %>%
  mutate ("Metric" = "FRic",
          "Taxon" = "Wrasses")


# Grunts
# Haemulidae
dat_Haemulidae_FRic <- fixef(env_haemulidae$model.ancova.FRic_sq,summary = F) %>%
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "sigma_Intercept",
                                                        "logSR",
                                                        "IlogSRE2",
                                                        "DatasetSimulatedBM",
                                                        "DatasetSimulatedEB",
                                                        "DatasetSimulatedOU",
                                                        "DatasetSimulatedRandom",
                                                        "logEffort",
                                                        "logSR.DatasetSimulatedBM",
                                                        "logSR.DatasetSimulatedEB",
                                                        "logSR.DatasetSimulatedOU",
                                                        "logSR.DatasetSimulatedRandom",
                                                        "IlogSRE2.DatasetSimulatedBM",
                                                        "IlogSRE2.DatasetSimulatedEB",
                                                        "IlogSRE2.DatasetSimulatedOU",
                                                        "IlogSRE2.DatasetSimulatedRandom",
                                                        "sigma_logSR"
                                                        ))) %>%
  mutate ("Metric" = "FRic",
          "Taxon" = "Grunts")


# bind all FRic data
dat_FRic <- rbind (dat_rodents_FRic,
       dat_Labridae_FRic,
       dat_Haemulidae_FRic)


# FEVE--------------------------
# Rodents
dat_rodents_FEve <- fixef(env_rodents$model.ancova.FEve_sq,summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "sigma_Intercept",
                                                        "logSR",
                                                        "IlogSRE2",
                                                        "DatasetSimulatedBM",
                                                        "DatasetSimulatedEB",
                                                        "DatasetSimulatedOU",
                                                        "DatasetSimulatedRandom",
                                                        "logEffort",
                                                        "logSR.DatasetSimulatedBM",
                                                        "logSR.DatasetSimulatedEB",
                                                        "logSR.DatasetSimulatedOU",
                                                        "logSR.DatasetSimulatedRandom",
                                                        "IlogSRE2.DatasetSimulatedBM",
                                                        "IlogSRE2.DatasetSimulatedEB",
                                                        "IlogSRE2.DatasetSimulatedOU",
                                                        "IlogSRE2.DatasetSimulatedRandom",
                                                        "sigma_logSR"
                                                        ))) %>%
  mutate ("Metric" = "FEve",
          "Taxon" = "Rodents")
  
# Wrasses

# Lbridae
dat_Labridae_FEve <- fixef(env_labridae$model.ancova.FEve,summary = F) %>%
   data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "sigma_Intercept",
                                                        "logSR",
                                                        "DatasetSimulatedBM",
                                                        "DatasetSimulatedEB",
                                                        "DatasetSimulatedOU",
                                                        "DatasetSimulatedRandom",
                                                        "logEffort",
                                                        "logSR.DatasetSimulatedBM",
                                                        "logSR.DatasetSimulatedEB",
                                                        "logSR.DatasetSimulatedOU",
                                                        "logSR.DatasetSimulatedRandom",
                                                        "sigma_logSR"
                                                        ))) %>%
  mutate ("Metric" = "FEve",
          "Taxon" = "Wrasses")


# Grunts
# Haemulidae
dat_Haemulidae_FEve <- fixef(env_haemulidae$model.ancova.FEve,summary = F) %>%
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "sigma_Intercept",
                                                        "logSR",
                                                        "DatasetSimulatedBM",
                                                        "DatasetSimulatedEB",
                                                        "DatasetSimulatedOU",
                                                        "DatasetSimulatedRandom",
                                                        "logEffort",
                                                        "logSR.DatasetSimulatedBM",
                                                        "logSR.DatasetSimulatedEB",
                                                        "logSR.DatasetSimulatedOU",
                                                        "logSR.DatasetSimulatedRandom",
                                                        "sigma_logSR"
                                                        ))) %>%
  mutate ("Metric" = "FEve",
          "Taxon" = "Grunts")


# bind all FRic data
dat_FEve <- rbind (dat_rodents_FEve,
       dat_Labridae_FEve,
       dat_Haemulidae_FEve)



# Bind data
dat_FD <- rbind(dat_FRic,
                dat_FEve)
dat_FD$CoefficientOld<-dat_FD$Coefficient
  
# organize factors
dat_FD$Coefficient <- dplyr::recode_factor(dat_FD$Coefficient,  
                                        "sigma_Intercept" = "eta[0]", 
                                        "sigma_logSR" = "eta[1]", 
                                        "Intercept" = "gamma[0]", 
                                        "DatasetSimulatedBM" = "gamma[0-BM]", 
                                        "DatasetSimulatedEB" ="gamma[0-EB]", 
                                        "DatasetSimulatedOU" = "gamma[0-OU]", 
                                        "DatasetSimulatedRandom" = "gamma[0-rdm]",
                                        "logSR" = "gamma[1]", 
                                        "IlogSRE2" = "gamma[2]", 
                                        "logEffort" = "gamma[3]", 
                                        "logSR.DatasetSimulatedBM" = "gamma[1-BM]",
                                        "logSR.DatasetSimulatedEB" = "gamma[1-EB]", 
                                        "logSR.DatasetSimulatedOU" = "gamma[1-OU]", 
                                        "logSR.DatasetSimulatedRandom" = "gamma[1-rdm]",
                                        "IlogSRE2.DatasetSimulatedBM" = "gamma[2-BM]",
                                        "IlogSRE2.DatasetSimulatedEB" = "gamma[2-EB]",
                                        "IlogSRE2.DatasetSimulatedOU" = "gamma[2-OU]",
                                        "IlogSRE2.DatasetSimulatedRandom" = "gamma[2-rdm]"
                                        
                                        )

# Taxa
dat_FD$Taxon <- factor(dat_FD$Taxon,levels=c("Rodents", "Wrasses", "Grunts"))
# Metric
dat_FD$Metric <- factor(dat_FD$Metric,levels=c("FRic", "FEve"))
# Model colors as per FIg. 2
dat_FD$col_ridges <-  dplyr::recode_factor(dat_FD$CoefficientOld,  
                                        "sigma_Intercept" = "#000000", 
                                        "sigma_logSR" = "#000000", 
                                        "Intercept" = "#000000", 
                                        "DatasetSimulatedBM" = "#0F00FF", 
                                        "DatasetSimulatedEB" ="#D98C00", 
                                        "DatasetSimulatedOU" = "#A4EBF3", 
                                        "DatasetSimulatedRandom" = "red4",
                                        "logSR" = "#000000", 
                                        "IlogSRE2" = "#000000", 
                                        "logEffort" = "#000000", 
                                        "logSR.DatasetSimulatedBM" = "#0F00FF",
                                        "logSR.DatasetSimulatedEB" = "#D98C00", 
                                        "logSR.DatasetSimulatedOU" = "#A4EBF3", 
                                        "logSR.DatasetSimulatedRandom" = "red4",
                                        "IlogSRE2.DatasetSimulatedBM" = "#0F00FF",
                                        "IlogSRE2.DatasetSimulatedEB" = "#D98C00",
                                        "IlogSRE2.DatasetSimulatedOU" = "#A4EBF3",
                                        "IlogSRE2.DatasetSimulatedRandom" = "red4"
                                        
                                        )



# plot
require(ggplot2)
dat_FD$col_lab <- NA
dat_FD$col_lab [grep ("EM", dat_FD$Coefficient)]<-1
dat_FD$col_lab [is.na(dat_FD$col_lab)]<-2


dat_FD %>%
  ggplot(aes(x=Estimate,y=Coefficient,fill=(col_ridges),col=(col_ridges))) +
  #geom_rect(aes(xmin=-6,xmax=3,ymin=1.5,ymax=5.5),fill="gray90",alpha=0.1) +
  #geom_rect(aes(xmin=-6,xmax=3,ymin=7.5,ymax=11.5),fill="gray90",alpha=0.1) +
  #geom_pointrange(aes(x=Estimate,y=Coefficient,xmin=Q2.5,xmax=Q97.5))+
  ggplot2::scale_fill_manual(values=levels(dat_FD$col_ridges))+
  ggplot2::scale_colour_manual(values=levels(dat_FD$col_ridges))+
  scale_y_discrete(labels=parse_format())+
  facet_grid(Taxon~Metric,scales="free_x")+
  geom_density_ridges(quantile_lines=TRUE,
                      quantile_fun=function(x,...)mean(x),
                      rel_min_height = 0.005,
                      alpha=0.5)+
  theme_bw(base_size = 14) +
  geom_vline(aes(xintercept=0),linetype=3)+
  theme(legend.position = "none")

ggsave(here("Output","Figures","Fig4.png"), width=7,height=8)


# -------------------------------------------------------

# Residuals
res_FRic <- rbind (
  # Rodents empirical
  fixef(env_rodents$mod_dat[[1]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionsouth",
                                                        "elevation",
                                                        "slope",
                                                        "forest",
                                                        "latitude"))) %>%
  mutate ("Metric" = "FRic",
          "Data" = "EM",
          "Taxon" = "Rodents")
  ,
  # Rodents BM
  fixef(env_rodents$mod_dat[[2]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionsouth",
                                                        "elevation",
                                                        "slope",
                                                        "forest",
                                                        "latitude"))) %>%
  mutate ("Metric" = "FRic",
          "Data" = "BM",
          "Taxon" = "Rodents")
  ,
  
  # Labridae empirical
  fixef(env_labridae$mod_dat[[1]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FRic",
          "Data" = "EM",
          "Taxon" = "Wrasses")
  
  ,
  # Labridae BM
  fixef(env_labridae$mod_dat[[2]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FRic",
          "Data" = "BM",
          "Taxon" = "Wrasses")
  
  ,

# Haemulidae empirical
  fixef(env_haemulidae$mod_dat[[1]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FRic",
          "Data" = "EM",
          "Taxon" = "Grunts")
  
  ,
  # Labridae BM
  fixef(env_haemulidae$mod_dat[[2]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FRic",
          "Data" = "BM",
          "Taxon" = "Grunts")

)

# organize factors
res_FRic$Coefficient <- dplyr::recode_factor(res_FRic$Coefficient,  
                                        "Intercept" = "Intercept", 
                                        "regionsouth" = "S region",
                                        "elevation" = "Elevation", 
                                        "slope" = "Slope", 
                                        "forest" = "Forest", 
                                        "latitude" = "Latitude",
                                        "regionoc_isl" = "Oceanic Islands",
                                        "regionse_reefs" = "S-SE region",
                                        "Reef.Area" = "Reef Area", 
                                        "BO_damean_lonlat" = "Turbidity", 
                                        "Present.Surface.Salinity.Mean" = "Salinity",
                                        "Present.Surface.Temperature.Mean" = "Temperature", 
                                        "decimalLatitude" = "Latitude")
# Taxa
res_FRic$Taxon <- factor(res_FRic$Taxon,levels=c("Rodents", "Wrasses", "Grunts"))
# data
res_FRic$Data <- factor(res_FRic$Data,levels=c("EM", "BM"))
res_FRic$facet <- "FRic"

# plot
res_FRic %>%
  ggplot(aes(x=Estimate,y=Coefficient)) +
  #geom_rect(aes(xmin=-6,xmax=3,ymin=1.5,ymax=5.5),fill="gray90",alpha=0.1) +
  #geom_rect(aes(xmin=-6,xmax=3,ymin=7.5,ymax=11.5),fill="gray90",alpha=0.1) +
  #geom_pointrange(aes(x=Estimate,y=Coefficient,xmin=Q2.5,xmax=Q97.5))+
  scale_fill_viridis_d(option="magma",begin=0.2,end=0.8)+
  #scale_y_discrete(labels=parse_format())+
  facet_grid(Taxon~facet+Data,scales="free_y")+
  #geom_density_ridges(quantile_lines=TRUE,
  #                    quantile_fun=function(x,...)mean(x),
  #                    rel_min_height = 0.005,
  #                    alpha=0.5)+
  stat_density_ridges(quantile_lines = TRUE,quantiles = c(0.025,0.5, 0.975),
                        rel_min_height = 0.005)+

  theme_bw(base_size = 14) +
  geom_vline(aes(xintercept=0),linetype=3)+
  theme(legend.position = "none")

ggsave(here("Output","Figures","FigS5.png"), width=7,height=7)


# FEve------------------------------------------------------------------

# Residuals
res_FEve <- rbind (
  # Rodents empirical
  fixef(env_rodents$mod_dat_FEve[[1]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                         "regionsouth",
                                                        "elevation",
                                                        "slope",
                                                        "forest",
                                                        "latitude"))) %>%
  mutate ("Metric" = "FEve",
          "Data" = "EM",
          "Taxon" = "Rodents")
  ,
  # Rodents BM
  fixef(env_rodents$mod_dat_FEve[[2]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                         "regionsouth",
                                                        "elevation",
                                                        "slope",
                                                        "forest",
                                                        "latitude"))) %>%
  mutate ("Metric" = "FEve",
          "Data" = "BM",
          "Taxon" = "Rodents")
  ,
  
  # Labridae empirical
  fixef(env_labridae$mod_dat_FEve[[1]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FEve",
          "Data" = "EM",
          "Taxon" = "Wrasses")
  
  ,
  # Labridae BM
  fixef(env_labridae$mod_dat_FEve[[2]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FEve",
          "Data" = "BM",
          "Taxon" = "Wrasses")
  
  ,

# Haemulidae empirical
  fixef(env_haemulidae$mod_dat_FEve[[1]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FEve",
          "Data" = "EM",
          "Taxon" = "Grunts")
  
  ,
  # Labridae BM
  fixef(env_haemulidae$mod_dat_FEve[[2]],summary = F) %>%
  
  data.frame()%>%
  melt() %>%
  mutate(Coefficient = variable,
         Estimate = value) %>%
  
  mutate (Coefficient = factor (Coefficient, levels = c("Intercept",
                                                        "regionoc_isl",
                                                        "regionse_reefs",
                                                        "Reef.Area",
                                                        "BO_damean_lonlat",
                                                        "Present.Surface.Salinity.Mean",
                                                        "Present.Surface.Temperature.Mean",
                                                        "decimalLatitude")))  %>%
  mutate ("Metric" = "FEve",
          "Data" = "BM",
          "Taxon" = "Grunts")

)

# organize factors
res_FEve$Coefficient <- dplyr::recode_factor(res_FEve$Coefficient,  
                                       "Intercept" = "Intercept", 
                                        "regionsouth" = "S region",
                                        "elevation" = "Elevation", 
                                        "slope" = "Slope", 
                                        "forest" = "Forest", 
                                        "latitude" = "Latitude",
                                        "regionoc_isl" = "Oceanic Islands",
                                        "regionse_reefs" = "S-SE region",
                                        "Reef.Area" = "Reef Area", 
                                        "BO_damean_lonlat" = "Turbidity", 
                                        "Present.Surface.Salinity.Mean" = "Salinity",
                                        "Present.Surface.Temperature.Mean" = "Temperature", 
                                        "decimalLatitude" = "Latitude")
# Taxa
res_FEve$Taxon <- factor(res_FEve$Taxon,levels=c("Rodents", "Wrasses", "Grunts"))
# data
res_FEve$Data <- factor(res_FEve$Data,levels=c("EM", "BM"))
res_FEve$facet <- "FEve"

# plot
res_FEve %>%
  ggplot(aes(x=Estimate,y=Coefficient)) +
  #geom_rect(aes(xmin=-6,xmax=3,ymin=1.5,ymax=5.5),fill="gray90",alpha=0.1) +
  #geom_rect(aes(xmin=-6,xmax=3,ymin=7.5,ymax=11.5),fill="gray90",alpha=0.1) +
  #geom_pointrange(aes(x=Estimate,y=Coefficient,xmin=Q2.5,xmax=Q97.5))+
  scale_fill_viridis_d(option="magma",begin=0.2,end=0.8)+
  #scale_y_discrete(labels=parse_format())+
  facet_grid(Taxon~facet+Data,scales="free_y")+

    stat_density_ridges(quantile_lines = TRUE,quantiles = c(0.025,0.5, 0.975),
                        rel_min_height = 0.005)+
#geom_density_ridges(quantile_lines=TRUE,
#                      quantile_fun=function(x,...)mean(x),
#                      rel_min_height = 0.005,
#                      alpha=0.5)+
  theme_bw(base_size = 14) +
  geom_vline(aes(xintercept=0),linetype=3)+
  theme(legend.position = "none")

ggsave(here("Output","Figures","FigS6.png"), width=7,height=7)


# summary of results
# FRic - empirical
# rodents
bayes_R2(env_rodents$mod_dat[[1]],summary = T)
bayes_R2(env_rodents$mod_dat[[2]],summary = T)# BM

# wrasses
bayes_R2(env_labridae$mod_dat[[1]],summary = T)
bayes_R2(env_labridae$mod_dat[[2]],summary = T)# BM

# grunts
bayes_R2(env_haemulidae$mod_dat[[1]],summary = T)
bayes_R2(env_haemulidae$mod_dat[[2]],summary = T)# BM


# FEve - empirical
# rodents
bayes_R2(env_rodents$mod_dat_FEve[[1]],summary = T)
bayes_R2(env_rodents$mod_dat_FEve[[2]],summary = T) #BM

# wrasses
bayes_R2(env_labridae$mod_dat_FEve[[1]],summary = T)
bayes_R2(env_labridae$mod_dat_FEve[[2]],summary = T)#BM

# grunts
bayes_R2(env_haemulidae$mod_dat_FEve[[1]],summary = T)
bayes_R2(env_haemulidae$mod_dat_FEve[[2]],summary = T)#BM

#end
