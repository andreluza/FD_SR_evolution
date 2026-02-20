
# ---------------------------------------------

# Run Models : Relationship between functional richness and evenness, considering also traits produced by the three models of evolution

# ----------------------------------------------
# Rodents
rm(list=ls())
source("R/packages.R")

# -----------------------------------------------------
require(here)
load (here("Output", "simulated_FD_BM_rodents.RData"))
load (here("Output", "simulated_FD_EB_rodents.RData"))
load (here("Output", "simulated_FD_OU_rodents.RData"))
load (here("Output", "empirical_FD_rodents.RData"))
#load (here("Output", "simulated_FD_random_rodents.RData"))

# load R data
load(here ( "Processed_data","image_rodents.RData"))
load(here ( "Processed_data","Spatial_covs_data.RData"))

# Create folder to store figures
dir.create (here("Output", "Figures"))

# density plot of macroevolutionary model parameters
# params
# BM

df_params<- lapply (seq (1,5), function (k){

  df_sigma_BM <- data.frame (Estimates = unlist(lapply (simul_param_BM, function (i) i[[k]]$opt$sigsq)),
                             Parameter = "sigma",
                             model = "BM",
                             Trait = names(simul_param_BM[[k]]))
  df_sigma_EB <- data.frame (Estimates = unlist(lapply (simul_param_EB, function (i) i[[k]]$opt$sigsq)),
                             Parameter = "sigma",
                             model = "EB",
                             Trait = names(simul_param_BM[[k]]))
  df_beta_EB <- data.frame (Estimates = unlist(lapply (simul_param_EB, function (i) i[[k]]$opt$a)),
                            Parameter = "beta", 
                            model = "EB",
                            Trait = names(simul_param_BM[[k]]))
  df_sigma_OU <- data.frame (Estimates = unlist(lapply (simul_param_OU, function (i) i[[k]]$opt$sigsq)),
                             Parameter = "sigma",
                             model = "OU",
                             Trait = names(simul_param_BM[[k]]))
  df_alpha_OU <- data.frame (Estimates = unlist(lapply (simul_param_OU, function (i) i[[k]]$opt$alpha)),
                             Parameter = "alpha",
                             model = "OU",
                             Trait = names(simul_param_BM[[k]]))
  # rbind
  df_density <- rbind (df_sigma_BM,
                       df_sigma_EB,
                       df_beta_EB,
                       df_sigma_OU,
                       df_alpha_OU)
})
# melt 
df_params<-do.call(rbind,df_params)
require(ggplot2)
# plot
# density plot
fig_params_rodents <- ggplot(df_params, 
                        aes(x=Estimates,
                            group=model,
                            color=model,
                            fill=model)) +
  geom_density(size=1,alpha=0.5)+
  #geom_histogram()+
  scale_fill_viridis_d(option = "magma",end=0.8,name="Model")+
  scale_colour_viridis_d(option = "magma",end=0.8,name="Model")+
  theme_classic()  + 
  facet_wrap(Trait~Parameter,scales = "free", labeller = label_parsed,ncol=3)+
  theme (legend.position = "right",
         axis.text.x = element_text(size=7))

fig_params_rodents
ggsave(here("Output", "Figures", "parameter_estimates_rodents.png"))


# ----------------------------------------------------------------

# empirical results
empirical_results <- data.frame (SR= apply(sapply(empirical_FD,"[[","nbsp"),1,mean),
                                 FRic= apply(sapply(empirical_FD,"[[","FRic"),1,mean),
                                 FEve=apply(sapply(empirical_FD,"[[","FEve"),1,mean),
                                 Dataset= "Empirical",
                                 Effort =  (spatial_data$Sampling_effort))

# random traits
#simulated_results_random <- data.frame (SR= apply(sapply(simulated_FD_random,"[[","nbsp"),1,mean),
#                                    FRic= apply(sapply(simulated_FD_random,"[[","FRic"),1,mean),
#                                    FEve=apply(sapply(simulated_FD_random,"[[","FEve"),1,mean),
#                                    Dataset= "SimulatedRandom",
#                                    Effort =  (spatial_data$Sampling_effort))

# average of simulated values (brownian motion)
simulated_results_BM <- data.frame (SR= apply(sapply(simulated_FD,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD,"[[","FEve"),1,mean),
                                    Dataset= "SimulatedBM",
                                    Effort =  (spatial_data$Sampling_effort))


# average of simulated values by EB
simulated_results_EB <- data.frame (SR= apply(sapply(simulated_FD_EB,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD_EB,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD_EB,"[[","FEve"),1,mean),
                                    Dataset = "SimulatedEB",
                                    Effort =  (spatial_data$Sampling_effort))

# average of simulated values by OU
simulated_results_OU <- data.frame (SR= apply(sapply(simulated_FD_OU,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD_OU,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD_OU,"[[","FEve"),1,mean),
                                    Dataset = "SimulatedOU",
                                    Effort =  (spatial_data$Sampling_effort))

# bind them
df_analyzes <- rbind(empirical_results,
                     #simulated_results_random,
                     simulated_results_BM,
                     simulated_results_EB,
                     simulated_results_OU)

# st transfo
spatial_data<-(st_transform(spatial_data,crs="EPSG:4326"))

# add covariates
df_analyzes <- cbind(df_analyzes,
      elevation = scale(spatial_data$elevation),
      slope = scale(spatial_data$slope),
      forest = scale(spatial_data$forest),
      latitude = scale(st_coordinates(spatial_data)[,"Y"]),
      #lat = st_coordinates(spatial_data)[,"Y"],
      region = ifelse (st_coordinates(spatial_data)[,"Y"] < -20, "south","north")
)

# RM sites with elevation == NA (Cerrado sites)
df_analyzes <- df_analyzes[!is.na(df_analyzes$elevation),]

# Number of sites used in the analyzes
nrow(df_analyzes[df_analyzes$Dataset == "Empirical",])

# Number of species used in the analyzes
nrow(empirical_FD[[1]]$x.axes)

# log responses and predictors - SR and effort
df_analyzes$logFRic <- log(df_analyzes$FRic)
df_analyzes$logFEve <- log(df_analyzes$FEve)
df_analyzes$logSR <- log(df_analyzes$SR)
df_analyzes$logEffort <- log(df_analyzes$Effort)

##---------------------------------------------------------
# analyses
# MCMC settings
nc<-3
ni<-10000
nb<-5000
nt<-10

# prior sigma
priors<-c(set_prior("normal(0,5)",class = "sigma"))

# run model (ancova)
model.ancova.FRic <-  brm (logFRic ~ (logSR*Dataset)+logEffort,
                           sigma ~ logSR,
                          data=df_analyzes,
                          prior=priors,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)

# add WAIC
model.ancova.FRic <- add_criterion(model.ancova.FRic, "loo", moment_match=T)

# run model (ancova)
model.ancova.FRic_sq <- brm (brmsformula(log(FRic) ~ ((logSR + I(logSR^2))*Dataset)+logEffort,
                          sigma ~ logSR),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic_sq)
tab_model(model.ancova.FRic_sq)

# add WAIC
model.ancova.FRic_sq <- add_criterion(model.ancova.FRic_sq, "loo", moment_match=T)

# compare
loo_compare(model.ancova.FRic, 
            model.ancova.FRic_sq)


# compare slopes (in main Ancova)
(m.lst.FRic <- emmeans::emtrends (model.ancova.FRic, "Dataset", var="logSR"))

# Analyse residuals in function of covariates
df_analyzes$Residuals_FRic <- residuals(model.ancova.FRic)[,"Estimate"]

# run one LM per data set
mod_dat <- lapply (unique(df_analyzes$Dataset)[1:2], function (i) {
  
  # fit the model
  dat <- df_analyzes [which(df_analyzes$Dataset == i),]
  # do LM
  mod_dat <- brm (formula = Residuals_FRic ~ region + elevation + slope + forest + latitude, 
                 data = dat,
                  family = gaussian (link="identity"),
                  chains=nc,
                  iter = ni,
                  warmup = nb,
                  thin=nt
                  )
  mod_dat
  
  }
)
tab_model(mod_dat[[1]])
tab_model(mod_dat[[2]])

# Functional  Evenness models ----------------------------------
#terms(model.ancova.FRic)
#priors<-c(set_prior("normal(0,5)",class = "sigma"))

# run model (ancova)
model.ancova.FEve <- brm (brmsformula(logFEve ~ (logSR*Dataset)+logEffort,
                          sigma ~ logSR),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)

# add WAIC
model.ancova.FEve <- add_criterion(model.ancova.FEve, "loo", moment_match=T)


#terms(model.ancova.FRic)
model.ancova.FEve_sq <- brm (brmsformula (log(FEve) ~  ((logSR + I(logSR^2))*Dataset)+logEffort,
                             sigma ~ logSR), 
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve_sq)
tab_model(model.ancova.FEve_sq)

# add WAIC
model.ancova.FEve_sq <- add_criterion(model.ancova.FEve_sq, "loo", moment_match=T)


# compare
loo_compare(model.ancova.FEve, 
            model.ancova.FEve_sq)


# compare slopes (in main Ancova)
(m.lst.FEve <- emmeans::emtrends (model.ancova.FEve_sq, "Dataset", var="logSR"))

# Analyse residuals in function of covariates
df_analyzes$Residuals_FEve <- residuals(model.ancova.FEve_sq)[,"Estimate"]

# run one LM per data set
mod_dat_FEve <- lapply (unique(df_analyzes$Dataset)[1:2], function (i) {
  
  # fit the model
  dat <- df_analyzes [which(df_analyzes$Dataset == i),]
  # do LM
  mod_dat <- brm ( formula = Residuals_FEve ~ region + elevation + slope + forest  + latitude, 
                 data = dat,
                  prior= priors,
                  family = gaussian (link="identity"),
                  chains=nc,
                  iter = ni,
                  warmup = nb,
                  thin=nt
                  )
  mod_dat
  
  
  }
)
tab_model(mod_dat_FEve[[1]])
tab_model(mod_dat_FEve[[2]])

# save results
save (model.ancova.FRic,
      model.ancova.FRic_sq,
      m.lst.FRic,
      model.ancova.FEve,
      model.ancova.FEve_sq,
      m.lst.FEve,
      mod_dat,
      mod_dat_FEve,
      file=here("Output", "GLM_test_rodents.RData"))
# end 