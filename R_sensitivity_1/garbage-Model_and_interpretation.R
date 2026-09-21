
# ---------------------------------------------

# Run Models : Relationship between functional richness and evenness, considering also traits produced by the three models of evolution

# ----------------------------------------------
# Rodents

source("R/packages.R")

# -----------------------------------------------------
require(here)
load (here("Output", "simulated_FD_BM_rodents.RData"))
load (here("Output", "simulated_FD_EB_rodents.RData"))
load (here("Output", "simulated_FD_OU_rodents.RData"))
load (here("Output", "empirical_FD_rodents.RData"))

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
                            group=Trait,
                            color=Trait,
                            fill=Parameter)) +
  geom_density(size=1,alpha=0.5)+
  #geom_histogram()+
  scale_fill_viridis_d(option = "magma")+
  scale_colour_viridis_d(option = "magma")+
  theme_classic()  + 
  facet_wrap(Parameter~model,scales = "free")+
  theme (legend.position = c(0.8,0.25),
         axis.text.x = element_text(size=7))
fig_params_rodents

# ----------------------------------------------------------------
# empirical results
empirical_results <- data.frame (SR= apply(sapply(empirical_FD,"[[","nbsp"),1,mean),
                                 FRic= apply(sapply(empirical_FD,"[[","FRic"),1,mean),
                                 FEve=apply(sapply(empirical_FD,"[[","FEve"),1,mean),
                                 Dataset= "Empirical",
                                 Effort =  (spatial_data$Sampling_effort))

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
                     simulated_results_BM,
                     simulated_results_EB,
                     simulated_results_OU)

# add covariates
df_analyzes <- cbind(df_analyzes,
      elevation = scale(spatial_data$elevation),
      slope = scale(spatial_data$slope),
      forest = scale(spatial_data$forest))

# RM sites with elevation == NA (Cerrado sites)
df_analyzes <- df_analyzes[!is.na(df_analyzes$elevation),]

# prepare data for fig 4
# neutral SES Rodents
SR_rod<-apply(sapply(empirical_FD,"[[","nbsp"),1,mean)
obsFEve_rod <- apply(sapply(empirical_FD,"[[","FEve"),1,mean) 
meanFEveSim_rod <- apply(sapply(simulated_FD,"[[","FEve"),1,mean)
sdFEveSim_rod <- apply(sapply(simulated_FD,"[[","FEve"),1,sd)

##---------------------------------------------------------
# analyses
# MCMC settings
nc<-3
ni<-10000
nb<-5000
nt<-10

# run model (ancova)
model.ancova.FRic <- brm (log(FRic) ~ (log(SR)*Dataset)+(elevation*Dataset)+(slope*Dataset)+(forest*Dataset)+log(Effort),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)


# run model 2 (ancova)
model.ancova.FRic_m2 <- brm (log(FRic) ~ (log(SR)*Dataset)+log(Effort),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic_m2)
tab_model(model.ancova.FRic_m2)

# add WAIC
model.ancova.FRic <- add_criterion(model.ancova.FRic, "loo", moment_match=T)
model.ancova.FRic_m2 <- add_criterion(model.ancova.FRic_m2, "loo", moment_match=T)

# Compare the models
# loo
loo(model.ancova.FRic,
    model.ancova.FRic_m2
    )

# plotting
p1<-plot(conditional_effects(model.ancova.FRic_m2,
                             method="posterior_epred",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = c(0.75,0.3)) ,
         points=T,
         point_args = list (width = 0.25,alpha=0.3)) [[1]] + 
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Richness (log(FRic))")

p1

# compare slopes (in main Ancova)
(m.lst.FRic <- emmeans::emtrends (model.ancova.FRic_m2, "Dataset", var="SR"))
# m.lst_tab.FRic <- summary(m.lst.FRic,point.est = mean)

# Analyse residuals in function of covariates
df_analyzes$Residuals_FRic <- residuals(model.ancova.FRic_m2)[,"Estimate"]

# run one LM per data set
i = "Empirical"
lapply (unique(df_analyzes$Dataset), function (i) {
  
  # fit the model
  dat <- df_analyzes [which(df_analyzes$Dataset == i),]
  # do LM
  mod_dat <- lm (data = dat, formula = Residuals_FRic ~  elevation + slope + forest)
  RsquareAdj(mod_dat)
  summary(mod_dat)
  
  }
)


# Functional  Evenness models ----------------------------------

#terms(model.ancova.FRic)
priors<-c(set_prior("normal(0,5)",class = "sigma")
          )

# run model (ancova)
model.ancova.FEve <- brm (log(FEve) ~ (log(SR)*Dataset)+(elevation*Dataset)+(slope*Dataset)+(forest*Dataset),
                          sigma ~ log(SR),
                          prior = priors,
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)

# run model (ancova)
model.ancova.FEve_m2 <- brm (log(FEve) ~ (log(SR)*Dataset)+log(Effort),
                             sigma ~ log(SR), 
                          data=df_analyzes,
                          prior = priors,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve_m2)
tab_model(model.ancova.FEve_m2)

# add WAIC
model.ancova.FEve <- add_criterion(model.ancova.FEve, "loo", moment_match=T)
model.ancova.FEve_m2 <- add_criterion(model.ancova.FEve_m2, "loo", moment_match=T)

# Compare the models
# loo
loo(model.ancova.FEve,
    model.ancova.FEve_m2
    )

# plotting
p2<-plot(conditional_effects(model.ancova.FEve,
                             method="posterior_epred",  # posterior_predict
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "none") ,
         points=T,
         point_args = list (width = 0.25,alpha=0.3)) [[1]] + 
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (log(FEve))")
p2

# Analyse residuals in function of covariates
df_analyzes$Residuals_FEve <- residuals(model.ancova.FEve)[,"Estimate"]

# run one LM per data set
lapply (unique(df_analyzes$Dataset), function (i) {
  
  # fit the model
  dat <- df_analyzes [which(df_analyzes$Dataset == i),]
  # do LM
  mod_dat <- lm (data = dat, formula = Residuals_FEve ~  elevation + slope + forest)
  RsquareAdj(mod_dat)
  summary(mod_dat)
  
  }
)


# organize plots
pdf(here("Output","Fig3_SM.pdf"), width=9,height=5)
  
  grid.arrange(p1,p2,nrow=1)

dev.off()

# 
cor(predict(model.ancova.FRic_m2)[,1],
    predict(model.ancova.FEve)[,1])

# compare slopes (in main Ancova)
(m.lst.FEve <- emmeans::emtrends (model.ancova.FEve, "Dataset", var="SR"))

# save results
save (model.ancova.FRic_m2,
      m.lst.FRic,
      model.ancova.FEve,
      m.lst.FEve,
      file=here("Output", "GLM_test_rodents.RData"))


# ---------------------------------------------------------------------
# fish - Labridae
load (here("Output", "simulated_FD_BM_labridae.RData"))
load (here("Output", "simulated_FD_EB_labridae.RData"))
load (here("Output", "simulated_FD_OU_labridae.RData"))
load (here("Output", "empirical_FD_labridae.RData"))

# load R data
load(here ( "Processed_data","image_labridae.RData"))
load(here ( "Processed_data","Spatial_covs_data_reefs.RData"))

# organize params
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

# plot
# density plot
fig_params_fish <- ggplot(df_params, 
                             aes(x=Estimates,
                                 group=Trait,
                                 color=Trait,
                                 fill=Parameter)) +
  geom_density(size=1,alpha=0.05)+
  scale_fill_viridis_d(option = "magma")+
  scale_colour_viridis_d(option = "magma")+
  theme_classic()  + 
  facet_wrap(Parameter~model,scales = "free")+
  theme (legend.position = c(0.81,0.23))
fig_params_fish

# =-----------------------------------------------------------------

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
simulated_FD_OU <- simulated_FD_OU [which(lapply (simulated_FD_OU,length) != 1)]
simulated_results_OU <- data.frame (SR= apply(sapply(simulated_FD_OU,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD_OU,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD_OU,"[[","FEve"),1,mean),
                                    Dataset = "SimulatedOU")

# bind them
df_analyzes <- rbind(empirical_results,
                     simulated_results_BM,
                     simulated_results_EB,
                     simulated_results_OU)

# add covariates
reef_covariates <- rbind(reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),])

# bind
df_analyzes <- cbind(df_analyzes,
      turbidity = scale(reef_covariates[,"BO_damean_lonlat"]),
      salinity = scale(reef_covariates[,"Present.Surface.Salinity.Mean"]),
      temperature = scale(reef_covariates[,"Present.Surface.Temperature.Mean"]),
      productivity = scale(reef_covariates[,"Present.Surface.Primary.productivity.Mean"]),
      area = scale(reef_covariates[,"Reef.Area"]))

# bind effort
df_analyzes <- cbind (df_analyzes,
                      effort_site[which(rowSums(subset_comm_data[[1]]) >0),])

# Remove NAs FRic
df_analyzes <- df_analyzes [!is.na(df_analyzes$FRic),] 

# correlations
cor(df_analyzes[,c(1:3,5:9,11)])

# prepare data for fig 4
# neutral SES fish
SR_fish<-apply(sapply(empirical_FD,"[[","nbsp"),1,mean)
obsFEve_fish <- apply(sapply(empirical_FD,"[[","FEve"),1,mean) 
meanFEveSim_fish <- apply(sapply(simulated_FD,"[[","FEve"),1,mean)
sdFEveSim_fish <- apply(sapply(simulated_FD,"[[","FEve"),1,sd)

##---------------------------------------------------------

# run model (ancova)
model.ancova.FRic <- brm (log(FRic) ~ (log(SR)*Dataset)+(Present.Surface.Temperature.Mean*Dataset)+(BO_damean_lonlat*Dataset)+(Reef.Area*Dataset)+log(effort),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)


# run model (ancova)
model.ancova.FRic_m2 <- brm (log(FRic) ~ (log(SR)*Dataset)+log(effort),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt
                          #control = list(max_treedepth = 15)
                          )

# summary of results
summary (model.ancova.FRic_m2)
tab_model(model.ancova.FRic_m2)

# add WAIC
model.ancova.FRic <- add_criterion(model.ancova.FRic, "loo", moment_match=T)
model.ancova.FRic_m2 <- add_criterion(model.ancova.FRic_m2, "loo", moment_match=T)

# Compare the models
# loo
loo(model.ancova.FRic,
    model.ancova.FRic_m2
    )

# plotting
p1<-plot(conditional_effects(model.ancova.FRic,
                             method="posterior_epred",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = c(0.75,0.3)) ,
         points=T) [[1]] + 
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Richness (FRic)")

p1

# compare slopes
m.lst.FRic <- emtrends (model.ancova.FRic, "Dataset", var="SR")
m.lst_tab.FRic <- summary(m.lst.FRic,point.est = mean)

# run model (ancova)
model.ancova.FEve <- brm (log(FEve) ~ (log(SR)*Dataset)+(Present.Surface.Temperature.Mean*Dataset)+(BO_damean_lonlat*Dataset)+(Reef.Area*Dataset)+log(effort),
                             sigma ~ log(SR), 
                          data=df_analyzes,
                          prior = priors,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)

# run model (ancova)
model.ancova.FEve_m2 <- brm (log(FEve) ~ (log(SR)*Dataset)+log(effort),
                             sigma ~ log(SR), 
                          data=df_analyzes,
                          prior = priors,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve_m2)
tab_model(model.ancova.FEve_m2)

# add WAIC
model.ancova.FEve <- add_criterion(model.ancova.FEve, "loo", moment_match=T)
model.ancova.FEve_m2 <- add_criterion(model.ancova.FEve_m2, "loo", moment_match=T)


# Compare the models
# loo
loo(model.ancova.FEve,
    model.ancova.FEve_m2
    )

# plotting
p2<-plot(conditional_effects(model.ancova.FEve_m2,
                             method="posterior_epred",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "none") ,
         points=T) [[1]] + 
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (FEve)")

p2

# organize plots
pdf(here("Output","Fig3_labridae.pdf"), width=9,height=5)
  grid.arrange(p1,p2,nrow=1)
dev.off()

# compare slopes
m.lst.FEve <- emtrends (model.ancova.FEve_m2, "Dataset", var="SR")
m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

# save results
save (model.ancova.FRic,
      m.lst.FRic,
      m.lst_tab.FRic,
      model.ancova.FEve_m2,
      m.lst.FEve,
      m.lst_tab.FEve,
      file=here("Output", "GLM_test_labridae.RData"))


# ----------------------------------------------------
# Haemulidae

load (here("Output", "simulated_FD_BM_haemulidae.RData"))
load (here("Output", "simulated_FD_EB_haemulidae.RData"))
load (here("Output", "simulated_FD_OU_haemulidae.RData"))
load (here("Output", "empirical_FD_haemulidae.RData"))

# load R data
load(here ( "Processed_data","image_haemulidae.RData"))
load(here ( "Processed_data","Spatial_covs_data_reefs.RData"))

# organize params
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

# plot
# density plot
fig_params_fish <- ggplot(df_params, 
                             aes(x=Estimates,
                                 group=Trait,
                                 color=Trait,
                                 fill=Parameter)) +
  geom_density(size=1,alpha=0.05)+
  scale_fill_viridis_d(option = "magma")+
  scale_colour_viridis_d(option = "magma")+
  theme_classic()  + 
  facet_wrap(Parameter~model,scales = "free")+
  theme (legend.position = c(0.81,0.23))
fig_params_fish

# =-----------------------------------------------------------------

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
simulated_FD_OU <- simulated_FD_OU [which(lapply (simulated_FD_OU,length) != 1)]
simulated_results_OU <- data.frame (SR= apply(sapply(simulated_FD_OU,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD_OU,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD_OU,"[[","FEve"),1,mean),
                                    Dataset = "SimulatedOU")

# bind them
df_analyzes <- rbind(empirical_results,
                     simulated_results_BM,
                     simulated_results_EB,
                     simulated_results_OU)

# add covariates
reef_covariates <- rbind(reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),])

# bind
df_analyzes <- cbind(df_analyzes,
      turbidity = scale(reef_covariates[,"BO_damean_lonlat"]),
      salinity = scale(reef_covariates[,"Present.Surface.Salinity.Mean"]),
      temperature = scale(reef_covariates[,"Present.Surface.Temperature.Mean"]),
      productivity = scale(reef_covariates[,"Present.Surface.Primary.productivity.Mean"]),
      area = scale(reef_covariates[,"Reef.Area"]))

# bind effort
df_analyzes <- cbind (df_analyzes,
                      effort_site[which(rowSums(subset_comm_data[[1]]) >0),])

# Remove NAs FRic
df_analyzes <- df_analyzes [!is.na(df_analyzes$FRic),] 

# correlations
cor(df_analyzes[,c(1:3,5:9,11)])

# prepare data for fig 4
# neutral SES fish
SR_fish<-apply(sapply(empirical_FD,"[[","nbsp"),1,mean)
obsFEve_fish <- apply(sapply(empirical_FD,"[[","FEve"),1,mean) 
meanFEveSim_fish <- apply(sapply(simulated_FD,"[[","FEve"),1,mean)
sdFEveSim_fish <- apply(sapply(simulated_FD,"[[","FEve"),1,sd)

##---------------------------------------------------------

# run model (ancova)
model.ancova.FRic <- brm (log(FRic) ~ (log(SR)*Dataset)+(Present.Surface.Temperature.Mean*Dataset)+(BO_damean_lonlat*Dataset)+(Reef.Area*Dataset)+log(effort),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)


# run model (ancova)
model.ancova.FRic_m2 <- brm (log(FRic) ~ (log(SR)*Dataset)+log(effort),
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt
                          #control = list(max_treedepth = 15)
                          )

# summary of results
summary (model.ancova.FRic_m2)
tab_model(model.ancova.FRic_m2)

# add WAIC
model.ancova.FRic <- add_criterion(model.ancova.FRic, "loo", moment_match=T)
model.ancova.FRic_m2 <- add_criterion(model.ancova.FRic_m2, "loo", moment_match=T)

# Compare the models
# loo
loo(model.ancova.FRic,
    model.ancova.FRic_m2
    )

# plotting
p1<-plot(conditional_effects(model.ancova.FRic,
                             method="posterior_epred",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = c(0.75,0.3)) ,
         points=T) [[1]] + 
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Richness (FRic)")

p1

# compare slopes
m.lst.FRic <- emtrends (model.ancova.FRic, "Dataset", var="SR")
m.lst_tab.FRic <- summary(m.lst.FRic,point.est = mean)

# run model (ancova)

#terms(model.ancova.FRic)
priors<-c(set_prior("normal(0,5)",class = "sigma")
          )
model.ancova.FEve <- brm (log(FEve) ~ (log(SR)*Dataset)+(Present.Surface.Temperature.Mean*Dataset)+(BO_damean_lonlat*Dataset)+(Reef.Area*Dataset)+log(effort),
                             sigma ~ log(SR), 
                          data=df_analyzes,
                          prior = priors,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)

# run model (ancova)
model.ancova.FEve_m2 <- brm (log(FEve) ~ (log(SR)*Dataset)+log(effort),
                             sigma ~ log(SR), 
                          data=df_analyzes,
                          prior = priors,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve_m2)
tab_model(model.ancova.FEve_m2)

# add WAIC
model.ancova.FEve <- add_criterion(model.ancova.FEve, "loo", moment_match=T)
model.ancova.FEve_m2 <- add_criterion(model.ancova.FEve_m2, "loo", moment_match=T)


# Compare the models
# loo
loo(model.ancova.FEve,
    model.ancova.FEve_m2
    )

# plotting
p2<-plot(conditional_effects(model.ancova.FEve_m2,
                             method="posterior_epred",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "none") ,
         points=T) [[1]] + 
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (FEve)")

p2

# organize plots
pdf(here("Output","Fig3_haemulidae.pdf"), width=9,height=5)
  grid.arrange(p1,p2,nrow=1)
dev.off()

# compare slopes
m.lst.FEve <- emtrends (model.ancova.FEve_m2, "Dataset", var="SR")
m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

# save results
save (model.ancova.FRic,
      m.lst.FRic,
      m.lst_tab.FRic,
      model.ancova.FEve_m2,
      m.lst.FEve,
      m.lst_tab.FEve,
      file=here("Output", "GLM_test_haemulidae.RData"))








# -------------------------------------------------

# Fig 4
# neutral SES 
neutral_ses_feve_rodent<-(obsFEve_rod-meanFEveSim_rod)/sdFEveSim_rod
neutral_ses_feve_fish<-(obsFEve_fish-meanFEveSim_fish)/sdFEveSim_fish

# dataframe for plotting
df_fig4<-rbind (
  data.frame (SES = neutral_ses_feve_rodent,
            SR = SR_rod,
            Organism = "Rodents"),
  data.frame (SES = neutral_ses_feve_fish,
              SR = SR_fish,
              Organism = "Reef fish"))
df_fig4$Significance <- ifelse (df_fig4$SES >= 1.96,
                       "Higher",
                       ifelse (df_fig4$SES <= -1.96,
                               "Lower",
                               "Equal"))
df_fig4$Significance<-factor(df_fig4$Significance,
                             levels = c("Higher","Equal", "Lower"))
# plot
require(ggplot2)
ggplot (data = df_fig4, aes (x=SR, y=SES))+
  geom_smooth(method="glm", formula = y~poly(x,2),col="black") + 
  theme_classic() + 
  facet_wrap(~Organism, scale="free") + 
  geom_point(data = df_fig4, 
             aes (x=SR, y=SES,col=Significance),
             size=2)  +
  scale_color_manual(values=c("#E69F00","#999999",  "#56B4E9"),
                     name = "Deviations from\na BM model",
                     labels = c("Positive",
                                "None",
                                "Negative"))+
  theme(legend.position = c(0.90,0.15),
        
        strip.text = element_text(size=15),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12)) + 
  ylab ("Standardized Effect Size") + 
  xlab ("Species richness gradient")
  
  

