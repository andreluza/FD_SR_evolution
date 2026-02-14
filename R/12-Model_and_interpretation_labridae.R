
# ---------------------------------------------

# Run Models : Relationship between functional richness and evenness, considering also traits produced by the three models of evolution
# Labridae

# ----------------------------------------------

rm(list=ls())
source("R/packages.R")

# -----------------------------------------------------
require(here)
##---------------------------------------------------------
# analyses
# MCMC settings
nc<-3
ni<-10000
nb<-5000
nt<-10

# prior sigma
priors<-c(set_prior("normal(0,5)",class = "sigma"))

# ---------------------------------------------------------------------
# fish - Labridae
load (here("Output", "simulated_FD_BM_labridae.RData"))
load (here("Output", "simulated_FD_EB_labridae.RData"))
load (here("Output", "simulated_FD_OU_labridae.RData"))
load (here("Output", "empirical_FD_labridae.RData"))
#load (here("Output", "simulated_FD_random_labridae.RData"))

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
fig_params_fish
ggsave(here("Output", "Figures", "parameter_estimates_labridae.png"))

# =-----------------------------------------------------------------

# empirical results
empirical_results <- data.frame (SR= apply(sapply(empirical_FD,"[[","nbsp"),1,mean),
                                 FRic= apply(sapply(empirical_FD,"[[","FRic"),1,mean),
                                 FEve=apply(sapply(empirical_FD,"[[","FEve"),1,mean),
                                 Dataset= "Empirical")

# random traits
#simulated_results_random <- data.frame (SR= apply(sapply(simulated_FD_random,"[[","nbsp"),1,mean),
#                                    FRic= apply(sapply(simulated_FD_random,"[[","FRic"),1,mean),
#                                    FEve=apply(sapply(simulated_FD_random,"[[","FEve"),1,mean),
#                                    Dataset= "SimulatedRandom")

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
                     #simulated_results_random,
                     simulated_results_BM,
                     simulated_results_EB,
                     simulated_results_OU)

# add covariates
reef_covariates <- rbind(reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         #reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),],
                         reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),])

# bind
df_analyzes <- cbind(df_analyzes,
                     region = reef_covariates$Region,
      turbidity = scale(reef_covariates[,"BO_damean_lonlat"]),
      salinity = scale(reef_covariates[,"Present.Surface.Salinity.Mean"]),
      temperature = scale(reef_covariates[,"Present.Surface.Temperature.Mean"]),
      productivity = scale(reef_covariates[,"Present.Surface.Primary.productivity.Mean"]),
      latitude = scale(reef_covariates[,"decimalLatitude"]),
      area = scale(reef_covariates[,"Reef.Area"]))

# bind effort
df_analyzes <- cbind (df_analyzes,
                      effort_site[which(rowSums(subset_comm_data[[1]]) >0),])

# Remove NAs FRic
df_analyzes <- df_analyzes [!is.na(df_analyzes$FRic),] 

# correlations
cor(df_analyzes[,c(1:3,6:9,11)])

# Number of sites used in the analyzes
nrow(df_analyzes[df_analyzes$Dataset == "Empirical",])

# Number of species used in the analyzes
nrow(empirical_FD[[1]]$x.axes)

# log responses and predictors - SR and effort
df_analyzes$logFRic <- log(df_analyzes$FRic)
df_analyzes$logFEve <- log(df_analyzes$FEve)
df_analyzes$logSR <- log(df_analyzes$SR)
df_analyzes$logEffort <- log(df_analyzes$effort)

##---------------------------------------------------------

# run model (ancova)
model.ancova.FRic <- brm (brmsformula(logFRic ~ (logSR*Dataset)+logEffort,
                          sigma ~ logSR),
                          data=df_analyzes,
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

# plotting
p1<-plot(conditional_effects(model.ancova.FRic_sq,
                             method="posterior_epred",
                             re_formula=NA,
                             robust=T,
                             effects = "logSR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = c(0.75,0.3)) ,
         points=T,
         point_args = list (width = 0.25,alpha=0.3)) [[1]]  + 
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness (ln)") + 
  
  ylab ("Functional Richness (ln)")

p1

# compare slopes
m.lst.FRic <- emtrends (model.ancova.FRic, "Dataset", var="logSR")
m.lst_tab.FRic <- summary(m.lst.FRic,point.est = mean)
m.lst_tab.FRic

# Analyse residuals in function of covariates
df_analyzes$Residuals_FRic <- residuals(model.ancova.FRic)[,"Estimate"]

# run one LM per data set
mod_dat<-lapply (unique(df_analyzes$Dataset)[1:2], function (i) {
  
  # fit the model
  dat <- df_analyzes [which(df_analyzes$Dataset == i),]
  # do LM
  mod_dat <- brm (formula = Residuals_FRic ~ region +  Reef.Area + 
                    BO_damean_lonlat + 
                    Present.Surface.Salinity.Mean+
                    Present.Surface.Temperature.Mean+
                    decimalLatitude,
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
tab_model(mod_dat[[1]])
tab_model(mod_dat[[2]])

# run model (ancova)
model.ancova.FEve <- brm (brmsformula(log(FEve) ~ (logSR*Dataset)+logEffort,
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

# plotting
p2<-plot(conditional_effects(model.ancova.FEve,
                             method="posterior_epred",
                             re_formula=NA,
                             robust=T,
                             effects = "logSR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "none") ,
         points=T,
         point_args = list (width = 0.25,alpha=0.3)) [[1]]  + 
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (FEve)")

p2

# organize plots
pdf(here("Output","Figures","Fig3_labridae.pdf"), width=9,height=5)
  grid.arrange(p1,p2,nrow=1)
dev.off()

# compare slopes
m.lst.FEve <- emtrends (model.ancova.FEve, "Dataset", var="logSR")
m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)
m.lst_tab.FEve


# Analyse residuals in function of covariates
df_analyzes$Residuals_FEve <- residuals(model.ancova.FEve)[,"Estimate"]

# run one LM per data set
mod_dat_FEve <- lapply (unique(df_analyzes$Dataset)[1:2], function (i) {
  
  # fit the model
  dat <- df_analyzes [which(df_analyzes$Dataset == i),]
  # do LM
  mod_dat <- brm (formula = Residuals_FEve ~ region +  Reef.Area + 
                    BO_damean_lonlat + 
                    Present.Surface.Salinity.Mean+
                    Present.Surface.Temperature.Mean+
                    decimalLatitude,
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
      m.lst_tab.FRic,
      model.ancova.FEve,
      model.ancova.FEve_sq,
      m.lst.FEve,
      m.lst_tab.FEve,
      mod_dat,
      mod_dat_FEve,
      file=here("Output", "GLM_test_labridae.RData"))
#load(here("Output", "GLM_test_labridae.RData"))
# end
