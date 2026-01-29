# output
source("R/packages.R")

# ----------------------------------------------
# labridae


# -----------------------------------------------------
require(here)
load (here("Output", "simulated_FD_BM_labridae.RData"))
load (here("Output", "simulated_FD_EB_labridae.RData"))
load (here("Output", "simulated_FD_OU_labridae.RData"))
load (here("Output", "empirical_FD_labridae.RData"))

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

# plot
require(ggplot2)
# density plot
fig_params_labridae <- ggplot(df_params, 
                        aes(x=Estimates,
                            group=Trait,
                            color=Trait,
                            fill=Parameter)) +
  geom_density(size=1,alpha=0.05)+
  scale_fill_viridis_d(option = "magma")+
  scale_colour_viridis_d(option = "magma")+
  theme_classic()  + 
  facet_wrap(~model+Parameter,scales = "free")+
  theme (legend.position = c(0.81,0.23),
         axis.text.x = element_text(size=7))
fig_params_labridae

# ----------------------------------------------------------------
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
simulated_results_OU <- data.frame (SR= apply(sapply(simulated_FD_OU,"[[","nbsp"),1,mean),
                                    FRic= apply(sapply(simulated_FD_OU,"[[","FRic"),1,mean),
                                    FEve=apply(sapply(simulated_FD_OU,"[[","FEve"),1,mean),
                                    Dataset = "SimulatedOU")

# bind them
df_analyzes <- rbind(empirical_results,
                     simulated_results_BM,
                     simulated_results_EB,
                     simulated_results_OU)

# prepare data for fig 4
# neutral SES labridae
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

model.ancova.FRic<-lm(FRic ~ poly(SR,2)*Dataset,
                    data=df_analyzes)
  
# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)
m.lst.FRic <- emmeans::emtrends (model.ancova.FRic, "Dataset", var="SR")


# FEve
model.ancova.FEve<-lm(FEve ~ poly(SR,2)*Dataset,
data=df_analyzes)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)
m.lst.FEve <- emmeans::emtrends (model.ancova.FEve, "Dataset", var="SR")


# fric
p1<-ggplot (df_analyzes, aes (x=SR, 
                          y=FRic, 
                          group=Dataset,
                          fill=Dataset,
                          colour=Dataset)) + 
  geom_point(position = position_jitter())+
  #stat_smooth(method = "lm", formula = y ~ x, colour = "red") +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  theme_bw()+
  scale_fill_viridis_d()+
  scale_colour_viridis_d()+
  xlim (c(3,8))+
  ylim(c(0,1))+ 
  theme_classic()+
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Richness (FRic)")

# feve
p2<-ggplot (df_analyzes, aes (x=SR, 
                          y=FEve, 
                          group=Dataset,
                          fill=Dataset,
                          colour=Dataset)) + 
  geom_point(position = position_jitter())+
  #stat_smooth(method = "lm", formula = y ~ x, colour = "red") +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2))+
  theme_bw()+
  scale_fill_viridis_d()+
  scale_colour_viridis_d()+
  xlim (c(3,8))+
  ylim(c(0,1))+ 
  theme_classic()+
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (FEve)")


# organize plots
pdf(here("Output","Fig3_labridae1.pdf"), width=9,height=5)
grid.arrange(p1+theme(legend.position = "top"),
             p2+theme(legend.position = "top"),nrow=1)
dev.off()

