# output
source("R/packages.R")

# ----------------------------------------------
# rodents


# -----------------------------------------------------
require(here)
load (here("Output", "simulated_FD_BM_rodents.RData"))
load (here("Output", "simulated_FD_EB_rodents.RData"))
load (here("Output", "simulated_FD_OU_rodents.RData"))
load (here("Output", "empirical_FD_rodents.RData"))

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
# density plot
fig_params_rodents <- ggplot(df_params, 
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
fig_params_rodents

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
model.ancova.FRic <- brm (FRic ~ poly(SR,2)*Dataset,
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)

# plotting
p1<-plot(conditional_effects(model.ancova.FRic,
                             method="fitted",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "top") ,
         points=T) [[1]] + 
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Richness (FRic)")

p1

# compare slopes

m.lst.FRic <- emtrends (model.ancova.FRic, "Dataset", var="SR")
# m.lst_tab.FRic <- summary(m.lst.FRic,point.est = mean)

# run model (ancova)
model.ancova.FEve <- brm (FEve ~ poly(SR,2)*Dataset,
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)

# plotting
p2<-plot(conditional_effects(model.ancova.FEve,
                             method="fitted",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "top") ,
         points=T) [[1]] + 
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (FEve)")

# organize plots
pdf(here("Output","Fig3_SM.pdf"), width=9,height=5)
grid.arrange(p1,p2,nrow=1)
dev.off()

# compare slopes
m.lst.FEve <- emtrends (model.ancova.FEve, "Dataset", var="SR")
# m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

# save results
save (model.ancova.FRic,
      m.lst.FRic,
      model.ancova.FEve,
      m.lst.FEve,
      file=here("Output", "GLM_test_rodents.RData"))


# ---------------------------------------------------------------------
# fish
load (here("Output", "simulated_FD_BM_fish.RData"))
load (here("Output", "simulated_FD_EB_fish.RData"))
load (here("Output", "simulated_FD_OU_fish.RData"))
load (here("Output", "empirical_FD_fish.RData"))

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
  facet_wrap(~model+Parameter,scales = "free")+
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

# prepare data for fig 4
# neutral SES fish
SR_fish<-apply(sapply(empirical_FD,"[[","nbsp"),1,mean)
obsFEve_fish <- apply(sapply(empirical_FD,"[[","FEve"),1,mean) 
meanFEveSim_fish <- apply(sapply(simulated_FD,"[[","FEve"),1,mean)
sdFEveSim_fish <- apply(sapply(simulated_FD,"[[","FEve"),1,sd)

##---------------------------------------------------------
# analyses
# MCMC settings
nc<-3
ni<-10000
nb<-5000
nt<-10

# run model (ancova)
model.ancova.FRic <- brm (FRic ~ poly(SR,2)*Dataset,
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FRic)
tab_model(model.ancova.FRic)

# plotting
p1<-plot(conditional_effects(model.ancova.FRic,
                             method="fitted",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "top") ,
         points=T) [[1]] + 
  
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Richness (FRic)")


# compare slopes

m.lst.FRic <- emtrends (model.ancova.FRic, "Dataset", var="SR")
m.lst_tab.FRic <- summary(m.lst.FRic,point.est = mean)

# run model (ancova)
model.ancova.FEve <- brm (FEve ~ poly(SR,2)*Dataset,
                          data=df_analyzes,
                          family = gaussian (link="identity"),
                          chains=nc,
                          iter = ni,
                          warmup = nb,
                          thin=nt)

# summary of results
summary (model.ancova.FEve)
tab_model(model.ancova.FEve)

# plotting
p2<-plot(conditional_effects(model.ancova.FEve,
                             method="fitted",
                             re_formula=NA,
                             robust=T,
                             effects = "SR:Dataset",
                             points=T,
                             prob = 0.95),
         
         theme = theme_classic() +
           
           theme (axis.title = element_text(size=15),
                  axis.text = element_text(size=12),
                  legend.position = "top") ,
         points=T) [[1]] + 
  
  scale_color_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  scale_fill_manual(values=c("#000000","#0F00FF","#D98C00","#A4EBF3")) + 
  
  xlab("Species richness gradient") + 
  
  ylab ("Functional Evenness (FEve)")

# organize plots
pdf(here("Output","Fig3_fish.pdf"), width=9,height=5)
grid.arrange(p1,p2,nrow=1)
dev.off()

# compare slopes
m.lst.FEve <- emtrends (model.ancova.FEve, "Dataset", var="SR")
m.lst_tab.FEve <- summary(m.lst.FEve,point.est = mean)

# save results
save (model.ancova.FRic,
      m.lst.FRic,
      m.lst_tab.FRic,
      model.ancova.FEve,
      m.lst.FEve,
      m.lst_tab.FEve,
      file=here("Output", "GLM_test.RData"))


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
             size=1.75)  +
  scale_color_manual(values=c("#E69F00","#999999",  "#56B4E9"),
                     name = "Deviations from\na BM model",
                     labels = c("Positive",
                                "None",
                                "Negative"))+
  theme(legend.position = c(0.90,0.15)) + 
  ylab ("Standardized Effect Size") + 
  xlab ("Species richness gradient")
  
  

