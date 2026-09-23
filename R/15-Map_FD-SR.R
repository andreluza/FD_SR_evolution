
# Maps and functional spaces - Fig. 1 and 4
# load functions
rm(list=ls())
source("R/packages.R")
source("R/functions.R")
require(dplyr)

# load fish data
load(here ("Processed_data","image_labridae.RData"))
load(here ( "Processed_data","Spatial_covs_data_reefs.RData"))
load(here ("Output","empirical_FD_labridae.RData"))

# coordinates
#coords_fish <- data.frame (Lon = aggregate(peixes$Lon, by =list(data=peixes$eventID_MOD), FUN=mean),
#                           Lat = aggregate(peixes$Lat, by =list(data=peixes$eventID_MOD), FUN=mean)[,2])
# add covariates
coords_fish <- reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),]

# Remove NAs FRic
#df_analyzes <- df_analyzes [!is.na(df_analyzes$FRic),] 

# axes
axes_fish <- empirical_FD[[1]]$x.axes
empirical_FD[[1]]$x.values/sum(empirical_FD[[1]]$x.values)

# bind FD

coords_fish$SR <- rowMeans(sapply (empirical_FD, "[[", "nbsp"))
coords_fish$FRic <- rowMeans(sapply (empirical_FD, "[[", "FRic"))
# bind FD
coords_fish$FEve <- rowMeans(sapply (empirical_FD, "[[", "FEve"))

# Value used to transform the data
coeff <- .065

plot_wrasses_FRic <- coords_fish %>%
  ggplot ()+
  
  geom_point(aes(x=decimalLatitude,y=SR),col="cyan4",shape=17) + 
  geom_smooth(aes(x=decimalLatitude,y=SR),se=F,col="cyan4") + 
  # islands
  geom_point(data = coords_fish %>%
               filter (Region == "oc_isl"),
               aes(x=decimalLatitude,y=SR),col="black",shape=3,size=3) + 

  geom_point(aes(x=decimalLatitude,y=FRic/coeff),col="red4") + # Divide by 10 to get the same range than the temperature
  geom_smooth(aes(x=decimalLatitude,y=FRic/coeff),col="red4",se=F) + 

  scale_y_continuous(
    
    # Features of the first axis
    name = "Species richness",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Functional richness")
  )+
  theme_bw(base_size = 14)+
  labs(x="Latitude",title="Wrasses")+
  theme(axis.title.y.right =  element_text(colour = "red4"),
        axis.text.y.right = element_text(colour = "red4"),
        axis.title.y.left =  element_text(colour = "cyan4"),
        axis.text.y.left = element_text(colour = "cyan4")) + 
  geom_vline(aes(xintercept= -18),linetype=3)
  
plot_wrasses_FRic

# FEVE
coeff <- .05
plot_wrasses_FEve <- coords_fish %>%
  ggplot ()+
  
  geom_point(aes(x=decimalLatitude,y=SR),col="cyan4",shape=17) + 
  geom_smooth(aes(x=decimalLatitude,y=SR),se=F,col="cyan4") + 
  ylim(c(0,15))+

  # islands
  geom_point(data = coords_fish %>%
               filter (Region == "oc_isl"),
               aes(x=decimalLatitude,y=SR),col="black",shape=3,size=3) + 

  geom_point(aes(x=decimalLatitude,y=FEve/coeff),col="red4") + # Divide by 10 to get the same range than the temperature
  geom_smooth(aes(x=decimalLatitude,y=FEve/coeff),col="red4",se=F) + 
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Species richness",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Functional evenness")
  )+
  
  theme_bw(base_size = 14)+
  labs(x="Latitude",title="")+
  theme(axis.title.y.right =  element_text(colour = "red4"),
        axis.text.y.right = element_text(colour = "red4"),
        axis.title.y.left =  element_text(colour = "cyan4"),
        axis.text.y.left = element_text(colour = "cyan4"))+ 
  geom_vline(aes(xintercept= -18),linetype=3)
  
plot_wrasses_FEve

# community
comm_fish <- match_comm_data[[1]]$comm

# Community that is not in the functional estimates
comm_fish <- comm_fish [as.numeric(names(empirical_FD[[1]]$FRic)),]

# closest to four coords
seq_search <- c(-5,-10,-17,-24.316838)#-24)
closest_fish <- lapply (as.list(seq_search), function (i) 
  closest (coords_fish$decimalLatitude,i))

# find the selected comms
sel_comms_fish <- lapply (closest_fish, function (i) 
  which(coords_fish$decimalLatitude == i[1])[1])
# which ones
# sites [unlist(sel_comms_fish)]

range_plot_fish <- range(comm_fish [unlist(sel_comms_fish),])

# find spp in each community and built the plot 

fish_space <- lapply (sel_comms_fish, function (i) {
  
        # community subset
        
        sel_comm_data <- comm_fish[i,]
        pres_fish <- sel_comm_data[which(sel_comm_data>0)]
          
        # complete trait space
        all <- cbind (axes_fish[,1:2],ext = F)
        a <- all [chull(all[,1:2], y = NULL),]
        
        # space occupied by the community
        setB<-cbind(all, ext1=ifelse(rownames(all) %in% names(pres_fish),
                                       F,
                                       T))
        pk <-setB[which(setB$ext1==F),]
        f <- pk [chull(pk, y = NULL),]
        # abundance on pk
        pk$abund <- as.numeric(pres_fish[match(rownames(pk), names(pres_fish))])
          
        # plot space
        plotA <- ggplot(a, aes(A1, A2)) + 
            geom_point() + theme_bw()+
            geom_polygon(data=a, aes (A1,A2),alpha=0.5,fill="gray") + 
            geom_polygon(data=f, aes (A1,A2,group=ext1, fill=ext1),alpha=0.5,
                         fill="black",size=3) +
            xlim(min (a$A1)-0.2,max (a$A1)+0.2) + 
          annotate("text",x=0.3,y=0.15,size=2.5,
                   label=paste ("SR=", empirical_FD[[1]]$nbsp[i],
                                "\nFRic=", round(empirical_FD[[1]]$FRic[i],2),
                                "\nFEve=", round(empirical_FD[[1]]$FEve[i],2))
                   ) + 
          geom_point(data=pk,aes (A1,A2,size=(abund)),
                     alpha=0.5,col="cyan") + 
          scale_size(name="CPUE",
                     limits=c(range_plot_fish[1],range_plot_fish[2]),
                     breaks=seq(range_plot_fish[1],range_plot_fish[2],2))+ 
          theme(axis.text = element_text(size=6),
                axis.title=element_text(size=8))
        ; # return
        plotA

})

array_wrasses <- grid.arrange(fish_space[[1]]+theme(legend.position=c(0.7,0.9),
                                                 legend.direction = "horizontal",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                            fish_space[[2]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[3]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[4]]+theme(legend.position="none",
                                                  axis.title.x = element_blank()),
                           ncol=1)


# Grunts ------------------------------------
# load fish data
load(here ("Processed_data","image_haemulidae.RData"))
load(here ("Output","empirical_FD_haemulidae.RData"))

# coordinates
#coords_fish <- data.frame (Lon = aggregate(peixes$Lon, by =list(data=peixes$eventID_MOD), FUN=mean),
#                           Lat = aggregate(peixes$Lat, by =list(data=peixes$eventID_MOD), FUN=mean)[,2])
# add covariates
coords_fish_grunts <- reef_covariates[which(rowSums(subset_comm_data[[1]]) >0),]

# Remove NAs FRic
#df_analyzes <- df_analyzes [!is.na(df_analyzes$FRic),] 

# axes
axes_fish_grunts <- empirical_FD[[1]]$x.axes
empirical_FD[[1]]$x.values/sum(empirical_FD[[1]]$x.values)

# bind FD

coords_fish_grunts$SR <- rowMeans(sapply (empirical_FD, "[[", "nbsp"))
coords_fish_grunts$FRic <- rowMeans(sapply (empirical_FD, "[[", "FRic"))
# bind FD
coords_fish_grunts$FEve <- rowMeans(sapply (empirical_FD, "[[", "FEve"))

# Value used to transform the data
coeff <- .05

plot_grunts_FRic <- coords_fish_grunts %>%
  ggplot ()+
  
  geom_point(aes(x=decimalLatitude,y=SR),col="cyan4",shape=17) + 
  geom_smooth(aes(x=decimalLatitude,y=SR),se=F,col="cyan4") + 
  # islands
  geom_point(data = coords_fish_grunts %>%
               filter (Region == "oc_isl"),
               aes(x=decimalLatitude,y=SR),col="black",shape=3,size=3) + 
    
  geom_point(aes(x=decimalLatitude,y=FRic/coeff),col="red4") + # Divide by 10 to get the same range than the temperature
  geom_smooth(aes(x=decimalLatitude,y=FRic/coeff),col="red4",se=F) + 

  scale_y_continuous(
    
    # Features of the first axis
    name = "Species richness",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Functional richness")
  )+
  theme_bw(base_size = 14)+
  labs(x="Latitude",title="Grunts")+
  theme(axis.title.y.right =  element_text(colour = "red4"),
        axis.text.y.right = element_text(colour = "red4"),
        axis.title.y.left =  element_text(colour = "cyan4"),
        axis.text.y.left = element_text(colour = "cyan4"))+
  geom_vline(aes(xintercept= -18),linetype=3)


plot_grunts_FRic

# FEVE
coeff <- 0.09
plot_grunts_FEve <- coords_fish_grunts %>%
  ggplot ()+
  
  geom_point(aes(x=decimalLatitude,y=SR),col="cyan4",shape=17) + 
  geom_smooth(aes(x=decimalLatitude,y=SR),se=F,col="cyan4") + 
  # islands
  geom_point(data = coords_fish_grunts %>%
               filter (Region == "oc_isl"),
               aes(x=decimalLatitude,y=SR),col="black",shape=3,size=3) + 

  geom_point(aes(x=decimalLatitude,y=FEve/coeff),col="red4") + # Divide by 10 to get the same range than the temperature
  geom_smooth(aes(x=decimalLatitude,y=FEve/coeff),col="red4",se=F) + 

  scale_y_continuous(
    
    # Features of the first axis
    name = "Species richness",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Functional evenness")
  )+
  theme_bw(base_size = 14)+
  labs(x="Latitude",title="")+
  theme(axis.title.y.right =  element_text(colour = "red4"),
        axis.text.y.right = element_text(colour = "red4"),
        axis.title.y.left =  element_text(colour = "cyan4"),
        axis.text.y.left = element_text(colour = "cyan4"))+
  geom_vline(aes(xintercept= -18),linetype=3)

plot_grunts_FEve + coord_flip()

# community
comm_fish <- match_comm_data[[1]]$comm

# Community that is not in the functional estimates
comm_fish <- comm_fish [as.numeric(names(empirical_FD[[1]]$FRic)),]

# closest to four coords
seq_search <- c(-5,-10,-17,-24.316838)#-24)
closest_fish <- lapply (as.list(seq_search), function (i) 
  closest (coords_fish_grunts$decimalLatitude,i))

# find the selected comms
sel_comms_fish <- lapply (closest_fish, function (i) 
  which(coords_fish_grunts$decimalLatitude == i[1])[1])

# which ones
# sites [unlist(sel_comms_fish)]
range_plot_fish <- range(comm_fish [unlist(sel_comms_fish),])

# find spp in each community and built the plot 
fish_space <- lapply (sel_comms_fish, function (i) {
  
        # community subset
        
        sel_comm_data <- comm_fish[i,]
        pres_fish <- sel_comm_data[which(sel_comm_data>0)]
          
        # complete trait space
        all <- cbind (axes_fish_grunts[,1:2],ext = F)
        a <- all [chull(all[,1:2], y = NULL),]
        
        # space occupied by the community
        setB<-cbind(all, ext1=ifelse(rownames(all) %in% names(pres_fish),
                                       F,
                                       T))
        pk <-setB[which(setB$ext1==F),]
        f <- pk [chull(pk, y = NULL),]
        # abundance on pk
        pk$abund <- as.numeric(pres_fish[match(rownames(pk), names(pres_fish))])
          
        # plot space
        plotA <- ggplot(a, aes(A1, A2)) + 
            geom_point() + theme_bw()+
            geom_polygon(data=a, aes (A1,A2),alpha=0.5,fill="gray") + 
            geom_polygon(data=f, aes (A1,A2,group=ext1, fill=ext1),alpha=0.5,
                         fill="black",size=3) +
            xlim(min (a$A1)-0.2,max (a$A1)+0.2) + 
          annotate("text",x=0.3,y=0.15,size=2.5,
                   label=paste ("SR=", empirical_FD[[1]]$nbsp[i],
                                "\nFRic=", round(empirical_FD[[1]]$FRic[i],2),
                                "\nFEve=", round(empirical_FD[[1]]$FEve[i],2))
                   ) + 
          geom_point(data=pk,aes (A1,A2,size=(abund)),
                     alpha=0.5,col="cyan") + 
          scale_size(name="CPUE",
                     limits=c(range_plot_fish[1],range_plot_fish[2]),
                     breaks=seq(range_plot_fish[1],range_plot_fish[2],6))+ 
          theme(axis.text = element_text(size=6),
                axis.title=element_text(size=8))
        ; # return
        plotA

})

array_grunts <- grid.arrange(fish_space[[1]]+theme(legend.position=c(0.7,0.9),
                                                 legend.direction = "horizontal",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                            fish_space[[2]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[3]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[4]]+theme(legend.position="none",
                                                  axis.title.x = element_blank()),
                           ncol=1)


# --------------------------------------------
## rodents
load(here ("Processed_data","image_rodents.RData"))
load(here ("Output", "empirical_FD_rodents.RData"))

# axes
axes_rodents <- empirical_FD[[1]]$x.axes
empirical_FD[[1]]$x.values/sum(empirical_FD[[1]]$x.values)

# bind FD
spatial_effort_data_LF$SR <- rowMeans(sapply (empirical_FD, "[[", "nbsp"))
spatial_effort_data_LF$FRic <- rowMeans(sapply (empirical_FD, "[[", "FRic"))
# bind FD
spatial_effort_data_LF$FEve <- rowMeans(sapply (empirical_FD, "[[", "FEve"))

# Value used to transform the data
coeff <- .065

plot_rodents_FRic <- spatial_effort_data_LF %>%
  ggplot ()+
  
  geom_point(aes(x=Latitude,y=SR),col="cyan4",shape=17) + 
  geom_smooth(aes(x=Latitude,y=SR),se=F,col="cyan4") + 

  geom_point(aes(x=Latitude,y=FRic/coeff),col="red4") + # Divide by 10 to get the same range than the temperature
  geom_smooth(aes(x=Latitude,y=FRic/coeff),col="red4",se=F) + 

  scale_y_continuous(
    
    # Features of the first axis
    name = "Species richness",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Functional richness")
  )+
  theme_bw(base_size = 14)+
  labs(x="Latitude",title="Rodents")+
  theme(axis.title.y.right =  element_text(colour = "red4"),
        axis.text.y.right = element_text(colour = "red4"),
        axis.title.y.left =  element_text(colour = "cyan4"),
        axis.text.y.left = element_text(colour = "cyan4"))+
  geom_vline(aes(xintercept=-20),linetype=3)

plot_rodents_FRic

# FEVE

# Value used to transform the data
coeff <- .075

plot_rodents_FEve <- spatial_effort_data_LF %>%
  ggplot ()+
  
  geom_point(aes(x=Latitude,y=SR),col="cyan4",shape=17) + 
  geom_smooth(aes(x=Latitude,y=SR),se=F,col="cyan4") + 

  geom_point(aes(x=Latitude,y=FEve/coeff),col="red4") + # Divide by 10 to get the same range than the temperature
  geom_smooth(aes(x=Latitude,y=FEve/coeff),col="red4",se=F) + 

  scale_y_continuous(
    
    # Features of the first axis
    name = "Species richness",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Functional evenness")
  )+
  theme_bw(base_size = 14)+
  labs(x="Latitude",title="")+
  theme(axis.title.y.right =  element_text(colour = "red4"),
        axis.text.y.right = element_text(colour = "red4"),
        axis.title.y.left =  element_text(colour = "cyan4"),
        axis.text.y.left = element_text(colour = "cyan4"))+
  geom_vline(aes(xintercept=-20),linetype=3)

plot_rodents_FEve

# community
comm_rodents <- match_comm_data[[1]]$comm

# closest to four coords
closest_rodents <- lapply (as.list(seq_search), function (i) 
  closest (spatial_effort_data_LF$Latitude,i))

# find the selected comms
sel_comms_rodents <- lapply (closest_rodents, function (i) 
  which(spatial_effort_data_LF$Latitude == i[1])[1])

# which ones
# sites [unlist(sel_comms_rodents)]
range_plot_rodents <- range(comm_rodents [unlist(sel_comms_rodents),])

# find spp in each community and built the plot 

rodents_space <- lapply (sel_comms_rodents, function (i) {
  
  # community subset
  sel_comm_data <- comm_rodents[i,]
  pres_rodents <- sel_comm_data[which(sel_comm_data>0)]
  
  # complete trait space
  all <- cbind (axes_rodents[,1:2],ext = F)
  a <- all [chull(all[,1:2], y = NULL),]
  
  # space occupied by the community
  setB<-cbind(all, ext1=ifelse(rownames(all) %in% names(pres_rodents),
                               F,
                               T))
  pk <-setB[which(setB$ext1==F),]
  f <- pk [chull(pk, y = NULL),]
  
  # abundance on pk
  pk$abund <- as.numeric(pres_rodents[match(rownames(pk), names(pres_rodents))])
  
  # plot space
  plotA <- ggplot(a, aes(A1, A2)) + 
    geom_point() + theme_bw()+
    geom_polygon(data=a, aes (A1,A2),alpha=0.5,fill="gray") + 
    geom_polygon(data=f, aes (A1,A2,group=ext1, fill=ext1),alpha=0.5,
                 fill="black",size=3) +
    xlim(min (a$A1)-0.2,max (a$A1)+0.2) + 
    annotate("text",x=2.5,y=2,size=2.5,
             label=paste ("SR=", empirical_FD$UNTITLED$nbsp[i],
                          "\nFRic=", round(empirical_FD$UNTITLED$FRic[i],2),
                          "\nFEve=", round(empirical_FD$UNTITLED$FEve[i],2))
    ) + 
    geom_point(data=pk,aes (A1,A2,size=(abund)),
               alpha=0.5,col="yellow") + 
    scale_size(name="CPUE",
               limits=c(range_plot_rodents[1],range_plot_rodents[2]),
               breaks=seq(range_plot_rodents[1],range_plot_rodents[2],0.002)) + 
    theme(axis.text = element_text(size=6),
          axis.title=element_text(size=8))
  ; # return
  plotA
  
})

array_rodents <- grid.arrange(rodents_space[[1]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                              rodents_space[[2]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                              rodents_space[[3]]+theme(legend.position=c(0.5,0.1),
                                                       legend.direction = "horizontal",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                              rodents_space[[4]]+theme(legend.position="none",
                                                       axis.title.x = element_blank()),
                              ncol=1)


# map of points
# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# cortar o mapa para ver a america do Sul e parte da central
wm <- ggplot() + 
  geom_sf (data=world, size = 0.1, 
           fill= "#aaaaaa",colour="#aaaaaa") +
  coord_sf (xlim = c(-20,-80),  ylim = c(-55, 10), expand = T,crs = st_crs(4326)) +
  theme_bw() + #xlab ("Longitude")  + ylab ("Latitude") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#f4f9f9",#darkslategray1
                                        colour = "#f4f9f9"),
        axis.text = element_text(size=5),
        axis.ticks=element_line(size=1),
        axis.title = element_text(size=6),
        title = element_blank(),
        plot.margin = unit(c(0,-0.8,0,0.3), "cm")) +
  xlab("Longitude") + ylab("Latitude")  
  
# parallel line
wm <- wm + annotate("segment", 
              x = -68, xend = -30, 
              y = seq_search,  yend = seq_search, 
              colour = "white",size=1,
              linetype = "dashed",alpha=0.5)
wm

# rodents
map_rodents <- wm +  geom_point (data=spatial_effort_data_LF,aes(x=Longitude, y=Latitude),
                                 size=2,
                                col="yellow")
# fish
map_fish_rodents <- map_rodents + 
  geom_point(data=coords_fish, aes(x=decimalLongitude, y=decimalLatitude),size=2,
             col="cyan")
map_fish_rodents

# arrange all
pdf(file=here("Output","Figures","Fig1.pdf"),height=7,width=10)

  grid.arrange(array_rodents,
               map_fish_rodents, 
               array_wrasses,
               array_grunts,
               
               ncol=5,nrow=8,
               layout_matrix = rbind (c(1,2,2,2,3,4),
                                      c(1,2,2,2,3,4),
                                      c(1,2,2,2,3,4),
                                      c(1,2,2,2,3,4),
                                      c(1,2,2,2,3,4)))

dev.off()

# arrange all
ggsave(file=here("Output","Figures","Fig2.pdf"),

  grid.arrange( plot_rodents_FRic+coord_flip(),
                plot_wrasses_FRic+coord_flip(),
                plot_grunts_FRic+coord_flip(),
                # FEve
                plot_rodents_FEve+coord_flip(),
                plot_wrasses_FEve+coord_flip(),
                plot_grunts_FEve+coord_flip(),
                
                ncol=3,nrow=2),
  height=7,width=12)
               

