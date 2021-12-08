
# maps
# load functions
source("R/packages.R")
source("R/functions.R")

# load fish data
load("image_fish.RData")

# coordinates
coords_fish <- data.frame (Lon = aggregate(peixes$Lon, by =list(data=peixes$eventID_MOD), FUN=mean),
                           Lat = aggregate(peixes$Lat, by =list(data=peixes$eventID_MOD), FUN=mean)[,2])

# axes
axes_fish <- empirical_FD[[1]]$x.axes

# community
comm_fish <- match_comm_data[[1]]$comm

# closest to four coords
seq_search <- c(-1,-15,-20,-27)
closest_fish <- lapply (as.list(seq_search), function (i) 
  closest (coords_fish$Lat,i))

# find the selected comms
sel_comms_fish <- lapply (closest_fish, function (i) 
  which(coords_fish$Lat == i[1])[1])
# which ones
# sites [unlist(sel_comms_fish)]

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
                     limits=c(0,10),
                     breaks=seq(0,10,2.5))+ 
          theme(axis.text = element_text(size=6),
                axis.title=element_text(size=8))
        ; # return
        plotA

})

array_fish <- grid.arrange(fish_space[[1]]+theme(legend.position="top",
                                                 axis.title.x = element_blank(),
                                                 axis.text.x = element_blank()),
                            fish_space[[2]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[3]]+theme(legend.position="none",
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank()),
                            fish_space[[4]]+theme(legend.position="none"),
                           ncol=1)

# --------------------------------------------
## rodents
load("image_rodents.RData")

# axes
axes_rodents <- empirical_FD$x.axes

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
             label=paste ("SR=", empirical_FD$nbsp[i],
                          "\nFRic=", round(empirical_FD$FRic[i],2),
                          "\nFEve=", round(empirical_FD$FEve[i],2))
    ) + 
    geom_point(data=pk,aes (A1,A2,size=(abund)),
               alpha=0.5,col="yellow") + 
    scale_size(name="CPUE",
               limits=c(0,0.3),
               breaks=seq(0,0.3,0.07)) + 
    theme(axis.text = element_text(size=6),
          axis.title=element_text(size=8))
  ; # return
  plotA
  
})

array_rodents <- grid.arrange(rodents_space[[1]]+theme(legend.position="top",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                              rodents_space[[2]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                              rodents_space[[3]]+theme(legend.position="none",
                                                    axis.title.x = element_blank(),
                                                    axis.text.x = element_blank()),
                              rodents_space[[4]]+theme(legend.position="none"),
                              ncol=1)


# map of points
# mapa mundi
world <- ne_countries(scale = "medium", returnclass = "sf")

# cortar o mapa para ver a america do Sul e parte da central
wm <- ggplot() + 
  geom_sf (data=world, size = 0.1, 
           fill= "#aaaaaa",colour="#aaaaaa") +
  coord_sf (xlim = c(-20,-65),  ylim = c(-35, 4), expand = T,crs = st_crs(4326)) +
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

# rodents
map_rodents <- wm +  geom_point (data=spatial_effort_data_LF,aes(x=Longitude, y=Latitude),
                                 size=2,
                                col="#FFF9B6")
# fish
map_fish_rodents <- map_rodents + 
  geom_point(data=coords_fish, aes(x=Lon.x, y=Lat),size=2,
             col="#0e49b5")

# arrange all
pdf(file=here("output","Fig2.pdf"),height=7,width=10)

grid.arrange(array_rodents,
             map_fish_rodents, 
             array_fish,
             ncol=5,nrow=5,
             layout_matrix = rbind (c(1,2,2,2,3),
                                    c(1,2,2,2,3),
                                    c(1,2,2,2,3),
                                    c(1,2,2,2,3),
                                    c(1,2,2,2,3)))

dev.off()

