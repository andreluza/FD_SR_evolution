# phylogenies Fig 1. (evotrait function)

# load functions
source("R/functions.R")
source("R/packages.R")

# simulate trait with varying signal
# highly conserved
niche_filling <-evotrait(nsp=10, phymin=3, phymax=+Inf, delta=TRUE)
# labile - convergent
labile  <-evotrait(nsp=10, phymin=-1, phymax=0.2,  delta=TRUE)
# random
brownian <-evotrait(nsp=10, phymin=0.95, phymax=1.05,  delta=TRUE)

# pdf
pdf(here("Output","trait_signal.pdf"),width=12,height=5)
par (mfrow=c(3,3))

# We want to confirm that the trait is  convergent. We can make a graph for that using the first phylogeny and the first vector of simulated trait values.
b<-contMap(niche_filling[[1]] [[1]], 
           niche_filling[[2]] [[1]],plot=F)
## change color scheme
b<-setMap(b,
          c("white","#F6F6F6","#9D9D9D",
            "#2D2424"))

plot(b,fsize=c(1,0.8),
     leg.txt="Trait value",
     lwd=7)


# We want to confirm that the trait is  highly conserved. We can make a graph for that using the first phylogeny and the first vector of simulated trait values.
d <- contMap (niche_filling[[1]] [[1]],
              brownian[[2]] [[1]],
              plot=F)
## change color scheme
d<-setMap(d,
          c("white","#F6F6F6","#AAA492","#9D9D9D",
            "#2D2424"))

plot(d,fsize=c(1,0.8),
     leg.txt="Trait value",
     lwd=7)

# We want to confirm that the trait is highly conserved. We can make a graph for that using the first phylogeny and the first vector of simulated trait values.
a<-contMap(niche_filling[[1]] [[1]], 
           labile[[2]] [[1]],plot=F)
## change color scheme
a<-setMap(a,
          c("white","#F6F6F6","#AAA492","#9D9D9D",
            "#2D2424"))

plot(a,fsize=c(1,0.8),
     leg.txt="Trait value",
     lwd=7)


dev.off()
