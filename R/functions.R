
### 

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}


## toupper for species names

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Here is the function (copy and paste it into R):
# frm here : https://rfunctions.blogspot.com/2016/08/simulating-phylogenetically-conserved.html

evotrait<-function(nsp=100,ntree=1,phymin=0.9,phymax=1.1,delta=FALSE,deltaval=0.01,noise=FALSE,meannoise=0,sdnoise="sdtrait",a=0,mu=0,sig2=1,bounds=c(-Inf,Inf),test=FALSE, nsim=999){
  listtree<-list()
  listtrait<-list()
  ksig<-numeric()
  pvals<-numeric()
  for(i in 1:ntree){
    physig<--Inf
    x<-1
    while(physig<phymin | physig>phymax | max(x)<=0){
      listtree[[i]]<-pbtree(n=nsp)
      if(delta==TRUE){deltatree<-rescale(listtree[[i]],"delta")}
      if(delta==TRUE){dtree<-deltatree(deltaval)}
      ifelse(delta==TRUE,x<-fastBM(dtree,a=a,mu=mu,sig2=sig2,bounds=bounds),x<-fastBM(listtree[[i]],a=a,mu=mu,sig2=sig2))
      if(noise==TRUE){x<-x+rnorm(nsp,mean = meannoise, sd = ifelse(sdnoise=="sdtrait",sd(x),sdnoise))}
      phys<-phylosig(listtree[[i]], x, method="K")
      physig<-phys[[1]]}
    if(test==TRUE){pval<-phylosig(listtree[[i]], x, method="K",test=test,nsim=nsim)[[2]]}
    ksig[i]<-physig
    if(test==TRUE){pvals[i]<-pval}
    listtrait[[i]]<-x}
  return(list(trees=listtree,traits=listtrait,ksig=ksig, pvals=pvals ))}



# Arguments:

# sp = pecies in simulated phylogenies. Default=100.
# phy= phylogenies to be simulated. Default=1.
# phymin: accepted phylogenetic signal in a simulated trait (Blomberg's K). Default=0.9.
# phymax: accepted phylogenetic signal in a simulated trait (Blomberg's K). Default=1.1.
# delta: if sim trees need to be delta transformed (Pagel's delta). If TRUE, trait evolution gets increasingly lower as it gets near the tips of a phylogeny. Default=FALSE.
# deltaval: value for rescaling the phylogeny under Pagel's delta transformation. Use low values to get highly conserved traits. Default=0.01.
# noise: if random normally distributed values should be added to the simulated trait. Increase the standard deviation of these values (sdnoise) to get highly labile traits. Default=FALSE.
# meannoise: mean of random normally distributed values. Default=0.
# sdnoise: standard deviation of the normally distributed values. Default="sdtrait" (standard deviation of the simulated trait).
# a: when simulating trait evolution, "a value for ancestral state at the root node" (see "fastBM" function; "phytools" package). Default=0.
# mu: when simulating trait evolution, "an optional value for the mean of random normal changes along branches of the tree - can be used to simulate a trend if mu!=0" (see "fastBM" function; "phytools" package). Default=0.
# sig2: when simulating trait evolution, "instantaneous variance of the BM process" (see "fastBM" function; "phytools" package). Default=1.
# bounds: when simulating trait evolution, "a vector with the lower and upper bounds (respectively) for bounded Brownian simulation - by default simulation is unbounded" (see "fastBM" function; "phytools" package). Default=c(-Inf,Inf).
# test: whether or not to calculate a p value for Blomberg's K (see "phylosig" function; "phytools" package). Default=FALSE.
# nsim: simulations when calculating p values (see "phylosig" function; "phytools" package). Default=999.#