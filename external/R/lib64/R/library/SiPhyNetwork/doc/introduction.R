## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.show='hold',setup,out.width='\\textwidth', fig.height = 4------------
library(ape)
library(SiPhyNetwork)
set.seed(82589933) ##set the seed for reproducibility. This is the exponent of the largest known Mersenne prime number

##First we need a function that describes how inheritance probabilities are drawn
inheritance.fxn <- make.beta.draw(10,10)
##We can see that this function makes draws from a beta(10,10) distribution
inheritance.fxn() 
inheritance.fxn()

##We also want to set the proportion of each type of hybrid event
hybrid_proportions <-c(0.5,  ##Lineage Generative
                       0.25, ##Lineage Degenerative
                       0.25) ##Lineage Neutral

##We can simulate to 7 extant tips with the SSA
ssa_nets<-sim.bdh.taxa.ssa(n=7,numbsim=20,
                    lambda=1,mu=0.2,
                    nu=0.20, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn)
ssa_net<-ssa_nets[[20]] ##The sim.bdh functions return a list of length numbsim. We get the 20th simulation
print(ssa_net)

##We can also simulate 7 extant taxa with the GSA. 
##We choose m=30 because it becomes very unlikely that at 30 tips we will ever return to 7
gsa_nets<-sim.bdh.taxa.gsa(m=30,n=7,numbsim=20,
                    lambda=1,mu=0.6,
                    nu=0.3, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn)
gsa_net<-gsa_nets[[19]]

##Simulate a network up to age 2
age_nets <-sim.bdh.age(age=2,numbsim=20,
                    lambda=1,mu=0.2,
                    nu=0.25, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn)
age_net<-age_nets[[8]]

## -----------------------------------------------------------------------------
age_net$inheritance ##This corresponds to the edges found in reticulation
age_net$reticulation

## -----------------------------------------------------------------------------
##We can simulate to 30 extant tips under the SSA. In this case the 30 acts as the m parameter of the GSA
ssa_nets<-sim.bdh.taxa.ssa(n=30,numbsim=10,
                    lambda=1,mu=0.2,
                    nu=0.20, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn)

my_net<-ssa_nets[[1]]
my_net<-network.gsa(net=my_net,ntaxa=5,complete=T,frac=1,stochsampling=F)

## -----------------------------------------------------------------------------

##Equal chance of all three types
hybprops1 <-c(1/3, ##Lineage Generative
              1/3, ##Lineage Degenerative
              1/3) ##Lineage Neutral

##Skewed chance of all three types
hybprops2 <-c(0.5, ##Lineage Generative
              0.2, ##Lineage Degenerative
              0.3) ##Lineage Neutral

##Only Lineage Generative Hybridization occurs
hybprops3 <-c(1, ##Lineage Generative
              0, ##Lineage Degenerative
              0) ##Lineage Neutral


##simulate where all 3 are equally likely
age_nets <-sim.bdh.age(age=2,numbsim=20,
                    lambda=1,mu=0.2,
                    nu=0.25, hybprops = hybprops1,
                    hyb.inher.fxn = inheritance.fxn)



## ---- fig.show='hold',fig.height = 4------------------------------------------
plot(ssa_net,main="SSA Network")
plot(gsa_net,main="GSA Network") 
plot(age_net,main="Age Network")

## ---- fig.show='hold',fig.height = 4------------------------------------------

ssa_pnet <-plottable.net(ssa_net)
gsa_pnet <-plottable.net(gsa_net)
age_pnet <-plottable.net(age_net)

plot(ssa_pnet,main="SSA Network")
plot(gsa_pnet,main="GSA Network")
plot(age_pnet,main="Age Network")

## ---- fig.show='hold',fig.height = 4------------------------------------------
isTreeChild(gsa_net)
isTreeBased(gsa_net)
isFUstable(gsa_net)
getNetworkLevel(gsa_net)

## -----------------------------------------------------------------------------

write.evonet(ssa_net,file='') ##we can see that inheritance probabilities aren't included here
my_newick<-write.net(ssa_net,file = '') ## if we include a file name the network will print to file instead of print on the console
print(my_newick)

my_net<-read.net(text=my_newick)
print(my_net)
str(my_net)##we can see that my_net has the inheritance element for inheritance probabilities


## ----fig.show='hold',fig.height = 3.5,fig.align='center',fig.width=7----------

pruned_gsa <- reconstructedNetwork(gsa_net)
plot(plottable.net(pruned_gsa),main='Reconstructed Phylogeny')
pruned_gsa <- incompleteSampling(pruned_gsa,rho=5/7,stochastic = F)
plot(plottable.net(pruned_gsa),main='Reconstructed Phylogeny with Incomplete Sampling')



## -----------------------------------------------------------------------------
##Here are some of the make functions
f1<-make.exp.decay(t=1,s=1)
f2<-make.linear.decay(threshold = 1)
f3<-make.stepwise(probs = c(1,0.5,0),distances = c(0.25,0.75,Inf))
f4<-make.polynomial.decay(threshold = 1,degree = 2)

##We can use any of these functions as the hyb.rate.fxn argument in a sim.bdh function
age_nets <-sim.bdh.age(age=2,numbsim=10,
                    lambda=1,mu=0.2,
                    nu=0.25, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn,
                    hyb.rate.fxn = f3)

## ----echo=F,fig.width=7,fig.height=7------------------------------------------
{
old_pars <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
x<-seq(0,1.2,by=0.001)
plot(x,f1(x),xlim = c(0,1.2),ylim = c(0,1),type = 'l',main = 'f1 Exponential Decay',xlab = 'Genetic Distance',ylab = 'Pr(Successful Hybridization)')
plot(x,f2(x),xlim = c(0,1.2),ylim = c(0,1),type = 'l',main = 'f2 Linear Decay',xlab = 'Genetic Distance',ylab = 'Pr(Successful Hybridization)')

y<-rep(NA,length(x))
for(i in 1:length(x)){
  y[i]<-f3(x[i])
}
plot(x,y,xlim = c(0,1.2),ylim = c(0,1),type = 'l',main = 'f3 Stepwise Function',xlab = 'Genetic Distance',ylab = 'Pr(Successful Hybridization)')
plot(x,f4(x),xlim = c(0,1.2),ylim = c(0,1),type = 'l',main = 'f4 Polynomial Decay',xlab = 'Genetic Distance',ylab = 'Pr(Successful Hybridization)')
}

##reset graphical parameters
par(old_pars)

## -----------------------------------------------------------------------------


initial_val<-2 ## The root starts off at 2N

###function for what happens at hybridization event
hyb_e_fxn <- function(parent_states,inheritance){
  ##For allopolyploidy we add the ploidy of both parents
  return(sum(parent_states)) 
}

##Function for determining whether hybridization occurs
hyb_c_fxn <-function(parent_states,hybrid_state){
  ##Hybridization occurs only when the ploidy is the same
  return(parent_states[1]==parent_states[2])
}


##Function for how the trait changes over time
t_fxn <- function(trait_states,timestep){
  ##We assume that autopolyploidy occur exponentially with rate lambda
  lambda<- 2 ##Rate of autopolyploidy
  
  ##The number of autopolyploidy events that occur on each lineage over the timestep
  nevents<-rpois(length(trait_states),timestep) 
  
  ##each event doubles the ploidy
  new_states<- trait_states * (2^nevents) 
  return(new_states)
}

##Function for how the trait changes at speciation events
s_fxn <-function(tip_state){
  ##Ploidy doesn't change at speciation events. 
  ##Both daughter lineages have the same ploidy as the parent
  return(c(tip_state,tip_state))
}

trait_model<-make.trait.model(initial_states = initial_val,
                              hyb.event.fxn = hyb_e_fxn,
                              hyb.compatibility.fxn = hyb_c_fxn,
                              time.fxn = t_fxn,
                              spec.fxn = s_fxn)


trait_nets <-sim.bdh.age(age=2,numbsim=10,
                    lambda=1,mu=0.2,
                    nu=0.25, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn,
                    trait.model = trait_model)


## -----------------------------------------------------------------------------
initial_val<-0 ## The root starts off at 0

###function for what happens at hybridization event
hyb_e_fxn <- function(parent_states,inheritance){
  ## Take a weighted average of the two traits
  return(sum(parent_states*c(inheritance,1-inheritance))) 
}

##Function for determining whether hybridization occurs
hyb_c_fxn <-function(parent_states,hybrid_state){
  #Make linear decay function that decreases hybrid success probability linearly
  #Hybridization cannot occur when the traits differ by more than 4
  decay.fxn <- make.linear.decay(4)
  
  #Trait dissimilarity
  diss<-abs(parent_states[1]-parent_states[2])
  
  success_prob<-decay.fxn(diss) 
  return(runif(1,0,1)<=success_prob)
}


##Function for how the trait changes over time
t_fxn <- function(trait_states,timestep){
  ##We assume brownian motion on the continuous trait
  sigma<- 2 ##sqrt(Rate of evolution)
  
  ##Make brownian motion draws for each lineage
  delta_x <- rnorm(length(trait_states),mean=0,sd=sigma*sqrt(timestep))
  
  return(delta_x+trait_states)
}

##Function for how the trait changes at speciation events
s_fxn <-function(tip_state){
  ##No change at speciation
  return(c(tip_state,tip_state))
}

trait_model<-make.trait.model(initial_states = initial_val,
                              hyb.event.fxn = hyb_e_fxn,
                              hyb.compatibility.fxn = hyb_c_fxn,
                              time.fxn = t_fxn,
                              spec.fxn = s_fxn)


trait_nets <-sim.bdh.age(age=2,numbsim=10,
                    lambda=1,mu=0.2,
                    nu=0.25, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn,
                    trait.model = trait_model)


## -----------------------------------------------------------------------------
initial_val<-0 ## The root starts off at 0

###function for what happens at hybridization event
hyb_e_fxn <- function(parent_states,inheritance){
  ## Take a weighted average of the two traits
  hyb_val<-sum(parent_states*c(inheritance,1-inheritance))
  
  ##Now allow transgressive evolution by making a random draw from a normal distribution
  transgression <- rnorm(1,0,4)
  
  return(hyb_val+transgression) 
}

##Function for determining whether hybridization occurs
hyb_c_fxn <-function(parent_states,hybrid_state){
  ##Lets say that anything outside 10% of the hybrid lineage trait value is a different nich
  niche_bound <-hybrid_state*0.1
  
  lower_bound <-hybrid_state-niche_bound
  upper_bound <-hybrid_state+niche_bound
  
  ##Check if the parental lienages are within the niche boundary
  if(all( (parent_states < lower_bound) | (parent_states > upper_bound) )){
    ##Parent lineages are outside the hybrid niche. Hybrid survives
    return(TRUE)
  }else{
    ##Parent lineage inside the hybrid niche. Hybrid breakdown
    return(FALSE)
  }
  
}


##Function for how the trait changes over time
t_fxn <- function(trait_states,timestep){
  ##We assume brownian motion on the continuous trait
  sigma<- 2 ##sqrt(Rate of evolution)
  
  ##Make brownian motion draws for each lineage
  delta_x <- rnorm(length(trait_states),mean=0,sd=sigma*sqrt(timestep))
  
  return(delta_x+trait_states)
}

##Function for how the trait changes at speciation events
s_fxn <-function(tip_state){
  ##No change at speciation
  return(c(tip_state,tip_state))
}

trait_model<-make.trait.model(initial_states = initial_val,
                              hyb.event.fxn = hyb_e_fxn,
                              hyb.compatibility.fxn = hyb_c_fxn,
                              time.fxn = t_fxn,
                              spec.fxn = s_fxn)


trait_nets <-sim.bdh.age(age=2,numbsim=10,
                    lambda=1,mu=0.2,
                    nu=0.25, hybprops = hybrid_proportions,
                    hyb.inher.fxn = inheritance.fxn,
                    trait.model = trait_model)


