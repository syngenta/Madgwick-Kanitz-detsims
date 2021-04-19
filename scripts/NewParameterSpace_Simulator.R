####################################################################################################################

#### Simulate new parameter space ####

# this is a general purpose simulator to be run on a server or PC
# the simulator is setup for two insecticides each with a resistance allele at separate loci, although other setups are possible
# the simulator provides 9 statistics per strategy per run relating to the time to resistance across allele frequency and population control
# the simulator runs one strategy at a time and for one combination of mitochondrial and nuclear loci for resistance alleles

# read in functions from package
library('detsims')

# packages and setup for parallelization
library(foreach)
library(doParallel)
registerDoParallel(2)  # use multicore, set to the number of our cores


####################################################################################################################

# make new parameter space
design = "random" # "factorial" or "random"
design_details = 1 # for factorial = number elements per set, random = power x of 10^x number of random samples

if(design=="factorial"){
  # parameter sets
  npset = 10^seq(npsr[1],npsr[2],length.out=design_details) # population size
  bkset = seq(bksr[1],bksr[2],length.out=design_details) # population growth rate
  gkset = seq(gksr[1],gksr[2],length.out=design_details) # adult death rate
  pkset = seq(pksr[1],pksr[2],length.out=design_details) # proportion of females exposed to insecticide
  ekset = seq(eksr[1],eksr[2],length.out=design_details) # proportion of males that have female exposure
  fAset = 10^-seq(fAsr[1],fAsr[2],length.out=design_details) # initial frequency of allele A
  fBset = 10^-seq(fBsr[1],fBsr[2],length.out=design_details) # initial frequency of allele B
  t1set = seq(t1sr[1],t1sr[2],length.out=design_details) # effectiveness of insecticide 1
  t2set = seq(t2sr[1],t2sr[2],length.out=design_details) # effectiveness of insecticide 2
  rtAset = seq(rtAsr[1],rtAsr[2],length.out=design_details) # resistance of allele A
  rtBset = seq(rtBsr[1],rtBsr[2],length.out=design_details) # resistance of allele B
  hrAset = seq(hrAsr[1],hrAsr[2],length.out=design_details) # dominance of resistance of allele A
  hrBset = seq(hrBsr[1],hrBsr[2],length.out=design_details) # dominance of resistance of allele B
  crAset = 10^-seq(crAsr[1],crAsr[2],length.out=design_details) # cost of resistance
  crBset = 10^-seq(crBsr[1],crBsr[2],length.out=design_details) # cost of resistance
  hcAset = seq(hcAsr[1],hcAsr[2],length.out=design_details) # dominance of resistance cost of allele A
  hcBset = seq(hcBsr[1],hcBsr[2],length.out=design_details) # dominance of resistance cost of allele B
  # the number of combinations
  var_comb = length(npset)*length(bkset)*length(gkset)*length(pkset)*length(ekset)*length(fAset)*length(fBset)*length(t1set)*length(t2set)*length(rtAset)*length(rtBset)*length(hrAset)*length(hrBset)*length(crAset)*length(crBset)*length(hcAset)*length(hcBset)
  # create matrix to store input combinations
  simdata = expand.grid(npset,bkset,gkset,pkset,ekset,fAset,t1set,rtAset,hrAset,crAset,hcAset,fBset,t2set,rtBset,hrBset,crBset,hcBset)
}

if(design=="random"){
  # the number of parameter combinations to sample
  var_comb = 10^design_details
  # create empty matrix to store parameter combinations
  simdata = matrix(0,var_comb,17)
  # populate simdata
  simdata[,1] = 10^runif(var_comb,npsr[1],npsr[2])
  simdata[,2] = rlnorm(var_comb,0,1) # runif(var_comb,bksr[1],bksr[2])
  simdata[,3] = runif(var_comb,gksr[1],gksr[2])
  simdata[,4] = runif(var_comb,pksr[1],pksr[2])
  simdata[,5] = runif(var_comb,eksr[1],eksr[2])
  simdata[,7] = runif(var_comb,t1sr[1],t1sr[2])
  simdata[,8] = runif(var_comb,rtAsr[1],rtAsr[2])
  simdata[,9] = runif(var_comb,hrAsr[1],hrAsr[2])
  simdata[,10] = 10^-runif(var_comb,crAsr[1],crAsr[2])
  simdata[,11] = runif(var_comb,hcAsr[1],hcAsr[2])
  simdata[,13] = runif(var_comb,t2sr[1],t2sr[2])
  simdata[,14] = runif(var_comb,rtBsr[1],rtBsr[2])
  simdata[,15] = runif(var_comb,hrBsr[1],hrBsr[2])
  simdata[,16] = 10^-runif(var_comb,crBsr[1],crBsr[2])
  simdata[,17] = runif(var_comb,hcBsr[1],hcBsr[2])
  for(i in 1:var_comb){ # frequencies are bounded between 1% and 1/N
    simdata[i,6] = 10^-runif(1,2,log10(simdata[i,1]))
    simdata[i,12] = 10^-runif(1,2,log10(simdata[i,1]))
  }
  # as data frame
  simdata = as.data.frame(simdata)
}

# provide column names for simdata
colnames(simdata) = c("Population Size","Intrinsic Birth Rate","Intrinsic Death Rate","Female Exposure","Male Exposure",
                      "Initial Frequency A","Effectiveness 1","Resistance Restoration A","Dominance of Resistance Restoration A","Resistance Cost A","Dominance of Resistance Cost A",
                      "Initial Frequency B","Effectiveness 2","Resistance Restoration B","Dominance of Resistance Restoration B","Resistance Cost B","Dominance of Resistance Cost B")

# save parameter space
write.csv(simdata,paste("PSData_",design,design_details,".csv",sep=""),row.names=F)


####################################################################################################################


# run all simulations

# strategies/scenarios to simulate
adaptive_inputs = rep(0,21)
inherit1_inputs = rep(c("m","m","n"),7)
inherit2_inputs = rep(c("m","n","n"),7)
revorder_inputs = c(rep("",15),rep("rev",3),rep("",3))
strategy_inputs = c(rep(1:5,each=3),rep(5,3),rep(4,3))
mixmaxps_inputs = c(rep(F,18),rep(T,3))

# generate output datafiles by simulation
for(l in 1:21){
  detsims_simulator(simdata,revorder_inputs[l],c(inherit1_inputs[l],inherit2_inputs[l]),mixmaxps_inputs[l],adaptive_inputs[l],strategy_inputs[l])
}


# END
