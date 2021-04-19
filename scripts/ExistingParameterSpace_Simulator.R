####################################################################################################################

#### Simulate pre-existing parameter space ####

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

# use existing parameter space

# read in data from file
simdata = read.csv("PSData_random6.csv")
simdata = as.data.frame(simdata)
var_comb = nrow(simdata)


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
