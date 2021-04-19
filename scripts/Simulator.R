####################################################################################################################

#### Strategy Comparison ####

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

#### Manual part of setup of parameter space for strategy comparison ####

# meta-parameters
loadps = "make" # "load" or "make" # NB: even if load must get manual entries correct
design = "random" # "factorial" or "random"
design_details = 2 # for factorial = number elements per set, random = power x of 10^x number of random samples
inheritance = c("n","n") # nn, mn, nm, mm
u = 1 # strategy number, picking element of strategyset: "SoloA","SoloB","Combination","Mixture","Rotation","Sequence"
v = 0 # strategy type, either: 0= as is and 1= adaptive (where solo for second insecticide after first has broken)
revin = "rev" # "rev" or "" for reverse inheritance or do not
mixmaxps = F # if true sets mixtures to maximum for a positive control, false else

# parameter set ranges as c(min,max)
npsr = c(2,9) # population size, NB: log scale
bksr = c(0,2) # population growth rate, NB: log-normal with mean=0 and sd=1 if random
gksr = c(0,1) # adult death rate
pksr = c(0,1) # proportion of females exposed to insecticide
eksr = c(0,1) # proportion of males that have female exposure
fAsr = c(2,9) # initial frequency of allele A, NB: log scale
fBsr = c(2,9) # initial frequency of allele B, NB: log scale
t1sr = c(0,1) # effectiveness of insecticide 1
t2sr = c(0,1) # effectiveness of insecticide 2
rtAsr = c(0,1) # resistance of allele A
rtBsr = c(0,1) # resistance of allele B
hrAsr = c(0,1) # dominance of resistance of allele A
hrBsr = c(0,1) # dominance of resistance of allele B
crAsr = c(0,1) # resistance cost of allele A
crBsr = c(0,1) # resistance cost of allele B
hcAsr = c(0,1) # dominance of resistance cost of allele A
hcBsr = c(0,1) # dominance of resistance cost of allele B
# simulation constants
strategy_start = 1 # the generation number when the strategy starts
gsf = strategy_start # here, a redundant function that allows for multiple different start times
gen = 500 # run simulation for max of gen generations
endsim = 0.5 # frequency threshold of resistance for the simulation to end
dk = 1000000  # constant specifying the decay in the toxicity of the insecticide, as the half-life of the treatment - if very large then decay is effectively ignored
rk = 36  # constant specifying after how many generations a bednet is to be replaced / retreated


#### Automated part of setup of parameter space for strategy comparison ####

# strategy set of names for the run
strategyset = c("SoloA","SoloB","Mosaic","Mixture","Rotation","Sequence") # NB: differs from strategy within the simulation code
strategytype = c("Fixed","Adaptive") # whether switch away from insecticide when 50% allele frequency is surpassed
vv = c(1,1) # initialize that adaptation of a strategy has not yet occured
# decode strategy number into strategy within simulation
if(u==1){
  strategy = "Combination" # Combination, Mixture, Rotation
  strategy_details = c(1,0) # the percentage of bed-nets that singly use an insecticide, numbers in order of insecticides_applied rows
}
if(u==2){
  strategy = "Combination" # Combination, Mixture, Rotation
  strategy_details = c(0,1) # the percentage of bed-nets that singly use an insecticide, numbers in order of insecticides_applied rows
}
if(u==3){
  strategy = "Combination" # Combination, Mixture, Rotation
  strategy_details = c(0.5,0.5) # the percentage of bed-nets that singly use an insecticide, numbers in order of insecticides_applied rows
}
if(u==4){
  strategy = "Mixture" # Combination, Mixture, Rotation
  strategy_details = c(0.5,0.5)
}
if(u==5){
  strategy = "Rotation" # Combination, Mixture, Rotation
  strategy_details = c(rk,rk) # the replacement schedule for a bed-net
}
if(u==6){
  strategy = "Sequential" # Combination, Mixture, Rotation
  strategy_details = c(endsim,endsim) # resistance threshold for switching strategy
}

if(loadps=="load"){
  simdata = read.csv("PSData_random2.csv")
  simdata = as.data.frame(simdata)
  var_comb = nrow(simdata)
}

if(loadps=="make"){
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
}

# switch order of A and B so that B goes first
if(revin=="rev"){
  inheritance = rev(inheritance)
  simdata = cbind(simdata[,1:5],simdata[,12:17],simdata[,6:11])
}

# remove dominance terms if for mitochondrial
if(inheritance[1]=="m"){
  simdata[,9] = 0
  simdata[,11] = 0
}
if(inheritance[2]=="m"){
  simdata[,15] = 0
  simdata[,17] = 0
}

# create matrix to store output data - NB: only save data from the control population of females
outdata = matrix(0,var_comb,10)
colnames(outdata) = c("A_t_50%","A_f_250","A_f_bar","B_t_50%","B_f_250","B_f_bar","nf_80%","nf_250","nf_bar","nf_ext")


####################################################################################################################

#### Run through simulations under variable parameter space ####

# start run loop
output = foreach(np=simdata[,1], bk=simdata[,2], gk=simdata[,3], pk=simdata[,4], ek=simdata[,5],
                 ifA=simdata[,6], t1=simdata[,7], rA=simdata[,8], hrA=simdata[,9], cA=simdata[,10], hcA=simdata[,11],
                 ifB=simdata[,12], t2=simdata[,13], rB=simdata[,14], hrB=simdata[,15], cB=simdata[,16], hcB=simdata[,17],.packages='detsims') %dopar% {

                   ## setup

                   sf = c(ifA,ifB) # starting frequencies

                   # insecticide matrix
                   #   Insecticide_Name :     use numbers 1-inf
                   #   Toxicity :             constant specifying the percentage toxicity to susceptible mosquitoes, within range 0 to 1
                   #   Assortment :           constant specifying the percentage negative assortment / repellency of the insecticide to susceptible mosquitoes, within range 0 to 1
                   #   Assortment_Exposure :  constant specifying the percentage exposure of repelled mosquitoes to an insecticide
                   insecticides_applied = matrix(0,2,4)
                   colnames(insecticides_applied) = c("Insecticide_Name","Toxicity","Assortment","Assorted_Exposure")
                   insecticides_applied[1,] = c(1,t1,0,0)
                   insecticides_applied[2,] = c(2,t2,0,0)
                   inumber = nrow(insecticides_applied) # short-hand definition, number of resistant alleles

                   # strategy setup
                   in_use = rep(0,inumber) # which insecticides are currently in use, start with none in generation 0
                   ins_dec = rep(0,inumber+1) # which insecticides are currently in decay, start with none in generation 0

                   # resistance matrix
                   #   Allele_Name :           use letters A-Z
                   #   Ploidy :                1 is haploid and 2 is diploid (no other values accepted)
                   #   Insecticide :           use insecticide name to describe what the resistance allele is resistant to, and use 0 if applies to all insecticides
                   #   Toxicity :              constant specifying the percentage reduction in toxicity of resistant mosquitoes, within range 0 to 1
                   #   Toxicity_Dominance :    constant specifying the dominance of the toxicity effect, within range 0 to 1
                   #   Assortment :            constant specifying the percentage increase in negative assortment / repellency of resistant mosquitoes, within range 0 to 1
                   #   Assortment_Dominance :  constant specifying the dominance of the assortment effect, within range 0 to 1
                   #   Cost :                  constant specifying the cost of resistance, within range 0 to 1
                   #   Cost_Dominance :        constant specifying the dominance of the cost of resistance effect, within range 0 to 1
                   resistant_alleles = matrix(0,2,9)
                   colnames(resistant_alleles) = c("Allele_Name","Ploidy","Insecticide","Toxicity","Toxicity_Dominance","Assortment","Assortment_Dominance","Cost","Cost_Dominance")
                   resistant_alleles[1,] = c("A",ifelse(inheritance[1]=="n",2,1),1,rA,hrA,0,0,cA,hcA) # haploid, decrease toxicity
                   resistant_alleles[2,] = c("B",ifelse(inheritance[2]=="n",2,1),2,rB,hrB,0,0,cB,hcB) # haploid, decrease toxicity
                   anumber = nrow(resistant_alleles) # short-hand definition, number of resistant alleles
                   anames = resistant_alleles[,1] # short-hand definition, names of resistant alleles
                   ploidy = as.numeric(resistant_alleles[,2]) # short-hand definition, ploidy of resistant loci

                   # with linkage disequilibrium: initialize allele frequency store
                   savedata = matrix(0,gen+1,2*anumber+2*2) # record resistant alleles only
                   colnames(savedata) = c(paste(resistant_alleles[,1],"_F",sep=""), paste(resistant_alleles[,1],"_M",sep=""), "Pop_Size_F", "Mean_Fitness_F", "Pop_Size_M", "Mean_Fitness_M")
                   savedata[1,] = c(rep(0,2*anumber),rep(c(0.5*np,np),times=2)) # starting condition
                   t50A = 0
                   t50B = 0
                   c80 = 0
                   ex1 = 0


                   ## start generations loop

                   for(t in 1:gen){


                     ## Generate population numbers, fitness, allele frequency and genotype frequency

                     # if first generation, then initialize fitness, allele frequency and genotype frequency
                     if(t==1){

                       # generate allele frequency
                       afset = rep(0,anumber) # all resistant alleles start at zero frequency when t=0
                       names(afset) = anames # transfer names
                       amset = rep(0,anumber) # all resistant alleles start at zero frequency when t=0
                       names(amset) = anames # transfer names

                       # generate genotype frequencies of males and females
                       ffset = initialize_freqset(afset,anames,ploidy)
                       fmset = initialize_freqset(amset,anames,ploidy)

                       # generate fitness of females and males assuming no insecticides
                       wfset = initialize_femalefitness(pk, 0, in_use, matrix(0,inumber,4), resistant_alleles, strategy, strategy_details, ins_dec, mixmaxps)
                       wmset = initialize_malefitness(wfset, 1, resistant_alleles)

                     }

                     # update in_use with the strategy for this generation
                     iu0 = in_use # in_use in previous generation
                     if(sum(vv)==2){in_use = update_strategy(t, strategy, strategy_start, strategy_details, afset, in_use, rk) } # update in use by strategy
                     if(v==1){if(sum(vv)==2){ # coerce so if adaptation of strategy has occured then consistent with adaptation not strategy
                       if(vv[1]==1&t50A>0){ # switch off if fA>50%
                         in_use[1] = 0
                         in_use[2] = 1
                         vv[1] = 0 # and record that switch insecticide 1 off
                       }
                       if(vv[2]==1&t50B>0){ # switch off if fB>50%
                         in_use[1] = 1
                         in_use[2] = 0
                         vv[2] = 0 # and record that switch insecticide 2 off
                       }
                     }}
                     ins_dec = insecticide_decay(t, iu0, in_use, dk, rk, ins_dec)

                     # if the first resistant allele enters this t then initialize fitness, allele frequency and genotype frequency
                     if(t==min(gsf)){

                       # generate allele frequency
                       if(sum(gsf==t)>0){
                         afset[gsf==t] = sf[gsf==t] # if a resistant allele enters at t=1 then input its starting frequency
                       }
                       if(sum(gsf==t)>0){
                         amset[gsf==t] = sf[gsf==t] # if a resistant allele enters at t=1 then input its starting frequency
                       }

                       # generate genotype frequencies
                       ffset = initialize_freqset(afset,anames,ploidy)
                       fmset = initialize_freqset(amset,anames,ploidy)

                       # save the mutation to the previous generation
                       savedata[t,seq(1,anumber)[gsf==t]] = sf[gsf==t]
                       savedata[t,seq(1+anumber,2*anumber)[gsf==t]] = sf[gsf==t]

                     }

                     # if a resistant allele enters this t but another has already exist then reinitialize fitness, allele frequency and genotype frequency
                     if(t>min(gsf)&&sum(gsf==t)>0){

                       # alleles that have already been added
                       aset_added = (gsf<(t-1))*0.5

                       # allele frequency vector including newly added alleles
                       aset_add = aset
                       aset_add[gsf==t] = sf[gsf==t]

                       # the element numbers of the loci to be added
                       add_loci = seq(1,anumber)[gsf==t]

                       # ... reinitialize fset with new alleles (from mutation)
                       ffset = reinitialize_freqset(aset_added,ffset,aset_add,anames,ploidy,add_loci)
                       fmset = reinitialize_freqset(aset_added,fmset,aset_add,anames,ploidy,add_loci)

                       # save the mutation to the previous generation
                       savedata[t,seq(1,anumber)[gsf==t]] = sf[gsf==t]
                       savedata[t,seq(1+anumber,2*anumber)[gsf==t]] = sf[gsf==t]

                     }

                     # generate fitness of all allele combinations under insecticide
                     if(t>=strategy_start){
                       wfset = initialize_femalefitness(pk, 0, in_use, insecticides_applied, resistant_alleles, strategy, strategy_details, ins_dec, mixmaxps)
                       wmset = initialize_malefitness(wfset, ek, resistant_alleles)
                     }

                     # number of females and males before selection
                     n0f = savedata[t,1+2*anumber] # female population size of previous generation
                     n0m = savedata[t,1+2*(1+anumber)] # female population size of previous generation


                     ## Selection and reproduction

                     # calculate frequencies post-selection (i.e. after insecticide mortality) of females and males
                     ffset_ps = ffset*wfset/sum(ffset*wfset)
                     fmset_ps = fmset*wmset/sum(fmset*wmset)

                     # calculate frequency of offspring genotypes
                     f1set = random_mating(ffset_ps, fmset_ps, resistant_alleles)

                     # calculate population
                     n1f = sum(n0f*ffset*wfset*(1-gk) + sum(n0f*ffset*wfset)*f1set*(1+bk*(1-(sum(n0f*ffset*wfset)/((bk*np)/(2*(1+bk-gk)))))))
                     n1m = sum(n0m*fmset*wmset*(1-gk) + sum(n0f*ffset*wfset)*f1set*(1+bk*(1-(sum(n0f*ffset*wfset)/((bk*np)/(2*(1+bk-gk)))))))

                     # calculate frequencies of males and females given overlapping generations - assume new generation has equal sex ratio
                     f1fset = (n0f*ffset*wfset*(1-gk) + sum(n0f*ffset*wfset)*f1set*(1+bk*(1-(sum(n0f*ffset*wfset)/((bk*np)/(2*(1+bk-gk)))))))/n1f
                     f1mset = (n0m*fmset*wmset*(1-gk) + sum(n0f*ffset*wfset)*f1set*(1+bk*(1-(sum(n0f*ffset*wfset)/((bk*np)/(2*(1+bk-gk)))))))/n1m

                     # calculate allele frequency in offspring genotypes
                     # females
                     gv = genotypic_values(anames,ploidy)
                     a1fset = rep(0,anumber)
                     for(j in 1:anumber){
                       a1fset[j] = sum(f1fset*gv[,j])
                     }
                     names(a1fset) = anames # transfer names
                     # males
                     a1mset = rep(0,anumber)
                     for(j in 1:anumber){
                       a1mset[j] = sum(f1mset*gv[,j])
                     }
                     names(a1mset) = anames # transfer names


                     ## Save and update

                     # save data
                     savedata[1+t,] = c(a1fset, a1mset, n1f, sum(ffset*wfset), n1m, sum(fmset*wmset))

                     # update matrices
                     afset = a1fset
                     amset = a1mset
                     ffset = f1fset
                     fmset = f1mset

                     # record when alleles have exceeded a resistance allele frequency threshold
                     if(t50A==0 & afset[1]>endsim){
                       t50A = t
                     }
                     if(t50B==0 & afset[2]>endsim){
                       t50B = t
                     }
                     if(c80==0 & n1f>(0.4*np)){
                       c80 = t
                     }
                     if(ex1==0 & n1f<1){
                       ex1 = t
                       break
                     }

                   }

                   # output : "A_t_50%","A_f_250","A_f_bar","B_t_50%","B_f_250","B_f_bar","n_80%","n_250","n_bar"
                   return(c(t50A, savedata[251,1], mean(savedata[1:251,1]),
                            t50B, savedata[251,2], mean(savedata[1:251,2]),
                            c80, savedata[251,5], mean(savedata[1:251,5]),ex1))

                 }

# reformat output
outdata[1:var_comb,1:10] = t(as.data.frame(output))

# save outdata
write.csv(outdata,paste("Data_",design,design_details,"_",ifelse(mixmaxps==T,"maxps",""),strategyset[u],"_",strategytype[v+1],"_",revin,paste(inheritance,collapse=""),".csv",sep=""),row.names=F)


# END
