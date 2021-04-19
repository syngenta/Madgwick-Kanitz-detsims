####################################################################################################################

#### Strategy Comparison ####

# this is a general purpose simulator to be run on a server or PC
# the simulator is setup for two insecticides each with a resistance allele at separate loci, although other setups are possible
# the simulator provides 9 statistics per strategy per run relating to the time to resistance across allele frequency and population control
# the simulator runs one strategy at a time and for one combination of mitochondrial and nuclear loci for resistance alleles

# read in functions from package
library('detsims')


####################################################################################################################

#### Manual part of setup of parameter space for strategy comparison ####

# meta-parameters
inheritance = c("m","n") # nn, mn, nm, mm
u = 4 # strategy number, picking element of strategyset: "SoloA","SoloB","Combination","Mixture","Rotation","Sequence"
v = 1 # strategy type, either: 0= as is and 1= adaptive (where solo for second insecticide after first has broken)
mixmaxps = F # if true sets mixtures to maximum for a positive control, false else

# parameter set
np = 10^9 # population size, NB: log scale
bk = 2 # population growth rate
gk = 1 # adult death rate
pk = 0.6 # proportion of females exposed to insecticide
ek = 0 # proportion of males that have female exposure
ifA = 10^-9 # initial frequency of allele A, NB: log scale
ifB = 10^-6 # initial frequency of allele B, NB: log scale
t1 = 0.9 # effectiveness of insecticide 1
t2 = 0.9 # effectiveness of insecticide 2
rA = 1 # resistance of allele A
rB = 1 # resistance of allele B
hrA = 0.5 # dominance of resistance of allele A
hrB = 0.5 # dominance of resistance of allele B
cA = 0.05 # resistance cost of allele A
cB = 0.05 # resistance cost of allele B
hcA = 0.5 # dominance of resistance cost of allele A
hcB = 0.5 # dominance of resistance cost of allele B
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

# remove dominance terms if for mitochondrial
if(inheritance[1]=="m"){
  hrA = 0
  hcA = 0
}
if(inheritance[2]=="m"){
  hrB = 0
  hcB = 0
}

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


####################################################################################################################

#### Run through simulations under variable parameter space ####

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
  if(sum(vv)==2){in_use = update_strategy(t, strategy, strategy_start, strategy_details, afset, in_use, rk) } # update in use by strategy   # TODO: This function is not declared anywhere.
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
  ins_dec = insecticide_decay(t, iu0, in_use, dk, rk, ins_dec)  # TODO: This function is not declared anywhere.

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

}

# output : "A_t_50%","A_f_250","A_f_bar","B_t_50%","B_f_250","B_f_bar","n_80%","n_250","n_bar"
outdata = matrix(c(t50A, savedata[251,1], mean(savedata[1:251,1]), t50B, savedata[251,2], mean(savedata[1:251,2]), c80, savedata[251,5], mean(savedata[1:251,5])),1,9)
colnames(outdata) = c("A_t_50%","A_f_250","A_f_bar","B_t_50%","B_f_250","B_f_bar","nf_80%","nf_250","nf_bar")
print(outdata)


####################################################################################################################

#### graphics ####

# ploidy assembly
ploidy_names = ploidy
ploidy_names[ploidy==1] = "Hap."
ploidy_names[ploidy==2] = "Dip."

# population size transformation
ps = savedata[1:(gen+1),1+2*anumber] + savedata[1:(gen+1),1+2*anumber]
psmax = max(ps)
psmin = min(ps)
logpsmax = ceiling(max(log10(ps)))
logpsmin = floor(min(log10(ps)))
logps = log10(ps)
tlogps = (logps-logpsmin)/(logpsmax-logpsmin)

# allele frequency transformation
logafmax = ceiling(max(log10(savedata[max(gsf):(gen+1),1:anumber]+savedata[max(gsf):(gen+1),(1+anumber):(2*anumber)])))
logafmin = floor(min(log10(savedata[max(gsf):(gen+1),1:anumber]+savedata[max(gsf):(gen+1),(1+anumber):(2*anumber)])))

# plot empty figure
par(mar=c(5,5,5,5))
plot(-1,-1,xlim=c(0,gen),ylim=c(-0.1,1.1),xlab="Time in Generations",ylab="",font.lab=2,las=1,xaxs="i",yaxs="i",yaxt="n")
# population size
points(seq(0,gen),tlogps,lwd=6,col="grey",type="l")
axis(4,at=seq(0,1,length.out=6),labels=seq(logpsmin,logpsmax,length.out=6),col="grey",col.axis="grey",las=1)
mtext("Log Population Size",side=4,col="grey",line=3,font=2)
# one individual population
points(c(0,gen),rep(1-(log10(np)/-logpsmin),2),type="l",lwd=3,lty=2,col="grey")
legend(0.85*gen,1-log10(np)/-logpsmin,"Min. Pop.",text.col="grey",bty="n",adj=0.3)

# allele frequency
palette(rainbow(anumber))
for(i in 1:anumber){
  points(seq(max(gsf)-1,gen),(log10(savedata[max(gsf):(gen+1),i]+savedata[max(gsf):(gen+1),i+anumber])-logafmin)/(logafmax-logafmin),type="l",lwd=4,col=i)
}
axis(2,at=seq(0,1,length.out=6),labels=seq(logafmin,logafmax,length.out=6),las=1)
mtext("Log Allele Frequency",side=2,line=3,font=2)

# legend and title
legend(0.6*gen,0.9,c("Pop. Size",paste(resistant_alleles[,1],resistant_alleles[,3],ploidy_names)),lwd=6,col=c("grey",rainbow(anumber)),bty="n")
title(paste("Strategy=",strategy,"; b=",bk,"; e=",ek,"; g=",gk,"; p=",pk,sep=""))
box()



# END
