#' A general purpose simulator function that takes a parameter space input and returns summary statistc outputs
#'
#'  @param simdata_input A matrix of parameter space, where columns are parameters and rows are random samples
#'  @param revin_input A string to determine whether the order that insecticide A and B are used is reversed (for sequences/rotations)
#'  @param inheritance_input A vector of the inheritance pattern of locus 1 and locus 2 (either "n"=nuclear or "m"=mitochondrial)
#'  @param mixmaxps_input A boolean of whether mixtures should be interpreted as k=variable (F) or k=1 (T)
#'  @param v_input A value of whether adaptive strategies are employed (yes=1; no=0)
#'  @param u_input A value describing which strategy should be used (between 1 and 6; see strategyset)
#'
#'  @return A vector of the genotype frequencies
#'  @export
#'  @examples
detsims_simulator <- function(simdata_input,revin_input,inheritance_input,mixmaxps_input,v_input,u_input){

  # simulation constants
  strategy_start = 1 # the generation number when the strategy starts
  gsf = strategy_start # here, a redundant function that allows for multiple different start times
  gen = 500 # run simulation for max of gen generations
  endsim = 0.5 # frequency threshold of resistance for the simulation to end
  dk = 1000000  # constant specifying the decay in the toxicity of the insecticide, as the half-life of the treatment - if very large then decay is effectively ignored
  rk = 36  # constant specifying after how many generations a bednet is to be replaced / retreated

  # strategy set of names for the run
  strategyset = c("SoloA","SoloB","Mosaic","Mixture","Rotation","Sequence") # NB: differs from strategy within the simulation code
  strategytype = c("Fixed","Adaptive") # whether switch away from insecticide when 50% allele frequency is surpassed
  vv = c(1,1) # initialize that adaptation of a strategy has not yet occured
  # decode strategy number into strategy within simulation
  if(u_input==1){
    strategy = "Combination" # Combination, Mixture, Rotation
    strategy_details = c(1,0) # the percentage of bed-nets that singly use an insecticide, numbers in order of insecticides_applied rows
  }
  if(u_input==2){
    strategy = "Combination" # Combination, Mixture, Rotation
    strategy_details = c(0,1) # the percentage of bed-nets that singly use an insecticide, numbers in order of insecticides_applied rows
  }
  if(u_input==3){
    strategy = "Combination" # Combination, Mixture, Rotation
    strategy_details = c(0.5,0.5) # the percentage of bed-nets that singly use an insecticide, numbers in order of insecticides_applied rows
  }
  if(u_input==4){
    strategy = "Mixture" # Combination, Mixture, Rotation
    strategy_details = c(0.5,0.5)
  }
  if(u_input==5){
    strategy = "Rotation" # Combination, Mixture, Rotation
    strategy_details = c(rk,rk) # the replacement schedule for a bed-net
  }
  if(u_input==6){
    strategy = "Sequential" # Combination, Mixture, Rotation
    strategy_details = c(endsim,endsim) # resistance threshold for switching strategy
  }

  # switch order of A and B so that B goes first
  if(revin_input=="rev"){
    inheritance_input = rev(inheritance_input)
    simdata_input = cbind(simdata_input[,1:5],simdata_input[,12:17],simdata_input[,6:11])
  }

  # remove dominance terms if for mitochondrial
  if(inheritance_input[1]=="m"){
    simdata_input[,9] = 0
    simdata_input[,11] = 0
  }
  if(inheritance_input[2]=="m"){
    simdata_input[,15] = 0
    simdata_input[,17] = 0
  }

  # create matrix to store output data - NB: only save data from the control population of females
  outdata = matrix(0,var_comb,10)
  colnames(outdata) = c("A_t_50%","A_f_250","A_f_bar","B_t_50%","B_f_250","B_f_bar","nf_80%","nf_250","nf_bar","nf_ext")


  ####################################################################################################################

  #### Run through simulations under variable parameter space ####

  # start run loop
  output = foreach(np=simdata_input[,1], bk=simdata_input[,2], gk=simdata_input[,3], pk=simdata_input[,4], ek=simdata_input[,5],
                   ifA=simdata_input[,6], t1=simdata_input[,7], rA=simdata_input[,8], hrA=simdata_input[,9], cA=simdata_input[,10], hcA=simdata_input[,11],
                   ifB=simdata_input[,12], t2=simdata_input[,13], rB=simdata_input[,14], hrB=simdata_input[,15], cB=simdata_input[,16], hcB=simdata_input[,17],.packages='detsims') %dopar% {

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
                     resistant_alleles[1,] = c("A",ifelse(inheritance_input[1]=="n",2,1),1,rA,hrA,0,0,cA,hcA) # haploid, decrease toxicity
                     resistant_alleles[2,] = c("B",ifelse(inheritance_input[2]=="n",2,1),2,rB,hrB,0,0,cB,hcB) # haploid, decrease toxicity
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
                         wfset = initialize_femalefitness(pk, 0, in_use, matrix(0,inumber,4), resistant_alleles, strategy, strategy_details, ins_dec, mixmaxps_input)
                         wmset = initialize_malefitness(wfset, 1, resistant_alleles)

                       }

                       # update in_use with the strategy for this generation
                       iu0 = in_use # in_use in previous generation
                       if(sum(vv)==2){in_use = update_strategy(t, strategy, strategy_start, strategy_details, afset, in_use, rk) } # update in use by strategy
                       if(v_input==1){if(sum(vv)==2){ # coerce so if adaptation of strategy has occured then consistent with adaptation not strategy
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
                         wfset = initialize_femalefitness(pk, 0, in_use, insecticides_applied, resistant_alleles, strategy, strategy_details, ins_dec, mixmaxps_input)
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
  write.csv(outdata,paste("Data_",ifelse(mixmaxps_input==T,"maxps",""),strategyset[u_input],"_",strategytype[v_input+1],"_",revin_input,paste(inheritance_input,collapse=""),".csv",sep=""),row.names=F)


}

# END
