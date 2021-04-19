#' Initialize female fitness of different genotypes
#'
#' Initialize fitness of bi-allele multi-locus genotypes of females from underlying selection coefficients frequencies
#'
#'  @param pk The constant for the percentage of the human population protected by the insecticide
#'  @param tk The constant for the toxicity cost of the insecticide e.g. mosquitio mortality
#'  @param ak The constant for the repellency of the insecticide as a (negative) assortment factor
#'  @param resistant_alleles A matrix that encodes the properties of the resistant alleles, as one allele per row
#'  The resistant alleles matrix takes the form of:
#'  name (i.e. a letter), 1=haploid or 2=diploid, T=toxicity or A=assortment, effect, dominance of effect, cost, dominance of cost
#'
#'  @return A vector of the fitness of different genotypes
#'  @export
#'  @examples
initialize_femalefitness <- function(pk, zk, in_use, insecticides_applied, resistant_alleles, strategy, strategy_details, ins_dec, mixmax){

  # errors in input constants
  #if(pk<0||pk>1){
  #  stop("The parameter pk must be in the range from 0 to 1.")
  #}
  #if(zk<0||zk>1){
  #  stop("The parameter zk must be in the range from 0 to 1.")
  #}

  # if there are errors in insecticides_applied, return matrix of its dimensions with full details of all the errors
  #ia_error = insecticides_applied
  # if outside of continuous range then error
  #ia_error[((insecticides_applied[,2]<0)+(insecticides_applied[,2]>1))==1,2] = "ERROR"
  #ia_error[((insecticides_applied[,3]<0)+(insecticides_applied[,3]>1))==1,3] = "ERROR"
  #ia_error[((insecticides_applied[,4]<0)+(insecticides_applied[,4]>1))==1,4] = "ERROR"

  # if there are errors in resistant_alleles, return matrix of its dimensions with full details of all the errors
  #ra_error = resistant_alleles
  # if lack categorical value then error
  #ra_error[((resistant_alleles[,2]==1)+(resistant_alleles[,2]==2))!=1,2] = "ERROR"
  # if outside of continuous range then error
  #ra_error[((resistant_alleles[,4]<0)+(resistant_alleles[,4]>1))==1,4] = "ERROR"
  #ra_error[((resistant_alleles[,5]<0)+(resistant_alleles[,5]>1))==1,5] = "ERROR"
  #ra_error[((resistant_alleles[,6]<0)+(resistant_alleles[,6]>1))==1,6] = "ERROR"
  #ra_error[((resistant_alleles[,7]<0)+(resistant_alleles[,7]>1))==1,7] = "ERROR"
  #ra_error[((resistant_alleles[,8]<0)+(resistant_alleles[,8]>1))==1,8] = "ERROR"
  #ra_error[((resistant_alleles[,9]<0)+(resistant_alleles[,9]>1))==1,9] = "ERROR"
  # collate errors and potentially return findings
  #if(sum(ra_error=="ERROR")>0){
  #  if(sum(((resistant_alleles[,2]==1)+(resistant_alleles[,2]==2))!=1)>0){
  #    print("Each resistant locus must be either haploid (=1) or diploid (=2).")}
  #  if(sum(((resistant_alleles[,4:9]<0)+(resistant_alleles[,4:9]>1))==1)>0){
  #    print("Each resistant allele parameter must be in the range from 0 to 1.")}
  #  print(ra_error)
  #  stop("As indicated in the print-out, there is one or more ERROR entries in the resistant_alleles matrix.")
  #}
  ## once resistant_alleles has no other errors, if dominance is given but is haploid then error
  #ra_error[((resistant_alleles[,2]==1)*(resistant_alleles[,5]!=0))>0,5] = "ERROR"
  #ra_error[((resistant_alleles[,2]==1)*(resistant_alleles[,7]!=0))>0,7] = "ERROR"
  #ra_error[((resistant_alleles[,2]==1)*(resistant_alleles[,9]!=0))>0,9] = "ERROR"
  #if(sum(ra_error=="ERROR")>0){
  #  print(ra_error)
  #  stop("As indicated in the print-out, there is one or more ERROR entries in the resistant_alleles matrix because a haploid locus has non-zero dominance.")
  #}

  ####

  # useful to calculate the total number of types, and the length of wfset
  tg = prod(1 + 2^(as.numeric(resistant_alleles[,2])-1)) # the total number of genotypes

  # genotype names and costs  only relate to the resistant_alleles matrix
  prevset = 1 # previous set size
  genotype_names = matrix(0,tg,nrow(resistant_alleles)) # labels for allele designations
  genotype_fitness_effect_cost = matrix(1,tg,nrow(resistant_alleles))
  for(j in 1:nrow(resistant_alleles)){ # take each resistant allele in turn
    # update short-hand definitions
    ck = as.numeric(resistant_alleles[j,8]) # cost
    hc = as.numeric(resistant_alleles[j,9]) # cost dominance
    # haploid
    if(resistant_alleles[j,2]==1){
      locus_alleles = c(toupper(resistant_alleles[j,1]),tolower(resistant_alleles[j,1])) # names set
      currset = 2*prevset # current set size, as there are two alleles of the resistant allele and wildtype
      costs = c(1-ck,1)
    }
    # diploid
    if(resistant_alleles[j,2]==2){
      locus_alleles = c(paste(toupper(resistant_alleles[j,1]),toupper(resistant_alleles[j,1]),sep=""),
                        paste(toupper(resistant_alleles[j,1]),tolower(resistant_alleles[j,1]),sep=""),
                        paste(tolower(resistant_alleles[j,1]),tolower(resistant_alleles[j,1]),sep="")) # names set
      currset = 3*prevset # current set size, as there are two alleles of the homozygote resistant allele, heterozygote and homozygote wildtype
      costs = c(1-ck,1-hc*ck,1)
    }
    # repeat elements to populate the matrices such that each row is an unique combination
    genotype_names[,j] = rep( rep(locus_alleles, each=prevset) , tg/currset )
    genotype_fitness_effect_cost[,j] = rep( rep(costs, each=prevset) , tg/currset )
    # update previous set size for next iteration of the loop
    prevset = prod(1 + 2^(as.numeric(resistant_alleles[,2])[1:j]-1))
  }

  # multiply effects across loci
  ce = apply(genotype_fitness_effect_cost,1,prod)

  ####

  # empty matrices for the cumulative effects across resistant alleles per insecticide
  tox_no_ass = matrix(0,tg,nrow(insecticides_applied))
  tox_ass = matrix(0,tg,nrow(insecticides_applied))
  no_ass = matrix(0,tg,nrow(insecticides_applied))

  # toxicity and assortment
  for(i in 1:nrow(insecticides_applied)){ # take each insecticide in turn

    # short-hand definitions, zero unless the insecticide is in use
    tk = 0 # toxicity
    ak = 0 # assortment
    xk = 0 # exposure with assortment

    # if the insecticide is in use then update short-hand definitions
    if(in_use[i]==1){
      # short-hand definitions
      tk = as.numeric(insecticides_applied[i,2]*ins_dec[i]) # toxicity
      if(strategy=="Mixture"&mixmax==F){
        mixk = (insecticides_applied[1,2]+insecticides_applied[2,2]-sqrt((((insecticides_applied[1,2]+insecticides_applied[2,2])^2)-2*insecticides_applied[1,2]*insecticides_applied[2,2]*(insecticides_applied[1,2]+insecticides_applied[2,2]))))/(2*insecticides_applied[1,2]*insecticides_applied[2,2])
        tk = as.numeric(mixk*insecticides_applied[i,2]*ins_dec[i])} # making mixtures have equal initial control to all other strategies
      ak = as.numeric(insecticides_applied[i,3]) # assortment
      xk = as.numeric(insecticides_applied[i,4]) # exposure with assortment
    }

    # genotype names and costs  only relate to the resistant_alleles matrix
    prevset = 1 # previous set size
    genotype_fitness_effect_toxicity = matrix(1,tg,nrow(resistant_alleles))
    genotype_fitness_effect_assortment = matrix(1,tg,nrow(resistant_alleles))

    for(j in 1:nrow(resistant_alleles)){ # take each resistant allele in turn

      # short-hand definitions, zero unless resistant allele interacts with the insecticide
      mtk = 0 # modifier of toxicity
      hmt = 0 # modifier of toxicity dominance
      mak = 0 # modifier of assortment
      hma = 0 # modifier of assortment dominance

      # if the resistant allele interacts with the insecticide : 0 generalist and >0 specialist
      if(resistant_alleles[j,3]==0||resistant_alleles[j,3]==insecticides_applied[i,1]){
        # short-hand definitions
        mtk = as.numeric(resistant_alleles[j,4]) # modifier of toxicity
        hmt = as.numeric(resistant_alleles[j,5]) # modifier of toxicity dominance
        mak = as.numeric(resistant_alleles[j,6]) # modifier of assortment
        hma = as.numeric(resistant_alleles[j,7]) # modifier of assortment dominance
      }

      # haploid
      if(resistant_alleles[j,2]==1){
        locus_alleles = c(toupper(resistant_alleles[j,1]),tolower(resistant_alleles[j,1])) # names set
        currset = 2*prevset # current set size, as there are two alleles of the resistant allele and wildtype
        toxicity = c(1-mtk,1)
        assortment = c(1-mak,1)
      }

      # diploid
      if(resistant_alleles[j,2]==2){
        locus_alleles = c(paste(toupper(resistant_alleles[j,1]),toupper(resistant_alleles[j,1]),sep=""),
                          paste(toupper(resistant_alleles[j,1]),tolower(resistant_alleles[j,1]),sep=""),
                          paste(tolower(resistant_alleles[j,1]),tolower(resistant_alleles[j,1]),sep="")) # names set
        currset = 3*prevset # current set size, as there are two alleles of the homozygote resistant allele, heterozygote and homozygote wildtype
        toxicity = c(1-mtk,1-hmt*mtk,1)
        assortment = c(1-mak,1-hma*mak,1)
      }

      # repeat elements to populate the matrices such that each row is an unique combination
      genotype_fitness_effect_toxicity[,j] = rep( rep(toxicity, each=prevset) , tg/currset )
      genotype_fitness_effect_assortment[,j] = rep( rep(assortment, each=prevset) , tg/currset )
      # update previous set size for next iteration of the loop
      prevset = prod(1 + 2^(as.numeric(resistant_alleles[,2])[1:j]-1))

    }

    # cumulative effects across resistant alleles per insecticide
    tox_no_ass[,i] = 1-apply(genotype_fitness_effect_toxicity,1,prod)*tk
    tox_ass[,i] = 1-apply(genotype_fitness_effect_toxicity,1,prod)*tk*xk
    no_ass[,i] = apply(genotype_fitness_effect_assortment,1,prod)*(1-ak)

  }

  # if mixture or rotation, multiply effects across insecticides to allow an insecticide to be exposed to both at once or one at a time
  if(strategy!="Combination"){
    tn = apply(tox_no_ass,1,prod)
    ta = apply(tox_ass,1,prod)
    na = apply(no_ass,1,prod)
  }

  # if combinations, add effects weighted by the probability (in strategy details) that experience them
  if(strategy=="Combination"){
    tn = apply(sweep(tox_no_ass,MARGIN=2,strategy_details,'*'),1,sum)
    ta = apply(sweep(tox_ass,MARGIN=2,strategy_details,'*'),1,sum)
    na = apply(sweep(no_ass,MARGIN=2,strategy_details,'*'),1,sum)
  }

  ####

  # assemble fitness from the composite of its parts and associate the genotype names
  wfset = ( (1-zk)*pk*na*tn + (1-zk)*pk*(1-na)*ta + (1-(1-zk)*pk) ) * ce
  # (1-zk)*pk*ra*rt : human & protected & not repelled (naturally or by resistance) = toxicity cost, which can be alleviated by resistance
  # (1-zk)*pk*(1-ra)*nt : human & protected & repelled (naturally or by resistance) = exposure to insecticide governed by xk
  # (1-(1-zk)*pk) : animal or unprotected human = baseline fitness 1
  # rc : costs are applied multiplicatively to all (even if have no cost, where the element will equal 1)
  names(wfset) = apply(genotype_names,1,paste,collapse="")

  return(wfset)

}


# END
