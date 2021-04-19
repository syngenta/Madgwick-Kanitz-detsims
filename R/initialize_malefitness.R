#' Initialize male fitness of different genotypes
#'
#' Initialize fitness of bi-allele multi-locus genotypes of males from underlying selection coefficients frequencies, assuming that a male has percentage similar fitness to a female but always pays the cost of resistance
#'
#'  @param wfset The vector of female fitness
#'  @param ek The constant for the degree to which males and females experience the same fitness
#'  @param resistant_alleles A matrix that encodes the properties of the resistant alleles, as one allele per row
#'  The resistant alleles matrix takes the form of:
#'  name (i.e. a letter), 1=haploid or 2=diploid, T=toxicity or A=assortment, effect, dominance of effect, cost, dominance of cost
#'
#'  @return A vector of the fitness of different genotypes
#'  @export
#'  @examples
initialize_malefitness <- function(wfset, ek, resistant_alleles){
  #if(ek<0||ek>1){
  #  stop("The parameter ek must be in the range from 0 to 1.")
  #}
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
  # once resistant_alleles has no other errors, if dominance is given but is haploid then error
  #ra_error[((resistant_alleles[,2]==1)*(resistant_alleles[,5]!=0))>0,5] = "ERROR"
  #ra_error[((resistant_alleles[,2]==1)*(resistant_alleles[,7]!=0))>0,7] = "ERROR"
  #ra_error[((resistant_alleles[,2]==1)*(resistant_alleles[,9]!=0))>0,9] = "ERROR"
  #if(sum(ra_error=="ERROR")>0){
  #  print(ra_error)
  #  stop("As indicated in the print-out, there is one or more ERROR entries in the resistant_alleles matrix because a haploid locus has non-zero dominance.")
  #}

  # generate empty genotype_names and genotype_fitness matrices
  l = nrow(resistant_alleles) # short-hand definition, equals the other lengths as well
  z = as.numeric(resistant_alleles[,2]) # short-hand definition, the ploidy of loci as a number
  tg = prod(1 + 2^(z-1)) # the total number of genotypes
  genotype_names = matrix(0,tg,l)
  genotype_fitness_effect_cost = matrix(0,tg,l)

  # populate genotype_names and genotype_fitness matrices
  prevset = 1 # previous set size
  for(j in 1:l){

    # short-hand definitions
    ck = as.numeric(resistant_alleles[j,8]) # cost
    hc = as.numeric(resistant_alleles[j,9]) # cost dominance

    # haploid
    if(resistant_alleles[j,2]==1){
      locus_alleles = c(toupper(resistant_alleles[j,1]),tolower(resistant_alleles[j,1])) # names set
      la = length(locus_alleles) # short-hand definition
      currset = la*prevset # current set size
      # default values of wildtypes assembling fitness
      effect_costs = rep(1,la)
      # if focal allele is present (upper-case) then attribute accordingly based on allele type
      effect_costs[1] = 1-ck
    }

    # diploid
    if(resistant_alleles[j,2]==2){
      locus_alleles = c(paste(toupper(resistant_alleles[j,1]),toupper(resistant_alleles[j,1]),sep=""),
                        paste(toupper(resistant_alleles[j,1]),tolower(resistant_alleles[j,1]),sep=""),
                        paste(tolower(resistant_alleles[j,1]),tolower(resistant_alleles[j,1]),sep="")) # names set
      la = length(locus_alleles) # short-hand definition, same for locus_alleles and alleles_freq
      currset = la*prevset # current set size
      # default values of wildtypes assembling fitness
      effect_costs = rep(1,la)
      # if focal allele is present (upper-case) then attribute accordingly based on allele type
      effect_costs[1:2] = c(1-ck, 1-hc*ck)
    }

    # repeat elements to populate the matrices such that each row is an unique combination
    genotype_names[,j] = rep( rep(locus_alleles, each=prevset) , tg/currset )
    genotype_fitness_effect_cost[,j] = rep( rep(effect_costs, each=prevset) , tg/currset )

    # redefine previous set size for next iteration of the loop
    prevset = prod(1 + 2^(z[1:j]-1))

  }

  # assemble fitness from the composite of its parts and associate the genotype names
  wmset = ek*wfset + (1-ek) * apply(genotype_fitness_effect_cost,1,prod)
  names(wmset) = apply(genotype_names,1,paste,collapse="")

  return(wmset)

}


# END
