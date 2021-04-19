#' Initialize genotypic value of offspring from ploidy with respect to the focal allele
#'
#' Initialize bi-allele multi-locus genotype frequency from underlying allele frequencies assuming no prior linkage disequilibrium
#'
#'  @param allele_names A vector of designated letters associated with each locus, with upper-case denoting the focal allele and lower-case denoting the non-focal (e.g. A vs a)
#'  @param ploidy A vector of the ploidies of each locus
#'
#'  @return A vector of the genotype frequencies
#'  @export
#'  @examples
genotypic_values <- function(allele_names,ploidy){
  if(sum(ploidy==1,ploidy==2)!=length(ploidy)){
    stop("Only haploid (ploidy=1) or diploid (ploidy=2) genetic architectures are permitted.")
  }

  if(length(allele_names)!=length(ploidy)){
    stop("Either too few names or ploidies have been specified.")
  }

  # generate empty genotype_names and genotype_freq matrices
  l = length(allele_names) # short-hand definition, equals the other lengths as well
  t = prod(1 + 2^(ploidy-1)) # the total number of genotypes
  genotype_names = matrix(0,t,l)
  genotype_value = matrix(0,t,l)

  # populate genotype_names and genotype_freq matrices
  prevset = 1 # previous set size
  for(j in 1:l){

    # haploid
    if(ploidy[j]==1){
      locus_alleles = c(toupper(allele_names[j]),tolower(allele_names[j])) # names set
      alleles_value = c(1, 0) # value set
      la = length(locus_alleles) # short-hand definition, same for locus_alleles and alleles_freq
      currset = la*prevset # current set size
    }

    # diploid
    if(ploidy[j]==2){
      locus_alleles = c(paste(toupper(allele_names[j]),toupper(allele_names[j]),sep=""),
                        paste(toupper(allele_names[j]),tolower(allele_names[j]),sep=""),
                        paste(tolower(allele_names[j]),tolower(allele_names[j]),sep="")) # names set
      alleles_value = c(1, 0.5, 0) # value set
      la = length(locus_alleles) # short-hand definition, same for locus_alleles and alleles_freq
      currset = la*prevset # current set size
    }

    # repeat elements to populate the matrices such that each row is an unique combination
    genotype_names[,j] = rep( rep(locus_alleles, each=prevset) , t/currset )
    genotype_value[,j] = rep( rep(alleles_value, each=prevset) , t/currset )

    # redefine previous set size for next iteration of the loop
    prevset = prod(1 + 2^(ploidy[1:j]-1))

  }

  # calculate the product and associate the names
  colnames(genotype_value) = toupper(allele_names)
  rownames(genotype_value) = apply(genotype_names,1,paste,collapse="")

  return(genotype_value)

}


# END
