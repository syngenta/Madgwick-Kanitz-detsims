#' Initialize genotype frequency from allele frequency
#'
#' Initialize bi-allele multi-locus genotype frequency from underlying allele frequencies assuming no prior linkage disequilibrium
#'
#'  @param allele_frequencies A vector of frequencies 'f' of the focal alleles (upper-case, e.g. A) across all loci, where non-focal alleles (lower-case, e.g. a) have frequency 1-f
#'  @param allele_names A vector of designated letters associated with each locus, with upper-case denoting the focal allele and lower-case denoting the non-focal (e.g. A vs a)
#'  @param ploidy A vector of the ploidies of each locus
#'
#'  @return A vector of the genotype frequencies
#'  @export
#'  @examples
initialize_freqset <- function(allele_frequencies,allele_names,ploidy){
  if(sum(ploidy==1,ploidy==2)!=length(ploidy)){
    stop("Only haploid (ploidy=1) or diploid (ploidy=2) genetic architectures are permitted.")
  }

  if(sum(allele_frequencies<0)==1||sum(allele_frequencies>1)==1){
    stop("Allele frequencies must be within the range between 0 and 1.")
  }

  if(length(allele_frequencies)!=length(allele_names)||length(allele_frequencies)!=length(ploidy)){
    stop("Either too few allele frequencies or names or ploidies have been specified.")
  }

  # generate empty genotype_names and genotype_freq matrices
  l = length(allele_frequencies) # short-hand definition, equals the other lengths as well
  t = prod(1 + 2^(ploidy-1)) # the total number of genotypes
  genotype_names = matrix(0,t,l)
  genotype_freqs = matrix(0,t,l)

  # populate genotype_names and genotype_freq matrices
  prevset = 1 # previous set size
  for(j in 1:l){

    # haploid
    if(ploidy[j]==1){
      locus_alleles = c(toupper(allele_names[j]),tolower(allele_names[j])) # names set
      alleles_freqs = c(allele_frequencies[j],1-allele_frequencies[j]) # frequency set
      la = length(locus_alleles) # short-hand definition, same for locus_alleles and alleles_freq
      currset = la*prevset # current set size
    }

    # diploid
    if(ploidy[j]==2){
      locus_alleles = c(paste(toupper(allele_names[j]),toupper(allele_names[j]),sep=""),
                        paste(toupper(allele_names[j]),tolower(allele_names[j]),sep=""),
                        paste(tolower(allele_names[j]),tolower(allele_names[j]),sep="")) # names set
      alleles_freqs = c(allele_frequencies[j]*allele_frequencies[j],
                       2*allele_frequencies[j]*(1-allele_frequencies[j]),
                       (1-allele_frequencies[j])*(1-allele_frequencies[j])) # frequency set
      la = length(locus_alleles) # short-hand definition, same for locus_alleles and alleles_freq
      currset = la*prevset # current set size
    }

    # repeat elements to populate the matrices such that each row is an unique combination
    genotype_names[,j] = rep( rep(locus_alleles, each=prevset) , t/currset )
    genotype_freqs[,j] = rep( rep(alleles_freqs, each=prevset) , t/currset )

    # redefine previous set size for next iteration of the loop
    prevset = prod(1 + 2^(ploidy[1:j]-1))

  }

  # calculate the product and associate the names
  freqset = apply(genotype_freqs,1,prod)
  names(freqset) = apply(genotype_names,1,paste,collapse="")

  return(freqset)

}


# END
