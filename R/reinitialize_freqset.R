#' Reinitialize genotype frequency from allele frequency
#'
#' Reinitialize bi-allele multi-locus genotype frequency from underlying allele frequencies assuming no prior linkage disequilibrium
#'
#'  @param aset_added A vector with the frequency of all alleles
#'  @param fset A vector with the frequency of all genotypes
#'  @param aset_add A vector of frequencies 'f' of the focal alleles (upper-case, e.g. A) across all loci, where non-focal alleles (lower-case, e.g. a) have frequency 1-f
#'  @param anames A vector of designated letters associated with each locus, with upper-case denoting the focal allele and lower-case denoting the non-focal (e.g. A vs a)
#'  @param ploidy A vector of the ploidies of each locus
#'  @param add_loci A vector of the number of loci to be added
#'
#'  @return A vector of the genotype frequencies
#'  @export
#'  @examples
reinitialize_freqset <- function(aset_added,fset,aset_add,anames,ploidy,add_loci){
  if(sum(ploidy==1,ploidy==2)!=length(ploidy)){
    stop("Only haploid (ploidy=1) or diploid (ploidy=2) genetic architectures are permitted.")
  }

  if(length(anames)!=length(ploidy)){
    stop("Either too few allele frequencies or names or ploidies have been specified.")
  }

  # generate empty genotype_names and genotype_freq matrices
  n = length(fset) # short-hand definition, number of genotypes
  la = length(add_loci) # short-hand definition, number of loci to add
  sf_al = aset_add[add_loci] # starting frequency of each locus to add

  # setup for loop
  fset_old = fset
  fset_new = fset*(1-sum(sf_al)) # decrease frequency of existing types (from mutation)
  aset_add_locus = aset_added

  for(j in la:1){ # reverse order (NB: just because easier to conceptualise, I think)

    # find new values associated with the addition of a mutant allele at this locus
    # generate iaal which carries indices associated with the fset_new that will generate
    aset_add_locus[add_loci[j]] = 0.5
    iaal = initialize_freqset(aset_add_locus,anames,ploidy)
    # find all distinct values where the mutatnt is being added
    fd = rep(0,n)
    fd[which(iaal>0)] = 1
    fd[which(fset_old>0)] = -1
    d = which(fd>0)

    # haploid
    if(ploidy[add_loci[j]]==1){
      fset_new[d] = sf_al[j]*fset_old[which(fset_old>0)] # divide mutants between the existing genotypes
    }

    # diploid
    if(ploidy[add_loci[j]]==2){
      homs = sf_al[j]*sf_al[j]*fset_old[which(fset_old>0)] # homozygotes
      hets = sf_al[j]*(1-sf_al[j])*fset_old[which(fset_old>0)] # heterozygotes
      fset_new[d] = c(homs, hets) # dhoms always before hets from the order of fset_new
    }

    # update fset_old and can check that still sums to 1
    fset_old = fset_old*(1-sf_al)
    fset_old[d] = fset_new[d]

  }

  names(fset_new) = names(fset)

  return(fset_new)

}


# END
