#' Sexual recombination that occurs freely with potential linkage disequilibrium
#'
#' Provides offspring bi-allele multi-locus frequencies under the assumption of panmixia, where female / male fitness can be different and selection can build linkage disequilibrium
#'
#'  @param ffset A vector with the frequency of all female genotypes
#'  @param fmset A vector with the frequency of all male genotypes
#'  @param resistant_alleles A matrix that encodes the properties of the resistant alleles, as one allele per row
#'
#'  @return A vector with the frequency of all genotypes
#'  @export
#'  @examples
random_mating <- function(ffset, fmset, resistant_alleles){

  if(length(ffset)!=length(fmset)){
    stop("Vectors ffset and fmset have different lengths.")
  }

  an = resistant_alleles[,1] # short-hand definition, names of alleles
  ploidy = as.numeric(resistant_alleles[,2]) # short-hand definition, ploidy of loci
  gv = genotypic_values(an,ploidy) # values of genotypes in terms of allelic replication
  l = length(ploidy) # short-hand definition, number of loci
  n = length(ffset) # short-hand definition, number of genotypes
  nf = 2^l # number of female genotypes
  nm = 2^sum(ploidy-1) # number of male genotypes
  gffset = rep(0,nf)
  gfmset = rep(0,nm)
  gfnames = matrix(0,nf,l)
  gmnames = matrix(0,nm,l)


  #################################################

  # NAMES

  fprevset = 1 # previous set size for females
  mprevset = 1 # previous set size for males
  for(j in 1:l){

    gla = c(toupper(an[j]), tolower(an[j])) # gametes names set
    la = length(gla) # short-hand definition
    fcurrset = la*fprevset # current set size
    gfnames[,j] = rep( rep(gla, each=fprevset) , nf/fcurrset )
    fprevset = 2^j

    # diploid, then consider males as well
    if(ploidy[j]==2){
      mcurrset = la*mprevset # current set size for males
      gmnames[,j] = rep( rep(gla, each=mprevset) , nm/mcurrset )
      mprevset = 2^sum(ploidy[1:j]==2)
    }

  }

  names(gffset) = apply(gfnames,1,paste,collapse="")
  names(gfmset) = apply(gmnames,1,paste,collapse="")


  #################################################

  # FREQUENCIES

  # setup female gamete frequencies
  gffset[1] = sum(apply(gv,1,prod)*ffset) # vector of gamete frequencies, entering first value
  fmat1 = gv # first matrix
  kf = 2 # seed value X for next fmatX definition

  # setup female gamete frequencies
  gvm = gv[,ploidy==2]
  if(is.vector(gv[,ploidy==2])==TRUE){ gvm = cbind(gvm, rep(1,length(gvm))) } # if there is only one male locus this ensures no errors going forward
  gfmset[1] = sum(apply(gvm,1,prod)*fmset)  # vector of gamete frequencies, entering first value
  mmat1 = gvm # first matrix
  km = 2 # seed value X for next mmatX definition
  jm = 1 # initial diploid locus number

  # calculate gamete frequencies
  for(j in 1:l){ # # go through all loci
    for(i in 1:(kf-1)){ # go through all previously created matrices
      eval(parse(text=(paste("fmat",kf," = fmat",i,sep="")))) # create new matrix from old
      eval(parse(text=(paste("fmat",kf,"[,",j,"] = (-1*gv[,",j,"]+1)",sep="")))) # replace one column with other allele
      eval(parse(text=(paste("gffset[",kf,"] = sum(apply(fmat",kf,",1,prod)*ffset)",sep="")))) # calculate frequency
      kf = kf + 1 # accumulative female locus number
    }
    if(ploidy[j]==2){
      for(i in 1:(km-1)){ # go through all previously created matrices
        eval(parse(text=(paste("mmat",km," = mmat",i,sep="")))) # create new matrix from old
        eval(parse(text=(paste("mmat",km,"[,",jm,"] = (-1*gvm[,",jm,"]+1)",sep="")))) # replace one column with other allele
        eval(parse(text=(paste("gfmset[",km,"] = sum(apply(mmat",km,",1,prod)*fmset)",sep="")))) # calculate frequency
        km = km + 1 # accumulative male locus number
      }
      jm = jm + 1 # accumulative diploid locus number
    }
  }


  #################################################

  # OFFSPRING GENOTYPE FREQUENCIES

  # start with fxm matrix which has the frequencies of female-male parent combinations
  # we want to know how to get from fxm matrix to f1set vector of genotype frequencies
  # so we calculate a set of coordinates in fxm fcoord and mcoord alongside the f1set element number
  # then at the end we go through fxm using the coordinates to attrbiute each matrix element to the correct element of f1set

  # we calculate the matching coordinates and element number in order of the loci as given to follow its native structure
  #   where each locus is given as A first then a, B first then b etc
  # the major complication is that a locus can be haploid or diploid
  # a haploid locus is female-only so we can create the new coordinates by
  #   duplicating the old mcoords
  #   adding the current max fcoord to the old fcoords to create new fcoords
  #   which finds the new f1set values in fxm for the alternate allele at this locus (a) below the values in fxm for the other locus (A)
  # a diploid locus is for both sexes so we can create the new coordinates by
  #   duplicating the old coords (for fcoords and mcoords) and adding the current max coord to the old coords to create new coords
  #   adding each set of new coords in the combinations of old coords such that:
  #     ignore old fcoords x old mcoords (which are already present in coords)
  #     add new fcoords x old mcoords and old fcoords x new mcoords (which are new but have identical f1set element numbers to one another)
  #     add new fcoords x new mcoords (which are new and have non-identical f1set element number)

  # tagtest is useful as the number of elements always increases, so it tells us how the algorithm is moving through fxm
  #   which reveals the hidden workings of the computational process of repsonding to a haploid/diploid locus

  fxm = gffset%o%gfmset # female by male reproduction matrix
  # enter starting values from first entry at coordinate c(1,1) in fxm
  fcoord = c(1) # first row/female coordinate
  mcoord = c(1) # first col/male coordiante
  gennum = c(1) # first element in f1set
  tagtest = c(1) # for testing
  kf = 1
  km = 1
  kg = 1
  kt = 1

  # generate mapping between fxm and f1set
  for(j in 1:l){
    # female potential coordinates
    fold = fcoord
    fnew = kf + fold
    kf = max(fnew)
    # male potential coordinates
    mold = mcoord
    # genotype potential number
    gold = kg+gennum
    # if haploid, simply 'go down' on previous
    if(ploidy[j]==1){
      # coordinate
      fcoord = append(fcoord,fnew)
      mcoord = append(mcoord,mold)
      # genotype number
      gennum = append(gennum,gold)
      kg = max(gennum)
      # test coordinate placement
      tagtest = append(tagtest,seq(kt+1,kt+length(fnew)))
      kt = max(tagtest)
    }
    # if diploid, go down + go right + go down&&right
    if(ploidy[j]==2){
      # coordinate
      mnew = km + mold
      km = max(mnew)
      fcoord = append(fcoord,c(fnew,fold,fnew))
      mcoord = append(mcoord,c(mold,mnew,mnew))
      # genotype number
      gnew = kg + gold
      gennum = append(gennum,c(gold,gold,gnew))
      kg = max(gennum)
      # test
      tagtest = append(tagtest,seq(kt+1,kt+length(c(fnew,fold,fnew))))
      kt = max(tagtest)
    }
  }

  # testing
  #fxm_test = fxm
  #for(i in 1:length(gennum)){
  #  fxm_test[fcoord[i],mcoord[i]] = tagtest[i] # prints fxm-like matrix with the order the gennum values for f1set were given in fxm
  #  #fxm_test[fcoord[i],mcoord[i]] = gennum[i] # prints fxm-like matrix with the element number in f1set
  #}
  #fxm_test

  # apply mapping of fcoord/mcoord/gennum to transfer frequencies from fxm to f1set
  f1set = rep(0,n) # empty offspring genotypes vector, i.e. frequency in next generations
  # populate empty vector from fxm
  for(i in 1:length(gennum)){
    f1set[gennum[i]] = f1set[gennum[i]] + fxm[fcoord[i],mcoord[i]]
  }


  #################################################

  # ASSEMBLE OUTPUT

  # transfer names
  names(f1set) = names(ffset)

  return(f1set)

}


# END
