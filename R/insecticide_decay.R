#' Attributes decay to insecticide toxicity
#'
#' Permits the effectiveness of an insecticide to decay over time according to its half-life and replacement / retreament schedule
#'
#'  @param tt Time in generations
#'  @param in_use0 The previous vector of in_use
#'  @param in_use1 The current vector of in_use
#'  @param dk Decay half-life in generations
#'  @param rk After how many generations does bednet replacement occur
#'  @param ins_dec0 The previous vector of ins_dec0
#'
#'  @return A vector of the fitness of decay and when last replaced bednet
#'  @export
#'  @examples
insecticide_decay <- function(tt, in_use0, in_use1, dk, rk, ins_dec0){

  # start with empty vector of insecticide decay
  id = rep(0,length(in_use1)+1)

  # if insecticides are in use and nothing has changed, attribute decay cycle
  if((sum(in_use1)>0)&&(sum(in_use0==in_use1)==length(in_use0))){
    id[in_use1>0] = 1/(1+exp(((tt-ins_dec0[length(in_use1)+1])%%rk)-dk))
    id[length(in_use1)+1] = ins_dec0[length(in_use1)+1]
    if(((tt-ins_dec0[length(in_use1)+1])%%rk)==0){id[length(in_use1)+1] = tt}
  }

  # if something has changed, restart decay
  if(sum(in_use0==in_use1)!=length(in_use0)){
    id[in_use1>0] = 1/(1+exp(-dk))
    id[length(in_use1)+1] = tt
  }

  return(id)

}

