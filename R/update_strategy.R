#' Update insecticides that are in use in line with the strategy
#'
#' Expresses strategy as a change in insecticide usage dependent upon conditions for when a strategy involves changing an insecticide
#'
#'  @param tt Current time t in generations
#'  @param strategy Either Combination, Mixture or Rotation
#'  @param strategy_start The time t in generations that the strategy starts
#'  @param strategy_details A vector informing the strategy
#'  @param res_freq The frequency of resistant alleles within the populaion
#'  @param in_use The previous vector of in_use
#'  @param rk THe replacement / retreatment schedule
#'
#'  @return A vector of the fitness of different genotypes
#'  @export
#'  @examples
update_strategy <- function(tt, strategy, strategy_start, strategy_details, res_freq, in_use, rk){

  # default position is for previously in-use insecticides to be in use
  iu = in_use

   # update strategy of in use insecticides if strategy is set to start or this is a bednet replacement generation
  if((tt>=strategy_start)&&(((tt-strategy_start)%%rk)==0)){

    #### fixed strategies ####

    # if Combination strategy, all insecticides are used simultaneously or not at all
    if(strategy=="Combination"){
      iu = rep(1,length(strategy_details))
      # unless assigned 0 value and not used
      if(sum(strategy_details==0)>0){
        iu[strategy_details==0] = 0
      }
    }

    # if Mixture strategy, all insecticides are used simultaneously or not at all
    if(strategy=="Mixture"){
      iu = rep(1,length(strategy_details))
    }

    # if Rotation strategy, switch between insecticides in the pattern of strategy details
    if(strategy=="Rotation"){
      iu = rep(0,length(strategy_details))
      # calculate the series of rotations
      tseries = c()
      for(pp in 1:length(strategy_details)){
        tseries = append(tseries, rep(pp, strategy_details[pp]))
      }
      # select one insecticide for use based on how far through tseries time tt currently is
      iu[tseries[1+((tt-strategy_start)%%sum(strategy_details))]] = 1
    }

    #### conditional strategies ####

    # if Sequential strategy, switch between insecticides when hit the threshold in strategy details
    if(strategy=="Sequential"){
      # default is to change nowt
      iu = in_use
      # if no previous insecticide in use then start with first one
      if(sum(in_use)==0){iu[1]=1}
      # else, move onto next insecticide if over resistance threshold
      if(sum(in_use)>0){
        ci = which(in_use==1) # current in-use insecticide
        if(res_freq[ci]>=strategy_details[ci]){ # is over resistance threshold? yes then:
          iu = rep(0,length(strategy_details)) # remove old strategy
          iu[(ci+1)%%(1+length(strategy_details))] = 1 # and move on in sequence
        }
      }
    }

    # if Adaptive strategy, start from mixture until both overcome thresholds in strategy details then switch to insecticide with rarer resistance
    if(strategy=="Adaptive"){
      # default is to change nowt
      iu = in_use
      # if no previous insecticide in use then start with mixture
      if(sum(in_use)==0){iu = rep(1,length(strategy_details)) }
      # else, only apply insecticides if not over resistance threshold
      if(sum(in_use)>0){
        ci = which.max(in_use==1) # current in-use insecticide
        if(sum(res_freq[ci]>=strategy_details[ci])>0){ # is over resistance threshold for one insecticide? yes then:
          iu = rep(1,length(strategy_details)) # remove old strategy
          iu[which(res_freq[ci]>=strategy_details[ci])] = 0 # turn off insecticides that resistance over threshold
        }
      }
    }

  }

  return(iu)

}


# END
