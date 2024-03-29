% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insecticide_decay.R
\name{insecticide_decay}
\alias{insecticide_decay}
\title{Attributes decay to insecticide toxicity}
\usage{
insecticide_decay(tt, in_use0, in_use1, dk, rk, ins_dec0)
}
\description{
Permits the effectiveness of an insecticide to decay over time according to its half-life and replacement / retreament schedule
}
\details{
@param tt Time in generations
 @param in_use0 The previous vector of in_use
 @param in_use1 The current vector of in_use
 @param dk Decay half-life in generations
 @param rk After how many generations does bednet replacement occur
 @param ins_dec0 The previous vector of ins_dec0

 @return A vector of the fitness of decay and when last replaced bednet
 @export
 @examples
}
