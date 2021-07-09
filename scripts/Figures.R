#' DEAR REVIEWER, download the data set on your browser and unzip to ./data/ with
#' this temporary Dryad link:
#' https://datadryad.org/stash/share/dlQTBmwAFmCpQMhunJvk2UFb4xJknABIoxAesfhZXa4


# Figures that involve special treatment of lost resistant alleles and population extinction

# packages
# library(detsims)
library(party)
library(partykit)
library(vioplot)
library(viridisLite)
library(matrixStats)
library(pracma)
# library(rdryad)

# parameter set ranges as c(min,max)
npsr = c(2,9) # population size, NB: log scale
#bksr = c(0,2) # population growth rate, NB: if random then log-normal with mean=0 and sd=1
gksr = c(0,1) # adult death rate
pksr = c(0,1) # proportion of females exposed to insecticide
eksr = c(0,1) # proportion of males that have female exposure
fAsr = -1*c(9,2) # initial frequency of allele A, NB: log scale
fBsr = -1*c(9,2) # initial frequency of allele B, NB: log scale
t1sr = c(0,1) # effectiveness of insecticide 1
t2sr = c(0,1) # effectiveness of insecticide 2
rtAsr = c(0,1) # resistance of allele A
rtBsr = c(0,1) # resistance of allele B
hrAsr = c(0,1) # dominance of resistance of allele A
hrBsr = c(0,1) # dominance of resistance of allele B
crAsr = -1*c(3,0.5) # resistance cost of allele A
crBsr = -1*c(3,0.5) # resistance cost of allele B
hcAsr = c(0,1) # dominance of resistance cost of allele A
hcBsr = c(0,1) # dominance of resistance cost of allele B
# simulation constants
strategyset = c("SoloA","SoloB","Mosaic","Mixture","Rotation","Sequence") # strategy that simulations have been run for
strategy_start = 1 # the generation number when the strategy starts
gsf = strategy_start # here, a redundant function that allows for multiple different start times
gen = 500 # run simulation for max of gen generations
endsim = 0.5 # frequency threshold of resistance for the simulation to end
dk = 1000000  # constant specifying the decay in the toxicity of the insecticide, as the half-life of the treatment - if very large then decay is effectively ignored
rk = 36  # constant specifying after how many generations a bednet is to be replaced / retreated
var_comb = 10^6 # number of random samples

# download data from Zenodo, if not existing locally, here: https://doi.org/10.5281/zenodo.5074770

# read in data
simdata = read.csv("data/PSData_random6.csv")
# nn
sdata_c_nn    = read.csv(paste("data/Data_random6_Mosaic_Fixed_nn.csv",sep=""))
sdata_x_nn    = read.csv(paste("data/Data_random6_Mixture_Fixed_nn.csv",sep=""))
sdata_psx_nn  = read.csv(paste("data/Data_random6_maxpsMixture_Fixed_nn.csv",sep=""))
sdata_r_nn    = read.csv(paste("data/Data_random6_Rotation_Fixed_nn.csv",sep=""))
sdata_r_revnn = read.csv(paste("data/Data_random6_Rotation_Fixed_revnn.csv",sep=""))
sdata_sA_nn   = read.csv(paste("data/Data_random6_SoloA_Fixed_nn.csv",sep=""))
sdata_sB_nn   = read.csv(paste("data/Data_random6_SoloB_Fixed_nn.csv",sep=""))
# mn
sdata_c_mn    = read.csv(paste("data/Data_random6_Mosaic_Fixed_mn.csv",sep=""))
sdata_x_mn    = read.csv(paste("data/Data_random6_Mixture_Fixed_mn.csv",sep=""))
sdata_psx_mn  = read.csv(paste("data/Data_random6_maxpsMixture_Fixed_mn.csv",sep=""))
sdata_r_mn    = read.csv(paste("data/Data_random6_Rotation_Fixed_mn.csv",sep=""))
sdata_r_revmn = read.csv(paste("data/Data_random6_Rotation_Fixed_revmn.csv",sep=""))
sdata_sA_mn   = read.csv(paste("data/Data_random6_SoloA_Fixed_mn.csv",sep=""))
sdata_sB_mn   = read.csv(paste("data/Data_random6_SoloB_Fixed_mn.csv",sep=""))
# mm
sdata_c_mm    = read.csv(paste("data/Data_random6_Mosaic_Fixed_mm.csv",sep=""))
sdata_x_mm    = read.csv(paste("data/Data_random6_Mixture_Fixed_mm.csv",sep=""))
sdata_psx_mm  = read.csv(paste("data/Data_random6_maxpsMixture_Fixed_mm.csv",sep=""))
sdata_r_mm    = read.csv(paste("data/Data_random6_Rotation_Fixed_mm.csv",sep=""))
sdata_r_revmm = read.csv(paste("data/Data_random6_Rotation_Fixed_revmm.csv",sep=""))
sdata_sA_mm   = read.csv(paste("data/Data_random6_SoloA_Fixed_mm.csv",sep=""))
sdata_sB_mm   = read.csv(paste("data/Data_random6_SoloB_Fixed_mm.csv",sep=""))




#### nn data processing ####

# NB: for sequence and rotation, weakest goes first

# raw data
ssdata_t50_nn_A = cbind(sdata_c_nn[,1],sdata_x_nn[,1],sdata_r_nn[,1],sdata_r_revnn[,4],sdata_sA_nn[,1],sdata_psx_nn[,1])
ssdata_t50_nn_B = cbind(sdata_c_nn[,4],sdata_x_nn[,4],sdata_r_nn[,4],sdata_r_revnn[,1],sdata_sB_nn[,4],sdata_psx_nn[,4])
ssdata_c80_nn = cbind(sdata_c_nn[,7],sdata_x_nn[,7],sdata_r_nn[,7],sdata_r_revnn[,7],sdata_sA_nn[,7],sdata_psx_nn[,7],sdata_sB_nn[,7])

# if no measurement (which comes to mean weak selection for resistance by elimination)
ssdata_t50_nn_A[ssdata_t50_nn_A==0] = 1000
ssdata_t50_nn_B[ssdata_t50_nn_B==0] = 1000
ssdata_c80_nn[ssdata_c80_nn==0] = 1000

# if no measurement because of selection against resistance
bsdata_A = simdata[,c(6,6,6,12,6,6)] # starting frequency
bsdata_B = simdata[,c(12,12,12,6,12,12)] # starting frequency
psdata_nn_A = cbind(sdata_c_nn[,2],sdata_x_nn[,2],sdata_r_nn[,2],sdata_r_revnn[,5],sdata_sA_nn[,2],sdata_psx_nn[,2]) # frequency at t=250
psdata_nn_B = cbind(sdata_c_nn[,5],sdata_x_nn[,5],sdata_r_nn[,5],sdata_r_revnn[,2],sdata_sB_nn[,5],sdata_psx_nn[,5]) # frequency at t=250
ssdata_t50_nn_A[psdata_nn_A<bsdata_A] = 1500
ssdata_t50_nn_B[psdata_nn_B<bsdata_B] = 1500

# if no measurement because of extinction set as 'weakest'
esdata_nn = cbind(sdata_c_nn[,10],sdata_x_nn[,10],sdata_r_nn[,10],sdata_r_revnn[,10],sdata_sA_nn[,10],sdata_psx_nn[,10])
ssdata_t50_nn_A[esdata_nn>0] = 0
ssdata_t50_nn_B[esdata_nn>0] = 0
ssdata_c80_nn[cbind(esdata_nn,sdata_sB_nn[,10])>0] = 0

# calculate statistics
# first to 50% frequency
mos = apply(cbind(ssdata_t50_nn_A[,1],ssdata_t50_nn_B[,1]),1,min) # mosaic
mix = apply(cbind(ssdata_t50_nn_A[,2],ssdata_t50_nn_B[,2]),1,min) # mixture
rot = cbind(ssdata_t50_nn_A[,3],ssdata_t50_nn_B[,4])[cbind(seq(1,var_comb),1+(ssdata_t50_nn_A[,5]>ssdata_t50_nn_B[,5]))] # rotation
seq = cbind(ssdata_t50_nn_A[,5],ssdata_t50_nn_B[,5])[cbind(seq(1,var_comb),1+(ssdata_t50_nn_A[,5]>ssdata_t50_nn_B[,5]))] # sequence = negative control
mixmax = apply(cbind(ssdata_t50_nn_A[,6],ssdata_t50_nn_B[,6]),1,min) # maximum mixture = positive control
ssdata_first_nn = cbind(mos,mix,rot,seq,mixmax)
ssdata_first_nn[ssdata_first_nn==0] = 2000 # mark extinction as best outcome

# second to 50% frequency
mos = apply(cbind(ssdata_t50_nn_A[,1],ssdata_t50_nn_B[,1]),1,max) # mosaic
mix = apply(cbind(ssdata_t50_nn_A[,2],ssdata_t50_nn_B[,2]),1,max) # mixture
rot = cbind(ssdata_t50_nn_B[,3],ssdata_t50_nn_A[,4])[cbind(seq(1,var_comb),1+(ssdata_t50_nn_A[,5]>ssdata_t50_nn_B[,5]))] # rotation
seq = ssdata_t50_nn_A[,5]+ssdata_t50_nn_B[,5] # sequence = negative control
seq[seq>500&seq!=1500&seq!=2000] = 1000 # make like others so only measure within 500 generations window
mixmax = apply(cbind(ssdata_t50_nn_A[,6],ssdata_t50_nn_B[,6]),1,max) # maximum mixture = positive control
ssdata_second_nn = cbind(mos,mix,rot,seq,mixmax)
ssdata_second_nn[ssdata_second_nn==0] = 2000 # mark extinction as best outcome

# 80% population control
mos = ssdata_c80_nn[,1] # mosaic
mix = ssdata_c80_nn[,2] # mixture
rot = ssdata_c80_nn[,3:4][cbind(seq(1,var_comb),1+(ssdata_t50_nn_A[,5]>ssdata_t50_nn_B[,5]))] # rotation
seq = cbind(ssdata_c80_nn[,5],ssdata_c80_nn[,7])[cbind(seq(1,var_comb),1+(ssdata_t50_nn_A[,5]>ssdata_t50_nn_B[,5]))] # sequence = negative control
#seq = ssdata_c80_nn[,5]+ssdata_c80_nn[,7] # alternative more conservative negative control
mixmax = ssdata_c80_nn[,6] # maximum mixture = positive control
ssdata_ctrl_nn = cbind(mos,mix,rot,seq,mixmax)
ssdata_ctrl_nn[ssdata_ctrl_nn==0] = 2000 # mark extinction as best outcome

# initial control
mos = (1-simdata[,4]) + simdata[,4]*(1-0.5*(simdata[,7]+simdata[,13]))
mixk = (simdata[,7]+simdata[,13]-sqrt((((simdata[,7]+simdata[,13])^2)-2*simdata[,7]*simdata[,13]*(simdata[,7]+simdata[,13]))))/(2*simdata[,7]*simdata[,13])
mix = (1-simdata[,4]) + simdata[,4]*(1-mixk*simdata[,7])*(1-mixk*simdata[,13])
mixmax = (1-simdata[,4]) + simdata[,4]*(1-simdata[,7])*(1-simdata[,13])
rot = ((1-simdata[,4]) + simdata[,4]*(1-simdata[,c(7,13)]))[cbind(seq(1,var_comb),1+(ssdata_t50_nn_A[,5]>ssdata_t50_nn_B[,5]))]
seq = ((1-simdata[,4]) + simdata[,4]*(1-simdata[,c(7,13)]))[cbind(seq(1,var_comb),1+(ssdata_t50_nn_A[,5]>ssdata_t50_nn_B[,5]))]
ssdata_init_nn = 100*(1-cbind(mos,mix,rot,seq,mixmax))




#### mn data processing ####

# NB: for sequence and rotation, weakest goes first

# raw data
ssdata_t50_mn_A = cbind(sdata_c_mn[,1],sdata_x_mn[,1],sdata_r_mn[,1],sdata_r_revmn[,4],sdata_sA_mn[,1],sdata_psx_mn[,1])
ssdata_t50_mn_B = cbind(sdata_c_mn[,4],sdata_x_mn[,4],sdata_r_mn[,4],sdata_r_revmn[,1],sdata_sB_mn[,4],sdata_psx_mn[,4])
ssdata_c80_mn = cbind(sdata_c_mn[,7],sdata_x_mn[,7],sdata_r_mn[,7],sdata_r_revmn[,7],sdata_sA_mn[,7],sdata_psx_mn[,7],sdata_sB_mn[,7])

# if no measurement (which comes to mean weak selection for resistance by elimination)
ssdata_t50_mn_A[ssdata_t50_mn_A==0] = 1000
ssdata_t50_mn_B[ssdata_t50_mn_B==0] = 1000
ssdata_c80_mn[ssdata_c80_mn==0] = 1000

# if no measurement because of selection against resistance
bsdata_A = simdata[,c(6,6,6,12,6,6)] # starting frequency
bsdata_B = simdata[,c(12,12,12,6,12,12)] # starting frequency
psdata_mn_A = cbind(sdata_c_mn[,2],sdata_x_mn[,2],sdata_r_mn[,2],sdata_r_revmn[,5],sdata_sA_mn[,2],sdata_psx_mn[,2]) # frequency at t=250
psdata_mn_B = cbind(sdata_c_mn[,5],sdata_x_mn[,5],sdata_r_mn[,5],sdata_r_revmn[,2],sdata_sB_mn[,5],sdata_psx_mn[,5]) # frequency at t=250
ssdata_t50_mn_A[psdata_mn_A<bsdata_A] = 1500
ssdata_t50_mn_B[psdata_mn_B<bsdata_B] = 1500

# if no measurement because of extinction set as 'weakest'
esdata_mn = cbind(sdata_c_mn[,10],sdata_x_mn[,10],sdata_r_mn[,10],sdata_r_revmn[,10],sdata_sA_mn[,10],sdata_psx_mn[,10])
ssdata_t50_mn_A[esdata_mn>0] = 0
ssdata_t50_mn_B[esdata_mn>0] = 0
ssdata_c80_mn[cbind(esdata_mn,sdata_sB_mn[,10])>0] = 0

# calculate statistics
# first to 50% frequency
mos = apply(cbind(ssdata_t50_mn_A[,1],ssdata_t50_mn_B[,1]),1,min) # mosaic
mix = apply(cbind(ssdata_t50_mn_A[,2],ssdata_t50_mn_B[,2]),1,min) # mixture
rot = cbind(ssdata_t50_mn_A[,3],ssdata_t50_mn_B[,4])[cbind(seq(1,var_comb),1+(ssdata_t50_mn_A[,5]>ssdata_t50_mn_B[,5]))] # rotation
seq = cbind(ssdata_t50_mn_A[,5],ssdata_t50_mn_B[,5])[cbind(seq(1,var_comb),1+(ssdata_t50_mn_A[,5]>ssdata_t50_mn_B[,5]))] # sequence = negative control
mixmax = apply(cbind(ssdata_t50_mn_A[,6],ssdata_t50_mn_B[,6]),1,min) # maximum mixture = positive control
ssdata_first_mn = cbind(mos,mix,rot,seq,mixmax)
ssdata_first_mn[ssdata_first_mn==0] = 2000 # mark extinction as best outcome

# second to 50% frequency
mos = apply(cbind(ssdata_t50_mn_A[,1],ssdata_t50_mn_B[,1]),1,max) # mosaic
mix = apply(cbind(ssdata_t50_mn_A[,2],ssdata_t50_mn_B[,2]),1,max) # mixture
rot = cbind(ssdata_t50_mn_B[,3],ssdata_t50_mn_A[,4])[cbind(seq(1,var_comb),1+(ssdata_t50_mn_A[,5]>ssdata_t50_mn_B[,5]))] # rotation
seq = ssdata_t50_mn_A[,5]+ssdata_t50_mn_B[,5] # sequence = negative control
seq[seq>500&seq!=1500&seq!=2000] = 1000 # make like others so only measure within 500 generations window
mixmax = apply(cbind(ssdata_t50_mn_A[,6],ssdata_t50_mn_B[,6]),1,max) # maximum mixture = positive control
ssdata_second_mn = cbind(mos,mix,rot,seq,mixmax)
ssdata_second_mn[ssdata_second_mn==0] = 2000 # mark extinction as best outcome

# 80% population control
mos = ssdata_c80_mn[,1] # mosaic
mix = ssdata_c80_mn[,2] # mixture
rot = ssdata_c80_mn[,3:4][cbind(seq(1,var_comb),1+(ssdata_t50_mn_A[,5]>ssdata_t50_mn_B[,5]))] # rotation
seq = cbind(ssdata_c80_mn[,5],ssdata_c80_mn[,7])[cbind(seq(1,var_comb),1+(ssdata_t50_mn_A[,5]>ssdata_t50_mn_B[,5]))] # sequence = negative control
#seq = ssdata_c80_mn[,5]+ssdata_c80_mn[,7] # alternative more conservative negative control
mixmax = ssdata_c80_mn[,6] # maximum mixture = positive control
ssdata_ctrl_mn = cbind(mos,mix,rot,seq,mixmax)
ssdata_ctrl_mn[ssdata_ctrl_mn==0] = 2000 # mark extinction as best outcome

# initial control
mos = (1-simdata[,4]) + simdata[,4]*(1-0.5*(simdata[,7]+simdata[,13]))
mixk = (simdata[,7]+simdata[,13]-sqrt((((simdata[,7]+simdata[,13])^2)-2*simdata[,7]*simdata[,13]*(simdata[,7]+simdata[,13]))))/(2*simdata[,7]*simdata[,13])
mix = (1-simdata[,4]) + simdata[,4]*(1-mixk*simdata[,7])*(1-mixk*simdata[,13])
mixmax = (1-simdata[,4]) + simdata[,4]*(1-simdata[,7])*(1-simdata[,13])
rot = ((1-simdata[,4]) + simdata[,4]*(1-simdata[,c(7,13)]))[cbind(seq(1,var_comb),1+(ssdata_t50_mn_A[,5]>ssdata_t50_mn_B[,5]))]
seq = ((1-simdata[,4]) + simdata[,4]*(1-simdata[,c(7,13)]))[cbind(seq(1,var_comb),1+(ssdata_t50_mn_A[,5]>ssdata_t50_mn_B[,5]))]
ssdata_init_mn = 100*(1-cbind(mos,mix,rot,seq,mixmax))




#### mm data processing ####

# NB: for sequence and rotation, weakest goes first

# raw data
ssdata_t50_mm_A = cbind(sdata_c_mm[,1],sdata_x_mm[,1],sdata_r_mm[,1],sdata_r_revmm[,4],sdata_sA_mm[,1],sdata_psx_mm[,1])
ssdata_t50_mm_B = cbind(sdata_c_mm[,4],sdata_x_mm[,4],sdata_r_mm[,4],sdata_r_revmm[,1],sdata_sB_mm[,4],sdata_psx_mm[,4])
ssdata_c80_mm = cbind(sdata_c_mm[,7],sdata_x_mm[,7],sdata_r_mm[,7],sdata_r_revmm[,7],sdata_sA_mm[,7],sdata_psx_mm[,7],sdata_sB_mm[,7])

# if no measurement (which comes to mean weak selection for resistance by elimination)
ssdata_t50_mm_A[ssdata_t50_mm_A==0] = 1000
ssdata_t50_mm_B[ssdata_t50_mm_B==0] = 1000
ssdata_c80_mm[ssdata_c80_mm==0] = 1000

# if no measurement because of selection against resistance
bsdata_A = simdata[,c(6,6,6,12,6,6)] # starting frequency
bsdata_B = simdata[,c(12,12,12,6,12,12)] # starting frequency
psdata_mm_A = cbind(sdata_c_mm[,2],sdata_x_mm[,2],sdata_r_mm[,2],sdata_r_revmm[,5],sdata_sA_mm[,2],sdata_psx_mm[,2]) # frequency at t=250
psdata_mm_B = cbind(sdata_c_mm[,5],sdata_x_mm[,5],sdata_r_mm[,5],sdata_r_revmm[,2],sdata_sB_mm[,5],sdata_psx_mm[,5]) # frequency at t=250
ssdata_t50_mm_A[psdata_mm_A<bsdata_A] = 1500
ssdata_t50_mm_B[psdata_mm_B<bsdata_B] = 1500

# if no measurement because of extinction set as 'weakest'
esdata_mm = cbind(sdata_c_mm[,10],sdata_x_mm[,10],sdata_r_mm[,10],sdata_r_revmm[,10],sdata_sA_mm[,10],sdata_psx_mm[,10])
ssdata_t50_mm_A[esdata_mm>0] = 0
ssdata_t50_mm_B[esdata_mm>0] = 0
ssdata_c80_mm[cbind(esdata_mm,sdata_sB_mm[,10])>0] = 0

# calculate statistics
# first to 50% frequency
mos = apply(cbind(ssdata_t50_mm_A[,1],ssdata_t50_mm_B[,1]),1,min) # mosaic
mix = apply(cbind(ssdata_t50_mm_A[,2],ssdata_t50_mm_B[,2]),1,min) # mixture
rot = cbind(ssdata_t50_mm_A[,3],ssdata_t50_mm_B[,4])[cbind(seq(1,var_comb),1+(ssdata_t50_mm_A[,5]>ssdata_t50_mm_B[,5]))] # rotation
seq = cbind(ssdata_t50_mm_A[,5],ssdata_t50_mm_B[,5])[cbind(seq(1,var_comb),1+(ssdata_t50_mm_A[,5]>ssdata_t50_mm_B[,5]))] # sequence = negative control
mixmax = apply(cbind(ssdata_t50_mm_A[,6],ssdata_t50_mm_B[,6]),1,min) # maximum mixture = positive control
ssdata_first_mm = cbind(mos,mix,rot,seq,mixmax)
ssdata_first_mm[ssdata_first_mm==0] = 2000 # mark extinction as best outcome

# second to 50% frequency
mos = apply(cbind(ssdata_t50_mm_A[,1],ssdata_t50_mm_B[,1]),1,max) # mosaic
mix = apply(cbind(ssdata_t50_mm_A[,2],ssdata_t50_mm_B[,2]),1,max) # mixture
rot = cbind(ssdata_t50_mm_B[,3],ssdata_t50_mm_A[,4])[cbind(seq(1,var_comb),1+(ssdata_t50_mm_A[,5]>ssdata_t50_mm_B[,5]))] # rotation
seq = ssdata_t50_mm_A[,5]+ssdata_t50_mm_B[,5] # sequence = negative control
seq[seq>500&seq!=1500&seq!=2000] = 1000 # make like others so only measure within 500 generations window
mixmax = apply(cbind(ssdata_t50_mm_A[,6],ssdata_t50_mm_B[,6]),1,max) # maximum mixture = positive control
ssdata_second_mm = cbind(mos,mix,rot,seq,mixmax)
ssdata_second_mm[ssdata_second_mm==0] = 2000 # mark extinction as best outcome

# 80% population control
mos = ssdata_c80_mm[,1] # mosaic
mix = ssdata_c80_mm[,2] # mixture
rot = ssdata_c80_mm[,3:4][cbind(seq(1,var_comb),1+(ssdata_t50_mm_A[,5]>ssdata_t50_mm_B[,5]))] # rotation
seq = cbind(ssdata_c80_mm[,5],ssdata_c80_mm[,7])[cbind(seq(1,var_comb),1+(ssdata_t50_mm_A[,5]>ssdata_t50_mm_B[,5]))] # sequence = negative control
#seq = ssdata_c80_mm[,5]+ssdata_c80_mm[,7] # alternative more conservative negative control
mixmax = ssdata_c80_mm[,6] # maximum mixture = positive control
ssdata_ctrl_mm = cbind(mos,mix,rot,seq,mixmax)
ssdata_ctrl_mm[ssdata_ctrl_mm==0] = 2000 # mark extinction as best outcome

# initial control
mos = (1-simdata[,4]) + simdata[,4]*(1-0.5*(simdata[,7]+simdata[,13]))
mixk = (simdata[,7]+simdata[,13]-sqrt((((simdata[,7]+simdata[,13])^2)-2*simdata[,7]*simdata[,13]*(simdata[,7]+simdata[,13]))))/(2*simdata[,7]*simdata[,13])
mix = (1-simdata[,4]) + simdata[,4]*(1-mixk*simdata[,7])*(1-mixk*simdata[,13])
mixmax = (1-simdata[,4]) + simdata[,4]*(1-simdata[,7])*(1-simdata[,13])
rot = ((1-simdata[,4]) + simdata[,4]*(1-simdata[,c(7,13)]))[cbind(seq(1,var_comb),1+(ssdata_t50_mm_A[,5]>ssdata_t50_mm_B[,5]))]
seq = ((1-simdata[,4]) + simdata[,4]*(1-simdata[,c(7,13)]))[cbind(seq(1,var_comb),1+(ssdata_t50_mm_A[,5]>ssdata_t50_mm_B[,5]))]
ssdata_init_mm = 100*(1-cbind(mos,mix,rot,seq,mixmax))




#### data divisions ####

# % that mt insecticide fails before nu
for(i in 1:5){
  print(paste("Mitochondrial resistance allele first-to-break for ",c("Mosaic","Mixture","Rotation","Sequence","MaxMixture")[i],": ",100*sum(ssdata_first_mn[,i]==ssdata_t50_mn_A[,c(1,2,3,5,6)][,i])/var_comb,"%",sep=""))
}



# % runs with measurement and if not why not
comp_type = matrix(0,5*9,8)
colnames(comp_type) = c("Inheritance","Statistic","Strategy",
                        "SimulationMeasurement","WeakSelectionForResistance","SelectionAgainstResistance","PopulationExtinction","Total")
comp_type[,1] = rep(rep(c("nn", "mn", "mm"),3),5)
comp_type[,2] = rep(rep(c("first","second","ctrl"),each=3),5)
# give table output
for(i in 1:9){for(j in 1:5){
  # data choice
  if(comp_type[i+(j-1)*9,1]=="nn"){
    if(comp_type[i+(j-1)*9,2]=="first"){ssdata = ssdata_first_nn} # time to first 50% frequency
    if(comp_type[i+(j-1)*9,2]=="second"){ssdata = ssdata_second_nn} # time to first 50% frequency
    if(comp_type[i+(j-1)*9,2]=="ctrl"){ssdata = ssdata_ctrl_nn} # time to 80% female population recovery
  }
  if(comp_type[i+(j-1)*9,1]=="mn"){
    if(comp_type[i+(j-1)*9,2]=="first"){ssdata = ssdata_first_mn} # time to first 50% frequency
    if(comp_type[i+(j-1)*9,2]=="second"){ssdata = ssdata_second_mn} # time to first 50% frequency
    if(comp_type[i+(j-1)*9,2]=="ctrl"){ssdata = ssdata_ctrl_mn} # time to 80% female population recovery
  }
  if(comp_type[i+(j-1)*9,1]=="mm"){
    if(comp_type[i+(j-1)*9,2]=="first"){ssdata = ssdata_first_mm} # time to first 50% frequency
    if(comp_type[i+(j-1)*9,2]=="second"){ssdata = ssdata_second_mm} # time to first 50% frequency
    if(comp_type[i+(j-1)*9,2]=="ctrl"){ssdata = ssdata_ctrl_mm} # time to 80% female population recovery
  }
  if(i<7){comp_type[i+(j-1)*9,3:8] = c(c("C","X","R","S","Z")[j],sum(ssdata[,j]<501),sum(ssdata[,j]==1000),sum(ssdata[,j]==1500),sum(ssdata[,j]==2000),1000000)}
  if(i>6){comp_type[i+(j-1)*9,3:8] = c(c("C","X","R","S","Z")[j],sum(ssdata[,j]<501)-sum(ssdata[,j]==1),sum(ssdata[,j]==1000),sum(ssdata[,j]==1),sum(ssdata[,j]==2000),1000000)}
}}
# save data
write.csv(comp_type,paste("Table_measurementtype.csv",sep=""),row.names=F)




#### barplot waterfall alternative metrics ####

# data
inher = rep(c("NN","MN","MM"),15)
stati = rep(rep(c("1st","2nd","ctrl"),each=3),5)
strat = rep(c("Mos","Mix","Rot","Seq","Max"),each=9)
xdata = inher
#for(i in 1:45){xdata[i] = paste(xdata[i],stati[i],sep="/")}
for(i in 1:45){xdata[i] = paste(xdata[i],strat[i],sep="/")}
ydata = matrix(as.numeric(comp_type[,4:7]),5*9,4)/(10^4)

# reorder data
x1data = c(xdata[stati=="1st"&inher=="NN"],xdata[stati=="1st"&inher=="MN"],xdata[stati=="1st"&inher=="MM"])
y1data = rbind(ydata[stati=="1st"&inher=="NN",],ydata[stati=="1st"&inher=="MN",],ydata[stati=="1st"&inher=="MM",])
x2data = c(xdata[stati=="2nd"&inher=="NN"],xdata[stati=="2nd"&inher=="MN"],xdata[stati=="2nd"&inher=="MM"])
y2data = rbind(ydata[stati=="2nd"&inher=="NN",],ydata[stati=="2nd"&inher=="MN",],ydata[stati=="2nd"&inher=="MM",])
xcdata = c(xdata[stati=="ctrl"&inher=="NN"],xdata[stati=="ctrl"&inher=="MN"],xdata[stati=="ctrl"&inher=="MM"])
ycdata = rbind(ydata[stati=="ctrl"&inher=="NN",],ydata[stati=="ctrl"&inher=="MN",],ydata[stati=="ctrl"&inher=="MM",])
# reorder rows to SRCXZ
x1data = x1data[c(4,3,1,2,5, 9,8,6,7,10, 14,13,11,12,15)]
y1data = y1data[c(4,3,1,2,5, 9,8,6,7,10, 14,13,11,12,15),]
x2data = x2data[c(4,3,1,2,5, 9,8,6,7,10, 14,13,11,12,15)]
y2data = y2data[c(4,3,1,2,5, 9,8,6,7,10, 14,13,11,12,15),]
xcdata = xcdata[c(4,3,1,2,5, 9,8,6,7,10, 14,13,11,12,15)]
ycdata = ycdata[c(4,3,1,2,5, 9,8,6,7,10, 14,13,11,12,15),]
# reverse rows for plotting in reverse
x1data = rev(x1data)
y1data = y1data[15:1,]
x2data = rev(x2data)
y2data = y2data[15:1,]
xcdata = rev(xcdata)
ycdata = ycdata[15:1,]

# plot
pdf("waterfall_BP_all.pdf",width=8*3.5,height=7.5)
par(mar=c(5,10,2.5,5),mfrow=c(1,4)) # set margins and layout
layout(matrix(c(1,1,2,2,3,3,4),1,7))

# 1st
y1data_w = y1data
y1data_w[c(1,5,6,10,11,15)] = NA
barplot(names.arg=x1data,horiz=T,t(y1data_w),las=1,xlab="Percentage of Simulated Runs (%)",main="First-to-Break",col=plasma(4),cex.axis=1.5,cex.names=1.5,cex.lab=2,cex.main=2)
barplot(names.arg=x1data,horiz=T,t(y1data),las=1,xlab="Percentage of Simulated Runs (%)",main="First-to-Break",col=plasma(4,alpha=0.5),cex.axis=1.5,cex.names=1.5,cex.lab=2,cex.main=2,add=T)
text(-18,9,"Inheritance / Strategy",xpd=NA,cex=2,srt=90)
text(-20,19,"A",font=2,xpd=NA,cex=4)
points(c(-14,100),c(6.1,6.1),type="l",lwd=2,xpd=T)
points(c(-14,100),c(12.1,12.1),type="l",lwd=2,xpd=T)

# 2nd
y2data_w = y2data
y2data_w[c(1,5,6,10,11,15)] = NA
barplot(names.arg=x2data,horiz=T,t(y2data_w),las=1,xlab="Percentage of Simulated Runs (%)",main="Second-to-Break",col=plasma(4),cex.axis=1.5,cex.names=1.5,cex.lab=2,cex.main=2)
barplot(names.arg=x2data,horiz=T,t(y2data),las=1,xlab="Percentage of Simulated Runs (%)",main="Second-to-Break",col=plasma(4,alpha=0.5),cex.axis=1.5,cex.names=1.5,cex.lab=2,cex.main=2,add=T)
text(-18,9,"Inheritance / Strategy",xpd=NA,cex=2,srt=90)
text(-20,19,"B",font=2,xpd=NA,cex=4)
points(c(-14,100),c(6.1,6.1),type="l",lwd=2,xpd=T)
points(c(-14,100),c(12.1,12.1),type="l",lwd=2,xpd=T)

# ctrl
ycdata_w = ycdata
ycdata_w[c(1,5,6,10,11,15)] = NA
barplot(names.arg=xcdata,horiz=T,t(ycdata_w),las=1,xlab="Percentage of Simulated Runs (%)",main="Control Failure",col=plasma(4),cex.axis=1.5,cex.names=1.5,cex.lab=2,cex.main=2)
barplot(names.arg=xcdata,horiz=T,t(ycdata),las=1,xlab="Percentage of Simulated Runs (%)",main="Control Failure",col=plasma(4,alpha=0.5),cex.axis=1.5,cex.names=1.5,cex.lab=2,cex.main=2,add=T)
text(-18,9,"Inheritance / Strategy",xpd=NA,cex=2,srt=90)
text(-20,19,"C",font=2,xpd=NA,cex=4)
points(c(-14,100),c(6.1,6.1),type="l",lwd=2,xpd=T)
points(c(-14,100),c(12.1,12.1),type="l",lwd=2,xpd=T)

# legend
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim=c(0,1),ylim=c(0,1))
legend(-0.3,0.5,c("Successful Measurement","Toward Threshold","Away from Threshold","Extinction"),title="1 million simulated runs:",pch=rep(15,4),col=plasma(4),bg="white",cex=2,yjust=0.5,xpd=NA)
legend(-0.3,0.1,c("NN/Seq/Away = 83","MN/Seq/Away = 49", "MM/Seq/Away = 25"),title="NB: Second-to-Break",pch=rep(15,3),col=plasma(4,alpha=0.5)[3],bg="white",cex=2,yjust=0.5,xpd=NA)
# "Mitochondrial Only","Mixed","Nuclear Only", & c(7,22,12,

dev.off()

# add texture
#barplot(names.arg=NA,c(rep(100,5),rep(0,10)),horiz=T,angle=90,density=17,col="black",add=T)
#barplot(names.arg=NA,c(rep(100,5),rep(0,10)),horiz=T,angle=0,density=17,col="black",add=T)
#barplot(names.arg=NA,c(rep(0,10),rep(100,5)),horiz=T,angle=45,density=17,col="black",add=T)
#barplot(names.arg=NA,c(rep(0,10),rep(100,5)),horiz=T,angle=135,density=17,col="black",add=T)




#### data visualisation : violin plot ####

# data type/choice, data filtering, statistic calculation and statistical comparison
id = c("nn", "mn", "mm")

for(i in 1:length(id)){

  # parameter
  dataset0 = id[i]

  # data choice
  if(dataset0=="nn"){
    ssdata_first = ssdata_first_nn # time to first 50% frequency
    ssdata_second = ssdata_second_nn # time to first 50% frequency
    ssdata_ctrl = ssdata_ctrl_nn # time to 80% female population recovery
    ssdata_init = ssdata_init_nn # initial control %
  }
  if(dataset0=="mn"){
    ssdata_first = ssdata_first_mn # time to first 50% frequency
    ssdata_second = ssdata_second_mn # time to first 50% frequency
    ssdata_ctrl = ssdata_ctrl_mn # time to 80% female population recovery
    ssdata_init = ssdata_init_mn # initial control %
  }
  if(dataset0=="mm"){
    ssdata_first = ssdata_first_mm # time to first 50% frequency
    ssdata_second = ssdata_second_mm # time to first 50% frequency
    ssdata_ctrl = ssdata_ctrl_mm # time to 80% female population recovery
    ssdata_init = ssdata_init_mm # initial control %
  }

  # set all data that is >500 to NA to exlude it from the plot NB: no need to do this for initial control
  ssdata_first[ssdata_first>500] = NA
  ssdata_second[ssdata_second>500] = NA
  ssdata_ctrl[ssdata_ctrl>500] = NA

  # rename columns for plotting
  colnames(ssdata_first) = c("Mosaic","Mixture","Rotation","-ve Control","+ve Control")
  colnames(ssdata_second) = c("Mosaic","Mixture","Rotation","-ve Control","+ve Control")
  colnames(ssdata_ctrl) = c("Mosaic","Mixture","Rotation","-ve Control","+ve Control")
  colnames(ssdata_init) = c("Mosaic","Mixture","Rotation","-ve Control","+ve Control")

  # plot
  pdf(paste("violin_DV_",dataset0,".pdf",sep=""),width=8*2,height=2*7.5)
  par(mar=c(2.5,5,5,5),mfrow=c(2,1)) # set margins and layout

  # plot t50
  vioplot(ssdata_first,ylab="",las=1,horizontal=F,col=rep(viridis(4),2),plotCentre="line",side="left",ylim=c(0,500))
  vioplot(ssdata_second,col=rep(viridis(4),2),plotCentre="line",side="right",add=T)
  text(-1.2,600,"A",font=2,xpd=NA,cex=2)
  mtext("Time to Resistance (generations)", side=2, line=3,cex=1.5)

  # plot c80
  vioplot(ssdata_ctrl,ylab="",las=1,horizontal=F,col=rep(viridis(4),2),plotCentre="line",side="left",ylim=c(0,500))
  vioplot(5*ssdata_init,col=rep(viridis(4),2),plotCentre="line",side="right",add=T)
  text(-1.2,600,"B",font=2,xpd=NA,cex=2)
  mtext("Time to Control Failure (generations)", side=2, line=3,cex=1.5)
  # add on second y-axis
  axis(side=4, at=seq(0,500,length.out=6), labels=seq(0,1,length.out=6), las=1)
  mtext("Initial Control (%)", side=4, line=3,cex=1.5)

  dev.off()

}



## violin panels

# data choice
ssdata_fnn = ssdata_first_nn # time to first 50% frequency
ssdata_snn = ssdata_second_nn # time to first 50% frequency
ssdata_fmn = ssdata_first_mn # time to first 50% frequency
ssdata_smn = ssdata_second_mn # time to first 50% frequency
ssdata_fmm = ssdata_first_mm # time to first 50% frequency
ssdata_smm = ssdata_second_mm # time to first 50% frequency

# set all data that is >500 to NA to exclude it from the plot
ssdata_fnn[ssdata_fnn>500] = NA
ssdata_snn[ssdata_snn>500] = NA
ssdata_fmn[ssdata_fmn>500] = NA
ssdata_smn[ssdata_smn>500] = NA
ssdata_fmm[ssdata_fmm>500] = NA
ssdata_smm[ssdata_smm>500] = NA

# rename columns for plotting
colnames(ssdata_fnn) = c("Mosaic","Mixture","Rotation","Sequence","Maximum")
colnames(ssdata_snn) = c("Mosaic","Mixture","Rotation","Sequence","Maximum")
colnames(ssdata_fmn) = c("Mosaic","Mixture","Rotation","Sequence","Maximum")
colnames(ssdata_smn) = c("Mosaic","Mixture","Rotation","Sequence","Maximum")
colnames(ssdata_fmm) = c("Mosaic","Mixture","Rotation","Sequence","Maximum")
colnames(ssdata_smm) = c("Mosaic","Mixture","Rotation","Sequence","Maximum")

# reorder data
ssdata_fnn = ssdata_fnn[,c(4,3,1,2,5)]
ssdata_snn = ssdata_snn[,c(4,3,1,2,5)]
ssdata_fmn = ssdata_fmn[,c(4,3,1,2,5)]
ssdata_smn = ssdata_smn[,c(4,3,1,2,5)]
ssdata_fmm = ssdata_fmm[,c(4,3,1,2,5)]
ssdata_smm = ssdata_smm[,c(4,3,1,2,5)]

# plot
pdf("violin_PV_all.pdf",width=8*2,height=3*7.5)
par(mar=c(5,10,8,2.5),mfrow=c(3,1)) # set margins and layout

# plot nn
vioplot(ssdata_fnn,ylab="",las=1,horizontal=F,col=c("white",viridis(5)[2:4],"white"),plotCentre="line",side="left",ylim=c(0,500),main="Nuclear Inheritance Only (NN)",cex.lab=1.5,cex.main=3,cex.names=2.5,cex.axis=2)
vioplot(ssdata_fnn,col=viridis(5,alpha=0.65),border=c("grey",rep("black",3),"grey"),rectCol=c("grey",rep("black",3),"grey"),lineCol=c("grey",rep("black",3),"grey"),plotCentre="line",side="left",add=T)
vioplot(ssdata_snn,col=c("white",viridis(5)[2:4],"white"),plotCentre="line",side="right",add=T)
vioplot(ssdata_snn,col=viridis(5,alpha=0.65),rectCol=c("grey",rep("black",3),"grey"),lineCol=c("grey",rep("black",3),"grey"),plotCentre="line",side="right",add=T)
text(-0.1,575,"A",font=2,xpd=NA,cex=4)
mtext("Time to Resistance (generations)", side=2, line=6,cex=2)
text(c(0.75,1.75,2.75,3.75,4.75),c(60,70,80,90,80),"1st",cex=2,col="white") # maybe white into violin itself?
text(c(1.25,2.25,3.25,4.25,5.25),c(140,180,170,200,180),"2nd",cex=2,col="white")
#axis(side=1,at=seq(1,5),labels=c("Mosaic","Mixture","Rotation","Sequence","Maximum"),cex=2)

# plot mn
vioplot(ssdata_fmn,ylab="",las=1,horizontal=F,col=c("white",viridis(5)[2:4],"white"),plotCentre="line",side="left",ylim=c(0,500),main="Mixed Inheritances (MN)",cex.lab=1.5,cex.main=3,cex.names=2.5,cex.axis=2)
vioplot(ssdata_fmn,col=viridis(5,alpha=0.5),plotCentre="line",side="left",add=T)
vioplot(ssdata_smn,col=c("white",viridis(5)[2:4],"white"),plotCentre="line",side="right",add=T)
vioplot(ssdata_smn,col=viridis(5,alpha=0.5),plotCentre="line",side="right",add=T)
text(-0.1,575,"B",font=2,xpd=NA,cex=4)
mtext("Time to Resistance (generations)", side=2, line=6,cex=2)
text(c(0.75,1.75,2.75,3.75,4.75),c(30,50,70,80,50),"1st",cex=2,col="white") # maybe white into violin itself?
text(c(1.25,2.25,3.25,4.25,5.25),c(110,150,140,150,130),"2nd",cex=2,col="white")

# plot mm
vioplot(ssdata_fmm,ylab="",las=1,horizontal=F,col=c("white",viridis(5)[2:4],"white"),plotCentre="line",side="left",ylim=c(0,500),main="Mitochondrial Inheritance Only (MM)",cex.lab=1.5,cex.main=3,cex.names=2.5,cex.axis=2)
vioplot(ssdata_fmm,col=viridis(5,alpha=0.5),plotCentre="line",side="left",add=T)
vioplot(ssdata_smm,col=c("white",viridis(5)[2:4],"white"),plotCentre="line",side="right",add=T)
vioplot(ssdata_smm,col=viridis(5,alpha=0.5),plotCentre="line",side="right",add=T)
text(-0.1,575,"C",font=2,xpd=NA,cex=4)
mtext("Time to Resistance (generations)", side=2, line=6,cex=2)
text(c(0.75,1.75,2.75,3.75,4.75),c(70,60,50,60,40),"1st",cex=2,col="white") # maybe white into violin itself?
text(c(1.25,2.25,3.25,4.25,5.25),c(80,120,100,100,80),"2nd",cex=2,col="white")

dev.off()




#### data visualisation : pairs (regression and measurement) ####

# data type/choice, data filtering, statistic calculation and statistical comparison
dataset0 = rep(c("nn", "mn", "mm"),4)
dataset3 = rep(c("first","second","ctrl","init"),each=3)
dataset4 = c("C","X","R","neg","pos") # NB: cycle through later
id = cbind(dataset0,dataset3)

for(i in 1:nrow(id)){

  # parameter
  dataset0 = id[i,1]
  dataset3 = id[i,2]

  # data choice
  if(dataset0=="nn"){
    if(dataset3=="first"){ssdata = ssdata_first_nn} # time to first 50% frequency
    if(dataset3=="second"){ssdata = ssdata_second_nn} # time to first 50% frequency
    if(dataset3=="ctrl"){ssdata = ssdata_ctrl_nn} # time to 80% female population recovery
    if(dataset3=="init"){ssdata = ssdata_init_nn} # initial control %
  }
  if(dataset0=="mn"){
    if(dataset3=="first"){ssdata = ssdata_first_mn} # time to first 50% frequency
    if(dataset3=="second"){ssdata = ssdata_second_mn} # time to first 50% frequency
    if(dataset3=="ctrl"){ssdata = ssdata_ctrl_mn} # time to 80% female population recovery
    if(dataset3=="init"){ssdata = ssdata_init_mn} # initial control %
  }
  if(dataset0=="mm"){
    if(dataset3=="first"){ssdata = ssdata_first_mm} # time to first 50% frequency
    if(dataset3=="second"){ssdata = ssdata_second_mm} # time to first 50% frequency
    if(dataset3=="ctrl"){ssdata = ssdata_ctrl_mm} # time to 80% female population recovery
    if(dataset3=="init"){ssdata = ssdata_init_mm} # initial control %
  }

  # set all data that is >500 to 500 to include it in a minimal way in this plot NB: no need to do this for initial control
  ssdata[ssdata>500] = 500

  # plot
  for(j in 1:5){
    # plot pair regression of parameter vs statistic
    png(paste("regression_DR_",dataset0,"_",dataset3,"_",dataset4[j],".png",sep=""),width=250*6,height=250*3)
    par(mfrow=c(3,6))
    for(k in 1:17){
      if(sum(k==c(1,2,6,10,12,16))>0){
        plot(log(simdata[ssdata[,j]>0,k],10),ssdata[ssdata[,j]>0,j],xlab=colnames(simdata)[k],ylab=ifelse(dataset3=="first","Time of First to 50% Frequency",ifelse(dataset3=="ctrl","Time to 80% Population Recovery",ifelse(dataset3=="second","Time of Second to 50% Frequency","Initial Control in First Generation"))),pch=".",font.lab=2)
        #points(log(simdata[ssdata[,j]==0,k],10),ssdata[ssdata[,j]==0,j],pch=".",col="red")
        abline(lm(ssdata[ssdata[,j]>0,j]~log(simdata[ssdata[,j]>0,k],10)),col="blue")
      }
      else{
        plot(simdata[ssdata[,j]>0,k],ssdata[ssdata[,j]>0,j],xlab=colnames(simdata)[k],ylab=ifelse(dataset3=="first","Time of First to 50% Frequency",ifelse(dataset3=="ctrl","Time to 80% Population Recovery",ifelse(dataset3=="second","Time of Second to 50% Frequency","Initial Control in First Generation"))),pch=".",font.lab=2)
        #points(simdata[ssdata[,j]==0,k],ssdata[ssdata[,j]==0,j],pch=".",col="red")
        abline(lm(ssdata[ssdata[,j]>0,j]~simdata[ssdata[,j]>0,k]),col="blue")
      }
      if(k==5){plot.new()} # blank
    }
    dev.off()

    # plot pair measurement of parameter vs parameter
    png(paste("regression_DM_",dataset0,"_",dataset3,"_",dataset4[j],".png",sep=""),width=250*17,height=250*17)
    par(mfrow=c(17,17))
    for(k in 1:17){
      if(sum(k==c(1,2,6,10,12,16))>0){ y = log(simdata[,k],10) }
      if(sum(k==c(1,2,6,10,12,16))==0){ y = simdata[,k] }
      for(l in 1:k){
        if(sum(l==c(1,2,6,10,12,16))>0){ x = log(simdata[,l],10) }
        if(sum(l==c(1,2,6,10,12,16))==0){ x = simdata[,l] }
        # subset data by valid or invalid values
        xb = x[ssdata[,j]!=500]
        xr = x[ssdata[,j]==500]
        yb = y[ssdata[,j]!=500]
        yr = y[ssdata[,j]==500]
        if(length(xb)<5000){sb = sample(seq(1,length(xb)),length(xb),replace=F)}
        if(length(xb)>=5000){sb = sample(seq(1,length(xb)),5000,replace=F)}
        if(length(xr)<5000){sr = sample(seq(1,length(xr)),length(xr),replace=F)}
        if(length(xr)>=5000){sr = sample(seq(1,length(xr)),5000,replace=F)}
        # plot
        plot(xb[sb],yb[sb],pch=".",xlab=colnames(simdata)[l],ylab=colnames(simdata)[k],font.lab=2,col=viridis(2,alpha=0.3)[1])
        points(xr[sr],yr[sr],col=viridis(2,alpha=0.3)[2],pch=".")
      }
      # blank plot
      if(k<17){ for(l in 1:(17-k)){plot.new()} }
    }
    dev.off()
  }
}




#### conditional inference decision tree ####

# tables to save
comp_mix_all = matrix(0,12,20)
colnames(comp_mix_all) = c("Inheritance","Statistic","Comparison","C","X","R","S","CX","CR","CS","XR","XS","RS","CXR","CXS","CRS","XRS","=","NA","Total")
comp_mix_other = matrix(0,48,8)
colnames(comp_mix_other) = c("Inheritance","Statistic","Comparison","Focal","O","J","=","Total")
comp_nomix_all = matrix(0,12,12)
colnames(comp_nomix_all) = c("Inheritance","Statistic","Comparison","C","R","S","CR","CS","RS","=","NA","Total")
comp_nomix_other = matrix(0,36,8)
colnames(comp_nomix_other) = c("Inheritance","Statistic","Comparison","Focal","O","J","=","Total")

# data type/choice, data filtering, statistic calculation and statistical comparison
# with mix
dataset0m = rep(rep(c("nn", "mn", "mm"),4),5)
dataset3m = rep(rep(c("first","second","ctrl","init"),each=3),5)
dataset4m = rep(c("all","C","X","R","S"),each=12)
dataset5m = rep("mix",length(dataset0m))
# without mix
dataset0n = rep(rep(c("nn", "mn", "mm"),4),4)
dataset3n = rep(rep(c("first","second","ctrl","init"),each=3),4)
dataset4n = rep(c("all","C","R","S"),each=12)
dataset5n = rep("nomix",length(dataset0n))
# assemble
id = cbind(c(dataset0m,dataset0n),c(dataset3m,dataset3n),c(dataset4m,dataset4n),c(dataset5m,dataset5n))

for(i in 1:nrow(id)){

  # parameter
  dataset0 = id[i,1]
  dataset3 = id[i,2]
  dataset4 = id[i,3]
  dataset5 = id[i,4]

  # data choice
  if(dataset0=="nn"){
    if(dataset3=="first"){ssdata = ssdata_first_nn[,1:4]} # time to first 50% frequency
    if(dataset3=="second"){ssdata = ssdata_second_nn[,1:4]} # time to first 50% frequency
    if(dataset3=="ctrl"){ssdata = ssdata_ctrl_nn[,1:4]} # time to 80% female population recovery
    if(dataset3=="init"){ssdata = ssdata_init_nn[,1:4]} # initial control %
  }
  if(dataset0=="mn"){
    if(dataset3=="first"){ssdata = ssdata_first_mn[,1:4]} # time to first 50% frequency
    if(dataset3=="second"){ssdata = ssdata_second_mn[,1:4]} # time to first 50% frequency
    if(dataset3=="ctrl"){ssdata = ssdata_ctrl_mn[,1:4]} # time to 80% female population recovery
    if(dataset3=="init"){ssdata = ssdata_init_mn[,1:4]} # initial control %
  }
  if(dataset0=="mm"){
    if(dataset3=="first"){ssdata = ssdata_first_mm[,1:4]} # time to first 50% frequency
    if(dataset3=="second"){ssdata = ssdata_second_mm[,1:4]} # time to first 50% frequency
    if(dataset3=="ctrl"){ssdata = ssdata_ctrl_mm[,1:4]} # time to 80% female population recovery
    if(dataset3=="init"){ssdata = ssdata_init_mm[,1:4]} # initial control %
  }

  # statistical comparison
  if(dataset4=="all"){ # plus and minus 10% difference clear for mixture/sequence to be favoured with near-equal category as well
    # null vector
    y = rep("NA",nrow(ssdata))
    # classification: sequence = S, mosaic = C, mixture = X, rotation = R
    if(dataset5=="mix"){
      y[ssdata[,1]>(1.1*apply(ssdata[,-1],1,max))&(0.9*ssdata[,1])>apply(ssdata[,-1],1,max)] = "C"
      y[ssdata[,2]>(1.1*apply(ssdata[,-2],1,max))&(0.9*ssdata[,2])>apply(ssdata[,-2],1,max)] = "X"
      y[ssdata[,3]>(1.1*apply(ssdata[,-3],1,max))&(0.9*ssdata[,3])>apply(ssdata[,-3],1,max)] = "R"
      y[ssdata[,4]>(1.1*apply(ssdata[,-4],1,max))&(0.9*ssdata[,4])>apply(ssdata[,-4],1,max)] = "S"

      y[!(ssdata[,2]>(1.1*ssdata[,1])&(0.9*ssdata[,2])>ssdata[,1]) & !(ssdata[,1]>(1.1*ssdata[,2])&(0.9*ssdata[,1])>ssdata[,2]) &
          ssdata[,1]>(1.1*apply(ssdata[,-c(1,2)],1,max))&(0.9*ssdata[,1])>apply(ssdata[,-c(1,2)],1,max) &
          ssdata[,2]>(1.1*apply(ssdata[,-c(1,2)],1,max))&(0.9*ssdata[,2])>apply(ssdata[,-c(1,2)],1,max)] = "CX"
      y[!(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) &
          ssdata[,1]>(1.1*apply(ssdata[,-c(1,3)],1,max))&(0.9*ssdata[,1])>apply(ssdata[,-c(1,3)],1,max) &
          ssdata[,3]>(1.1*apply(ssdata[,-c(1,3)],1,max))&(0.9*ssdata[,3])>apply(ssdata[,-c(1,3)],1,max)] = "CR"
      y[!(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) &
          ssdata[,1]>(1.1*apply(ssdata[,-c(1,4)],1,max))&(0.9*ssdata[,1])>apply(ssdata[,-c(1,4)],1,max) &
          ssdata[,4]>(1.1*apply(ssdata[,-c(1,4)],1,max))&(0.9*ssdata[,4])>apply(ssdata[,-c(1,4)],1,max)] = "CS"
      y[!(ssdata[,3]>(1.1*ssdata[,2])&(0.9*ssdata[,3])>ssdata[,2]) & !(ssdata[,2]>(1.1*ssdata[,3])&(0.9*ssdata[,2])>ssdata[,3]) &
          ssdata[,2]>(1.1*apply(ssdata[,-c(2,3)],1,max))&(0.9*ssdata[,2])>apply(ssdata[,-c(2,3)],1,max) &
          ssdata[,3]>(1.1*apply(ssdata[,-c(2,3)],1,max))&(0.9*ssdata[,3])>apply(ssdata[,-c(2,3)],1,max)] = "XR"
      y[!(ssdata[,4]>(1.1*ssdata[,2])&(0.9*ssdata[,4])>ssdata[,2]) & !(ssdata[,2]>(1.1*ssdata[,4])&(0.9*ssdata[,2])>ssdata[,4]) &
          ssdata[,2]>(1.1*apply(ssdata[,-c(2,4)],1,max))&(0.9*ssdata[,2])>apply(ssdata[,-c(2,4)],1,max) &
          ssdata[,4]>(1.1*apply(ssdata[,-c(2,4)],1,max))&(0.9*ssdata[,4])>apply(ssdata[,-c(2,4)],1,max)] = "XS"
      y[!(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3]) & !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) &
          ssdata[,3]>(1.1*apply(ssdata[,-c(3,4)],1,max))&(0.9*ssdata[,3])>apply(ssdata[,-c(3,4)],1,max) &
          ssdata[,4]>(1.1*apply(ssdata[,-c(3,4)],1,max))&(0.9*ssdata[,4])>apply(ssdata[,-c(3,4)],1,max)] = "RS"

      y[!(ssdata[,2]>(1.1*ssdata[,1])&(0.9*ssdata[,2])>ssdata[,1]) & !(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,2]>(1.1*ssdata[,3])&(0.9*ssdata[,2])>ssdata[,3]) &
          !(ssdata[,1]>(1.1*ssdata[,2])&(0.9*ssdata[,1])>ssdata[,2]) & !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) & !(ssdata[,3]>(1.1*ssdata[,2])&(0.9*ssdata[,3])>ssdata[,2]) &
          ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4] & ssdata[,2]>(1.1*ssdata[,4])&(0.9*ssdata[,2])>ssdata[,4] & ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]] = "CXR"
      y[!(ssdata[,2]>(1.1*ssdata[,1])&(0.9*ssdata[,2])>ssdata[,1]) & !(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) & !(ssdata[,2]>(1.1*ssdata[,4])&(0.9*ssdata[,2])>ssdata[,4]) &
          !(ssdata[,1]>(1.1*ssdata[,2])&(0.9*ssdata[,1])>ssdata[,2]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) & !(ssdata[,4]>(1.1*ssdata[,2])&(0.9*ssdata[,4])>ssdata[,2]) &
          ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3] & ssdata[,2]>(1.1*ssdata[,3])&(0.9*ssdata[,2])>ssdata[,3] & ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3]] = "CXS"
      y[!(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) & !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) &
          !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) & !(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3]) &
          ssdata[,1]>(1.1*ssdata[,2])&(0.9*ssdata[,1])>ssdata[,2] & ssdata[,3]>(1.1*ssdata[,2])&(0.9*ssdata[,3])>ssdata[,2] & ssdata[,4]>(1.1*ssdata[,2])&(0.9*ssdata[,4])>ssdata[,2]] = "CRS"
      y[!(ssdata[,3]>(1.1*ssdata[,2])&(0.9*ssdata[,3])>ssdata[,2]) & !(ssdata[,4]>(1.1*ssdata[,2])&(0.9*ssdata[,4])>ssdata[,2]) & !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) &
          !(ssdata[,2]>(1.1*ssdata[,3])&(0.9*ssdata[,2])>ssdata[,3]) & !(ssdata[,2]>(1.1*ssdata[,4])&(0.9*ssdata[,2])>ssdata[,4]) & !(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3]) &
          ssdata[,2]>(1.1*ssdata[,1])&(0.9*ssdata[,2])>ssdata[,1] & ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1] & ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]] = "XRS"
      y[!(ssdata[,2]>(1.1*ssdata[,1])&(0.9*ssdata[,2])>ssdata[,1]) & !(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) &
          !(ssdata[,1]>(1.1*ssdata[,2])&(0.9*ssdata[,1])>ssdata[,2]) & !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) &
          !(ssdata[,2]>(1.1*ssdata[,3])&(0.9*ssdata[,2])>ssdata[,3]) & !(ssdata[,2]>(1.1*ssdata[,4])&(0.9*ssdata[,2])>ssdata[,4]) & !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) &
          !(ssdata[,3]>(1.1*ssdata[,2])&(0.9*ssdata[,3])>ssdata[,2]) & !(ssdata[,4]>(1.1*ssdata[,2])&(0.9*ssdata[,4])>ssdata[,2]) & !(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3])] = "="
      # factor levels
      fl = c("C","X","R","S","CX","CR","CS","XR","XS","RS","CXR","CXS","CRS","XRS","=","NA")
    }
    if(dataset5=="nomix"){
      y[ssdata[,1]>(1.1*apply(ssdata[,c(-1,-2)],1,max))&(0.9*ssdata[,1])>apply(ssdata[,c(-1,-2)],1,max)] = "C"
      y[ssdata[,3]>(1.1*apply(ssdata[,c(-3,-2)],1,max))&(0.9*ssdata[,3])>apply(ssdata[,c(-3,-2)],1,max)] = "R"
      y[ssdata[,4]>(1.1*apply(ssdata[,c(-4,-2)],1,max))&(0.9*ssdata[,4])>apply(ssdata[,c(-4,-2)],1,max)] = "S"

      y[!(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) &
          ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4] &
          ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]] = "CR"
      y[!(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) &
          ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3] &
          ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3]] = "CS"
      y[!(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3]) & !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) &
          ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1] &
          ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]] = "RS"

      y[!(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) &
          !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) &
          !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) & !(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3])] = "="
      # factor levels
      fl = c("C","R","S","CR","CS","RS","=","NA")
    }
  }
  if(dataset4!="all"){ # plus and minus 10% difference clear for mixture/sequence to be favoured with near-equal category as well
    if(dataset5=="mix"){
      j = which(dataset4==c("C","X","R","S"))
      # null vector
      y = rep("O",nrow(ssdata)) # some other strategy is best
      y[ssdata[,j]>(1.1*apply(ssdata[,-j],1,min))&(0.9*ssdata[,j])>apply(ssdata[,-j],1,min)] = "J" # the focal strategy is jointly best
      y[!(ssdata[,2]>(1.1*ssdata[,1])&(0.9*ssdata[,2])>ssdata[,1]) & !(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) &
          !(ssdata[,1]>(1.1*ssdata[,2])&(0.9*ssdata[,1])>ssdata[,2]) & !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) &
          !(ssdata[,2]>(1.1*ssdata[,3])&(0.9*ssdata[,2])>ssdata[,3]) & !(ssdata[,2]>(1.1*ssdata[,4])&(0.9*ssdata[,2])>ssdata[,4]) & !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) &
          !(ssdata[,3]>(1.1*ssdata[,2])&(0.9*ssdata[,3])>ssdata[,2]) & !(ssdata[,4]>(1.1*ssdata[,2])&(0.9*ssdata[,4])>ssdata[,2]) & !(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3])] = "=" # the focal strategy is equally best
      y[ssdata[,j]>(1.1*apply(ssdata[,-j],1,max))&(0.9*ssdata[,j])>apply(ssdata[,-j],1,max)] = dataset4 # the focal strategy is outright the best
      # factor levels
      fl = c(c("C","X","R","S")[j],"O","J","=")
    }
    if(dataset5=="nomix"){
      j = which(dataset4==c("C","X","R","S"))
      # null vector
      y = rep("O",nrow(ssdata)) # some other strategy is best
      y[ssdata[,j]>(1.1*apply(ssdata[,c(-2,-j)],1,min))&(0.9*ssdata[,j])>apply(ssdata[,c(-2,-j)],1,min)] = "J" # the focal strategy is jointly best
      y[!(ssdata[,3]>(1.1*ssdata[,1])&(0.9*ssdata[,3])>ssdata[,1]) & !(ssdata[,4]>(1.1*ssdata[,1])&(0.9*ssdata[,4])>ssdata[,1]) &
          !(ssdata[,1]>(1.1*ssdata[,3])&(0.9*ssdata[,1])>ssdata[,3]) & !(ssdata[,1]>(1.1*ssdata[,4])&(0.9*ssdata[,1])>ssdata[,4]) &
          !(ssdata[,3]>(1.1*ssdata[,4])&(0.9*ssdata[,3])>ssdata[,4]) & !(ssdata[,4]>(1.1*ssdata[,3])&(0.9*ssdata[,4])>ssdata[,3])] = "=" # the focal strategy is equally best
      y[ssdata[,j]>(1.1*apply(ssdata[,c(-2,-j)],1,max))&(0.9*ssdata[,j])>apply(ssdata[,c(-2,-j)],1,max)] = dataset4 # the focal strategy is outright the best
      # factor levels
      fl = c(c("C","X","R","S")[j],"O","J","=")
    }
  }

  # reformat data
  plotdata = cbind(simdata,y)
  colnames(plotdata) = c("Population_Size","Intrinsic_Birth","Intrinsic_Death","Female_Exposure","Male_Exposure",
                         "Initial_Frequency_A","Effectiveness_1","Resistance_Restoration_A","Restoration_Dominance_A","Resistance_Cost_A","Cost_Dominance_A",
                         "Initial_Frequency_B","Effectiveness_2","Resistance_Restoration_B","Restoration_Dominance_B","Resistance_Cost_B","Cost_Dominance_B","y")

  # give table output
  if(dataset4=="all"&dataset5=="mix"){comp_mix_all[i,] = c(dataset0,dataset3,dataset4,table(factor(plotdata$y,levels=fl)),1000000)}
  if(dataset4!="all"&dataset5=="mix"){comp_mix_other[i-12,] = c(dataset0,dataset3,dataset4,table(factor(plotdata$y,levels=fl)),1000000)}
  if(dataset4=="all"&dataset5=="nomix"){comp_nomix_all[i-60,] = c(dataset0,dataset3,dataset4,table(factor(plotdata$y,levels=fl)),1000000)}
  if(dataset4!="all"&dataset5=="nomix"){comp_nomix_other[i-72,] = c(dataset0,dataset3,dataset4,table(factor(plotdata$y,levels=fl)),1000000)}

  # filter to remove = and all infrequent combinations for plotting
  rf = table(factor(plotdata$y,levels=fl))
  if(dataset4!="all"){for(j in 1:length(fl)){if((rf[j]/sum(rf))<0.01){plotdata = plotdata[plotdata[,18]!=fl[j],]}}}
  if(dataset4=="all"){for(j in 1:length(fl)){if((rf[j]/sum(rf))<0.05){plotdata = plotdata[plotdata[,18]!=fl[j],]}}}
  if(sum(table(factor(plotdata$y,levels=fl))>0)>1){
    # build model of y using parameter variables of plotdata
    simdata_tree = partykit::ctree(y ~ Population_Size+Intrinsic_Birth+Intrinsic_Death+Female_Exposure+Male_Exposure+Initial_Frequency_A+Effectiveness_1+Resistance_Restoration_A+Restoration_Dominance_A+Resistance_Cost_A+Cost_Dominance_A+Initial_Frequency_B+Effectiveness_2+Resistance_Restoration_B+Restoration_Dominance_B+Resistance_Cost_B+Cost_Dominance_B, data=as.data.frame(plotdata),control=ctree_control(minbucket=0.05*nrow(plotdata)))
    # plot classification tree
    pdf(paste("DT_",dataset5,"/DT_",dataset0,"_",dataset3,"_",dataset4,"_",dataset5,".pdf",sep=""),width=40,height=9)
    if(dataset4!="all"){plot(simdata_tree)}
    if(dataset4=="all"){plot(simdata_tree,terminal_panel=node_barplot(simdata_tree,rot=90,just="center"))}
    dev.off()
  }

}

# save tables
write.csv(comp_mix_all,paste("DT_mix_STable_mix_all.csv",sep=""),row.names=F)
write.csv(comp_mix_other,paste("DT_mix_STable_mix_other.csv",sep=""),row.names=F)
write.csv(comp_nomix_all,paste("DT_nomix_STable_nomix_all.csv",sep=""),row.names=F)
write.csv(comp_nomix_other,paste("DT_nomix_STable_nomix_other.csv",sep=""),row.names=F)




#### mn by effectiveness, partner effectiveness, exposure ####

# data
ssdata_fnn = ssdata_first_nn # time to first 50% frequency
ssdata_fmn = ssdata_first_mn # time to first 50% frequency
ssdata_fmm = ssdata_first_mm # time to first 50% frequency

# reorder data
ssdata_fnn = ssdata_fnn[,c(4,3,1,2,5)]
ssdata_fmn = ssdata_fmn[,c(4,3,1,2,5)]
ssdata_fmm = ssdata_fmm[,c(4,3,1,2,5)]

# plot
pdf("waterfall_PTbyParam_MN.pdf",width=8*2.5,height=7.5*1.5)
par(mar=c(5,8,5,5),mfrow=c(2,4)) # set margins and layout
layout(matrix(c(1,4,1,4,2,5,2,5,3,6,3,6,7,7),2,7))

## probability of resistance

# focal insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Focal Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fmn[round(simdata[,7],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
  #if(i==5){points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=0.25)[i])}
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"A",font=2,xpd=NA,cex=3)

# partner insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Partner Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fmn[simdata[,7]>0.8&round(simdata[,13],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"B",font=2,xpd=NA,cex=3)

# female exposure
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Female Exposure",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fmn[simdata[,7]>0.8&round(simdata[,4],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"C",font=2,xpd=NA,cex=3)

## time to resistance

# focal insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,250),xlab="Focal Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fmn[ssdata_fmn[,i]<501&round(simdata[,7],2)==(j/100),i]
    tres_mean[j+1] = mean(xd)
    tres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,275,"D",font=2,xpd=NA,cex=3)


# partner insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,300),xlab="Partner Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fmn[ssdata_fmn[,i]<501&simdata[,7]>0.8&round(simdata[,13],2)==(j/100),i]
    if(length(xd)<2){
      tres_mean[j+1] = NA
      tres_ci[,j+1] = NA
    }
    if(length(xd)>1){
      tres_mean[j+1] = mean(xd)
      tres_se[j+1] = sd(xd)/sqrt(length(xd))
      #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
      #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
    }
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,330,"E",font=2,xpd=NA,cex=3)

# female exposure
plot(-10,-10,xlim=c(0,1),ylim=c(0,500),xlab="Female Exposure",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fmn[ssdata_fmn[,i]<501&simdata[,7]>0.8&round(simdata[,4],2)==(j/100),i]
    if(length(xd)<2){
      tres_mean[j+1] = NA
      tres_ci[,j+1] = NA
    }
    if(length(xd)>1){
      tres_mean[j+1] = mean(xd)
      tres_se[j+1] = sd(xd)/sqrt(length(xd))
      #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
      #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
    }
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,550,"F",font=2,xpd=NA,cex=3)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim=c(0,1),ylim=c(0,1))
legend(-1,0.5,c("Sequences","Rotations","Mosaics","Mixtures","Maximum"),pch=rep(16,5),col=viridis(5),bg="white",cex=2.5,yjust=0.5,xpd=NA)

dev.off()




#### nn by effectiveness, partner effectiveness, exposure ####

# plot
pdf("waterfall_PTbyParam_NN.pdf",width=8*2.5,height=7.5*1.5)
par(mar=c(5,8,5,5),mfrow=c(2,4)) # set margins and layout
layout(matrix(c(1,4,1,4,2,5,2,5,3,6,3,6,7,7),2,7))

## probability of resistance

# focal insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Focal Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fnn[round(simdata[,7],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
  #if(i==5){points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=0.25)[i])}
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"A",font=2,xpd=NA,cex=3)

# partner insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Partner Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fnn[simdata[,7]>0.8&round(simdata[,13],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"B",font=2,xpd=NA,cex=3)

# female exposure
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Female Exposure",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fnn[simdata[,7]>0.8&round(simdata[,4],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"C",font=2,xpd=NA,cex=3)

## time to resistance

# focal insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,250),xlab="Focal Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fnn[ssdata_fnn[,i]<501&round(simdata[,7],2)==(j/100),i]
    tres_mean[j+1] = mean(xd)
    tres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,275,"D",font=2,xpd=NA,cex=3)


# partner insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,300),xlab="Partner Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fnn[ssdata_fnn[,i]<501&simdata[,7]>0.8&round(simdata[,13],2)==(j/100),i]
    if(length(xd)<2){
      tres_mean[j+1] = NA
      tres_ci[,j+1] = NA
    }
    if(length(xd)>1){
      tres_mean[j+1] = mean(xd)
      tres_se[j+1] = sd(xd)/sqrt(length(xd))
      #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
      #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
    }
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,330,"E",font=2,xpd=NA,cex=3)

# female exposure
plot(-10,-10,xlim=c(0,1),ylim=c(0,500),xlab="Female Exposure",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fnn[ssdata_fnn[,i]<501&simdata[,7]>0.8&round(simdata[,4],2)==(j/100),i]
    if(length(xd)<2){
      tres_mean[j+1] = NA
      tres_ci[,j+1] = NA
    }
    if(length(xd)>1){
      tres_mean[j+1] = mean(xd)
      tres_se[j+1] = sd(xd)/sqrt(length(xd))
      #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
      #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
    }
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,550,"F",font=2,xpd=NA,cex=3)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim=c(0,1),ylim=c(0,1))
legend(-1,0.5,c("Sequences","Rotations","Mosaics","Mixtures","Maximum"),pch=rep(16,5),col=viridis(5),bg="white",cex=2.5,yjust=0.5,xpd=NA)

dev.off()




#### mm by effectiveness, partner effectiveness, exposure ####

# plot
pdf("waterfall_PTbyParam_MM.pdf",width=8*2.5,height=7.5*1.5)
par(mar=c(5,8,5,5),mfrow=c(2,4)) # set margins and layout
layout(matrix(c(1,4,1,4,2,5,2,5,3,6,3,6,7,7),2,7))

## probability of resistance

# focal insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Focal Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fmm[round(simdata[,7],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
  #if(i==5){points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=0.25)[i])}
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"A",font=2,xpd=NA,cex=3)

# partner insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Partner Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fmm[simdata[,7]>0.8&round(simdata[,13],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"B",font=2,xpd=NA,cex=3)

# female exposure
plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xlab="Female Exposure",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  pres_mean = rep(0,101)
  pres_se = rep(0,101)
  pres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = 1*(ssdata_fmm[simdata[,7]>0.8&round(simdata[,4],2)==(j/100),i]<1001)
    pres_mean[j+1] = mean(xd)
    pres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #pres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  pres_ci[1,] = pres_mean + 1.96*pres_se
  pres_ci[2,] = pres_mean - 1.96*pres_se
  # smooth
  pres_mean = movavg(pres_mean,5,"s")
  pres_ci[1,] = movavg(pres_ci[1,],5,"s")
  pres_ci[2,] = movavg(pres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),pres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(pres_ci[1,],rev(pres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Probability of Resistance",side=2,line=4,cex=1.5,font=2)
text(-0.2,1.1,"C",font=2,xpd=NA,cex=3)

## time to resistance

# focal insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,250),xlab="Focal Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fmm[ssdata_fmm[,i]<501&round(simdata[,7],2)==(j/100),i]
    tres_mean[j+1] = mean(xd)
    tres_se[j+1] = sd(xd)/sqrt(length(xd))
    #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
    #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,275,"D",font=2,xpd=NA,cex=3)


# partner insecticide effectiveness
plot(-10,-10,xlim=c(0,1),ylim=c(0,300),xlab="Partner Insecticide Effectiveness",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fmm[ssdata_fmm[,i]<501&simdata[,7]>0.8&round(simdata[,13],2)==(j/100),i]
    if(length(xd)<2){
      tres_mean[j+1] = NA
      tres_ci[,j+1] = NA
    }
    if(length(xd)>1){
      tres_mean[j+1] = mean(xd)
      tres_se[j+1] = sd(xd)/sqrt(length(xd))
      #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
      #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
    }
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,330,"E",font=2,xpd=NA,cex=3)

# female exposure
plot(-10,-10,xlim=c(0,1),ylim=c(0,500),xlab="Female Exposure",ylab="",cex.axis=1.5,cex.lab=2,font.lab=2,xaxs="i",yaxs="i",las=1)
for(i in c(1,5,2,3,4)){
  # mean and se
  tres_mean = rep(0,101)
  tres_se = rep(0,101)
  tres_ci = matrix(0,2,101)
  for(j in 0:100){
    xd = ssdata_fmm[ssdata_fmm[,i]<501&simdata[,7]>0.8&round(simdata[,4],2)==(j/100),i]
    if(length(xd)<2){
      tres_mean[j+1] = NA
      tres_ci[,j+1] = NA
    }
    if(length(xd)>1){
      tres_mean[j+1] = mean(xd)
      tres_se[j+1] = sd(xd)/sqrt(length(xd))
      #bp = boot(xd,function(x,i) mean(x[i]),R=1000)
      #tres_ci[,j+1] = as.vector(boot.ci(bp,conf=0.95,type="perc")[[4]])[4:5]
    }
  }
  # confidence interval
  tres_ci[1,] = tres_mean + 1.96*tres_se
  tres_ci[2,] = tres_mean - 1.96*tres_se
  # smooth
  tres_mean = movavg(tres_mean,5,"s")
  tres_ci[1,] = movavg(tres_ci[1,],5,"s")
  tres_ci[2,] = movavg(tres_ci[2,],5,"s")
  # plot
  points(seq(0,1,0.01),tres_mean,type="l",lwd=2,col=viridis(5,alpha=1)[i])
  polygon(c(seq(0,1,0.01),seq(1,0,-0.01)),c(tres_ci[1,],rev(tres_ci[2,])),col=viridis(5,alpha=0.5)[i],border=viridis(5,alpha=0.5)[i])
}
box()
mtext("Time to First-to-Break",side=2,line=4,cex=1.5,font=2)
text(-0.2,550,"F",font=2,xpd=NA,cex=3)

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim=c(0,1),ylim=c(0,1))
legend(-1,0.5,c("Sequences","Rotations","Mosaics","Mixtures","Maximum"),pch=rep(16,5),col=viridis(5),bg="white",cex=2.5,yjust=0.5,xpd=NA)

dev.off()




# END
