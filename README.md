## Description 
Package detsims is a framework and model for running deterministic simulations 
of resistance evolution in mosquitoes. It is designed to fully reproduce the 
results in Madgwick & Kanitz 2021 "Comparison  of resistance-management 
strategies with new insecticides for mosquito control using bed-nets", from 
model definitions, simulation executions and visualization of results, including
published figures.

The final version of the publication can be found [here](https://malariajournal.biomedcentral.com/articles/10.1186/s12936-022-04083-z)


## Installation
```R
# You will need devtools to install evolveR. If you do not have it yet, run:
install.packages("devtools")

# Install the Syngenta development version from GitLab:
devtools::install_github("syngenta/Madgwick-Kanitz-detsims", ref='main')
```

## Package structure
This package follows a traditional R-package setup, with its functions defined 
under folder 'R', their associated help files in 'man'. Importantly though, we 
used also used a folder called 'scripts' here: the files contained in this 
folder are pure R scripts designed to execute an example run or the simulations 
in 'OneShot.R'; and the full range of results into ./data/ with 'Simulator.R'.
Please, use the latter with caution as it would execute 1 million simulations
on the machine it's been called from, using the parameter definitions given at 
the beginning of the script under section 'meta-parameters'.


Note also that there are 4 inheritance definitions (Nuclear, Mitochondrial):
* NN - c(n,n)
* MM - c(m,m)
* MN - c(m,n)
* NM - c(n,m)


And 6 strategy numbers as follows:
  1. soloA
  2. soloB
  3. Combination (Mosaic in the main text)
  4. Mixture
  5. Rotation
  6. Sequence
  

These can be combined in all possible ways, but in the paper, we report a subset
of 15 of these inheritance-strategy combinations. See the paper for more details.

We have executed these simulations in a high-performance computer cluster.
Additional information and assistance on this can be provided on demand (see 
contact information below in section 'Issues')

The script for generating all figures in the paper and supplement (plus extras)
is also under 'vignettes': 'Figures.R'. Use with caution as this is compilation
of all visualization explored in the development of this work leading to a fair
amount of computing time, as well.


## Issues
Issues in this package can be reported using GitHub's Issues above. Other
inquiries can be made by e-mail to Ricardo(dot)Kanitz(at)syngenta(dot)com and 
Philip(dot)Madgwick(at)syngenta(dot)com.
