####################################################################
# Drift/Isolation model
# Jordan Hoffman and Janna Willoughby
#
# Usage: set parameters for species of interest and create data frame
# with genotype data. (Data should be organized with two columns per 
# locus cannot include a column with sample labels. Missing data code
# can be input.) Working directory must also be set. Simulations require 
# use of strataG, which links to Structure.  See strataG documentation
# for instalation instructions.
# 
#
####################################################################
setwd("~/model1.5")
directory = getwd()
outdir    = paste(directory,"/output/",sep="")                    #directory to save model output  
source(paste(directory, "/source/FunctionSourcer.R", sep = ''))   #source functions and set source directory

#input files and info
genotypes = GENOTYPES_DATAFRAME                                          # input data frame with genotypes
Anames    = c("A112", "A130", "C4", "C114", "D102", "R9", "Loc7")        # input locus names

# simulation parameters
popsize  = 50        # population size to simulate 
simyears = 501       # total years to run simulation
survival = 0.96      # survival rate
agecap   = 62        # lifespan (maximum allowable age)
reps     = 100       # number of replicate simulations
structK  = 3         # number of K to use in Structure analyses

# genotypic data parameters
numloci     = 7      # number of loci in input file
missingdata = -9     # missing data code for genotypic data


RunSims(Anames, popsize, simyears, survival, agecap, reps, structK, numloci, missingdata)