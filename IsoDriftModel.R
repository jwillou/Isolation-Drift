####################################################################
# Drift/Isolation model
# Jordan Hoffman and Janna Willoughby
#
# Usage: Set parameters for species of interest and intput alleles
# and associated allele frequencies. (Data should be organized such that 
# allele and the corresponding frequency are input into their respective 
# variables in the same order.) A working directory must also be set, 
# and must contain this file (IsoDriftModel.R), a directory called 
# 'output' and a directory called 'source' that contains FunctionSourcer.R,
# Advance.R, and RunSims.R. Simulations require use of strataG, which links 
# to Structure, as well as diveRsity and adegenet packages. See package 
# documentation for instalation instructions. In some situations, the 
# library path variable may need to be modified. This can be done in the
# FunctionSourcer.R file.
# 
#
####################################################################
setwd("")

directory = getwd()
outdir    = paste(directory,"/output/",sep="")                    #directory to save model output  
source(paste(directory, "/source/FunctionSourcer.R", sep = ''))   #source functions and set source directory

#input information

#alleles
alleles = list(A112 = c(146,148,150,152,154,157,159,160,161,168,172,191,193),
               A130 = c(234,236,250,252,254,256,260,262,264,266,270,272,274,276,278,280),
               C4   = c(199,203,207,211,215,219,223,227,231,239,259,263),
               C114 = c(204,216,220,224,228,232,236,248),
               D102 = c(202,204,208,216,220,224,248),
               R9   = c(205,208,218,220,222,225) )

#allele frequencies
allelefreqs = list(A112freq = c(0.024,0.167,0.135,0.415,0.021,0.011,0.196,0.008,0.008,0.008,0.003,0.003,0.003),
                   A130freq = c(0.017,0.020,0.006,0.003,0.078,0.075,0.003,0.009,0.003,0.049,0.003,0.017,0.095,0.391,0.213,0.020),
                   C4freq   = c(0.018,0.416,0.018,0.009,0.186,0.021,0.105,0.129,0.045,0.003,0.006,0.045),
                   C114freq = c(0.003,0.008,0.492,0.011,0.370,0.109,0.005,0.003),
                   D102freq = c(0.006,0.012,0.003,0.871,0.089,0.006,0.012),
                   R9freq   = c(0.005,0.011,0.014,0.957,0.008,0.005) )

# simulation parameters
popsize  = c(50,100,200,350,500)        # population sizes to simulate 
simyears = 501                          # total years to run the isolation portion of the simulation (does not include delay)
survival = 0.96                         # survival rate
agecap   = seq(2, 102, 20)              # maximum age
reps     = 100                          # replicates
structK  = 3                            # number of K for structure analyses
levels   = seq(0, 500, 25)              # years to run structure
delay    = 75                           # number of years between initiation of large pop and isolation second pop

RunSims(alleles, allelefreqs, popsize, simyears, survival, agecap, reps, structK, levels, delay)
