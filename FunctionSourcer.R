#Set working directory, import packages, source functions, 
setwd(paste(directory,"/source/", sep = ''))    # set temp working directory 

#import packages

#install.packages("diveRsity", lib="/scratch/snyder/j/jwillou/isodrift/Rlibs/", repos='http://cran.us.r-project.org')
#install.packages("psych", lib="/scratch/snyder/j/jwillou/isodrift/Rlibs/", repos='http://cran.us.r-project.org')
.libPaths("/scratch/snyder/j/jwillou/isodrift/Rlibs") 
library(diveRsity)
library(strataG)
library(adegenet)

#source functions
source(paste(getwd(), "/RunSims.R", sep = ''))
source(paste(getwd(), "/Advance.R", sep = ''))

