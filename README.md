# Isolation-Drift
R code that calls simulates drift and isoaltion and calls Structure to evaluate the impact

Use IsoDriftModel.R to control: 
  -input of genotype object (data frame with 1 row per individual and 2 columns per locus)
  -simulation parameters
  -running of program

IsoDrifModel.R will call FunctionSourcer.R (stored in source directory), which calls RunSims.R (also in source directory).
Output written to output directory. See zip file for example setup.

