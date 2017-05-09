# Isolation-Drift
R code that simulates drift and isolation and calls Structure to evaluate the impact. 
Analysis also includes measurement of heterozygosity, Fst, and DAPC (discriminant analysis of
principle components) to futher quantify the effects of isolation in the two populations.

Use IsoDriftModel.R to control: 
  -input of genotypes (alleles and allele frequencies)
  -simulation parameters
  -running of program

IsoDriftModel.R will call FunctionSourcer.R (stored in source directory), which calls RunSims.R (also in source directory).
Output written to a directory called 'output'. See isodriftmodel2.9.zip file for the setup of directories and subdirectories.
