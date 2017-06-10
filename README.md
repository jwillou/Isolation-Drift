# Isolation-Drift
R code that simulates drift and isolation and calls Structure to evaluate the impact. 
Analysis also includes measurement of heterozygosity, Fst, and DAPC (discriminant analysis of
principle components) to futher quantify the effects of isolation in the two populations.

Use IsoDriftModel.R to control: 
  1. input of genotypes (alleles and allele frequencies)
  2. simulation parameters
  3. running of program

IsoDriftModel.R will call FunctionSourcer.R (stored in source directory), sources calls 
RunSims.R and Advance.R (also in source directory). Output written to a directory called 'output'. 
See isodriftmodel2.9.zip file for the setup of directories and subdirectories.

Program will output 4 files for each of the population sizes simulated. Files written will include:
  1. dpc.csv: results of dicriminant analysis of principle components. Columns contents (in order): 
  population size, replicate number (set in IsoDriftModel.R as reps), year, maximum age allowed in 
  the simulation, proportion of individuals in source population that assigned to the source population, 
  proporiton of individuals in sink (isolated) population that assigned to the larger source population.
  2. eva.csv: results of structure analysis. Columns contents (in order): population size, replicate 
  number, year, maximum age allowed in the simulation, penalized log likelihood for K = 1, penalized log 
  likelihood for K = 2, penalized log likelihood for K = 3.
  3. fst.csv: estimates of Fst between source and sink populations. Columns contents (in order): population 
  size, replicate number, year, maximum age allowed in the simulation, and Fst between the two simulated 
  populaitons.
  4. het.csv: heterozygosity estimated each year in the sink population. Columns contents (in order): 
  population size, replicate number, year, maximum age allowed in the simulation, and heterozygosity estimate.
  
