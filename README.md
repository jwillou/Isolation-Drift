# Isolation-Drift
R code that simulates drift and isolation and calls Structure to evaluate the impact. 
Analysis also includes measurement of heterozygosity, Fst, and DAPC (discriminant analysis of
principle components) to futher quantify the effects of isolation in the two populations.

Use IsoDriftModel.R to control: 
  1. input of genotypes (alleles and allele frequencies)
  2. simulation parameters
  3. running of program

IsoDriftModel.R will call FunctionSourcer.R (stored in source directory), which calls RunSims.R (also in source directory).
Output written to a directory called 'output'. See isodriftmodel2.9.zip file for the setup of directories and subdirectories.

Program will output 4 files for each of the population sizes simulated. Assuming the simulated population size is 100 individuals 
(i.e. popsize = 100 in IsoDriftModel.R), files written will include:
  1. dpc100.csv: results of dicriminant analysis of principle components. Columns are, in order, replicate number (set in 
  IsoDriftModel.R as reps), year, maximum age allowed in the simulation, proportion of individuals in source population that 
  assigned to the source population, proporiton of individuals in sink (isolated) population that assigned to the larger source 
  population.
  2. eva100.csv: results of structure analysis. Columns are, in order, replicate number, year, maximum age allowed in the simulation, 
  penalized log likelihood for K = 1, penalized log likelihood for K = 2, penalized log likelihood for K = 3.
  3. fst100.csv: estimates of Fst between source and sink populations. Columns are, in order, replicate number, year, maximum age 
  allowed in the simulation, and Fst between the two simulated populaitons.
  4. het100.csv: heterozygosity estimated each year in the sink population. Columns are, in order, replicate number, year, maximum 
  age allowed in the simulation, and heterozygosity estimate.
  
