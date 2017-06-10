Advance = function(pop, numloci, a, lcols){
  #keep previous year, assign new genotypes to new individuals (replacing dead) based on frequencies in previous year
  #age individuals and determine mortality
  pop[,(numloci*2)+1] = pop[,(numloci*2)+1] + 1
  
  #determine 'random' mortality, mark those selected as dead (0)
  tokill = sample(seq(1, nrow(pop), 1), size = round((1-survival)*nrow(pop)), replace=FALSE)
  pop[tokill,(numloci*2)+2] = 0
  
  #mark individuals older than agecap as dead (0)
  pop[pop[,(numloci*2)+1]>agecap[a],(numloci*2)+2] = 0
  
  #generate list of dead
  toreplace = which(pop[,(numloci*2)+2]==0)
  pop[toreplace,] = NA
  
  #retain copy of last year
  ppop = pop
  
  #iterate over loci to set genotypes
  locus = 1
  for(l in lcols){
    #generate list of available alleles
    alleles[[locus]]     = sort(as.numeric(unique(c(ppop[,l], ppop[,l+1]))))
    
    #generate list of allele frequencies
    allelefreqs[[locus]] = as.numeric(sort(table(c(ppop[,l], ppop[,l+1]))) / (nrow(ppop[!is.na(ppop[,1]),,drop=FALSE])*2))
    
    #replace rows of 'dead' indivuals with newly sampled allele frequencies
    for(r in 1:length(toreplace)){
      if(length(alleles[[locus]])==1){
        pop[toreplace[r], l]   = alleles[[locus]]
        pop[toreplace[r], l+1] = alleles[[locus]]
      }else{
        pop[toreplace[r], l]   = sample(x=alleles[[locus]], size=1, prob=allelefreqs[[locus]], replace=TRUE)
        pop[toreplace[r], l+1] = sample(x=alleles[[locus]], size=1, prob=allelefreqs[[locus]], replace=TRUE)
      }
    }
    locus = locus + 1  
    
    #mark age = 0
    pop[toreplace, (numloci*2)+1]  = 0
    
    #reset to alive
    pop[toreplace, (numloci*2)+2]  = 1 
  }
  pops = list(ppop=ppop, pop=pop)
  return(pops)
}
