RunSims = function(Anames, popsize, simyears, survival, agecap, reps, structK, numloci, missingdata){
  #convert new variable names to old names
  Nc = popsize     # population sizes to simulate c(50,100,200,350,500)
  gt = simyears    # total years to run simulation (must be larger than popsize)
  SR = survival    # survival rate
  agecap = agecap  # maximum age
  x = reps         # replicates
  K = structK      # number of K for structure analyses
  L = numloci      # number of loci in input file

  # determine allele frequencies from input genotype data frame
  # check input data
  if(dim(genotypes)[2] %% 2 != 0){
    "genotype data not formatted properly"
    break
  }
  
  #iterate over pairs of columns
  names = 1
  for(l in 1:(L*2)){
    #skip every other column
    if(l %% 2 == 0){
      names = names + 1
      next
    }
    #get all alelles, remove missing values
    genos = c(genotypes[,l], genotypes[,l+1])
    genos = as.numeric(gsub(missingdata, NA, genos))
    genos = genos[!is.na(genos)]
    ugenos = sort(unique(genos))
    
    #alleles
    name = paste(Anames[names], sep="")  
    assign(name, ugenos) 
    
    #allele frequencies
    freqs = NULL
    for(a in 1:length(ugenos)){
      freqs = c(freqs, length(grep(ugenos[a], genos, fixed=TRUE))/length(genos))
    }
    name = paste(Anames[names], "freq", sep="") 
    assign(name, freqs)
  }
  #example starting data when genotype data frame is unavailable
  #All alleles for each loci
  #A112 = c(179,181,189,191,195)
  #A130 = c(190,194,204)
  #C4   = c(152,154,158,160,162,172)
  #C114 = c(200,203,209)
  #D102 = c(126,128,132,134,135,136,138,140,152)
  #R9   = c(189,191,193,197,199)
  #Loc7 = c(77,85,89,93)
  #Allele Freq
  #A112freq = c(0.082,0.021,0.598,0.284,0.015)
  #A130freq = c(0.896,0.047,0.057)
  #C4freq   = c(0.005,0.146,0.245,0.547,0.047,0.010)
  #C114freq = c(0.036,0.567,0.397)
  #D102freq = c(0.120,0.234,0.151,0.036,0.172,0.130,0.036,0.089,0.031)
  #R9freq   = c(0.484,0.057,0.130,0.224,0.104)
  #Loc7freq = c(0.577,0.057,0.314,0.052)
    
  #HWE p val result array
  HWEp     = array(NA,c(gt,L,x))
  agerange = c(1:(agecap - 1))
  Y = length(Nc)                #End format array
  m = 0                         # migration rate
  
  #set up various matrices for storing input and results
  for(a in 1:length(Anames)){
    name = paste("Amatrix", Anames[a], sep="")  #Alleles matrix (g,x)
    assign(name, matrix(0,gt,x))                
    name = paste(Anames[a], "matrix", sep="")   #Alleles at different gens for each iteration
    assign(name, matrix(NA,gt,x))              
    name = paste("avg", Anames[a], sep="")      #Number of alleles averages across iterations
    assign(name, as.matrix(data.frame(c1 = seq(1,gt), c2 = rep(NA, gt))))               
  }
  
  #begin iterating over Nc values 
  for (n in Nc){
    age = array(0,c(n,gt+1,x))
    #heterozygosity matrix (g,x)
    Hmatrix  = matrix(0,gt,x)
    Hmatrix2 = matrix(0,gt,x)
    
    for (j in 1:x){
      levels = seq(0, 500, 25)
      #output array for sims == A
      A = array(NA,c(n,2,L,gt+1))
      for(i in 1:n){
        for(a in 1:length(Anames)){
          A[i,,a,1]=sample(get(paste(Anames[a])),2,replace=T,prob=get(paste(Anames[a], "freq", sep="")))
          #calculate new allele frequency (g=0)
          name = paste("g0", Anames[a], "freq", sep="")
          assign(name, table(factor(A[,,a,1],levels=get(Anames[a])))/(2*n))
        }
      }
      #loop array generation considering mortality survival rate of 0 means non-overlapping generations
      for(g in 1:gt){
        A[,,,g+1]<-A[,,,g]
        for(a in 1:length(Anames)){
          name = paste("g1", Anames[a], "freq", sep="")
          assign(name, table(factor(A[,,a,g],levels=get(Anames[a])))/(2*n))
        }
        
        for (i in 1:n){
          age[i,1,j] = sample(agerange,1,replace=TRUE,prob=NULL)
          if(runif(1)<SR){
            #survives
            A[i,,,g+1] = A[i,,,g]
            # +1 year to age
            age[i,g+1,j] = age[i,g,j]+1
          }else{
            if(runif(1)<m){
              for(a in 1:length(Anames)){
                A[i,,a,g+1] = sample(get(paste(Anames[a])),2,replace=TRUE,prob=get(paste(Anames[a], "freq", sep="")))
              }
            }else{
              for(a in 1:length(Anames)){
                vars = get(paste("g1",Anames[a], "freq", sep=""))
                A[i,,a,g+1] = sample(get(paste(Anames[a])),2,replace=TRUE,prob=c(vars))
              }
            }
            #age reset
            age[i,g+1,j] = 1
          }
          if(age[i,g,j]==agecap){
            #dies, sample alleles from allele freq of prev timestep population
            for(a in 1:length(Anames)){
              A[i,,a,g+1] = sample(get(paste(Anames[a])),2,replace=TRUE,prob=c(get(paste("g1",Anames[a], "freq", sep=""))))
            }
            #age reset
            age[i,g+1,j] = 1
          }
        }#i in 1:n
        
        for(a in 1:length(Anames)){
          #
          values = get(paste(Anames[a], "matrix", sep=""))
          values[g,j] = nrow(table(A[,,a,g]))
          name = paste(Anames[a], "matrix", sep="")
          assign(name, values)
          
          #Number of Alleles
          values = get(paste("Amatrix", Anames[a], sep=""))
          values[g,j] = length(table(factor(A[,,a,g])))
          name = paste(paste("Amatrix", Anames[a], sep=""))
          assign(name, values)  
        }
        #heterozygosity
        for (h in 1:n){
          for(a in 1:length(Anames)){
            if (A[h,1,a,g]==A[h,2,a,g]){ Hmatrix[g,j]=Hmatrix[g,j]+0
            }else{                       Hmatrix[g,j]=Hmatrix[g,j]+1}
          }
        }#h
        b=L*n
        
        Hmatrix2[g,j] = Hmatrix[g,j]/b
      }#gt  
      
      ##need to format appropriately for STRUCTURE in strataG, then run structure
      Eout = NULL
      for(s in 1:length(levels)){
        structure = array(NA,c(n,(2+(2*L)),1))
        gs = matrix(NA,100,(2+(2*L)))
        if(s==1){
          structure[,1,1]=c(1:50)
          structure[,2,1]="A"
          index = levels[s]+1
        }else{
          structure[,1,1]=c(51:100)
          structure[,2,1]="B"
          index = levels[s]
        }
        spot = 3
        for(a in 1:length(Anames)){
          structure[,spot,1]  = A[,1,a,index]
          spot = spot + 1
          structure[,spot,1]  = A[,2,a,index]
          spot = spot + 1
        }
        if(s==1){
          hold = structure
          next
        }else{
          gs[1:50,]   = hold[1:50,,1]
          gs[51:100,] = structure[1:50,,1]
          write.table(gs, file=paste(outdir, "gs",n,"_",j,".txt", sep=""), row.names=FALSE, col.names=FALSE)
          s  = read.table(paste(outdir, "gs",n,"_",j,".txt", sep=""))
          g  = gtypes(gen.data=s,id.col=1,strata.col=2,locus.col=3,dna.seq=NULL,description=NULL,delete.missing.strata=TRUE,code.start=1)
          st = structure.run(g,k.range=1:K,num.k.rep=5,numreps=100000,burnin=100000,noadmix=TRUE,freqscorr=TRUE,num.cores=2)
          E  = structure.evanno(st,plot=FALSE)
          Eout = rbind(Eout, E$mean.ln.k)
          Evanno_out = cbind(rep(j, nrow(Eout)), rep(n, nrow(Eout)), Eout)
          write.csv(Evanno_out,paste(outdir, "Evanno",n,"_",j,".csv", sep = ""),sep=",",col.names=FALSE, append=FALSE)
          
          #write Heterozygosity
          write.csv(Hmatrix2,paste(outdir, "Heterozygosity",n,"_",j,".csv", sep = ""),sep=",",col.names=FALSE, append=FALSE)
        }
      }
    }#n  
  }
}