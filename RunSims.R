RunSims = function(alleles, allelefreqs, popsize, simyears, survival, agecap, reps, structK, levels, delay){
  #set up some needed variables
  lnames  = names(alleles)
  numloci = as.numeric(length(lnames))
  
  #keep copy of alleles and frequencies
  originalalleles     = alleles
  originalallelefreqs = allelefreqs
  
  #iterate over popsize values 
  for (n in 1:length(popsize)){
    #iterate over maximum ages
    for(a in 1:length(agecap)){
      #iterate over replicates
      for (j in 1:reps){
        #het matrix
        hetmatrix = matrix(nrow = (simyears + 1), ncol = 5) 
        colnames(hetmatrix) = c("popsize", "rep", "year", "maxage", "H")
        
        #fst matrix 
        fstmatrix = matrix(nrow = (simyears + 1), ncol = 5)
        colnames(fstmatrix) = c("popsize", "rep", "year", "maxage", "fst")
        
        #dapc matrix 
        dpcmatrix = matrix(nrow = (simyears + 1), ncol = 6)
        colnames(dpcmatrix) = c("popsize", "rep", "year", "maxage", "OintoO", "CintoO")
        
        #evanno matrix 
        evamatrix = NULL #first three cols will be rep, year, and max age, then penalized log likelihood for K=1, K=2, ... K=structK 
        
        #simulate popualtion over years
        for(g in 1:(simyears + delay)){
          #####assign alleles based on initial starting frequencies in year 1####
          if(g==1){
            #set alleles/allele frequencies to starting values
            alleles = originalalleles
            allelefreqs = originalallelefreqs
            
            #initiate population, assign genotypes
            current = matrix(nrow=1000, ncol=numloci*2)
            lcols = seq(1,numloci*2, 2)
            locus = 1
            for(l in lcols){
              #if only one allele, set allele for all alleles/individuals
              if(length(alleles[[locus]])==1){
                current[, l]   = alleles[[locus]]
                current[, l+1] = alleles[[locus]]
              #otherwise, sample alleles randomly, proportional to their frequencies
              }else{
                current[,l]   = sample(x=alleles[[locus]], size=nrow(current), prob=allelefreqs[[locus]], replace=TRUE)
                current[,l+1] = sample(x=alleles[[locus]], size=nrow(current), prob=allelefreqs[[locus]], replace=TRUE)
              }
              locus = locus + 1
            }
            
            #assign age randomly, from 0 to max age allowed
            age = sample(seq(0, agecap[a], 1), size = nrow(current), replace=TRUE)
            current = cbind(current, age)
            
            #set all individuals as alive (==1)
            alive   = rep(1, nrow(current))
            current = cbind(current, alive)
            
            #retain original source popualtion
            original = current
            remove(current)
          }
          
          ####allow large pop to evolve####
          neworiginal = Advance(original, numloci, a, lcols)
          original  = neworiginal[[2]]
          poriginal = neworiginal[[1]]
          remove(neworiginal)
          
          #####create subsetted population#### 
          if(g==delay){
            current = original[sample(seq(1, nrow(original), 1), popsize),,drop=FALSE]
          }

          #### allow subsetted pop to age; measure population parameters of interest ####
          if(g>=delay){
            #age population
            newcurrent = Advance(current, numloci, a, lcols)
            previous = newcurrent[[1]]
            current  = newcurrent[[2]]
            remove(newcurrent)
            
            #first, create variables to hold current year and original population data
            holdoriginal = original
            holdcurrent  = current
            
            #sub-sample 50 individuals from each population (to represent logistical realities)
            if(nrow(current) >= 50){
              tokeep = sample(seq(1,nrow(current), 1), size=50, replace=FALSE)
              current = current[tokeep,]
            }
            if(nrow(original) >= 50){
              tokeep = sample(seq(1,nrow(original), 1), size=50, replace=FALSE)
              original = original[tokeep,]
            }
            
            ####heterozygosity####
            lcols = seq(1,numloci*2, 2)
            hets = matrix(nrow=nrow(current), ncol=1)
            hets[,1] = 0
            for(ll in lcols){
              hrows = which(current[,ll] != current[,ll+1])
              hets[hrows,1] = hets[hrows,1] + 1
            }
            hets[,1] = hets[,1]/(numloci)
            
            #record H, year, and replicate number
            hetmatrix[(g-(delay-1)),1] = popsize[n]     #simulated population size
            hetmatrix[(g-(delay-1)),2] = j              #replicate
            hetmatrix[(g-(delay-1)),3] = g-(delay-1)    #year
            hetmatrix[(g-(delay-1)),4] = agecap[a]      #max age
            hetmatrix[(g-(delay-1)),5] = mean(hets[,1]) #het
            
            ####fst####
            #begin setting up genepop file
            write.table("#no comment", "../output/tempgen.gen", col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")

            #iterate over loci and combine alleles from 2 columns (1 locus) into 1 column
            tempfst  = NULL
            firstfst = NULL
            locus    = 1
            lcols    = seq(1,numloci*2, 2)
            for(h in lcols){
              write.table(lnames[locus], "../output/tempgen.gen", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="\t")
              #pop 1
              firstfst = cbind(firstfst, paste(original[,h], original[,h+1], sep=""))
              #pop 2
              tempfst  = cbind(tempfst, paste(current[,h], current[,h+1], sep=""))
              locus = locus + 1
            }

            #fix sample names, add ,
            names    = paste(seq(1,nrow(firstfst), 1), rep(",", nrow(firstfst)), sep="")
            firstfst = cbind(names, firstfst)
            names    = paste(seq(1,nrow(tempfst), 1), rep(",", nrow(tempfst)), sep="")
            tempfst  = cbind(names, tempfst)

            #write POP, then pop1 data, POP, then pop2 data
            write.table("POP", "../output/tempgen.gen", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="\t")
            write.table(firstfst, "../output/tempgen.gen", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="\t")
            write.table("POP", "../output/tempgen.gen", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="\t")
            write.table(tempfst, "../output/tempgen.gen", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE, sep="\t")

            #generate fst estimates
            fst = diffCalc("../output/tempgen.gen", "../output/tempgenout.gen", fst=TRUE, pairwise=FALSE, bs_locus=FALSE, bs_pairwise=FALSE, boots=0, ci_type="loci", para=FALSE)

            #record fst, year, and replicate number
            fstmatrix[(g-(delay-1)),1] = popsize[n]                           #simulated population size
            fstmatrix[(g-(delay-1)),2] = j                                    #replicate
            fstmatrix[(g-(delay-1)),3] = g-(delay-1)                          #year
            fstmatrix[(g-(delay-1)),4] = agecap[a]                            #max age
            fstmatrix[(g-(delay-1)),5] = fst$std_stats$Fst[length(lnames)+1]  #fst

            ####dapc####
            #prepare data for dapc procedure
            todapc = rbind(original[,1:(numloci*2)], current[,1:(numloci*2)])
            groups = c(rep("O", nrow(original)), rep("C", nrow(current)))

            #run dapc
            dapc.out = dapc(x=todapc, grp=groups, n.pca=30, n.da=3)

            #extract and record probability of original pop individuals assigned to original pop and current indv. to original
            dpcmatrix[(g-(delay-1)),1] = popsize[n]                           #simulated population size
            dpcmatrix[(g-(delay-1)),2] = j                                    #replicate
            dpcmatrix[(g-(delay-1)),3] = g-(delay-1)                          #year
            dpcmatrix[(g-(delay-1)),4] = agecap[a]                            #max age
            dpcmatrix[(g-(delay-1)),5] = mean(dapc.out$posterior[1:nrow(original), 2])                                   #Original into original
            dpcmatrix[(g-(delay-1)),6] = mean(dapc.out$posterior[(nrow(original)+1):(nrow(current)+nrow(original)), 2])  #Current into original

            ####structure####
            #only run durring years specified in levels
            if(!is.na(match(g, levels))){
              #generate file to putinto structure (borrowing from above)
              original.structure = cbind(seq(1, nrow(original), 1), rep("O", nrow(original)), original[,1:(numloci*2)])
              current.structure  = cbind(seq(nrow(original)+1, nrow(original)+nrow(current), 1), rep("C", nrow(current)), current[,1:(numloci*2)])
              into.structure = rbind(original.structure, current.structure)

              #generate file to update column names of into.structure and update the names
              into.structure.colnames = c("ID", "loc")
              for(al in 1:length(lnames)){
                into.structure.colnames = c(into.structure.colnames, lnames[al], lnames[al])
              }
              colnames(into.structure) = into.structure.colnames

              #generate file for strucutre
              structure.file = gtypes(gen.data=into.structure,id.col=1,strata.col=2,locus.col=3,dna.seq=NULL,description=NULL,delete.missing.strata=TRUE,code.start=1)

              #run structure
              struct.out = structure.run(structure.file, k.range=c(1:structK), num.k.rep=2, numreps=100000, burnin=100000, noadmix=TRUE, freqscorr=TRUE, num.cores=2)
              E = structure.evanno(struct.out,plot=FALSE)

              #record mean penealized log likelihood, year, and replicate number
              torecord = rbind(as.matrix(popsize[n]), as.matrix(j), as.matrix(g-(delay-1)), as.matrix(agecap[a]), as.matrix(E$mean.ln.k))
              torecord = t(torecord)
              evamatrix = rbind(evamatrix, torecord)
            }

            #replace sub-setted variables with full datasets saved above
            original = holdoriginal
            current  = holdcurrent
          }
        }#g
        current = previous = original = NULL
        write.table(hetmatrix,paste(outdir, "het.csv", sep = ""),sep=",",col.names=FALSE, row.names=FALSE, append=TRUE)
        write.table(fstmatrix,paste(outdir, "fst.csv", sep = ""),sep=",",col.names=FALSE, row.names=FALSE, append=TRUE)
        write.table(dpcmatrix,paste(outdir, "dpc.csv", sep = ""),sep=",",col.names=FALSE, row.names=FALSE, append=TRUE)
        write.table(evamatrix,paste(outdir, "eva.csv", sep = ""),sep=",",col.names=FALSE, row.names=FALSE, append=TRUE)
      }#j
    }#a
  }#n
}#RunSims
