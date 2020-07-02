###Community assembly processes 
##Neutral model 
#Adam Burns - 2/10/2015
#aburns2@uoregon.edu
#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
#comun: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.

sncm.fit<- function(comun, stats=TRUE){
	require(minpack.lm)
	require(Hmisc)
	require(stats4)
	
	options(warn=-1)

	#Calculate the number of individuals per community
	N <- round(mean(apply(comun, 1, sum)))
	
	#Calculate the average relative abundance of each taxa across communities
	p.m <- apply(comun, 2, mean)
	p.m <- p.m[p.m != 0]
	p <- p.m/N
	
	#Calculate the occurrence frequency of each taxa across communities
	comun.bi <- 1*(comun>0)
	freq <- apply(comun.bi, 2, mean)
	freq <- freq[freq != 0]

	#Combine
	C <- merge(p, freq, by=0)
	C <- C[order(C[,2]),]
	C <- as.data.frame(C)
	C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
	p <- C.0[,2]
	freq <- C.0[,3]
	names(p) <- C.0[,1]
	names(freq) <- C.0[,1]

	#Calculate the limit of detection
	d = 1/N

	##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
	m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
	m.ci <- confint(m.fit, 'm', level=0.95)
	
	##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
	sncm.LL <- function(m, sigma){
		R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
		R = dnorm(R, 0, sigma)
		-sum(log(R))
	}
	m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
	
	##Calculate Akaike's Information Criterion (AIC)
	aic.fit <- AIC(m.mle, k=2)
	bic.fit <- BIC(m.mle)

	##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
	freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
	Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
	
	pred.ci <- binconf(freq.pred*nrow(comun), nrow(comun), alpha=0.05, method="wilson", return.df=TRUE)
	
	##Calculate AIC for binomial model
	bino.LL <- function(mu, sigma){
		R = freq - pbinom(d, N, p, lower.tail=FALSE)
		R = dnorm(R, mu, sigma)
		-sum(log(R))
	}
	bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
	
	aic.bino <- AIC(bino.mle, k=2)
	bic.bino <- BIC(bino.mle)
	
	##Goodness of fit for binomial model
	bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
	Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
	RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))

	bino.pred.ci <- binconf(bino.pred*nrow(comun), nrow(comun), alpha=0.05, method="wilson", return.df=TRUE)
	
		##Results
	if(stats==TRUE){
		fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), RMSE=numeric(), RMSE.bino=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
		fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, Rsqr, Rsqr.bino, RMSE, RMSE.bino, aic.fit, bic.fit, aic.bino, bic.bino, N, nrow(comun), length(p), d)
		return(fitstats)
	} else {
		A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
		A <- as.data.frame(A)
		colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
		B <- A[order(A[,1]),]
		return(B)
	}
}

a<-sncm.fit(core_otu_absolute,stats=F)#stats=F
a
#model parameters 
a<-sncm.fit(core_otu_absolute,stats=TRUE)#stats=T

##Null model 
phylo: Phylogenetic tree of each OTU
comun: A community table with samples as rows and OTUs as columns. 

##Beta_NTI
Beta_NTI<-function(phylo,comun,beta.reps=999){
   require(picante)

   comun=t(comun)
   match.phylo.comun = match.phylo.data(phylo, t(comun))
   beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.comun$data),cophenetic(match.phylo.comun$phy),abundance.weighted=T))

   rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.comun$data),ncol(match.phylo.comun$data),beta.reps))
   for (rep in 1:beta.reps) {
       rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.comun$data),taxaShuffle(cophenetic(match.phylo.comun$phy)),abundance.weighted=T,exclude.conspecifics = F))
       print(c(date(),rep))
   }

   weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.comun$data),ncol=ncol(match.phylo.comun$data))
   for(columns in 1:(ncol(match.phylo.comun$data)-1)) {
          for(rows in (columns+1):ncol(match.phylo.comun$data)) {
                  rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
                  weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
                  rm("rand.vals")
          }
   }
   
  rownames(weighted.bNTI) = colnames(match.phylo.comun$data);
  colnames(weighted.bNTI) = colnames(match.phylo.comun$data);
  return(as.dist(weighted.bNTI))
}

#data
phylo=read.tree("OTUs.tre")
#screen for core otu tree 
phylo<-prune.sample(core_otu_absolute,phylo)
comun=t(core_otu_absolute)
core_index_beta_NTI=Beta_NTI(phylo,comun,beta.reps=999)

##RC_bray
raup_crick= function(comun, reps=999){
  require(ecodist) 
   
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(comun)
  gamma<-ncol(comun)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(comun), row.names(comun)))
  ##make the comun matrix into a new, pres/abs. matrix:
  ceiling(comun/max(comun))->comun.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(comun.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(comun, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(comun)-1)){
    for(null.two in (null.one+1):nrow(comun)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(comun.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(comun[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(comun.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(comun[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.comun = rbind(com1,com2); # null.comun;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.comun,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(comun[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop

  results<-as.dist(results)
  return(results)
}

#data
comun=core_otu_absolute 
index_RC=raup_crick(comun, reps=999)

##Niche breadth 
#comun: A community table with samples as rows and taxa as columns.
comun=core_otu_absolute
Com.niche <- function(comun, stats=TRUE){
    require(spaa)
    comun<-comun[,colSums(comun)>0]
    B<-niche.width(comun,method="levins")
    B_com<-1:nrow(comun)
    for(i in 1:nrow(comun)){
       a<-comun[i,]
       a<-a[a>0]
       B_com[i]<-mean(as.numeric(B[,names(a)]))
    }
    return(B_com)
}

habitat_niche_breadth=as.matrix(Com.niche(comun))
rownames(habitat_niche_breadth)=rownames(comun)
colnames(habitat_niche_breadth)="habitat_niche_breadth"
core_otu_bcom=habitat_niche_breadth