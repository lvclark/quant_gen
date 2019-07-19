# Lindsay V. Clark
# November 15, 2016
# Feel free to modify and distribute.  This R function comes with no warranty.

# Generalized function to estimate genetic correlation between traits measured on the same field
# plots (the same plants).
# from equation 9 of Howe et al. (2000; doi:10.1007/s001220051525)
# see also http://www.knowledgebank.irri.org/ricebreedingcourse/Lesson_8_Correlations_among_traits__implications_for_screening.htm

# formula: Should start with "Response ~" and then include variables to be used as predictors.
#          Formatted for lme4.  An example would be "Response ~ (1|Genotype) + (1|Loc/Rep) + (1|Genotype:Loc)"
# data:    Data frame containing genotype IDs, traits, and any predictor variables.
# traitNames: Names of the columns in data containing traits between which we want to estimate genetic correlation.
# genIDName:  Name of column in data containing genotype ID (for estimation of genetic components of
#             variance and covariance).  In the above example this would be "Genotype".
# subset:  Subset of data to analyze. (Numeric, boolean, or character vector indicating rows to keep.)

genCorrTraits <- function(formula, data, traitNames, genIDName, subset = 1:dim(data)[1]){
  require(lme4)
  nTrait <- length(traitNames)
  if(nTrait < 2){
    stop("At least two traitNames needed.")
  }
  if(!all(traitNames %in% names(data))){
    stop("All traitNames need to be column names of data.")
  }
  if(length(genIDName) != 1){
    stop("Only one genIDName allowed.")
  }
  if(!genIDName %in% names(data)){
    stop("genIDName not found in column names of data.")
  }
  data <- data[subset,]
  nObs <- dim(data)[1]
  
  # fill in variances along the diagonal only if there is no missing data at all
  dovar <- !any(is.na(data[,traitNames]))
  
  # manipulate formula to build new ones for linear models
  OriginalTerms <- attr(terms(formula), "term.labels")
  for(i in grep("|", OriginalTerms, fixed = TRUE)){ # fix formatting of random effects
    OriginalTerms[i] <- paste("(", OriginalTerms[i], ")", sep = "")
  }
  
  # get genetic covariance and correlation for pairs of traits
  gencov <- gencor <- matrix(NA, nrow = nTrait, ncol = nTrait,
                             dimnames = list(traitNames, traitNames))
  for(tr1 in traitNames[-nTrait]){
    for(tr2 in traitNames[(match(tr1, traitNames)+1):nTrait]){
      message(paste("Evaluating", tr1, "with", tr2))
      # get the genetic variances, just using plots where both traits were measured
      obsBothTrait <- !is.na(data[[tr1]]) & !is.na(data[[tr2]])
      if(sum(obsBothTrait) == 0) next
      if(dovar && !is.na(gencov[tr1, tr1])){
        var1 <- gencov[tr1, tr1] # don't recalculate variance if already done
      } else {
        thismodel <- lme4::lmer(reformulate(termlabels = OriginalTerms, response = tr1), data = data,
                                   subset = obsBothTrait)
        thisVC <- as.data.frame(VarCorr(thismodel))
        var1 <- thisVC$vcov[match(genIDName, thisVC$grp)]
        if(dovar) gencov[tr1, tr1] <- var1 # add variance to diagonal
      }
      if(dovar && !is.na(gencov[tr2, tr2])){
        var2 <- gencov[tr2, tr2]
      } else {
        thismodel <- lme4::lmer(reformulate(termlabels = OriginalTerms, response = tr2), data = data,
                                   subset = obsBothTrait)
        thisVC <- as.data.frame(VarCorr(thismodel))
        var2 <- thisVC$vcov[match(genIDName, thisVC$grp)]
        if(dovar) gencov[tr2, tr2] <- var2
      }
      
      # now add the two variables together and get the genetic variance for the new trait
      data$Response <- data[[tr1]] + data[[tr2]]
      thismodel <- lme4::lmer(reformulate(termlabels = OriginalTerms, response = "Response"), 
                                 data = data)
      thisVC <- as.data.frame(VarCorr(thismodel))
      var12 <- thisVC$vcov[match(genIDName, thisVC$grp)]
      # get the genetic covariance from these three variances
      thisGenCov <- (var12 - var1 - var2)/2
      gencov[tr1, tr2] <- thisGenCov
      gencov[tr2, tr1] <- thisGenCov
      # get genetic correlation
      thisGenCor <- thisGenCov/sqrt(var1*var2)
      gencor[tr1, tr2] <- thisGenCor
      gencor[tr2, tr1] <- thisGenCor
    }
  }
  
  return(list(genetic.covariance = gencov,
              genetic.correlation = gencor))
}
