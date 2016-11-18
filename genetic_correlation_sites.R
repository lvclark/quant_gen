# Lindsay V. Clark
# November 18, 2016
# Feel free to modify and distribute.  This R function comes with no warranty.

# Function to find genetic correlation between traits measured on different plants,
# e.g. the same trait measured in different environments.
# See equation 8 of Howe et al. (2000 doi:10.1007/s001220051525)
# see also http://www.knowledgebank.irri.org/ricebreedingcourse/Lesson_8_Correlations_among_traits__implications_for_screening.htm

# formula: The formula for the linear model that will be used to get variance components for
#          estimating heritability within environment.  Should start with "Response" and then include
#          terms from the column names of data, formatted for lme4.  It is important to include
#          genotype (genIDName) and replicates in the model.  Environment (siteName) should not
#          be included in the model.
#          Example: Response ~ (1|Genotype) + (1|Rep)
#          Example: Response ~ (1|Genotype) + (1|Year/Rep)
#          Example: Response ~ (1|Genotype) + (1|Year/Rep) + (1|Genotype:Year)
# data: A data frame with columns for genotype or family IDs, environment IDs, replicates,
#       any other grouping factors, and phenotypes.
# genIDName: Name of column in data containing genotype ID (for estimation of genetic component of
#            variance).  In the above example this would be "Genotype".
# siteName: The name of the column listing the environment for each observation.  These are the
#           environments between which genetic correlations will be estimated.
# traitNames: Names of columns containing phenotype data.  For each phenotype listed in traitNames,
#             a matrix of genetic correlations between environments will be output.
# subset:  Subset of data to analyze. (Numeric, boolean, or character vector indicating rows to keep.)

genCorrSites <- function(formula, data, genIDName, siteName, traitNames, subset = 1:dim(data)[1]){
  require(lme4)
  if(!all(traitNames %in% names(data))){
    stop("All traitNames need to be column names of data.")
  }
  if(length(genIDName) != 1){
    stop("Only one genIDName allowed.")
  }
  if(!genIDName %in% names(data)){
    stop("genIDName not found in column names of data.")
  }
  if(length(siteName) != 1){
    stop("Only one sitName allowed.")
  }
  if(!siteName %in% names(data)){
    stop("siteName not found in column names of data.")
  }
  
  data <- data[subset,]
  nTrait <- length(traitNames)
  nObs <- dim(data)[1]
  # manipulate formula to build new ones for linear models
  OriginalTerms <- attr(terms(formula), "term.labels")
  for(i in grep("|", OriginalTerms, fixed = TRUE)){ # fix formatting of random effects
    OriginalTerms[i] <- paste("(", OriginalTerms[i], ")", sep = "")
  }
  # get site names
  sites <- as.character(unique(data[[siteName]]))
  sites <- sites[!is.na(sites)]
  nSites <- length(sites)
  
  # set up lists to contain output
  gencorrByTrait <- list()
  length(gencorrByTrait) <- nTrait
  names(gencorrByTrait) <- traitNames
  linemeansByTrait <- nindByTrait <- gencorrByTrait
  
  # loop through traits to do analysis
  for(tr in traitNames){
    # matrix to contain results
    gencorrByTrait[[tr]] <- nindByTrait[[tr]] <- matrix(NA, nrow = nSites, ncol = nSites,
                                  dimnames = list(sites, sites))
    # get line means
    linemeansByTrait[[tr]] <- tapply(data[[tr]], list(data[[genIDName]], data[[siteName]]),
                                     mean, na.rm = TRUE)
    # loop through pairs of sites
    for(s1 in sites[-nSites]){
      for(s2 in sites[(match(s1, sites) + 1):nSites]){
        # genotypes for which there is data for both sites for this trait
        theseGen <- dimnames(linemeansByTrait[[tr]])[[1]][!is.na(linemeansByTrait[[tr]][,s1]) &
                                                            !is.na(linemeansByTrait[[tr]][,s2])]
        nindByTrait[[tr]][s1,s2] <- nindByTrait[[tr]][s2,s1] <- length(theseGen)
        if(length(theseGen) == 0) next
        # subsets of data frame to use
        subset1 <- which(data[[genIDName]] %in% theseGen & data[[siteName]] == s1)
        subset2 <- which(data[[genIDName]] %in% theseGen & data[[siteName]] == s2)
        # skip if there are too few observations
        if(sum(!is.na(data[[tr]][subset1])) <= length(theseGen) |
           sum(!is.na(data[[tr]][subset2])) <= length(theseGen)){
          next
        }
        # run linear model to get variances
        tryCatch(
        model1 <- lmer(reformulate(termlabels = OriginalTerms, response = tr),
                       data = data, subset = subset1),
        warning = function(w){
          print(paste("Warning when estimating heritability for ", tr, s1, "with", s2))
          print(w)
        })
        tryCatch(
        model2 <- lmer(reformulate(termlabels = OriginalTerms, response = tr),
                       data = data, subset = subset2),
        warning = function(w){
          print(paste("Warning when estimating heritability for ", tr, s2, "with", s1))
          print(w)
        })
        grp1 <- ngrps(model1) # number of groups
        grp1 <- grp1[-grep(genIDName, names(grp1))] # eliminate genotype groups since we don't need those
        vc1 <- as.data.frame(VarCorr(model1)) # variances
        grp2 <- ngrps(model2)
        grp2 <- grp2[-grep(genIDName, names(grp2))]
        vc2 <- as.data.frame(VarCorr(model2))
        # genetic variance
        gvar1 <- vc1$vcov[match(genIDName, vc1$grp)]
        gvar2 <- vc2$vcov[match(genIDName, vc2$grp)]
        # residual variance divided by total number of replicates
        resvar1 <- vc1$vcov[match("Residual", vc1$grp)]/max(grp1)
        resvar2 <- vc2$vcov[match("Residual", vc2$grp)]/max(grp2)
        # any genotype interaction terms, divided by the appropriate coefficient
        intvar1 <- 0
        intvar2 <- 0
        for(trm in grep(paste(genIDName, ":", sep = ""), vc1$grp)){
          thisvar1 <- vc1$vcov[trm]
          thisvar2 <- vc2$vcov[trm]
          trmName <- substring(vc1$grp[trm], nchar(genIDName)+2, nchar(vc1$grp[trm]))
          if(!trmName %in% names(grp1)){
            stop("Error matching terms for genotype interaction effects.")
          }
          thiscoef1 <- grp1[termName]
          thiscoef2 <- grp2[termName]
          intvar1 <- intvar1 + thisvar1/thiscoef1
          intvar2 <- intvar2 + thisvar2/thiscoef2
        }
        # estimate heritability
        herit1 <- gvar1/(gvar1 + resvar1 + intvar1)
        herit2 <- gvar2/(gvar2 + resvar2 + intvar2)
        # estimate correlation between line means
        phencor <- cor(linemeansByTrait[[tr]][theseGen,s1],
                       linemeansByTrait[[tr]][theseGen,s2])
        # estimate genetic correlation
        gencor <- phencor/sqrt(herit1 * herit2)
        # add to output
        gencorrByTrait[[tr]][s1,s2] <- gencorrByTrait[[tr]][s2,s1] <- gencor
      }
    }
  }
  return(list(line.means = linemeansByTrait, genetic.correlations = gencorrByTrait,
              number.of.genotypes.used = nindByTrait))
}
