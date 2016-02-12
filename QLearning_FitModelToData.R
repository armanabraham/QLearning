

datToFit <- droplevels(subset(glmData, SubjectID=='s011' & SessionID %in% c(1) ))
datToFit <- droplevels(subset(glmData, Condition==1)) # Condition 1 of the experiment
learningCoefs <- QLearning_FitModelToData(datToFit, constrainedFit=FALSE, fittingMethod='BFGS')

pp <- subset(learningCoefs, (beta < 10) & !(Condition %in% c(9,10,11,12)))
#pp <- learningCoefs
ggplot(pp) + 
  #scale_y_log10() + 
  #scale_x_log10() +
  geom_point(aes(x=alpha, y=beta), alpha=0.7) + 
  theme_few() + 
  facet_grid(~Condition, labeller=ConditionLabels)

## To plot means
ppM <- ddply(pp, .(SubjectID, Condition), summarise, 
             seAlpha=sd(alpha)/sqrt(length(alpha)),
             seBeta=sd(beta)/sqrt(length(beta)),
             alpha=mean(alpha), 
             beta=mean(beta))
ppM$seAlpha[is.na(ppM$seAlpha)] <- 0
ppM$seBeta[is.na(ppM$seBeta)] <- 0

ggplot(ppM, aes(x=alpha, y=beta)) + 
  theme_few() + 
  #scale_y_log10() + 
  #scale_x_log10() +
  geom_errorbarh(aes(xmin=alpha-seAlpha, xmax=alpha+seAlpha), color="grey70", size=0.2) + 
  geom_errorbar(aes(ymin=beta-seBeta, ymax=beta+seBeta), color="grey70", size=0.2) +
  geom_point(aes(x=alpha, y=beta, fill=SubjectID),size=3.5, shape=21, color="white") + 
  scale_fill_manual(values=colorRampPalette(brewer.pal(12, "Set1"))(15)) + 
  facet_grid(~Condition, labeller=ConditionLabels)

###################################################################
## Fitting the model to the data. The fit is done using grid search
## on a range of alphas and betas. The alphas and betas that produce the 
## lowest log-likelihood are selected as model parameters
QLearning_FitModelToData <- function(inputData, ## glmData
                                     constrainedFit=FALSE, 
                                     fittingMethod='BFGS')  # When using constrained fit, the default fitting method is L-BFGS-B to work. 
{
  # To store parameter estimates of the model, together with information about 
  # each subject, condition and run
  learningCoefs <- data.frame(SubjectID=character(), SessionID=character(), Condition=as.numeric(), alpha=double(), beta=double(), logLhood=double())
  
  # Using grid search. Can be slower than using random start points
  alphas <- seq(0.1, 1, 0.1)
  betas <- seq(0, 10)
  initValues <- expand.grid(alphas=alphas, betas=betas)
  print(paste('Performing grid search with', nrow(initValues), 'initial values per fit'))
  
  for (ixSubject in levels(droplevels(inputData$SubjectID))) {
    oneSubjectData <- subset(inputData, inputData$SubjectID == ixSubject)
    for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
      oneSessionData <- subset(oneSubjectData, oneSubjectData$SessionID == ixSession)
      print(paste("Subject", ixSubject, "session", ixSession, "condition", unique(oneSessionData$Condition)))
      # Use unconstrained minimization which is done on a grid of values specified by initValues
      # This is without limiting the search range. Seems to work well with grid search, but standalone can sometimes produce spurious results
      
      if (!constrainedFit) 
        # Default is unconstrained fit using grid search
        res <- apply(initValues, 1, function(par) optim(fn=QLearning_Model1, par, control=list(maxit=5000), method=fittingMethod, oneRunData=oneSessionData))
      else
        # Constrained fitting using grid search
        res <- apply(initValues, 1, function(par) optim(fn=QLearning_Model1, par, lower=c(0.001,0.001), upper=c(1,10), control=list(maxit=5000), oneRunData=oneSessionData, hessian=T, method="L-BFGS-B"))

      # Extract parameters of the fitted models
      fitVals <- t(sapply(1:nrow(initValues), function(x) res[[x]]$par))
      fitValsAndLogLike <- t(sapply(1:nrow(initValues), function(x) c(res[[x]]$par, LogLike=res[[x]]$value)))
      # Show the best fitted values of alpha and beta based on smallest log-likelihood fit
      bestParams <- as.numeric(fitValsAndLogLike[which.min(fitValsAndLogLike[,3]),])
      learningCoefs <- rbind(learningCoefs, data.frame(SubjectID=ixSubject, SessionID=ixSession, Condition=unique(oneSessionData$Condition), alpha=bestParams[1], beta=bestParams[2], logLhood=bestParams[3]))
      #browser()
    }
  }
  return(learningCoefs)
}

###################################################################
## Fitting Model3 to the data. The fit is done using grid search
## on a range of corrAlphas, failAlphas, and betas. 
## The alphas and betas that produce the 
## lowest log-likelihood are selected as model parameters
QLearning_FitModel3ToData <- function(inputData, ## glmData
                                     constrainedFit=FALSE, 
                                     fittingMethod='BFGS')  # When using constrained fit, the default fitting method is L-BFGS-B to work. 
{
  # To store parameter estimates of the model, together with information about 
  # each subject, condition and run
  learningCoefs <- data.frame(SubjectID=character(), SessionID=character(), Condition=as.numeric(), alpha=double(), beta=double(), logLhood=double())
  
  
  
  ### TESTING MODEL 3: Fitting corrAlpha, failAlpha, and beta
  # Using grid search. Can be slower than using random start points
  corrAlphas <- seq(0.1, 1.6, 0.2)
  failAlphas <- seq(0.1, 1.6, 0.2)
  betas <- seq(0, 10, 2)
  initValues <- expand.grid(corrAlphas=corrAlphas, failAlphas=failAlphas, betas=betas)
  print(paste('Performing grid search with', nrow(initValues), 'initial values per fit'))
  
  for (ixSubject in levels(droplevels(inputData$SubjectID))) {
    oneSubjectData <- subset(inputData, inputData$SubjectID == ixSubject)
    for (ixSession in levels(droplevels(oneSubjectData$SessionID))) {
      oneSessionData <- subset(oneSubjectData, oneSubjectData$SessionID == ixSession)
      print(paste("Subject", ixSubject, "session", ixSession, "condition", unique(oneSessionData$Condition)))
      # Use unconstrained minimization which is done on a grid of values specified by initValues
      # This is without limiting the search range. Seems to work well with grid search, but standalone can sometimes produce spurious results
      
      if (!constrainedFit) 
        # Default is unconstrained fit using grid search
        res <- apply(initValues, 1, function(par) optim(fn=QLearning_Model3, par, control=list(maxit=5000), method=fittingMethod, oneRunData=oneSessionData))
      else
        # Constrained fitting using grid search
        res <- apply(initValues, 1, function(par) optim(fn=QLearning_Model3, par, lower=c(0.001,0.001, 0.001), upper=c(1.5, 1.5, 15), control=list(maxit=5000), oneRunData=oneSessionData, hessian=T, method="L-BFGS-B"))
      
      # Extract parameters of the fitted models
      fitVals <- t(sapply(1:nrow(initValues), function(x) res[[x]]$par))
      fitValsAndLogLike <- t(sapply(1:nrow(initValues), function(x) c(res[[x]]$par, LogLike=res[[x]]$value)))
      # Show the best fitted values of alpha and beta based on smallest log-likelihood fit
      bestParams <- as.numeric(fitValsAndLogLike[which.min(fitValsAndLogLike[,dim(fitValsAndLogLike)[2]]),])
      learningCoefs <- rbind(learningCoefs, data.frame(SubjectID=ixSubject, SessionID=ixSession, Condition=unique(oneSessionData$Condition), CorrAlpha=bestParams[1], FailAlpha=bestParams[2], Beta=bestParams[3], logLhood=bestParams[4]))
      #browser()
    }
  }
  return(learningCoefs)
}


