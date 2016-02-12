library(ggplot2)
library(ggthemes)
library(Hmisc)

source('QLearning_Functions.R')
run <- QLearning_GenerateTrials(alpha=0.5, beta=1, nTrials=300, rightwardBias = 0.5)
run <- droplevels(subset(glmData, SubjectID=='s001' & SessionID==1))

# Quick fit. Don't trust this much, use more thorough fit below
optim(fn=QLearning_Model1, c(0.1, 0.1), control=list(maxit=5000), oneRunData=run)

system.time({print(optim(fn=QLearning_Model2, c(0.1, 0.1, 0.5, 0.5), oneRunData=run, method="BFGS"))})

nlm(f = QLearning_Model2, p=c(0.1, 0.1, 0.5, 0.5), oneRunData=run, steptol=1e-4,gradtol=1e-4)

nlm(f = QLearning_Model2, p=c(0.5, 1, 1, -1), oneRunData=run, steptol=1e-4,gradtol=1e-4)

optifix(c(0.1, 0.1, 1,0), fn=QLearning_Model2, oneRunData=run, fixed=c(FALSE, FALSE, TRUE, TRUE), method='BFGS')

#######################################################
# Fit model2 with corrReward and failReward fixed to 
# values obtained from the probabilistic choice model
pcmBiases <- subset(allWeights, Parameter %in% c('PrevCorr1', 'PrevFail1') & Regularized=='yes')
run <- droplevels(subset(glmData, SubjectID=='s001' & SessionID==6))
alphas <- seq(-1.1, 1.6, 0.8)
betas <- seq(0.1, 10, 1)
initValues <- expand.grid(alphas=alphas, betas=betas)
res <- apply(initValues, 1, function(par) optifix(c(0.1, 0.1, 0.49814377,-0.26585219), fn=QLearning_Model2, oneRunData=run, 
                                                  fixed=c(FALSE, FALSE, TRUE, TRUE), method='BFGS'))
fitVals <- t(sapply(1:nrow(initValues), function(x) res[[x]]$par))
fitValsAndLogLike <- t(sapply(1:nrow(initValues), function(x) c(res[[x]]$par, LogLike=res[[x]]$value)))
# Show the best fitted values of alpha and beta based on smallest log-likelihood fit
fitValsAndLogLike[which.min(fitValsAndLogLike[,3]),]




# Repeat fitting multiple times with multiple simulations with same learnign params
# To see the range of spread of fitted values
fitVals <- data.frame(alpha=as.numeric(), beta=as.numeric())
for (ixSim in 1:20) {
  run <- QLearning_GenerateTrials(alpha=0.1, beta=4, nTrials=500, rightwardBias = 0.5)
  fit <- optim(fn=QLearning_Model1, c(0.3, 0.5), control=list(maxit=5000), oneRunData=run)
  fitVals <- rbind(fitVals, data.frame(alpha=fit$par[1], beta=fit$par[2]))
}

#Qlearning_AlphaBetaRunningEstimate(run, trueAlpha=0.2, trueBeta=2)
#QLearning_Model1(c(0.2, 4),run, plotProbabilities=T)


#QLearning_ValidateAlphaBeta(nValidations=20, nRunTrials=300)

## Using ranodom alphas and betas
nTestVals <- 20
upperBeta <- 10
alphas <- runif(nTestVals)
betas <- runif(nTestVals) * upperBeta
initValues <- cbind(alphas, betas)

# Using grid search. Can be slower than using random start points
alphas <- seq(0.1, 1, 0.1)
betas <- seq(0, 10)
initValues <- expand.grid(alphas=alphas, betas=betas)

# This is by limiting the search range. Sometimes doesn't produce good results. 
res <- apply(initValues, 1, function(par) optim(fn=QLearning_OneRunLhood, par, lower=c(0.001,0.001), upper=c(1,10), control=list(maxit=5000), oneRunData=run, hessian=T, method="L-BFGS-B"))
# This is without limiting the search range. Seems to work well with grid search, but standalone can sometimes produce spurious results
res <- apply(initValues, 1, function(par) optim(fn=QLearning_OneRunLhood, par, control=list(maxit=5000), oneRunData=run))

fitVals <- t(sapply(1:nrow(initValues), function(x) res[[x]]$par))
fitValsAndLogLike <- t(sapply(1:nrow(initValues), function(x) c(res[[x]]$par, LogLike=res[[x]]$value)))
# Show the best fitted values of alpha and beta based on smallest log-likelihood fit
fitValsAndLogLike[which.min(fitValsAndLogLike[,3]),]

QLearning_Model1(fitVals[344,],run, plotProbabilities=T)
QLearning_Model1(c(0.06424065,1.51000304),run, plotProbabilities=T)


########
### TESTING MODEL 2: Fitting alpha, beta, corrRewards, failRewards
# Using grid search. Can be slower than using random start points
alphas <- seq(0.1, 1.2, 0.8)
betas <- seq(0.1, 3, 1)
corrRewards <- seq(-10, 10, 7)
failRewards <- seq(-10, 10, 7)
initValues <- expand.grid(alphas=alphas, betas=betas, corrRewards=corrRewards, failRewards=failRewards)
print(paste('Performing grid search with', nrow(initValues), 'initial values per fit'))

res <- apply(initValues, 1, function(par) optim(fn=QLearning_Model2, par, control=list(maxit=5000), oneRunData=run, method='BFGS'))
fitVals <- t(sapply(1:nrow(initValues), function(x) res[[x]]$par))
fitValsAndLogLike <- t(sapply(1:nrow(initValues), function(x) c(res[[x]]$par, LogLike=res[[x]]$value)))
bestParams <- fitValsAndLogLike[which.min(fitValsAndLogLike[,5]),]

########
### TESTING MODEL 3: Fitting corrAlpha, failAlpha, and beta
# Using grid search. Can be slower than using random start points
corrAlphas <- seq(0.1, 1.2, 0.3)
failAlphas <- seq(0.1, 1.2, 0.3)
betas <- seq(0, 10, 5)
initValues <- expand.grid(corrAlphas=corrAlphas, failAlphas=failAlphas, betas=betas)
print(paste('Performing grid search with', nrow(initValues), 'initial values per fit'))
res <- apply(initValues, 1, function(par) optim(fn=QLearning_Model3, par, control=list(maxit=5000), oneRunData=run, method='BFGS'))
fitVals <- t(sapply(1:nrow(initValues), function(x) res[[x]]$par))
fitValsAndLogLike <- t(sapply(1:nrow(initValues), function(x) c(res[[x]]$par, LogLike=res[[x]]$value)))
bestParams <- fitValsAndLogLike[which.min(fitValsAndLogLike[,dim(fitValsAndLogLike)[2]]),]



#optim(fn=QLearning_OneRunLhood, params, lower=c(0,0), upper=c(1,upperBeta), control=list(maxit=5000), oneSessionData=run, hessian=T, method="L-BFGS-B")


############################################################################
## Generate trials with a bunch of randomly selected alpha and beta values
## Then fit the model to estimate alpha and beta values from those generated 
## trials. Then plot true and estimated alpha/beta values as a scatterplot 
## and connect relevant pairs together using geom_segment to see how far or close 
## are true and estimated values based on values of alpha and beta
QLearning_ValidateAlphaBeta <- function(nValidations,		## Sample of random alpha and beta pairs to choose 
										nRunTrials=100			## Number of trials in a run that will be generated using sampled alpha/beta
										) {
											
	upperBeta <- 10	## The upper value of beta
	alphas <- runif(nValidations)	## alpha ranges from 0 to 1
	betas <- runif(nValidations) * upperBeta
	params <- data.frame(cbind(TrueAlpha=alphas, TrueBeta=betas))
	params$TrueLikelihood <- NaN
	params$EstimatedLikelihood <- params$EstimatedBeta <- params$EstimatedAlpha <- NaN
	
	## For each alpha/beta pair, run a simulation 
	for (ixValidation in 1:nValidations) {
		print(ixValidation)
		## Generate trials using a pair of alpha and beta
		simulatedTrials <- QLearning_GenerateTrials(alpha= params$TrueAlpha[ixValidation], beta=params$TrueBeta[ixValidation], nTrials=nRunTrials)
		params$TrueLikelihood[ixValidation] <- QLearning_Model1(c(params$TrueAlpha[ixValidation], params$TrueBeta[ixValidation]), oneSessionData=simulatedTrials)
		## Unconstrained with initial params set at 0.5 and 1
		res <- optim(fn=QLearning_Model1, c(0.5, 1), oneSessionData=simulatedTrials)
		## Unconstrained with initial params set at true alpha and beta
		#res <- optim(fn=QLearning_OneRunLhood, c(params$TrueAlpha[ixValidation], beta=params$TrueBeta[ixValidation]), oneSessionData=simulatedTrials)
		## Constrained search with initial params set at true alpha and beta
		#res <- optim(fn=QLearning_OneRunLhood, c(params$TrueAlpha[ixValidation], params$TrueBeta[ixValidation]), lower=c(0,0), upper=c(1,upperBeta), control=list(maxit=5000), oneSessionData=simulatedTrials, hessian=T, method="L-BFGS-B")
		
		params$EstimatedAlpha[ixValidation] <- res$par[1]
		params$EstimatedBeta[ixValidation] <- res$par[2]
		params$EstimatedLikelihood[ixValidation] <- res$value
	}
	ggplot(data=params) + 
			theme_few() + 
			scale_y_log10() +
		    xlab("Alpha") + 
		    ylab("Beta") +
			geom_segment(aes(x=TrueAlpha, y=TrueBeta, xend=EstimatedAlpha, yend=EstimatedBeta), color="grey90") + 
			geom_point(aes(x=TrueAlpha, y=TrueBeta), color="#2b83ba", size=3) + 
			geom_point(aes(x=EstimatedAlpha, y=EstimatedBeta), color="#fdae61", size=3) 
	
	axisMax <- max(c(params$TrueLikelihood, params$EstimatedLikelihood))
	ggplot(data=params) + 
			theme_few() + 
			#scale_y_log10() +
			#geom_segment(aes(x=TrueAlpha, y=TrueBeta, xend=EstimatedAlpha, yend=EstimatedBeta), color="grey90") + 
			geom_point(aes(x=TrueLikelihood, y= EstimatedLikelihood), color="blue", size=3) +
			ylim(limits=c(0, axisMax)) + 
			xlim(limits=c(0, axisMax)) + 
			coord_equal()
			#+ 
			#geom_point(aes(x=EstimatedAlpha, y=EstimatedBeta), color="red", size=3) 
			
	browser()
	# res <- apply(randInitialValues, 1, 
			# function(par) optim(fn=QLearning_OneRunLhood, par, lower=c(0,0), upper=c(5,50), control=list(maxit=5000), oneSessionData=run, hessian=T, method="L-BFGS-B"))
}



############################################################################
## Estimates learning rate on each trial by running 
## optimisation on a subset of trials 
Qlearning_AlphaBetaRunningEstimate <- function(run, 
									trueAlpha=NaN, 
									trueBeta=NaN) {
	nTrials <- nrow(run)
	runningEstimate <- data.frame(alpha=rep(NaN, nTrials), beta=rep(NaN, nTrials))
	for (ixTrial in 2:nTrials) {
		res <- optim(fn=QLearning_OneRunLhood, c(0.01, 2), oneSessionData=run[1:ixTrial,])
		runningEstimate$alpha[ixTrial] <- res$par[1]
		runningEstimate$beta[ixTrial] <- res$par[2]
		print(ixTrial)
	}
	gAlpha <- qplot(1:ixTrial, runningEstimate$alpha, geom="line") + theme_few()
	if (!is.nan(trueAlpha)) gAlpha <- gAlpha + geom_hline(yintercept=trueAlpha, linetype = "dashed")
	gBeta <- qplot(1:ixTrial, runningEstimate$beta, geom="line") + scale_y_log10() + theme_few()
	if (!is.nan(trueBeta)) gBeta <- gBeta + geom_hline(yintercept=trueBeta, linetype = "dashed")
	browser()
} 


