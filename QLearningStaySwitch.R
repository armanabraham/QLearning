## Natural bias
oneSessionData <- subset(glmData, SubjectID=="s001" & SessionID==1)

## Success/switch bias induction
## Very poor fit in case of this data
oneSessionData <- subset(glmData, SubjectID=="s001" & SessionID==34)

## Success/stay bias induction
oneSessionData <- subset(glmData, SubjectID=="s001" & SessionID==30)


nTestVals <- 10
upperBeta <- 10
alphas <- runif(nTestVals)
betas <- runif(nTestVals) * upperBeta
randInitialValues <- cbind(alphas, betas)
res <- apply(randInitialValues, 1, function(par) optim(fn= QLearningStaySwitch, par, lower=c(0,0), upper=c(1,upperBeta), control=list(maxit=5000, fnscale=-1), oneSessionData=oneSessionData, hessian=T, method="L-BFGS-B"))
fitVals <- t(sapply(1:nTestVals, function(x) res[[x]]$par))
QLearningStaySwitch(fitVals[1,],oneSessionData, plotProbabilities=T)


#################################################################################################
## *** FUNCTION to fit model parameters using regulirized logistic regression (glmnet)
QLearningStaySwitch <- function (params,				# alpha and beta params to be fitted 
								 oneSessionData, 		# subject responses
								 plotProbabilities=F	# Plot pRight and pLeft 	
									) {
	# Unpack params
	alpha <- params[1]
	beta <- params[2]
	
	## Initial Q (value) function is set to 0 for both left and right choices
	nTrials <- nrow(oneSessionData)
	initValues <- t(rep(0, nTrials))
	Qvalue <- data.frame(WSLS=rep(0, nTrials), NonWSLS=rep(0, nTrials))
	## Array showing if on that trial subject used WSLS. 1 did, -1 did not 
	isTrialWSLS <- as.integer(rep(0, nTrials))
	
	## Starting from trial 2, update the value function Q
	## using the Q-Learning rule: Q(t+1) = Q(t) + alpha*(R(T) - Q(t))
	for (ixTrial in 2:nTrials) {
		## Was it success stay
		if ((oneSessionData$CorrIncorr[ixTrial-1]==1) && (oneSessionData$VisualField[ixTrial-1] == oneSessionData$VisualField[ixTrial]) ) {
			Qvalue$WSLS[ixTrial] <- Qvalue$WSLS[ixTrial-1] + alpha * (1 - Qvalue$WSLS[ixTrial-1])
			Qvalue$NonWSLS[ixTrial] <- Qvalue$NonWSLS[ixTrial-1]
			isTrialWSLS[ixTrial] <- 1
			#print("Win Stay")
		} else if ((oneSessionData$CorrIncorr[ixTrial-1]==0) && (oneSessionData$VisualField[ixTrial-1] != oneSessionData$VisualField[ixTrial])) {
			Qvalue$WSLS[ixTrial] <- Qvalue$WSLS[ixTrial-1] + alpha * (1 - Qvalue$WSLS[ixTrial-1])
			Qvalue$NonWSLS[ixTrial] <- Qvalue$NonWSLS[ixTrial-1]
			isTrialWSLS[ixTrial] <- 1
			#print("Fail Switch")
		} else if ((oneSessionData$CorrIncorr[ixTrial-1]==1) && (oneSessionData$VisualField[ixTrial-1] != oneSessionData$VisualField[ixTrial]) ) {
			## Win switch
			Qvalue$NonWSLS[ixTrial] <- Qvalue$NonWSLS[ixTrial-1] + alpha * (1 - Qvalue$NonWSLS[ixTrial-1])
			Qvalue$WSLS[ixTrial] <- Qvalue$WSLS[ixTrial-1]
			isTrialWSLS[ixTrial] <- -1
			#print("** Win Switch")
		} else if ((oneSessionData$CorrIncorr[ixTrial-1]==0) && (oneSessionData$VisualField[ixTrial-1] == oneSessionData$VisualField[ixTrial])) {
			Qvalue$NonWSLS[ixTrial] <- Qvalue$NonWSLS[ixTrial-1] + alpha * (1 - Qvalue$NonWSLS[ixTrial-1])
			Qvalue$WSLS[ixTrial] <- Qvalue$WSLS[ixTrial-1]
			isTrialWSLS[ixTrial] <- -1
			#print("** Fail Stay")				
		} else {
			Qvalue[ixTrial,] <- Qvalue[ixTrial-1,]
			#print("**** NONE REGISTERED ****")
		}
	}
	## Compute probabilities using Softmax rule based on reward values

	probs <- data.frame(WSLS=rep(NaN, nTrials), NonWSLS=rep(NaN, nTrials))
	probs$WSLS <- exp(Qvalue$WSLS*beta) / (exp(Qvalue$NonWSLS*beta) + exp(Qvalue$WSLS*beta))
	probs$NonWSLS <- exp(Qvalue$NonWSLS*beta) / (exp(Qvalue$NonWSLS*beta) + exp(Qvalue$WSLS*beta))
	browser()
	## Plot data
	#browser()
	if (plotProbabilities) { 
		#mProbs <- melt(probs)
		#g1 <- ggplot(data=mProbs, aes(x=1:nTrials, y=value, color=variable, group=variable)) + geom_line()
		g1 <- ggplot(data=probs, aes(x=1:length(WSLS))) + theme_few() + geom_line(aes(y=WSLS), color="Blue") + geom_line(aes(y=NonWSLS), color="Red") + ylab("p") + xlab("Trial")
		g2 <- ggplot(data=Qvalue, aes(x=1:length(WSLS))) + theme_few() + geom_line(aes(y= WSLS), color="Blue") + geom_line(aes(y=NonWSLS), color="Red") + ylab("V") + xlab("Trial")
		g3 <- grid.arrange(g2,g1)
		print(g3)
	}
	#browser()
	## Compute from data if trial 
	
	
	lhood <- LikelihoodWSLS(isTrialWSLS, probs$WSLS)
	#browser()
	print(paste(alpha, beta, lhood))
	return(lhood)
}


###########################################
## Returns log likelihood 
LikelihoodWSLS <- function(trialType, 		# Is trial WSLS (1) or not (-1)
					   WSLS	# Probability of WSLS
						) {
	-sum(log(ifelse(trialType ==1, WSLS, 1 - WSLS)))
}
