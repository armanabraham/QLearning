


## @knitr DefineFunctions
library(reshape)
library(ggplot2)
library(ggthemes)
library(gridExtra)

#source('~/proj/arman/Venus/DataCollection/Venus_DataAnalysisFunctions.R')

## Define colors for Left and Right choices
colChooseRight <- '#2b83ba'
colChooseLeft <- '#d7191c'

#' Model1: QLearning with alpha and beta
#' 
#' Fit QLearning model for given alpha and beta and return
#' -2*log-likelihood of the model 
#'
#' @export 
QLearning_Model1 <- function (params,				# alpha and beta params to be fitted 
                              oneRunData, 		# subject responses
                              plotProbabilities=F,	# Plot pRight and pLeft 
                              displayParams=F)	
{
  # Unpack params
  alpha <- as.numeric(params[1])
  beta <- as.numeric(params[2])
  
  # ## When using Rmalschains package for optimization
  # ## we need to receive oneRunData from environment 
  # ## because Rmalschains uses RCPP package
  # if (missing(oneRunData)) {
  # if (!is.null(env[["oneRunData"]])) {
  # oneRunData <- env[["oneRunData"]]
  # browser()
  # }
  # }
  
  ## Initial Q (value) function is set to 0 for both left and right choices
  nTrials <- nrow(oneRunData)
  initValues <- t(rep(0, nTrials))
  Qvalue <- data.frame(Left=rep(0, nTrials), Right=rep(0, nTrials))
  
  # Get Q value of the first trial
  if (oneRunData$CorrIncorr[1]) {
    # Because initial Q value was 0, we only assign 'alpha' to Q on first trial
    if (oneRunData$Response[1]==1) Qvalue$Left[1] <- alpha else Qvalue$Right[1] = alpha # dividing by 5 to avoid giving large alpha in the beginning
  }
  
  ## Starting from trial 2, update the value function Q
  ## using the Q-Learning rule: Q(t+1) = Q(t) + alpha*(R(T) - Q(t))
  for (ixTrial in 2:nTrials) {
    if (oneRunData$CorrIncorr[ixTrial]==1) {  ## If correct response
      if (oneRunData$Response[ixTrial]==1) { ## on left
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + alpha * (1 - Qvalue$Left[ixTrial-1])
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
        
      } else  {		## or was it on right
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + alpha * (1 - Qvalue$Right[ixTrial-1])
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
    } else {	## If response was not correct
      if (oneRunData$Response[ixTrial]==1) {
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + alpha * (0 - Qvalue$Left[ixTrial-1])
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
        
      } else {
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + alpha * (0 - Qvalue$Right[ixTrial-1])
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
      
    } 
  }
  
  #browser()
  ## Compute probabilities using Softmax rule based on reward values
  probs <- data.frame(pLeft=rep(NaN, nTrials), pRight=rep(NaN, nTrials))
  #probs$Left <- exp(beta * Qvalue$Left) / (exp(beta*Qvalue$Left) + exp(beta*Qvalue$Right))
  #probs$Right <- exp(Qvalue$Right*beta) / (exp(Qvalue$Left*beta) + exp(Qvalue$Right*beta))
  
  probs$pRight <- QLearning_Softmax(data.frame(Qvalue, beta))$pRight
  probs$pLeft <- 1-probs$pRight
  
  ## Plot data
  #browser()
  if (plotProbabilities) { 
    #browser()
    dev.new(width=10.851064, height=5.851064)
    mProbs <- cbind(probs, x=1:nTrials)
    mProbs <- melt(mProbs, id="x")
    ixFailTrials <- which(oneRunData$CorrIncorr==0)
    g1 <- ggplot(data=mProbs) + theme_few() + geom_line(aes(x=x, y=value, color=variable)) + 
      scale_color_manual(values=c(colChooseLeft, colChooseRight)) + 
      theme(legend.title=element_blank()) + 
      geom_vline(xintercept = ixFailTrials, size = 0.3, colour = "grey80", linetype = "dashed") + 
      #scale_y_continuous(expand=c(0,0)) + 
      scale_x_continuous(expand=c(0.005,0.01)) + 
      ylab("p") + xlab("Trial") 
    mQval <- cbind(Qvalue, x=1:nTrials)
    mQval <- melt(mQval, id="x")
    g2 <- ggplot(data=mQval) + theme_few() + geom_line(aes(x=x, y=value, color=variable)) + 
      scale_color_manual(values=c(colChooseLeft, colChooseRight)) + 
      theme(legend.title=element_blank()) + 
      geom_vline(xintercept = ixFailTrials, size = 0.3, colour = "grey80", linetype = "dashed") + 
      scale_x_continuous(expand=c(0.005,0.01)) + 
      ylab("V") + xlab("") 
    g3 <- QLearning_PlotSparklineOfRun(oneRunData) + coord_fixed(ratio = 15)
    grid.arrange(g3,g2,g1)
    #print(g4)
  }
  
  lhood <- QLearning_LogLikelihood(oneRunData$y, probs$pRight)
  if (displayParams) print(paste(alpha, beta, lhood))
  return( ifelse((is.nan(lhood) | is.infinite(lhood)), 100000, lhood) )
  ## Return large log likelihood, if NaN or Inf was computed. 
  ## NaN indicates there was likely division by 0. Inf means log likelihood was 
  ## too big to be computed. 
}

# Fit Q-Learning model to data given parameters of alpha, beta, CorrReward
## and failReward. The function returns logLikelihood of the model. 
QLearning_Model2 <- function (params,  			# alpha and beta params to be fitted 
                              oneRunData, 		# subject responses
                              plotProbabilities=F,	# Plot pRight and pLeft 
                              displayParams=F)	
{
  # Unpack params
  alpha <- as.numeric(params[1])
  beta <- as.numeric(params[2])
  corrReward <- as.numeric(params[3])
  failReward <- as.numeric(params[4])
  
  ## Initial Q (value) function is set to 0 for both left and right choices
  nTrials <- nrow(oneRunData)
  initValues <- t(rep(0, nTrials))
  Qvalue <- data.frame(Left=rep(0, nTrials), Right=rep(0, nTrials))
  
  # Get Q value of the first trial
  if (oneRunData$CorrIncorr[1]) {
    # Because initial Q value was 0, we only assign 'alpha' to Q on first trial
    if (oneRunData$Response[1]==1) Qvalue$Left[1] <- alpha else Qvalue$Right[1] = alpha # dividing by 5 to avoid giving large alpha in the beginning
  }
  
  ## Starting from trial 2, update the value function Q
  ## using the Q-Learning rule: Q(t+1) = Q(t) + alpha*(R(T) - Q(t))
  for (ixTrial in 2:nTrials) {
    if (oneRunData$CorrIncorr[ixTrial]==1) {  ## If correct response
      if (oneRunData$Response[ixTrial]==1) { ## on left
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + alpha * (corrReward - Qvalue$Left[ixTrial-1])
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
        
      } else  {		## or was it on right
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + alpha * (corrReward - Qvalue$Right[ixTrial-1])
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
    } else {	## If response was not correct
      if (oneRunData$Response[ixTrial]==1) {
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + alpha * (failReward - Qvalue$Left[ixTrial-1])
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
        
      } else {
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + alpha * (failReward - Qvalue$Right[ixTrial-1])
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
    } 
  }
  
  #browser()
  ## Compute probabilities using Softmax rule based on reward values
  probs <- data.frame(pLeft=rep(NaN, nTrials), pRight=rep(NaN, nTrials))
  #probs$Left <- exp(beta * Qvalue$Left) / (exp(beta*Qvalue$Left) + exp(beta*Qvalue$Right))
  #probs$Right <- exp(Qvalue$Right*beta) / (exp(Qvalue$Left*beta) + exp(Qvalue$Right*beta))
  
  probs$pRight <- QLearning_Softmax(data.frame(Qvalue, beta))$pRight
  probs$pLeft <- 1-probs$pRight
  
  #   ## Plot data
  #   #browser()
  #   if (plotProbabilities) { 
  #     #browser()
  #     dev.new(width=10.851064, height=5.851064)
  #     mProbs <- cbind(probs, x=1:nTrials)
  #     mProbs <- melt(mProbs, id="x")
  #     ixFailTrials <- which(oneRunData$CorrIncorr==0)
  #     g1 <- ggplot(data=mProbs) + theme_few() + geom_line(aes(x=x, y=value, color=variable)) + 
  #       scale_color_manual(values=c(colChooseLeft, colChooseRight)) + 
  #       theme(legend.title=element_blank()) + 
  #       geom_vline(xintercept = ixFailTrials, size = 0.3, colour = "grey80", linetype = "dashed") + 
  #       #scale_y_continuous(expand=c(0,0)) + 
  #       scale_x_continuous(expand=c(0.005,0.01)) + 
  #       ylab("p") + xlab("Trial") 
  #     mQval <- cbind(Qvalue, x=1:nTrials)
  #     mQval <- melt(mQval, id="x")
  #     g2 <- ggplot(data=mQval) + theme_few() + geom_line(aes(x=x, y=value, color=variable)) + 
  #       scale_color_manual(values=c(colChooseLeft, colChooseRight)) + 
  #       theme(legend.title=element_blank()) + 
  #       geom_vline(xintercept = ixFailTrials, size = 0.3, colour = "grey80", linetype = "dashed") + 
  #       scale_x_continuous(expand=c(0.005,0.01)) + 
  #       ylab("V") + xlab("") 
  #     g3 <- QLearning_PlotSparklineOfRun(oneRunData) + coord_fixed(ratio = 15)
  #     grid.arrange(g3,g2,g1)
  #     #print(g4)
  #   }
  
  lhood <- QLearning_LogLikelihood(oneRunData$y, probs$pRight)
  if (displayParams) print(paste(alpha, beta, lhood))
  return( ifelse((is.nan(lhood) | is.infinite(lhood)), 100000, lhood) )
  ## Return large log likelihood, if NaN or Inf was computed. 
  ## NaN indicates there was likely division by 0. Inf means log likelihood was 
  ## too big to be computed. 
}


#' Model3: QLearning with three parameters
#' 
#' Fitting parameters: corrAlpha, failAlpha, beta
#' Return: -2*log-likelihood of the model 
#' 
#' This model tries to find if subjects have different learning rates
#' for success and failure 
#' @export
#######################################################################
## Fitting parameters: corrAlpha, failAlpha, beta
## Fitting learning rate parameters seprately for success and fail trials 
QLearning_Model3 <- function (params,    		# alpha and beta params to be fitted 
                              oneRunData, 		# subject responses
                              plotProbabilities=F,	# Plot pRight and pLeft 
                              displayParams=F)	
{
  # Unpack params
  corrAlpha <- as.numeric(params[1])
  failAlpha <- as.numeric(params[2])
  beta <- as.numeric(params[3])
  
  ## Initial Q (value) function is set to 0 for both left and right choices
  nTrials <- nrow(oneRunData)
  initValues <- t(rep(0, nTrials))
  Qvalue <- data.frame(Left=rep(0, nTrials), Right=rep(0, nTrials))
  
  # Get Q value of the first trial
  if (oneRunData$CorrIncorr[1]) {
    # Because initial Q value was 0, we only assign 'alpha' to Q on first trial
    if (oneRunData$Response[1]==1) Qvalue$Left[1] <- corrAlpha else Qvalue$Right[1] = corrAlpha # dividing by 5 to avoid giving large alpha in the beginning
  }
  
  ## Starting from trial 2, update the value function Q
  ## using the Q-Learning rule: Q(t+1) = Q(t) + alpha*(R(T) - Q(t))
  for (ixTrial in 2:nTrials) {
    if (oneRunData$CorrIncorr[ixTrial]==1) {  ## If correct response
      if (oneRunData$Response[ixTrial]==1) { ## on left
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + corrAlpha * (1 - Qvalue$Left[ixTrial-1])
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
        
      } else  {		## or was it on right
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + corrAlpha * (1 - Qvalue$Right[ixTrial-1])
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
    } else {	## If response was not correct
      if (oneRunData$Response[ixTrial]==1) {
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + failAlpha * (0 - Qvalue$Left[ixTrial-1])
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
        
      } else {
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + failAlpha * (0 - Qvalue$Right[ixTrial-1])
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
    } 
  }
  
  #browser()
  ## Compute probabilities using Softmax rule based on reward values
  probs <- data.frame(pLeft=rep(NaN, nTrials), pRight=rep(NaN, nTrials))
  #probs$Left <- exp(beta * Qvalue$Left) / (exp(beta*Qvalue$Left) + exp(beta*Qvalue$Right))
  #probs$Right <- exp(Qvalue$Right*beta) / (exp(Qvalue$Left*beta) + exp(Qvalue$Right*beta))
  
  probs$pRight <- QLearning_Softmax(data.frame(Qvalue, beta))$pRight
  probs$pLeft <- 1-probs$pRight
  
  #   ## Plot data
  #   #browser()
  #   if (plotProbabilities) { 
  #     #browser()
  #     dev.new(width=10.851064, height=5.851064)
  #     mProbs <- cbind(probs, x=1:nTrials)
  #     mProbs <- melt(mProbs, id="x")
  #     ixFailTrials <- which(oneRunData$CorrIncorr==0)
  #     g1 <- ggplot(data=mProbs) + theme_few() + geom_line(aes(x=x, y=value, color=variable)) + 
  #       scale_color_manual(values=c(colChooseLeft, colChooseRight)) + 
  #       theme(legend.title=element_blank()) + 
  #       geom_vline(xintercept = ixFailTrials, size = 0.3, colour = "grey80", linetype = "dashed") + 
  #       #scale_y_continuous(expand=c(0,0)) + 
  #       scale_x_continuous(expand=c(0.005,0.01)) + 
  #       ylab("p") + xlab("Trial") 
  #     mQval <- cbind(Qvalue, x=1:nTrials)
  #     mQval <- melt(mQval, id="x")
  #     g2 <- ggplot(data=mQval) + theme_few() + geom_line(aes(x=x, y=value, color=variable)) + 
  #       scale_color_manual(values=c(colChooseLeft, colChooseRight)) + 
  #       theme(legend.title=element_blank()) + 
  #       geom_vline(xintercept = ixFailTrials, size = 0.3, colour = "grey80", linetype = "dashed") + 
  #       scale_x_continuous(expand=c(0.005,0.01)) + 
  #       ylab("V") + xlab("") 
  #     g3 <- QLearning_PlotSparklineOfRun(oneRunData) + coord_fixed(ratio = 15)
  #     grid.arrange(g3,g2,g1)
  #     #print(g4)
  #   }
  
  lhood <- QLearning_LogLikelihood(oneRunData$y, probs$pRight)
  if (displayParams) print(paste(corrAlpha, failAlpha, beta, lhood))
  return( ifelse((is.nan(lhood) | is.infinite(lhood)), 100000, lhood) )
  ## Return large log likelihood, if NaN or Inf was computed. 
  ## NaN indicates there was likely division by 0. Inf means log likelihood was 
  ## too big to be computed. 
}


############################################################################
## Generates trials using Q-Learning model based on provided 
## alpha (learning rate) and beta (inverse temperature)
## The function can also generate biased trials, such that 
## more correct responses are generated on one side using rightwardBias parameter
QLearning_GenerateTrials <- function(alpha, 				# Learning rate
                                     beta, 					# beta (inverse temperature)
                                     rightwardBias=0.5, 	# probability of rightward placed trials. When 0.5, no bias is present. 
                                     nTrials)				# Number of trials
{
  
  ## Setup variable to store values
  Qvalue <- data.frame(Left=rep(0.0, nTrials), Right=rep(0.0, nTrials))
  run <- data.frame(Stim=rep(0, nTrials), Response=rep(0, nTrials), CorrIncorr=rep(0, nTrials))
  # Generate stimuli ensuring that they have a bias to be presented on the right
  run$Stim <- sample(c(1, 2), nTrials, prob=c(1-rightwardBias, rightwardBias), replace=T)
  ## On first trial, we choose our response with 50/50 probability
  run$Response[1] <- sample(c(1,2), 1)
  run$CorrIncorr[1] <- ifelse(run$Response[1]==run$Stim[1], 1, 0)
  run$y <- ifelse(run$Response[1]==1, -1, 1)
  
  ## Array to store probabilities of responding on each trial
  probs <- data.frame(pLeft=rep(NaN, nTrials), pRight=rep(NaN, nTrials))
  ## On the first trial, the response was made with 50/50 probability
  # based on whether random guess response was correct, we will initialize 
  # Q and p values
  if (run$CorrIncorr[1]) {
    # Because initial Q value was 0, we only assign 'alpha' to Q on first trial
    if (run$Response[1]==1) Qvalue$Left[1] <- alpha else Qvalue$Right[1] = alpha # dividing by 5 to avoid giving large alpha in the beginning
  }
  
  probs$pRight[1] <- as.numeric(QLearning_Softmax(data.frame(Qvalue[1,], beta)))
  probs$pLeft[1] <- 1-probs$pRight[1]
  
  #browser()
  for (ixTrial in 2:nTrials) {
    ## Here is how the differential learning works
    ## 1. Compute Qvalue for current trial based on the outcome of previous trial
    ##	Q(t+1) = Q(t) + alpha*(R(t) - Q(t))
    ## 2. Compute probability of responding right or left using Q(t)
    ## p(R) = Softmax(Q(t+1)$Left, Q(t+1)$Right, beta)
    ## Make a new response using p(R) and fill out if that response was correct or incorrect
    
    #browser()
    # First, we make a response based on probabilities of choosing left or right 
    # that were computed on the previous trial
    run$Response[ixTrial] <- sample(c(2,1), 1, prob=c(probs$pRight[ixTrial-1], 1-probs$pRight[ixTrial-1]))
    ## Check if the response was correct 
    run$CorrIncorr[ixTrial] <- ifelse(run$Response[ixTrial]==run$Stim[ixTrial], 1, 0) 
    
    if (run$CorrIncorr[ixTrial]) {
      # Update Q value based on above response
      if (run$Response[ixTrial]==1) {
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + alpha * (1 - Qvalue$Left[ixTrial-1]) 
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
      }
      else {
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + alpha * (1 - Qvalue$Right[ixTrial-1]) 
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
    } 
    # If response was incorect, we need to decrease the Q value
    if (!run$CorrIncorr[ixTrial]) {
      # Update Q value based on above response
      if (run$Response[ixTrial]==1) {
        Qvalue$Left[ixTrial] <- Qvalue$Left[ixTrial-1] + alpha * (0 - Qvalue$Left[ixTrial-1])
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1]
      }
      
      else { 
        Qvalue$Right[ixTrial] <- Qvalue$Right[ixTrial-1] + alpha * (0 - Qvalue$Right[ixTrial-1]) 
        Qvalue$Left[ixTrial] <-  Qvalue$Left[ixTrial-1]
      }
    } 
    
    # Based on updated Q value, let's get probabilities of choosing left or right 
    # using Softmax ruls
    probs$pRight[ixTrial] <- as.numeric(QLearning_Softmax(data.frame(Qvalue[ixTrial,], beta)))
    probs$pLeft[ixTrial] <- 1-probs$pRight[ixTrial]
  }
  ## Copy of responses coded as -1 (Left) or 1 (Right)
  run$y <- ifelse(run$Response==1, -1, 1)
  
  # Bind data
  run <- cbind(run, Qvalue, probs)
  colnames(run)[names(run) %in% c('Left', 'Right')] <- c('QLeft', 'QRight')
  return(run)	
}


#############  How does the beta value change the slope of the psychometric function ##############
## Compute Softmax for a set of Left/Right choice values and betas (temperature) 
## Returns a plot of Softmax functions for requsted beta values. 
QLearning_ImpactOfBeta <- function(betasToTest=c(0.1, 0.5, 1, 2, 3, 5, 10), ## Beta values for which Softmax will be computed 
                                   QRightMax=1,	# Max value for Right choice 
                                   QLeftMax=1 	# Max value for Right choice
) {
  
  valueIncrement <- 0.01 ## How fine should values be computed. If too slow, make the step a bit bigger. 
  ## Matrix of choosing left and right values by beta 
  valuesByBeta <-expand.grid(QRight=seq(0, QRightMax, valueIncrement), QLeft=seq(0, QLeftMax, valueIncrement), beta=betasToTest) 
  ## Add a column with probabilities of choosing rightwar for all values and betas. 
  ## The pRight is computed using the Softmax rule (general form of logistic regression)
  datToPlot <- ddply(valuesByBeta, .(beta), transform, pRight=QLearning_Softmax(data.frame(QLeft, QRight, beta)))
  ## Get relative value as diff btw right and left values to plot as x axis. 
  datToPlot <- within(datToPlot, {QRelative=QRight-QLeft})
  ## Same as above line
  ##datToPlot$QRelative <- dd$QRight - dd$QLeft
  
  #browser()
  g <- ggplot(data=datToPlot) + 
    theme_few() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) + 
    theme(axis.line = element_line(colour = "#a9a9a9", size = 0.3)) + 
    theme(axis.ticks.x = element_line(colour = "#a9a9a9", size = 0.3)) + 
    geom_line(aes(x=QRelative, y=pRight, color=as.factor(beta)), size=1) + 
    scale_color_brewer(palette="Paired") + 
    xlab("Relative value (Right-Left)") + 
    ylab("Proportion rightward choice") + 
    guides(color=guide_legend(title="beta"))  
  return(g)
}

############### Convert action values into probabilitiies ###################
## Q1 and Q2 are values related to choosing 1 or 2
## Beta is a parameter which is called inverse temperature (1/T). 
## Higher values of beta make the sigmoid to be steeper (more sensitivity)
## Smaller values of beta make choices very random
QLearning_Softmax <- function(vals      ## Vector containing value for option1, value for option2 and beta
)	{
  ## unpack input params
  Q1 <- vals[[1]]		## Value associated with choosing 1 (Left)
  Q2 <- vals[[2]]		## Value associated with choosing 2 (Right)
  beta <- unique(vals[[3]])		## Inverse temperature (beta)
  
  ## Compute and return probability of choosing 1 using Softmax rule 
  pRight <- exp(Q2*beta) / (exp(Q1*beta) + exp(Q2*beta))
  return(as.data.frame(pRight))
}


##########################################################
## Intuitive plot of stimulus presentation side over time.
## Helps to discover unusual patterns such as long streaks 
## when inducing bias in success/stay condition. 
QLearning_PlotSparklineOfRun <- function(runData,   	 	# Input data in the format of "glmData"
                                         plotResults=F, 	# If False, returns list of ggplot objects. If True, plots those objects, each subject in separate window
                                         plotInColor=T		# Success and failure are marked by two different colors. Otherwise, black and white
){
  ## Add column with trial numbers
  #runData <- droplevels(ddply(runData, mutate, TrialID=1:length(Response)))
  runData <- cbind(runData, TrialID=1:length(runData$Response))
  
  ## Generate ggplot graphs as binary sparklines 
  g <- ggplot(data=runData, aes(x=TrialID, y=Response)) + theme_few() 
  g <- g + geom_linerange(subset=.(Response==1), aes(ymin=1.5, ymax=Response, color="Left"), size=0.3, ) + 
    geom_linerange(subset=.(Response==2), aes(ymin=1.5, ymax=Response, color="Right"), size=0.3) +
    geom_linerange(subset=.(CorrIncorr==0), aes(ymin=1.5, ymax=Response, color="Fail"), size=0.3) + 
    scale_color_manual(values=c("black", colChooseLeft, colChooseRight)) + 
    scale_x_continuous(expand=c(0.01,0.01)) + 
    scale_y_continuous(breaks=c(1,2), labels=c("L","R")) + 
    theme(axis.ticks.x=element_blank()) + 
    theme(legend.title=element_blank()) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank()) + 
    ylab("Response") +
    xlab("")
  
  if (plotResults) print(g)
  return(g)
  
} 

###########################################
## Returns log likelihood 
QLearning_LogLikelihood <- function(y, 		# Participant responses 
                                    pRight)	# Probability of responding right
{
  # Before computing log-lhood, some adjustments need to be made
  # First, p that correspond to responses, are one trial before the response
  # So we can shift pRight down by one trial and leave out of computation the 1st trial 
  # because we don't know the true probability of responding to the first 
  # trial
  
  # To use Lag, load library Hmisc
  #browser()
  
  pRight <- Lag(pRight, 1)
  pRight <- pRight[2:length(pRight)]
  y <- y[2:length(y)]
  -sum(log(ifelse(y==1, pRight, 1 - pRight)))
}

