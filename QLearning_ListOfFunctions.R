# List of all function I wrote for RL. 
# More functions were added later and those need to be added

# ----------- File: GenerateTrialsUsingQLearning.R -----------

QLearning_ValidateAlphaBeta
## Generate trials with a bunch of randomly selected alpha and beta values
## Then fit the model to estimate alpha and beta values from those generated 
## trials. Then plot true and estimated alpha/beta values as a scatterplot 
## and connect relevant pairs together using geom_segment to see how far or close 
## are true and estimated values based on values of alpha and beta

Qlearning_AlphaBetaRunningEstimate
############################################################################
## Estimates learning rate on each trial by running 
## optimisation on a subset of trials 

QLearning_GenerateTrials
############################################################################
## Generates random trials with a given rightward or leftward bias 
## Then run the Q-Learning model on those trials and generates responses to those 
## trials using provided alpha and beta parameters. 

# ------------- File: QLearningLikelihood_NEW.R --------------

QLearning_Model1
## Need to figure out what this function does. This is an important function

QLearning_LogLikelihood
## Returns log likelihood of choosing right given subject responses

PlotSparklineOfRun
## Plot of subject responses on left or right and 
## whether those responses were rewarded or not 
