load('20150211_cond1_learningCoefs_model1.RData')
learningCoefsCond1 <- learningCoefs
load('20150214_cond2and3_learningCoefs_model1.RData')
learningCoefsCond2and3 <- learningCoefs
load('20150214_cond7and8_learningCoefs_model1.RData')
learningCoefsCond7and8 <- learningCoefs
learnCoefsAll <- rbind(learningCoefsCond1, learningCoefsCond2and3, learningCoefsCond7and8)



meanLrnCoefs <- ddply(learnCoefsAll, .(SubjectID, Condition), summarise, alpha=mean(alpha), beta=mean(beta))
meanLrnCoefs <- subset(meanLrnCoefs, SubjectID %in% c('s001', 's005', 's007', 's009', 's010', 's012', '014'))

ggplot(data=meanLrnCoefs, aes(x=alpha, y=beta, color=SubjectID)) + 
  theme_few() + 
  scale_color_brewer(palette='Set1') + 
  facet_grid(~Condition, labeller=ConditionLabels) + 
  geom_point(size=2)

meanByCond <- ddply(learnCoefsAll, .(Condition), summarise, alpha=mean(alpha), beta=mean(beta))
ggplot(data=meanByCond, aes(x=alpha, y=beta, color=as.factor(Condition))) + 
  theme_few() + 
  scale_color_brewer(palette='Set1') + 
  #facet_grid(~Condition, labeller=ConditionLabels) + 
  geom_point(size=3)
