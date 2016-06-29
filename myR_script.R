# Re-run everything!!!
# first run this: rm(list = ls())
# this deletes all previous objects

## Pakages we need for our project
library(mvabund)
library(vegan)
library(coda)
library(rjags)
library(boral)
#library(mgcv)
library(effects)
##################################
##################################
##### 1- Isolation success analysis:
## Isolation success data input

# Get the input data from here first:
# https://figshare.com/s/fbf753647e6e86467707

MyAbund = read.csv(file="morphotype_matrix_incubator.csv", 
                   header = T, row.names = 1)
MetaData = read.csv(file="metadata.csv", header = T, row.names = 1)
# you write out tables with the write.csv command. Check the helpfile with ?write.csv
# are the rownames matching?
# rownames(MyAbund) == rownames(MetaData)
# Isolation success
IsolSucc = apply(MyAbund,1, sum)


# Combine with the rest of the metadata
MetaData = cbind(MetaData, success = IsolSucc)

# Re-arrange the order of the factor predictors. branches and dust are compared to leaves.
MetaData$source <- factor(MetaData$source, levels = c("Leaf", "Branch", "Dust"))
MetaData$time <- factor(MetaData$time, levels = c("1", "2", "3"))

# Sample names:
# B1T1L1:
# B: locality of collection (B = Bisoton)
# 1: sampling time 1
# T1: the code of the tree at that site 
# L1: leaf 1 of that tree (can be B=branch, D=dust)

# Some data overview
# table(MetaData$locality, MetaData$time)
# table(MetaData$locality, MetaData$source)
# table(MetaData$time, MetaData$source)
# table(MetaData$tree, MetaData$source)
# table(MetaData$tree, MetaData$time)
# table(MetaData$tree, MetaData$source)
# table(MetaData$success, MetaData$source)

# # Continuous variables
# par(mfrow = c(2,1), mar=c(2,2,1,1))
# hist(MetaData$success, nclass=20)
# summary(MetaData$success)
# hist(MetaData$dust, xlab="Dust (mg/cm2)")

# Some simple plots
# Time?
# boxplot(MetaData$success ~ MetaData$time)
# # Locality and dust quantity
# boxplot(MetaData$dust ~ MetaData$locality)
# # Dust and sampling time
# boxplot(MetaData$dust ~ MetaData$time)
# # Dust quantity at localities at different times
# boxplot(MetaData$dust[MetaData$time == 1] ~ 
#           MetaData$locality[MetaData$time == 1])
# boxplot(MetaData$dust[MetaData$time == 2] ~ 
#           MetaData$locality[MetaData$time == 2])
# boxplot(MetaData$dust[MetaData$time == 3] ~ 
#           MetaData$locality[MetaData$time == 3])

# Models of culturing success
# success.m1 = glm(success ~ locality + time + source + dust, data=MetaData, 
#                  family = poisson(link = "log"))
success.m = glm(success ~ locality + time + source*dust, data=MetaData, 
               family = poisson(link = "log"))

# Compare model fit: succes.m2 seems better
# AIC(success.m1)
# AIC(success.m2)
# 
# par(mfrow = c(2,2))
# plot(success.m1)
# plot(success.m2)

### Culturing success model results
success.sammary=summary(success.m)
success.anova=anova(success.m, test="LRT")

# Visualize the culturing success model results
# Q1: Dust deposition and isolation success
# Use the model predictions here
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(MetaData$dust, fitted(success.m2), xlab=c("Dust (mg/cm2)"), 
     ylab=c("Isolation success"))
# Try to add a trendline + confidence interval for the trendline here
abline(success.m2$coefficients["(Intercept)"], success.m2$coefficients["dust"])

# Q2: Locality and isolation success
boxplot(fitted(success.m2) ~ MetaData$locality, 
        ylab="Isolation success")

# Q4: Source of isolation
boxplot(fitted(success.m2) ~ MetaData$source,
        ylab="Source of isolation")

# Dust - source interactions
par(mar = c(12,4,3,1))
boxplot(MetaData$success ~ MetaData$locality*MetaData$source, 
        las=2)

## plot the interaction
success.effect = effect("source:dust",success.m2 ,  multiline=TRUE)
plot(success.effect)
# Effect summaries
succ.eff.sum = summary(success.effect)

# Coefficients and 95% confidence intervals
Effs = cbind(Coefficients = coefficients(success.m2), confint(success.m2))

###########################
###########################
#### 2- Diversity analysis:
### Diversity Data input

# factors
#my.metadata$time=factor(my.metadata$time) 
locality=factor(MetaData$locality)

# I changed source to Source, as source() seems to be an R function
Source = factor(MetaData$source, levels = c("Leaf", "Branch", "Dust"))
MetaData$time<-factor(my.metadata$time, levels = c("1","2","3"))
#dust=(my.metadata$dust)

#### Richness (Species number) model
# Richness is the nmber of culture observations in the samples. 
# Some samples had more, than one observed species.
hist(IsolSucc)

# For the evaluation of richness we need to remove the samples with zero observations.
NotZero = IsolSucc > 0 #filter for zero-observation samples

# Keep only the samples with at least one observed species
AbundNotZero=MyAbund[NotZero,]

# Richness in the samples
Richness = specnumber(AbundNotZero)
hist(Richness)
hist(log(Richness))

# Remove the samples with zero observation from the metadata
MetaRich = MetaData[NotZero,]

## fitting the models for Richness
# Richness.m1=lm(Richness ~ locality+ time + source*dust, data = MetaRich)
# Richness.m2=glm(Richness ~ locality+ time + source+dust, data = MetaRich, 
#                 family=poisson(link = "log"))
Richness.m= glm(Richness ~ locality + time + source*dust, data = MetaRich, 
                family=poisson(link = "log"))
 
# model diagnostic plots
# Change plotting parameters, so all diagnostic plots are on the same sheet
# par(mfrow=c(2,2))
# plot(Richness.m1)
# plot(Richness.m2)
# plot(Richness.m3)
# summary(Richness.m1)
# summary(Richness.m2)
# summary(Richness.m3)
# Look for the AICs: Akaike Information Criteria and model selection
# AIC(Richness.m1)
# AIC(Richness.m2)
# AIC(Richness.m3)
# The AIC says that the simple linear model (m1) is the best fit 
# (although still crappy according to the diagnostic plot)

### Richness Model Evaluations
Richness.anova=anova(Richness.m, test="LRT")
Richness.summary=summary(Richness.m)

# Visualize Richness changes
# What does it mean? With predicted values
par(mfrow=c(1,1))
boxplot(fitted(Richness.m) ~ MetaRich$source, xlab="Source", ylab="Richness")
# Richness in branch and dust similar, but richness stat. sign. lower in leaf

# Interaction plots:
Richness.effect = effect("source:dust",Richness.m,multiline=TRUE,  ylim=c(-10,10))
plot(Richness.effect, ylab = "Richness")
Rich.eff.sum = summary(Richness.effect)

# par(mfrow=c(1,3), mar = c(5,3,2,1))
# for (i in levels(MetaOne$source)) {
#   plot(c(min(MetaOne$dust), max(MetaOne$dust)),
#        c(min(Rich.eff.sum$lower[i,]), max(Rich.eff.sum$upper[i,])), 
#        type="n", xlab = paste(i), ylab = "Richness")
#   lines(c(min(MetaOne$dust), max(MetaOne$dust)),
#         c(min(Rich.eff.sum$effect[i,]), max(Rich.eff.sum$effect[i,])))
#   lines(c(min(MetaOne$dust), max(MetaOne$dust)),
#         c(min(Rich.eff.sum$lower[i,]), max(Rich.eff.sum$lower[i,])), 
#         lty="dashed")
#   lines(c(min(MetaOne$dust), max(MetaOne$dust)),
#         c(min(Rich.eff.sum$upper[i,]), max(Rich.eff.sum$upper[i,])),
#         lty="dashed")}

### Shannon and Simpson models.
# Keep only samples with at least two OTUs
RichNotOne = Richness > 1
AbundNotOne=AbundNotZero[RichNotOne,]

# This keeps observations with at least two OTUs
MetaNotOne = MetaRich[RichNotOne,] 

# Calculate diversity indices
shannon = diversity(AbundNotOne,index = "shannon")
simpson = diversity(AbundNotOne,index = "simpson")
par(mfrow = c(2,2))
hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

## Fishers alpha log-series
# Check what Fisher's alpha means. 
# You get very high values with only a few species observations.
#Fisheralpha= fisher.alpha(fungl.abun.not.one)
#hist(Fisheralpha)
## Evenness
#RichnessMinTwo = Richness[Richness > 1]
#Evenness= shannon/log(RichnessMinTwo)
#hist(Evenness)

## Shannon models: let's fit all of them first
shannon.m=lm(shannon ~ locality+time+source*dust, data=MetaNotOne)

# # plot all diagnostics
 # par(mfrow=c(2,2))
 # plot(shannon.m)
# plot(shannon.m2)
# plot(shannon.m)
# 
# # AICs - the third is the best fit + it has the best diagnostic plots
# AIC(shannon.m1)
# AIC(shannon.m2)
# AIC(shannon.m4)

# These are probably not necessary, as these models are worse fits
# summary(shannon.m1) 
# anova(shannon.m1, test="Chisq")
# summary(shannon.m2) 
# anova(shannon.m2, test = "Chisq")

# Best Shannon model statistics
shannon.summary=summary(shannon.m) 
shannon.anova=anova(shannon.m, test = "LRT")

# Interaction plots:
shannon.effect = effect("source:dust",shannon.m,multiline=TRUE,  ylim=c(-10,10))

# Put here other boxplots if you wnat to show other effects, e.g. time or locality

# Plot the interaction effect
# Either plot them on three different plots, 
# or scale the ylim values for each of them.
plot(shannon.effect, ylim = c(-1,2))

# Effect summary
shan.eff.sum = summary(shannon.effect)

## Simpson models:
simpson.m = lm(simpson ~ locality + time + source*dust, data=MetaNotOne)
# simpson.m2=glm(simpson~time+locality+source+dust, data=MetaOne,family=Gamma(link = "log"))
# simpson.m3=glm(simpson~locality+time+source*dust, data = MetaRich,family = Gamma(link="log"))
# 
# # Diagnostics
# par(mfrow = c(2,2))
# plot(simpson.m)
# plot(simpson.m2)
# plot(simpson.m3)
# 
# # AICs - again everything says the m3 is the best.
# AIC(simpson.m1)
# AIC(simpson.m2)
# AIC(simpson.m3)
# 
# # Model statistics
# # summary(simpson.m1)
# # anova(simpson.m1, test= "Chisq")
# # summary(simpson.m2)
# # anova(simpson.m2, test= "Chisq")

simpson.summary=summary(simpson.m)
simpson.anova=anova(simpson.m, test= "Chisq")

# Other boxplot if needed here:

# Interaction plots:
simpson.effect = effect("source:dust",simpson.m,multiline=TRUE,  ylim=c(-10,10))

# plot interaction effects
plot(simpson.effect, ylim = c(0,1))

# Effect summary
Simp.eff.sum = summary(simpson.effect)


# # Actually, let's plot the model predictions.
# #plot(MetaOne$dust, fitted(shannon.m3), xlab="Dust concentration")
# par(mfrow =c(1,3))
# boxplot(fitted(Richness.m1) ~ MetaRich$source, xlab="Source of isolation", ylab="Richness")
# ## why the ylab doesn't work here?
# boxplot(fitted(shannon.m3) ~ MetaOne$source, xlab="Shannon")
# boxplot(fitted(simpson.m3) ~ MetaOne$source, xlab="Simpson")

#### Fishers alpha log-series: do we really need this?
## Family for this model?
#hist(Fisheralpha)
#Fisheralpha.m=glm(Fisheralpha~time+source*dust+locality, data= MetaOne)

#par(mfrow=c(2,2))
#plot(Fisheralpha.m)

#Fisher.summary=summary(Fisheralpha.m)
#fisher.anova=anova(Fisheralpha.m, test="Chisq")

###Final models:
### I changed the order of the factors in the model and it did not really change the results and these three
### models are final
# Richness.m1=lm(RichnessNotZero~time+locality+source*dust, data = MetaRich)
# shannon.m3=glm(shannon~locality+time+source*dust, data=MetaObs, family=Gamma(link="log"))
# simpson.m3=glm(simpson~locality+time+source*dust, data = MetaObs,family = Gamma(link="log"))

### 3- Community composition models
# Put the model first

#### Model-based ordination

# Remove species present in only one sample
CountInSample = apply(fungl.abun.not.zero,2,sum)

# I needed to remove the species that were seen in a few samples.
fung.abun.reduced = fungl.abun.not.zero[,CountInSample > 3]
fung.abun.reduced = fung.abun.reduced[apply(fung.abun.reduced,1,sum) > 0,]

MetaOrd = MetaRich[rownames(MetaRich) %in% rownames(fung.abun.reduced),]

# not working
ModelOrd <- boral(fung.abun.reduced, family = "negative.binomial", num.lv = 2, 
                  n.burnin = 10, n.iteration = 100, n.thin = 1)

## model diagnostics
par(mfrow = c(2,2))
plot(ModelOrd, ask = FALSE )

# Ordination plot. Please make it similar as the NMDS.
dev.off()
par(mar=c(4,4,1,1))
ordibora= ordiplot(ModelOrd$lv.median, choices = c(1,2), type = "none", cex =0.5 
         ,display = "sites", xlim = c(-0.3,0.3))
points(ordibora,"sites", pch=20 ,col=as.numeric(MetaOrd$source))
# ordispider(ordibora,MetaOrd$source, col= "gray" )
mylegend = legend(0.7, 0.9, c("Leaf","Branch","Dust"), 
                  fill=c(3,1,2), border="white", bty="n")
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col=3,
                          alpha=200,kind="se",conf=0.95, 
                          show.groups=(c("Leaf"))))
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col=1,
                          alpha=150,kind="se",conf=0.95,
                          show.groups=(c("Branch")))) 
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col=2,
                          alpha=200,kind="se",conf=0.95,
                          show.groups=(c("Dust"))))
# ordisurf of dust: show dust + source interaction on the whole community

### Why this plot is different from NMDS?
# Because the NMDS plot likely has overdispersion effects.
# dev.off()
# par(mar=c(4,4,1,1))
# plot(ModelOrd$lv.median, col=as.numeric(MetaOrd$source), 
#             pch=19, main="Latent variable model", las=1.5)

# The NMDS will not be necessary anymore. If you run the two lines below, 
# you will see how locality and dispersal are mixed up because of the overdispersion.
# # This is very strong for the green samples (leaves?).
# NMDS1<-metaMDS(fung.abun.reduced)
# plot(NMDS1$points, ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2", 
#      col=as.numeric(MetaOrd$source), pch=19)
# 
# 



# 
# png(file="mean-variance-plot.jpeg", units="mm", height=90, width=90, 
#     pointsize=10, bg="white", res=1200)
# par(mar = c(4,4,1,1))
# plot(0.1,0.1, type="n", xlim=c(0.1,100), 
#      ylim=c(0.1, max(apply(Corfun0,2,var))),
#      xlab="Mean (log scale)", ylab="Variance (log scale)", log="xy", xaxt="n", yaxt="n")
# for (i in 1:length(colnames(Corfun0))) {
#   points(mean(Corfun0[,i]), var(Corfun0[,i]), pch=20, cex=0.7)}
# axis(1, at=c(0.1, 1,5,10), labels=c(0,1,5,10))
# ticks.2 = seq(1,9, by=2)
# labels.2 <- sapply(ticks.2, function(i) as.expression(bquote(10^ .(i))))
# axis(2, at=c(0.1, 10, 1000, 100000, 10000000), 
#      labels=labels.2)


funcor.mvabund = mvabund(Corfun0)
plot(funcor.mvabund)

## Model selection
# ### if I remember correctly you said that I don't have to do this part and I should just go with the most complicated 
# ## model from diversity models
# 
# # funcor.m1 = manyglm(funcor.mvabund ~ locality+time+source+dust, data= metacor0,
# #                  family="negative.binomial", show.residuals=T)
# # summary(funcor.m1)
# 
# funcor.m2 = manyglm(funcor.mvabund ~ locality+time+source*dust, data= metacor0,
#                     family="negative.binomial", show.residuals=T)
# 
# # funcor.m3 = manyglm(funcor.mvabund ~ locality*time+source*dust , data= metacor0, 
# #                     family="negative.binomial", show.residuals=T) 
# # summary(funcor.m3)
# # anova(funcor.m2, funcor.m3, nBoot=50)
# # 
# # funcor.m4= manyglm(funcor.mvabund ~ locality*dust+time+source*dust , data= metacor0, 
# #                    family="negative.binomial", show.residuals=T)
# # summary(funcor.m4)
# # anova(funcor.m3,funcor.m4, nBoot = 50)
# 
# 
# ### model 3 is the better option comparing to other models. but I can't find any explanation for the 
# ## interaction effect of time and locality! aparently somthing happend in those times in those localities but 
# ### we don't know what!!
# ### so I want to continue with model 2
# ## factors:
# metacor0$time=factor(metacor0$time, levels = c("1","2","3"))
# metacor0$source = factor(metacor0$source, levels = c("Leaf", "Branch", "Dust"))
# metacor0$locality=factor(metacor0$locality)

## Selected model: funcor.m2
funcor.m2 = manyglm(funcor.mvabund ~ locality+time+source*dust, data= metacor0,
                    family="negative.binomial", show.residuals=T)

plot.manyglm(funcor.m2)

### ask about the error

Model.m2.summary= summary.manyglm(funcor.m2, nBoot=300, test="LR",p.uni="adjusted", 
                                  resamp="montecarlo")

## Analysis of variance explained by the predictors
funcor.anova.m2 = anova.manyglm(funcor.m2, nBoot=300, test="LR", p.uni="adjusted", 
                     resamp="montecarlo")

# ## OTUs significantly affected by the source?? (at p<=0.001???)
m2.p.anova <- as.data.frame(funcor.anova.m2$uni.p)

fun.m2.source = colnames(m2.p.anova)[m2.p.anova["source",]<=0.1]
colnames(m2.p.anova)[m2.p.anova["source",]<=0.001]
colnames(m2.p.anova)[m2.p.anova["locality",]<=0.1]
colnames(m2.p.anova)[m2.p.anova["time",]<=0.1]
# try also putting the interaction here


## Visualization of source*dus interactions 
## Coefficients
funcor.m2.coef = as.data.frame(funcor.m2$coefficients)

## How do I know if i have an outlier or not?


## mean-centering the contrasts
fun.coef.mean.contrast = funcor.m2.coef - apply(funcor.m2.coef,2,mean)


### plot effects

?predict.glm
predict(funcor.m2, type="response")

plot(metacor0$dust, predict(funcor.m2, type="response")[,"Gnomoniaceae_sp_66"])

plot(metacor0$dust[metacor0$source == "Leaf"], predict(funcor.m2, type="response")
     [metacor0$source == "Leaf","Gnomoniaceae_sp_66"])
plot(metacor0$dust[metacor0$source == "Branch"], predict(funcor.m2, type="response")
     [metacor0$source == "Branch","Gnomoniaceae_sp_66"])
plot(metacor0$dust[metacor0$source == "Dust"], predict(funcor.m2, type="response")
     [metacor0$source == "Dust","Gnomoniaceae_sp_66"])
plot(metacor0$dust[metacor0$source == "Leaf"],
     Corfun0[metacor0$source == "Leaf","Gnomoniaceae_sp_66"])
plot(metacor0$dust[metacor0$source == "Leaf"], predict(funcor.m2, type="response")
     [metacor0$source == "Leaf","Gnomoniaceae_sp_66"])


###source*dust interactions plot
### Do a simple glm for the two species affected by the source* dust interaction
library(MASS)
Microsphaeriopsis.model= glm.nb(Corfun0$Microsphaeriopsis_olivacea ~ locality+time+source*dust, data= metacor0)
Micros.anova= anova(Microsphaeriopsis.model)                             

Micros.effec= effect("source:dust",Microsphaeriopsis.model ,multiline=TRUE,  ylim=c(-10,10))
Micros.effec.sum = summary(Micros.effec)
par(mfrow=c(1,3), mar = c(5,3,2,1))
for (i in levels(metacor0$source)) {
  plot(c(min(metacor0$dust), max(metacor0$dust)),
       c(min(Micros.effec.sum$lower[i,]), max(Micros.effec.sum$upper[i,])), 
       type="n", xlab = paste(i), ylab = "")
  lines(c(min(metacor0$dust), max(metacor0$dust)),
        c(min(Micros.effec.sum$effect[i,]), max(Micros.effec.sum$effect[i,])))
  lines(c(min(metacor0$dust), max(metacor0$dust)),
        c(min(Rich.eff.sum$lower[i,]), max(Rich.eff.sum$lower[i,])), 
        lty="dashed")
  lines(c(min(metacor0$dust), max(metacor0$dust)),
        c(min(Micros.effec.sum$upper[i,]), max(Micros.effec.sum$upper[i,])),
        lty="dashed")}












Aureobasidium.model= glm.nb (Corfun0$Aureobasidium_sp_A30 ~ locality+time+source*dust, data= metacor0)

plot(effect("source:dust",Aureobasidium.model ,multiline=TRUE,confidence.level = 0.95))

anova(Aureobasidium.model)


### similarity analysis using anosim and adonis functions (corfun0=fungal abundance matrix for core OTUs and metacor0 
## is the metadata)
braysimi= anosim(Corfun0,metacor0$source,permutations = 999,distance = "bray")

jaccardsimi= anosim(Corfun0,metacor0$source, permutations = 999, distance = "jaccard")


corfunAdon= adonis(Corfun0~locality+time+source*dust,data= metacor0,permutations = 999, method = "bray")

### permanova for all of the species

permanova.total=adonis(fungl.abun2~locality+time+source*dust,data=MetaObs,permutations = 999,method = "bray")



