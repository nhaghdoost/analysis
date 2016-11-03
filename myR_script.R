
## Effects of dust deposition on endophytic fungal diversity and community composition of Persian Oak##
## Niloufar Hagh Doust, Moslem Akbarinia, Naser Safaie, Hamed Yousefzadeh, Miklos Balint

## Pakages we need for our project
library(mvabund)
library(vegan)
library(coda)
library(rjags)
library(boral)
#library(mgcv)
library(effects)
library(MASS)
library(reshape)
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
MetaData$tree <- factor(MetaData$tree)
MetaData$locality <- factor(MetaData$locality)
# Sample names:
# B1T1L1:
# B: locality of collection (B = Bisoton)
# 1: sampling time 1
# T1: the code of the tree at that site 
# L1: leaf 1 of that tree (can be B=branch, D=dust)

#### Isolation success model
success.m = glm (success ~ locality + time + source * dust, data=MetaData, 
               family = poisson(link = "log"))

### Culturing success model results
success.sammary = summary(success.m)
success.anova = anova(success.m, test="LRT")

# Visualize the culturing success model results
# Q1: Dust deposition and isolation success
# Use the model predictions here
# par(mfrow=c(1,1), mar=c(5,5,2,2))
# plot(MetaData$dust, fitted(success.m), xlab=c("Dust (mg/cm2)"), 
#      ylab=c("Isolation success"))
# lines(MetaData$dust, predict.lm(success.m), col= " blue")
# # Try to add a trendline + confidence interval for the trendline here. need mor etime for this!
# abline(coef = coef(success.m))
# abline(coef= success.m$coefficients["(Intercept)"], success.m
#       $coefficients["dust"], col= "blue")
# abline (col="blue", h= coef(success.m), )
# plot again
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(MetaData$success ~ MetaData$dust, xlab=c("Dust Deposition (mg/cm2)"), 
      ylab=c("Isolation success"))
lines(MetaData$dust, predict.lm(success.m), col= " blue")

# Q2: Locality and isolation success
boxplot(fitted(success.m) ~ MetaData$locality, 
        ylab="Isolation success")

# Q4: Source of isolation
boxplot(fitted(success.m) ~ MetaData$source,
        ylab="Source of isolation")

# Dust - source interactions
par(mar = c(12,4,3,1))
boxplot(MetaData$success ~ MetaData$locality*MetaData$source, 
        las=2)

## plot the interaction
success.effect = effect("source:dust",success.m , multiline=TRUE)

# suc.eff <- plot(success.effect, xlab = "Dust Deposition (mg/cm2)", ylab = "Isolation success", 
#     main = NULL, ylim = c(-2.5,1.5))
# 
# suc.eff # plot original object
# suc.eff$x.scales$cex <- c(1.3, 1.3)
# suc.eff$xlab$cex <- 1.3
# suc.eff$y.scales$cex <- c(1.3, 1.3)
# suc.eff$ylab$cex <- 1.3
# suc.eff # plot edited object

# Effect summaries
succ.eff.sum = summary(success.effect)

# Coefficients and 95% confidence intervals
Effs = cbind(Coefficients = coefficients(success.m), confint(success.m))

#######################################
#######################################
#### 2- Diversity analysis:
### Diversity Data input

# factors
#my.metadata$time=factor(my.metadata$time) 
locality=factor(MetaData$locality)

# I changed source to Source, as source() seems to be an R function
Source = factor(MetaData$source, levels = c("Leaf", "Branch", "Dust"))
MetaData$time <- factor(MetaData$time, levels = c("1","2","3"))
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

## fitting the model for Richness
# Richness.m1=lm(Richness ~ locality+ time + source*dust, data = MetaRich)
# Richness.m2=glm(Richness ~ locality+ time + source+dust, data = MetaRich, 
#                 family=poisson(link = "log"))
Richness.m= glm(Richness ~ locality + time + source*dust, data = MetaRich, 
                family=poisson(link = "log"))

### Richness Model Evaluations
Richness.anova = anova (Richness.m, test="LRT")
Richness.summary = summary (Richness.m)

# Visualize Richness changes

par(mfrow=c(1,1))
boxplot(fitted(Richness.m) ~ MetaRich$source, xlab="Source", ylab="Richness")
# Richness in branch and dust similar, but richness stat. sign. lower in leaf

# Interaction plots:
Richness.effect = effect("source:dust",Richness.m, multiline=TRUE, confidence.level = 0.95)
# Rich.eff<- plot(Richness.effect, ylab = "Richness", xlab = "Dust Deposition (mg/cm2)", 
#                 main= NULL, ylim = c(-1.5,2))
# Rich.eff # plot original object
# Rich.eff$x.scales$cex <- c(1.3, 1.3)
# Rich.eff$xlab$cex <- 1.3
# Rich.eff$y.scales$cex <- c(1.3, 1.3)
# Rich.eff$ylab$cex <- 1.3
# Rich.eff # plot edited object
# Rich.eff.sum = summary(Richness.effect)


### Shannon and Simpson models.
# Keep only samples with at least two OTUs
RichNotOne = Richness > 1
AbundNotOne=AbundNotZero[RichNotOne,]

# This keeps observations with at least two OTUs
MetaNotOne = MetaRich[RichNotOne,] 

# Calculate diversity indices
shannon = diversity(AbundNotOne,index = "shannon")
simpson = diversity(AbundNotOne,index = "simpson")

hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))
dev.off()
## Fishers alpha log-series
# Check what Fisher's alpha means. 
# You get very high values with only a few species observations.
#Fisheralpha= fisher.alpha(fungl.abun.not.one)
#hist(Fisheralpha)
## Evenness
#RichnessMinTwo = Richness[Richness > 1]
#Evenness= shannon/log(RichnessMinTwo)
#hist(Evenness)

## Shannon model
shannon.m = lm (shannon ~ locality + time + source*dust, data = MetaNotOne)

## Shannon model statistics
shannon.summary = summary (shannon.m) 
shannon.anova = anova (shannon.m, test = "LRT")

# Interaction plots:
shannon.effect = effect("source:dust", shannon.m, multiline=TRUE)
shan.eff.sum = summary (shannon.effect)

# shan.eff<-plot(shannon.effect,ylab = "Shannon Diversity", xlab = "Dust Deposition (mg/cm2)",
#      main= NULL, ylim = c(-0.1,1.5) )
# shan.eff # plot original object
# shan.eff$x.scales$cex <- c(1, 1)
# shan.eff$xlab$cex <- 1.3
# shan.eff$y.scales$cex <- c(1.3, 1.3)
# shan.eff$ylab$cex <- 1.3
# shan.eff # plot edited object

## Simpson model:
simpson.m = lm (simpson ~ locality + time + source*dust, data=MetaNotOne)
simpson.summary = summary (simpson.m)
simpson.anova = anova (simpson.m, test= "Chisq")
# Interaction plots:
simpson.effect = effect("source:dust", simpson.m, multiline=TRUE)
# 
# 
# # plot interaction effects
# simp.eff<- plot(simpson.effect,ylab = "Simpson Diversity", xlab = "Dust Deposition (mg/cm2)", 
#      main= NULL, ylim = c(-0.1,1))
# simp.eff # plot original object
# simp.eff$x.scales$cex <- c(1, 1)
# simp.eff$xlab$cex <- 1.3
# simp.eff$y.scales$cex <- c(1.3, 1.3)
# simp.eff$ylab$cex <- 1.3
# simp.eff # plot edited object
# # Effect summary
# Simp.eff.sum = summary (simpson.effect)

#############################
#############################
### MERG ALL THE INTERACTION PLOTS 
success.effect = effect("source:dust",success.m , multiline=TRUE)
Richness.effect = effect("source:dust",Richness.m, multiline=TRUE, confidence.level = 0.95)
shannon.effect = effect("source:dust", shannon.m, multiline=TRUE)
simpson.effect = effect("source:dust", simpson.m, multiline=TRUE)
### The plot
dev.off()
par(mfrow=c(1,3))

#Leaf
plot(c(0.01,0.06), c(0.3,1.3), type = "n", xlab="Dust deposition (mg/cm2)",
     main = "Leaf", ylab = "",cex.axis = 1.5,
     cex.lab = 1.5, cex.main= 1.5)
points(c(0.01,0.06), c(0.3756226,0.2955457), type = "l",lwd=2, col= "black")#isolation success
points(c(0.01,0.06), c(1.245071,0.7484556), type = "l",lwd=2, col="red")#Richness
points(c(0.01,0.05), c(0.7573420,0.4898882), type = "l",lwd=2, col="blue")#Shannon
points(c(0.01,0.05), c(0.5374573,0.5828385), type = "l",lwd=2, col="green")#Simpson

#Branch
plot(c(0.01,0.06), c(0.5,2.1), type = "n", xlab="Dust deposition (mg/cm2)",
     main = "Branch", ylab = "",cex.axis = 1.5,
     cex.lab = 1.5,cex.main= 1.5)
points(c(0.01,0.06), c(2.0987844,0.8687591), type = "l",lwd=2, col= "black")#isolation success
points(c(0.01,0.06), c(2.015726,1.3100015), type = "l",lwd=2, col="red")#Richness
points(c(0.01,0.05), c(0.8670190,0.6722335), type = "l",lwd=2, col="blue")#Shannon
points(c(0.01,0.05), c(0.5550835,0.4743116), type = "l",lwd=2, col="green")#Simpson

#Dust
plot(c(0.01,0.06), c(0.5,3.2), type = "n", xlab="Dust deposition (mg/cm2)",  
     main = "Deposited dust", ylab = "",cex.axis = 1.5,
     cex.lab = 1.5,cex.main= 1.5)
points(c(0.01,0.06), c(2.1772373,2.5623075), type = "l",lwd=2, col= "black")#isolation success
points(c(0.01,0.06), c(2.144628,3.1278729), type = "l",lwd=2, col="red")#Richness
points(c(0.01,0.05), c(0.9417183,1.1881632), type = "l",lwd=2, col="blue")#Shannon
points(c(0.01,0.05), c(0.5854507,0.6784688), type = "l",lwd=2, col="green")#Simpson
legend("topleft",legend = c("Isolation success","Richness","Shannon diversity",
                             "Simpson diversity"),lwd=2,
      col= c("black","red","blue","green"),  border="white", bty="n", cex = 1.5)


# ## check to see if we have the same patterns in each locality and each time
# ###############
# ############### Locality
# ##### Isolation success
# ### I think it is better to have only one interaction model instead of several models.
# success.local= glm(MetaData$success~source*dust*locality, data = MetaData,
#                    family = poisson(link = "log") )
# anova(success.local)
# summary(success.local)
# plot(effect("source:dust:locality",success.local ))
# #library(sjPlot)
# #library(ggplot2)
# sjp.int(success.local)
# ## have no idea how to plot this
# 
# ##Locality models
# success.B = glm (MetaData$success[MetaData$locality == "BISOTON"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "BISOTON",],
#                  family = poisson(link = "log"))
# suc.B.anova= anova(success.B)
# suc.B.sum= summary(success.B)
# 
# 
# 
# plot(effect("dust",success.B
# effect("dust",success.B
# 
# success.H = glm (MetaData$success[MetaData$locality == "HASAN ABAD"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "HASAN ABAD",],
#                  family = poisson(link = "log"))
# suc.H.anova= anova(success.H)
# suc.H.sum= summary(success.H)
# 
# success.KH = glm (MetaData$success[MetaData$locality == "KHOSRO ABAD"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "KHOSRO ABAD",],
#                   family = poisson(link = "log"))
# suc.KH.anova= anova(success.KH)
# suc.KH.sum= summary(success.KH)
# 
# success.K = glm (MetaData$success[MetaData$locality == "KEREND"] ~  
#                     source * dust, data = MetaData[MetaData$locality == "KEREND",],
#                  family = poisson(link = "log"))
# suc.K.anova= anova(success.K)
# suc.K.sum= summary(success.K)
# 
# success.S = glm (MetaData$success[MetaData$locality == "SORKHE DIZE"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "SORKHE DIZE",],
#                  family = poisson(link = "log"))
# suc.S.anova= anova(success.S)
# suc.S.sum= summary(success.S)
# 
# ##plotting the locality  models
# ## effect objects from each model
# succ.b.eff = effect("source:dust", success.B, multiline=TRUE)#black
# succ.H.eff = effect("source:dust", success.H, multiline=TRUE)#red
# succ.KH.eff = effect("source:dust", success.KH, multiline=TRUE)#blue
# succ.k.eff = effect("source:dust", success.K, multiline=TRUE)#green
# succ.S.eff = effect("source:dust", success.S, multiline=TRUE)#pink
# rich.b.eff=effect("source:dust", Richness.B, multiline=TRUE)#black
# rich.h.eff=effect("source:dust", Richness.H, multiline=TRUE)#red
# rich.kh.eff=effect("source:dust", Richness.KH, multiline=TRUE)#blue
# rich.K.eff=effect("source:dust", Richness.K, multiline=TRUE)#Green
# rich.s.eff=effect("source:dust", Richness.S, multiline=TRUE)#pink
# shann.b.eff=effect("source:dust", shannon.B, multiline=TRUE)#black
# shann.h.eff=effect("source:dust", shannon.H, multiline=TRUE)#red
# shann.kh.eff=effect("source:dust", shannon.KH, multiline=TRUE)#blue
# shann.K.eff=effect("source:dust", shannon.K, multiline=TRUE)#Green
# shann.s.eff=effect("source:dust", shannon.S, multiline=TRUE)#pink
# simp.b.eff=effect("source:dust", simpson.B, multiline=TRUE)#black
# simp.h.eff=effect("source:dust", simpson.H, multiline=TRUE)#red
# simp.kh.eff=effect("source:dust", simpson.KH, multiline=TRUE)#blue
# simp.K.eff=effect("source:dust", simpson.K, multiline=TRUE)#Green
# simp.s.eff=effect("source:dust", simpson.S, multiline=TRUE)#pink
# ### The plot
# dev.off()
# par(mfrow=c(4,3))
# 
# ##isolation success
# #Leaf
# plot(c(0,0.06), c(0,2.5), type = "n", xlab="Dust deposition (mg/cm2)",ylab="Isolation success", 
#      main = "Leaf")
# points(c(0.01, 0.05), c(0.2748741,0.1686716), type = "l", col= "black")
# points(c(0.004, 0.014), c(0.2718915,0.0592935), type = "l", col="red")
# points(c(0.01, 0.06), c(0.7055354,1.044323), type = "l", col="blue")
# points(c(0.002, 0.014), c(0.7593056,0.3068659), type = "l", col="green")
# points(c(0.01, 0.04), c(0.01939612,2.340067), type = "l", col= "pink")
# #Branch
# plot(c(0,0.06), c(0,8.5), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Isolation success",
#      main = "Branch")
# points(c(0.01, 0.05), c(2.0790332,1.9561498), type = "l", col= "black")
# points(c(0.004, 0.014), c(2.8895416,0.7379235), type = "l", col="red")
# points(c(0.01, 0.06), c(2.6565222,1.592049), type = "l", col="blue")
# points(c(0.002, 0.014), c(2.1588502,1.3147342), type = "l", col="green")
# points(c(0.01, 0.04), c(1.93057642,8.454086), type = "l", col= "pink")
# #Dust
# plot(c(0,0.06), c(1.5,4), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Isolation success", 
#      main = "Dust")
# points(c(0.01, 0.05), c(1.9263297, 2.8308144), type = "l", col= "black")
# points(c(0.004, 0.014), c(1.8297366, 2.1506237), type = "l", col="red")
# points(c(0.01, 0.06), c(2.3279795, 3.652412), type = "l", col="blue")
# points(c(0.002, 0.014), c(1.6220277, 3.4907928), type = "l", col="green")
# points(c(0.01, 0.04), c(2.15173242, 3.140973), type = "l", col= "pink")
# ###Richness
# #Leaf
# plot(c(0,0.06), c(0,2.5), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Richness" 
#      )
# points(c(0.01, 0.05), c(1.191061,0.9498106), type = "l", col= "black")
# points(c(0.004, 0.014), c(1.974127,0.371652), type = "l", col="red")
# points(c(0.01, 0.06), c(1.508286,1.332691), type = "l", col="blue")
# points(c(0.002, 0.014), c(2.102010,0.8823088), type = "l", col="green")
# points(c(0.01, 0.04), c(1.000000,1.000000), type = "l", col= "pink")
# #Branch
# plot(c(0,0.06), c(1.5,3), type = "n", xlab="Dust deposition (mg/cm2)",ylab="Richness" 
#     )
# points(c(0.01, 0.05), c(2.033527, 2.1161365), type = "l", col= "black")
# points(c(0.004, 0.014), c(2.027248, 1.542201), type = "l", col="red")
# points(c(0.01, 0.06), c(2.426941, 2.819684), type = "l", col="blue")
# points(c(0.002, 0.014), c(1.863773, 2.5688436), type = "l", col="green")
# points(c(0.01, 0.04), c(2.095453, 2.043207), type = "l", col= "pink")
# #Dust
# plot(c(0,0.06), c(1.5,3.5), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Richness"
#     )
# points(c(0.01, 0.05), c(1.906941, 3.0415398), type = "l", col= "black")
# points(c(0.004, 0.014), c(2.211765, 2.017266), type = "l", col="red")
# points(c(0.01, 0.06), c(1.780617, 2.164073), type = "l", col="blue")
# points(c(0.002, 0.014), c(1.489016, 2.3935918), type = "l", col="green")
# points(c(0.01, 0.04), c(2.119341, 2.406066), type = "l", col= "pink")
# ###Shannon diversity
# #Leaf
# plot(c(0,0.06), c(-0.1,1.7), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Shannon diversity", 
#      )
# points(c(0.01, 0.05), c(0.7529783,1.0985422), type = "l", col= "black")
# points(c(0.004, 0.014), c(0.6931472,0.6931472), type = "l", col="red")
# points(c(0.01, 0.06), c(0.7089145,1.0335732), type = "l", col="blue")
# points(c(0.002, 0.014), c(1.3483858,1.684787), type = "l", col="green")
# points(c(0.01, 0.04), c(0,0), type = "l", col= "pink")
# #Branch
# plot(c(0,0.06), c(0.5,1.1), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Shannon diversity", 
#      )
# points(c(0.01, 0.05), c(0.7564696,0.7729018), type = "l", col= "black")
# points(c(0.004, 0.014), c(0.9195578,0.6375401), type = "l", col="red")
# points(c(0.01, 0.06), c(0.9228958,0.8483512), type = "l", col="blue")
# points(c(0.002, 0.014), c(1.0869643,0.638904), type = "l", col="green")
# points(c(0.01, 0.04), c(0.9357071,0.5588031), type = "l", col= "pink")
# #Dust
# plot(c(0,0.06), c(0.5,1.5), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Shannon diversity", 
#      )
# points(c(0.01, 0.05), c(0.9145406,1.2601045), type = "l", col= "black")
# points(c(0.004, 0.014), c(0.9591446,0.6159573), type = "l", col="red")
# points(c(0.01, 0.06), c(0.8851114,1.2226282), type = "l", col="blue")
# points(c(0.002, 0.014), c(0.8729864,1.209388), type = "l", col="green")
# points(c(0.01, 0.04), c(0.9382430,1.1192351), type = "l", col= "pink")
# ##Simpson diversity
# #Leaf
# plot(c(0,0.06), c(-0.1,1), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Simpson diversity", 
#      )
# points(c(0.01, 0.05), c(0.5263859,0.6787818), type = "l", col= "black")
# points(c(0.004, 0.014), c(0.5000000,0.5000000), type = "l", col="red")
# points(c(0.01, 0.06), c(0.5154674,0.8339494), type = "l", col="blue")
# points(c(0.002, 0.014), c(0.7356800,0.8627564), type = "l", col="green")
# points(c(0.01, 0.04), c(0,0), type = "l", col= "pink")
# #Branch
# plot(c(0,0.06), c(0.3,0.7), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Simpson diversity", 
#      )
# points(c(0.01, 0.05), c(0.5080873,0.5174149), type = "l", col= "black")
# points(c(0.004, 0.014), c(0.5782333,0.4667959), type = "l", col="red")
# points(c(0.01, 0.06), c(0.5743985,0.5394737), type = "l", col="blue")
# points(c(0.002, 0.014), c(0.6528082,0.4788129), type = "l", col="green")
# points(c(0.01, 0.04), c(0.5741159,0.3950454), type = "l", col= "pink")
# #Dust
# plot(c(0,0.06), c(0.3,0.8), type = "n", xlab="Dust deposition (mg/cm2)", ylab="Simpson diversity", 
#      )
# points(c(0.01, 0.05), c(0.5621166,0.7145125), type = "l", col= "black")
# points(c(0.004, 0.014), c(0.5992121,0.4655260), type = "l", col="red")
# points(c(0.01, 0.06), c(0.5626564,0.6665441), type = "l", col="blue")
# points(c(0.002, 0.014), c(0.5639695,0.6910460), type = "l", col="green")
# points(c(0.01, 0.04), c(0.5906979,0.6560817), type = "l", col= "pink")
# 
# 
# ####### Richness
# 
# Richness.B = glm (Richness[MetaData$locality == "BISOTON"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "BISOTON",],
#                  family = poisson(link = "log"))
# rich.B.anova= anova(Richness.B)
# rich.B.sum= summary(Richness.B)
# 
# Richness.H = glm (Richness[MetaData$locality == "HASAN ABAD"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "HASAN ABAD",],
#                  family = poisson(link = "log"))
# rich.H.anova= anova(Richness.H)
# rich.H.sum= summary(Richness.H)
# 
# Richness.KH = glm (Richness[MetaData$locality == "KHOSRO ABAD"] ~  
#                     source * dust, data = MetaData[MetaData$locality == "KHOSRO ABAD",],
#                   family = poisson(link = "log"))
# rich.KH.anova= anova(Richness.KH)
# rich.KH.sum= summary(Richness.KH)
# 
# Richness.K = glm (Richness[MetaData$locality == "KEREND"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "KEREND",],
#                  family = poisson(link = "log"))
# rich.K.anova= anova(Richness.K)
# rich.K.sum= summary(Richness.K)
# 
# Richness.S = glm (Richness[MetaData$locality == "SORKHE DIZE"] ~  
#                    source * dust, data = MetaData[MetaData$locality == "SORKHE DIZE",],
#                  family = poisson(link = "log"))
# rich.S.anova= anova(Richness.S)
# rich.S.sum= summary(Richness.S)
# 
# ####### Shannon Diversity
# shannon.B = lm (shannon[MetaNotOne$locality == "BISOTON"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "BISOTON",])
# shann.B.anova= anova(shannon.B)
# shann.B.sum= summary(shannon.B)
# 
# shannon.H = lm (shannon[MetaNotOne$locality == "HASAN ABAD"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "HASAN ABAD",])
# shann.H.anova= anova(shannon.H)
# shann.H.sum= summary(shannon.H)
# 
# shannon.KH = lm (shannon[MetaNotOne$locality == "KHOSRO ABAD"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "KHOSRO ABAD",])
# shann.KH.anova= anova(shannon.KH)
# shann.KH.sum= summary(shannon.KH)
# 
# shannon.K = lm (shannon[MetaNotOne$locality == "KEREND"] ~  
#                     source * dust, data = MetaNotOne[MetaNotOne$locality == "KEREND",])
# shann.K.anova= anova(shannon.K)
# shann.K.sum= summary(shannon.K)
# 
# shannon.S = lm (shannon[MetaNotOne$locality == "SORKHE DIZE"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "SORKHE DIZE",])
# shann.S.anova= anova(shannon.S)
# shann.S.sum= summary(shannon.S)
# 
# ####### Simpson Diversity
# simpson.B = lm (simpson[MetaNotOne$locality == "BISOTON"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "BISOTON",])
# simp.B.anova= anova(simpson.B)
# simp.B.sum= summary(simpson.B)
# 
# simpson.H = lm (simpson[MetaNotOne$locality == "HASAN ABAD"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "HASAN ABAD",])
# simp.H.anova= anova(simpson.H)
# simp.H.sum= summary(simpson.H)
# 
# simpson.KH = lm (simpson[MetaNotOne$locality == "KHOSRO ABAD"] ~  
#                    source * dust, data = MetaNotOne[MetaNotOne$locality == "KHOSRO ABAD",])
# simp.KH.anova= anova(simpson.KH)
# simp.KH.sum= summary(simpson.KH)
# 
# simpson.K = lm (simpson[MetaNotOne$locality == "KEREND"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "KEREND",])
# simp.K.anova= anova(simpson.K)
# simp.K.sum= summary(simpson.K)
# 
# simpson.S = lm (simpson[MetaNotOne$locality == "SORKHE DIZE"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$locality == "SORKHE DIZE",])
# simp.S.anova= anova(simpson.S)
# simp.S.sum= summary(simpson.S)
# 
# 
# 
# #########
# ######### TIME
# ### ISOLATION SUCCESS
# success.1 = glm (MetaData$success[MetaData$time == "1"] ~  
#                    source * dust, data = MetaData[MetaData$time == "1",],
#                  family = poisson(link = "log"))
# suc.1.anova= anova(success.1)
# suc.1.sum= summary(success.1)
# 
# success.2 = glm (MetaData$success[MetaData$time == "2"] ~  
#                    source * dust, data = MetaData[MetaData$time == "2",],
#                  family = poisson(link = "log"))
# suc.2.anova= anova(success.2)
# suc.2.sum= summary(success.2)
# 
# success.3 = glm (MetaData$success[MetaData$time == "3"] ~  
#                     source * dust, data = MetaData[MetaData$time == "3",],
#                   family = poisson(link = "log"))
# suc.3.anova= anova(success.3)
# suc.3.sum= summary(success.3)
# 
# ####### Richness
# Richness.1 = glm (Richness[MetaData$time == "1"] ~  
#                     source * dust, data = MetaData[MetaData$time == "1",],
#                   family = poisson(link = "log"))
# rich.1.anova= anova(Richness.1)
# rich.1.sum= summary(Richness.1)
# 
# Richness.2 = glm (Richness[MetaData$time == "2"] ~  
#                     source * dust, data = MetaData[MetaData$time == "2",],
#                   family = poisson(link = "log"))
# rich.2.anova= anova(Richness.2)
# rich.2.sum= summary(Richness.2)
# 
# Richness.3 = glm (Richness[MetaData$time == "3"] ~  
#                      source * dust, data = MetaData[MetaData$time == "3",],
#                    family = poisson(link = "log"))
# rich.3.anova= anova(Richness.3)
# rich.3.sum= summary(Richness.3)
# 
# ####### Shannon Diversity
# shannon.1 = lm (shannon[MetaNotOne$time == "1"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$time == "1",])
# shann.1.anova= anova(shannon.1)
# shann.1.sum= summary(shannon.1)
# 
# shannon.2 = lm (shannon[MetaNotOne$time == "2"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$time == "2",])
# shann.2.anova= anova(shannon.2)
# shann.2.sum= summary(shannon.2)
# 
# shannon.3 = lm (shannon[MetaNotOne$time == "3"] ~  
#                    source * dust, data = MetaNotOne[MetaNotOne$time == "3",])
# shann.3.anova= anova(shannon.3)
# shann.3.sum= summary(shannon.3)
# 
# ####### Simpson Diversity
# simpson.1 = lm (simpson[MetaNotOne$time == "1"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$time == "1",])
# simp.1.anova= anova(simpson.1)
# simp.1.sum= summary(simpson.1)
# 
# simpson.2 = lm (simpson[MetaNotOne$time == "2"] ~  
#                   source * dust, data = MetaNotOne[MetaNotOne$time == "2",])
# simp.2.anova= anova(simpson.2)
# simp.2.sum= summary(simpson.2)
# 
# simpson.3 = lm (simpson[MetaNotOne$time == "3"] ~  
#                    source * dust, data = MetaNotOne[MetaNotOne$time == "3",])
# simp.3.anova= anova(simpson.3)
# simp.3.sum= summary(simpson.3)
#######################################################
#######################################################
### 3- Community composition model
fun.Mvabund = mvabund(MyAbund)
plot(fun.Mvabund)

fun.Mvabund.m = manyglm (fun.Mvabund ~ locality + time + source*dust, data= MetaData,
                    family="negative.binomial", show.residuals=T)
plot.manyglm(fun.Mvabund.m)
fun.Mvabund.m.sum = summary.manyglm (fun.Mvabund.m, nBoot=100, test="LR",p.uni="adjusted", 
                                  resamp="montecarlo")
## Analysis of variance explained by the predictors
fun.Mvabund.m.anova = anova.manyglm (fun.Mvabund.m, nBoot=300, test="LR", p.uni="adjusted", 
                                resamp="montecarlo")
fun.Mvabund.m.anova100 = anova.manyglm (fun.Mvabund.m, nBoot=100, test="LR", p.uni="adjusted", 
                                     resamp="montecarlo")
## ploting the coef not working. 
coef.manyglm=coef(fun.Mvabund.m)
plot(coef.manyglm)
## Coefficients
fun.Mvabund.m.coef = as.data.frame(fun.Mvabund.m$coefficients)
## mean-centering the contrasts
fun.Mvabund.m.coef.contrast = fun.Mvabund.m.coef - apply (fun.Mvabund.m.coef,2,mean)

# ## OTUs significantly affected by the source??
mvabund.m.anova <- as.data.frame(fun.Mvabund.m.anova$uni.p)
fun.Mvabund.m.source = colnames(mvabund.m.anova)[mvabund.m.anova ["source",] <= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["source",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["locality",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["time",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["dust",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["source:dust",]<= 0.05]
### Rownames of OTUs affected by the source
row.names(MyAbund)[MyAbund[ ,"Alternaria_sp_A9"] >0]
row.names(MyAbund)[MyAbund[ ,"Aspergillus_sp_A20"] >0]
row.names(MyAbund)[MyAbund[ ,"Aureobasidium_sp_A17"] >0]
row.names(MyAbund)[MyAbund[ ,"Aureobasidium_sp_A30"] >0]
row.names(MyAbund)[MyAbund[ ,"Byssochlamys_spectabilis_"] >0]
row.names(MyAbund)[MyAbund[ ,"Cladosporium_herbarum_A8"] >0]
row.names(MyAbund)[MyAbund[ ,"Dothideomycetes_sp_A1"] >0]
row.names(MyAbund)[MyAbund[ ,"Gnomoniaceae_sp_66"] >0]
row.names(MyAbund)[MyAbund[ ,"Microsphaeriopsis_olivacea"] >0]
row.names(MyAbund)[MyAbund[ ,"Preussia_sp_A31"] >0]
row.names(MyAbund)[MyAbund[ ,"Quambalaria_cyanescens"] >0]
row.names(MyAbund)[MyAbund[ ,"Ustilago_A14"] >0]
row.names(MyAbund)[MyAbund[ ,"Ustilago_A15"] >0]
row.names(MyAbund)[MyAbund[ ,"Ustilago_A16"] >0]
### source*dust interactions plot
### Negative binomial glm for the species affected by the source* dust interaction
Microsphaeriopsis.model= glm.nb(MyAbund$Microsphaeriopsis_olivacea ~ locality +
                                  time + source*dust,data= MetaData)
Micros.anova = anova (Microsphaeriopsis.model, test = "Chisq")                             
Micros.summary= summary(Microsphaeriopsis.model)
# ploting the effect
Micros.effec= effect("source:dust", Microsphaeriopsis.model, multiline = FALSE)
Micros.effec.sum = summary (Micros.effec)
Mic.eff<-plot(Micros.effec, ylab = "Microsphaeropsis olivacea abundance",
              xlab = "Dust Deposition (mg/cm2)",main= NULL, ylim = c(-30,5),
              se=FALSE,ci.style="none")
Mic.eff # plot original object
Mic.eff$x.scales$cex <- c(1.3, 1.3)
Mic.eff$xlab$cex <- 1.3
Mic.eff$ylab$cex <- 1.3
Mic.eff # plot edited object

#### visualize how other OTUs are affected by the Source*Dust interaction
model.list = list()
for( i in names(MyAbund)){
  i.model = glm.nb(MyAbund[,i] ~ locality +
                                  time + source*dust,data= MetaData)
  # plot(effect("source:dust", i.model))
  model.list[[i]] <- i.model
  
}

pdf(file = "all species reaction.pdf")
for( i in names(model.list)){
  plot(effect("source:dust", model.list[[i]]))
  print(model.list[[i]])
}
dev.off()

pdf(file = "all species reaction.pdf", paper = "USr")
plot(effect("source:dust", model.list[["Alternaria_sp_A22"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Alternaria_sp_A22 abundance",main= NULL)
plot(effect("source:dust", model.list[["Alternaria_sp_A25"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Alternaria_sp_A25 abundance",main= NULL)
plot(effect("source:dust", model.list[["Alternaria_sp_A76"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Alternaria_sp_A76 abundance",main= NULL)
plot(effect("source:dust", model.list[["Alternaria_sp_A9"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Alternaria_sp_A9 abundance",main= NULL)
plot(effect("source:dust", model.list[["Arthrinium_marii"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Arthrinium_marii abundance",main= NULL)
plot(effect("source:dust", model.list[["Aspergillus_sp_A20"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Aspergillus_sp_A20 abundance",main= NULL)
plot(effect("source:dust", model.list[["Aspergillus_sp_B20"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Aspergillus_sp_B20 abundance",main= NULL)
plot(effect("source:dust", model.list[["Aspergillus_sp_C7"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Aspergillus_sp_C7 abundance",main= NULL)
plot(effect("source:dust", model.list[["Aureobasidium_sp_A17"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Aureobasidium_sp_A17 abundance",main= NULL)
plot(effect("source:dust", model.list[["Aureobasidium_sp_A30"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Aureobasidium_sp_A30 abundance",main= NULL)
plot(effect("source:dust", model.list[["B25"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "B25 abundance",main= NULL)
plot(effect("source:dust", model.list[["Beauveria_sp_B15"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Beauveria_sp_B15 abundance",main= NULL)
plot(effect("source:dust", model.list[["Biscogniauxia_mediterranea"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Biscogniauxia_mediterranea abundance",main= NULL)
plot(effect("source:dust", model.list[["Byssochlamys_spectabilis_"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Byssochlamys_spectabilis abundance",main= NULL)
plot(effect("source:dust", model.list[["Chaetomiaceae_sp_A37"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Chaetomiaceae_sp_A37 abundance",main= NULL)
plot(effect("source:dust", model.list[["Cladosporium_herbarum_A8"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Cladosporium_herbarum_A8 abundance",main= NULL)
plot(effect("source:dust", model.list[["Comoclathris_sedi"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Comoclathris_sedi abundance",main= NULL)
plot(effect("source:dust", model.list[["Coniochaeta_sp_A85"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Coniochaeta_sp_A85 abundance",main= NULL)
plot(effect("source:dust", model.list[["Coniothyrium_sp_A41"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Coniothyrium_sp_A41 abundance",main= NULL)
plot(effect("source:dust", model.list[["Coniothyrium_sp_A43"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Coniothyrium_sp_A43 abundance",main= NULL)
plot(effect("source:dust", model.list[["Cytospora_sp_AC35"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Cytospora_sp_AC35 abundance",main= NULL)
plot(effect("source:dust", model.list[["Cytospora_sp_C2"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Cytospora_sp_C2 abundance",main= NULL)
plot(effect("source:dust", model.list[["Diatrype_sp_C1"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Diatrype_sp_C1 abundance",main= NULL)
plot(effect("source:dust", model.list[["Diatrypella_sp_C6"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Diatrypella_sp_C6 abundance",main= NULL)
plot(effect("source:dust", model.list[["Dikarya_sp_A38"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Dikarya_sp_A38 abundance",main= NULL)
plot(effect("source:dust", model.list[["Dothideomycetes_sp_A1"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Dothideomycetes_sp_A1 abundance",main= NULL)
plot(effect("source:dust", model.list[["Dothideomycetes_sp_A48"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Dothideomycetes_sp_A48 abundance",main= NULL)
plot(effect("source:dust", model.list[["Dothideomycetes_sp_A79"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Dothideomycetes_sp_A79 abundance",main= NULL)
plot(effect("source:dust", model.list[["Endoconidioma_populi_A39"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Endoconidioma_populi_A39",main= NULL)
plot(effect("source:dust", model.list[["Fusarium_sp_46"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Fusarium_sp_46 abundance",main= NULL)
plot(effect("source:dust", model.list[["Gnomoniaceae_sp_70"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Gnomoniaceae_sp_70 abundance",main= NULL)
plot(effect("source:dust", model.list[["Gnomoniaceae_sp_66"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Gnomoniaceae_sp_66 abundance",main= NULL)
plot(effect("source:dust", model.list[["Helotiales_sp_A68"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Helotiales_sp_A68 abundance",main= NULL)
plot(effect("source:dust", model.list[["Humicola_sp_A52"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Humicola_sp_A52 abundance",main= NULL)
plot(effect("source:dust", model.list[["leotiomyceta_sp_A42"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "leotiomyceta_sp_A42 abundance",main= NULL)
plot(effect("source:dust", model.list[["Paraphaeosphaeria_sp_B10"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Paraphaeosphaeria_sp_B10 abundance",main= NULL)
plot(effect("source:dust", model.list[["Penicillium_sp_A21"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Penicillium_sp_A21 abundance",main= NULL)
plot(effect("source:dust", model.list[["Penicillium_sp_A3"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Penicillium_sp_A3 abundance",main= NULL)
plot(effect("source:dust", model.list[["Pleosporaceae_sp_A5"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Pleosporaceae_sp_A5 abundance",main= NULL)
plot(effect("source:dust", model.list[["Pleosporaceae_sp_A78"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Pleosporaceae_sp_A78 abundance",main= NULL)
plot(effect("source:dust", model.list[["Pleosporaceae_sp_B27"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Pleosporaceae_sp_B27 abundance",main= NULL)
plot(effect("source:dust", model.list[["Pleosporales_sp_A63"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Pleosporales_sp_A63 abundance",main= NULL)
plot(effect("source:dust", model.list[["Pleosporineae_sp_A23"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Pleosporineae_sp_A23 abundance",main= NULL)
plot(effect("source:dust", model.list[["Preussia_africana"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Preussia_africana abundance",main= NULL)
plot(effect("source:dust", model.list[["Preussia_australis"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Preussia_australis abundance",main= NULL)
plot(effect("source:dust", model.list[["Preussia_complex_sp_A36"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Preussia_complex_sp_A36 abundance",main= NULL)
plot(effect("source:dust", model.list[["Preussia_funiculata"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Preussia_funiculata abundance",main= NULL)
plot(effect("source:dust", model.list[["Preussia_intermedia"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Preussia_intermedia abundance",main= NULL)
plot(effect("source:dust", model.list[["Preussia_sp_A31"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Preussia_sp_A31 abundance",main= NULL)
plot(effect("source:dust", model.list[["Quambalaria_cyanescens"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Quambalaria_cyanescens abundance",main= NULL)
plot(effect("source:dust", model.list[["Schizothecium_sp_B14"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Schizothecium_sp_B14 abundance",main= NULL)
plot(effect("source:dust", model.list[["Scytalidium_thermophilum"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Scytalidium_thermophilum abundance",main= NULL)
plot(effect("source:dust", model.list[["Seimatosporium_pezizoides"]]), xlab = "Dust Deposition (mg/cm2)", ylab = "Seimatosporium_pezizoides abundance",main= NULL)
dev.off()
#### ploting all the OTUs that shows a reaction to the interaction
Alternaria_sp_A25.model.effec= effect("source:dust",Alternaria_sp_A25.model)
Alternaria_sp_A9.model.effec= effect("source:dust", Alternaria_sp_A9.model)
Aureobasidium_sp_A30.model.effec= effect("source:dust", Aureobasidium_sp_A30.model)
Byssochlamys.effec= effect("source:dust", Byssochlamys_spectabilis_.model)
Cladosporium_herbarum_A8.model.effec= effect("source:dust", Cladosporium_herbarum_A8.model)
Penicillium_sp_A21.model.effec= effect("source:dust", Penicillium_sp_A21.model)
Pleosporaceae_sp_A5.model.effec= effect("source:dust", Pleosporaceae_sp_A5.model)
### The reaction plot
dev.off()
par(mfrow=c(1,3))
#Leaf
plot(c(0.01,0.06), c(0.0,2.220446e-2), type = "n", xlab="Dust deposition (mg/cm2)",
     main = "Leaf", ylab = "",cex.axis = 1.5,
     cex.lab = 1.5, cex.main= 1.5)
points(c(0.01,0.06), c(0.00119931,7.373837e-14), type = "l", lty=6, lwd=2, col="black")#A25
points(c(0.01,0.06), c(0.006563875,1.436627e-07), type = "l",lty=5, lwd=2, col="blue")#A9
points(c(0.01,0.06), c(2.205382e-06,2.220446e-16), type = "l",lty=3, lwd=2, col= "red")#A30
points(c(0.01,0.06), c(0.004881679,0.001257212), type = "l",lty=1,lwd=2, col= "green")#Byssochlamys
points(c(0.01,0.06), c(0.002904384,0.01811648), type = "l",lwd=2,lty=6, col="red")#A8
points(c(0.01,0.06), c(0.002465486,2.867857e-13), type = "l",lwd=2,lty=4, col= "blue")#A21
points(c(0.01,0.06), c(0.0009805355,5.550239e-09), type = "l",lwd=2,lty=3, col= "green")#A5
#Branch
plot(c(0.01,0.06), c(0.0,4.224453e-2), type = "n", xlab="Dust deposition (mg/cm2)",
     main = "Branch", ylab = "",cex.axis = 1.5,
     cex.lab = 1.5,cex.main= 1.5)
points(c(0.01,0.06), c(0.03120315,1.575250e-05), type = "l",lty=6, lwd=2, col="black")#25
points(c(0.01,0.06), c(0.039682698,5.224453e-04), type = "l",lty=5, lwd=2, col="blue")#A9
points(c(0.01,0.06), c(3.372862e-03,3.891702e-01), type = "l",lty=3, lwd=2, col= "red")#A30
points(c(0.01,0.06), c(0.311248787,0.273451531), type = "l",lty=1,lwd=2, col= "green")#Byssochlamys this is not shown
points(c(0.01,0.06), c(0.032479263,0.02415204), type = "l",lwd=2,lty=6, col="red")#A8
points(c(0.01,0.06), c(0.010226362,3.445531e-08), type = "l",lwd=2,lty=4, col= "blue")#A21
points(c(0.01,0.06), c(0.0144565757,7.301888e-04), type = "l",lwd=2,lty=3, col= "green")#A5
#Dust
plot(c(0.01,0.06), c(0.00,10e-1), type = "n", xlab="Dust deposition (mg/cm2)",
     main = "Deposited dust", ylab = "",cex.axis = 1.5,
     cex.lab = 1.5,cex.main= 1.5)
points(c(0.01,0.06), c(0.05675319,4.689172e-02), type = "l",lty=6, lwd=2, col="black")#A25
points(c(0.01,0.06), c(0.154575206,9.479314e-02), type = "l",lty=5, lwd=2, col="blue")#A9
points(c(0.01,0.06), c(1.537366e-01,1.470763e-01), type = "l",lty=3, lwd=2, col= "red")#A30
points(c(0.01,0.06), c(0.052085705,0.102723379), type = "l",lty=1, lwd=2, col= "green")#Byssochlamys
points(c(0.01,0.06), c(0.277383476,0.15360726), type = "l",lwd=2,lty=6, col="red")#A8
points(c(0.01,0.06), c(0.009849907,4.371290e-05), type = "l", lwd=2,lty=4, col= "blue")#A21
points(c(0.01,0.06), c(0.0238803385,8.354105e-01), type = "l",lwd=2,lty=3, col= "green")#A5

legend(0.04, 1.03, c("A25", "A9", "A30","Bys","A8","A21","A5"), col = c("black", "blue", 
                 "red","green","red","blue","green"),box.col="white",
       lty = c(6,5,3,1,6,4,3), cex = 1.5)
#### Model-based ordination
# Remove species present in only one sample
CountInSample = apply(AbundNotZero,2,sum)
# I needed to remove the species that were seen in a few samples.
fung.abun.reduced = AbundNotZero[,CountInSample > 3]
fung.abun.reduced = fung.abun.reduced[apply(fung.abun.reduced,1,sum) > 0,]
MetaOrd = MetaRich[rownames(MetaRich) %in% rownames(fung.abun.reduced),]
ModelOrd <- boral(fung.abun.reduced, family = "negative.binomial", num.lv = 2, 
                  n.burnin = 10, n.iteration = 100, n.thin = 1)
## model diagnostics
par(mfrow = c(2,2))
plot(ModelOrd, ask = FALSE )
# Ordination plot.
dev.off()
par(mar=c(4,4,1,1))
ordibora= ordiplot(ModelOrd$lv.median, choices = c(1,2),display = "sites",
                   type = "none", cex = 1.5 , xlim = c(-0.3,0.3), cex.axis = 1.5,
         cex.lab = 1.2)
points(ordibora,"sites", pch=20 ,col=as.numeric(MetaOrd$source))
#ordispider(ordibora,MetaOrd$source, col= "gray" )
mylegend = legend(0.7, 0.9, c("Leaf","Branch","Dust"), 
                  c("black", "red", "green"), border="white", bty="n", cex = 1.3)
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col="black",
                          alpha=200,kind="se",conf=0.95, 
                          show.groups=(c("Leaf"))))
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col="red",
                          alpha=150,kind="se",conf=0.95,
                          show.groups=(c("Branch")))) 
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col="green",
                          alpha=200,kind="se",conf=0.95,
                          show.groups=(c("Dust"))))
# ordisurf of dust: show dust + source interaction on the whole community
#ordisurf(ModelOrd$lv.median, MetaOrd$dust, add=T, col = "blue", labcex= 1.3)
# Insert the source into the rownames of the point coordinates
ModelToDist = ModelOrd$lv.median
rownames(ModelToDist) = paste(MetaOrd$source, rownames(ModelToDist), sep = ".")
# Calculate the distance among the sites
ModelDist = dist(ModelToDist, method = "euclidean")

BranchLeaf_Dist = melt( as.matrix(ModelDist)[grep("Branch", rownames(ModelToDist)), 
                                             grep("Leaf", rownames(ModelToDist))] )

BranchDust_Dist = melt( as.matrix(ModelDist)[grep("Branch", rownames(ModelToDist)), 
                                             grep("Dust", rownames(ModelToDist))] )

DustLeaf_Dist = melt( as.matrix(ModelDist)[grep("Dust", rownames(ModelToDist)), 
                                           grep("Leaf", rownames(ModelToDist))] )

GroupedDist = rbind(BranchLeaf_Dist, BranchDust_Dist, DustLeaf_Dist)
GroupedDist = data.frame(GroupedDist, 
                         groups = c(rep("Branch2Leaf", nrow(BranchLeaf_Dist)), 
                                    rep("Branch2Dust", nrow(BranchDust_Dist)),
                                    rep("Dust2Leaf", nrow(DustLeaf_Dist))))

# Post-hoc testing of group distances
Tukey.dist= TukeyHSD (aov(value ~ groups, data = GroupedDist))

# It makes sense if you did a good surface sterilization: probably it is easier 
# for a fungus to go into a leaf than into a branch
boxplot(GroupedDist$value ~GroupedDist$groups)


