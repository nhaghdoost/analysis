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

# Sample names:
# B1T1L1:
# B: locality of collection (B = Bisoton)
# 1: sampling time 1
# T1: the code of the tree at that site 
# L1: leaf 1 of that tree (can be B=branch, D=dust)

# Models of culturing success
# success.m1 = glm(success ~ locality + time + source + dust, data=MetaData, 
#                  family = poisson(link = "log"))
success.m = glm (success ~ locality + time + source * dust, data=MetaData, 
               family = poisson(link = "log"))

# Compare model fit: succes.m2 seems better
# AIC(success.m1)
# AIC(success.m2)
# 
# par(mfrow = c(2,2))
# plot(success.m1)
# plot(success.m2)

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
success.effect = effect("source:dust",success.m ,  multiline=TRUE)
plot(success.effect, xlab = "Dust Deposition (mg/cm2)", ylab = "Isolation success", ylim = c(-2.5,1.5))
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
# What does it mean? With predicted values
par(mfrow=c(1,1))
boxplot(fitted(Richness.m) ~ MetaRich$source, xlab="Source", ylab="Richness")
# Richness in branch and dust similar, but richness stat. sign. lower in leaf

# Interaction plots:
Richness.effect = effect("source:dust",Richness.m, multiline=TRUE, confidence.level = 0.95)
plot(Richness.effect, ylab = "Richness", xlab = "Dust Deposition (mg/cm2)", ylim = c(-1.5,2))
Rich.eff.sum = summary(Richness.effect)

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
shannon.m = lm (shannon ~ locality + time + source*dust, data = MetaNotOne)

## Shannon model statistics
shannon.summary = summary (shannon.m) 
shannon.anova = anova (shannon.m, test = "LRT")

# Interaction plots:
shannon.effect = effect("source:dust", shannon.m, multiline=TRUE)
shan.eff.sum = summary (shannon.effect)
plot(shannon.effect,ylab = "Shannon Diversity", xlab = "Dust Deposition (mg/cm2)", ylim = c(-0.1,1.5) )

# Put here other boxplots if you wnat to show other effects, e.g. time or locality
# Plot the interaction effect
# Either plot them on three different plots, 
# or scale the ylim values for each of them.
# plot(shannon.effect, ylim = c(-1,2))


## Simpson model:
simpson.m = lm (simpson ~ locality + time + source*dust, data=MetaNotOne)
simpson.summary = summary (simpson.m)
simpson.anova = anova (simpson.m, test= "Chisq")

# Other boxplot if needed here:

# Interaction plots:
simpson.effect = effect("source:dust", simpson.m, multiline=TRUE)

# plot interaction effects
plot(simpson.effect,ylab = "Simpson Diversity", xlab = "Dust Deposition (mg/cm2)", ylim = c(-0.1,1))

# Effect summary
Simp.eff.sum = summary (simpson.effect)

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
fun.Mvabund.m.anova = anova.manyglm (fun.Mvabund.m, nBoot=100, test="LR", p.uni="adjusted", 
                                resamp="montecarlo")

# ## OTUs significantly affected by the source?? (at p<=0.001???)
mvabund.m.anova <- as.data.frame(fun.Mvabund.m.anova$uni.p)

fun.Mvabund.m.source = colnames(mvabund.m.anova)[mvabund.m.anova ["source",] <= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["source",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["locality",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["time",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["dust",]<= 0.05]
colnames(mvabund.m.anova)[mvabund.m.anova["source:dust",]<= 0.01]

###source*dust interactions plot
### Negative binomial glm for the species affected by the source* dust interaction

Microsphaeriopsis.model= glm.nb(MyAbund$Microsphaeriopsis_olivacea ~ locality + time + source*dust,
                                data= MetaData)
Micros.anova = anova (Microsphaeriopsis.model, test = "Chisq")                             

Micros.effec= effect("source:dust", Microsphaeriopsis.model)
Micros.effec.sum = summary (Micros.effec)
plot(Micros.effec, ylab = "Microsphaeriopsis_olivacea", xlab = "Dust Deposition (mg/cm2)", ylim = c(-35,3))

## Coefficients
fun.Mvabund.m.coef = as.data.frame(fun.Mvabund.m$coefficients)

## mean-centering the contrasts
fun.Mvabund.m.coef.contrast = fun.Mvabund.m.coef - apply (fun.Mvabund.m.coef,2,mean)

# ### plot effects
# 
# predict(funcor.m2, type="response")
# 
# plot(metacor0$dust, predict(funcor.m2, type="response")[,"Gnomoniaceae_sp_66"])
# 
# plot(metacor0$dust[metacor0$source == "Leaf"], predict(funcor.m2, type="response")
#      [metacor0$source == "Leaf","Gnomoniaceae_sp_66"])
# plot(metacor0$dust[metacor0$source == "Branch"], predict(funcor.m2, type="response")
#      [metacor0$source == "Branch","Gnomoniaceae_sp_66"])
# plot(metacor0$dust[metacor0$source == "Dust"], predict(funcor.m2, type="response")
#      [metacor0$source == "Dust","Gnomoniaceae_sp_66"])
# plot(metacor0$dust[metacor0$source == "Leaf"],
#      Corfun0[metacor0$source == "Leaf","Gnomoniaceae_sp_66"])
# plot(metacor0$dust[metacor0$source == "Leaf"], predict(funcor.m2, type="response")
#      [metacor0$source == "Leaf","Gnomoniaceae_sp_66"])


###source*dust interactions plot
### Do a simple glm for the two species affected by the source* dust interaction
# 
# Microsphaeriopsis.model= glm.nb(Corfun0$Microsphaeriopsis_olivacea ~ locality+time+source*dust, data= metacor0)
# Micros.anova= anova(Microsphaeriopsis.model)                             
# 
# Micros.effec= effect("source:dust",Microsphaeriopsis.model ,multiline=TRUE,  ylim=c(-10,10))
# Micros.effec.sum = summary(Micros.effec)
# 
# 
# Aureobasidium.model= glm.nb (Corfun0$Aureobasidium_sp_A30 ~ locality+time+source*dust, data= metacor0)
# 
# plot(effect("source:dust",Aureobasidium.model ,multiline=TRUE,confidence.level = 0.95))
# 
# anova(Aureobasidium.model)
### similarity analysis using anosim and adonis functions (corfun0=fungal abundance matrix for core OTUs and metacor0 
# ## is the metadata)
# braysimi= anosim(Corfun0,metacor0$source,permutations = 999,distance = "bray")
# 
# jaccardsimi= anosim(Corfun0,metacor0$source, permutations = 999, distance = "jaccard")
# 
# 
# corfunAdon= adonis(Corfun0~locality+time+source*dust,data= metacor0,permutations = 999, method = "bray")
# 
# ### permanova for all of the species
# 
permanova.total = adonis (MyAbund ~ locality + time + source*dust, data= MetaData 
                       ,permutations = 999, method = "bray")

#### Model-based ordination

# Remove species present in only one sample
CountInSample = apply(AbundNotZero,2,sum)

# I needed to remove the species that were seen in a few samples.
fung.abun.reduced = AbundNotZero[,CountInSample > 3]
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
#ordispider(ordibora,MetaOrd$source, col= "gray" )
mylegend = legend(0.7, 0.9, c("Leaf","Branch","Dust"), 
                  fill=c("green", "red", "gray"), border="white", bty="n")
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col="green",
                          alpha=200,kind="se",conf=0.95, 
                          show.groups=(c("Leaf"))))
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col="red",
                          alpha=150,kind="se",conf=0.95,
                          show.groups=(c("Branch")))) 
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col="gray",
                          alpha=200,kind="se",conf=0.95,
                          show.groups=(c("Dust"))))
# ordisurf of dust: show dust + source interaction on the whole community
ordisurf(ModelOrd$lv.median, MetaOrd$dust, add=T, col = "blue")


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




