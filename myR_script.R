



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

## ploting the coef
## Coefficients
fun.Mvabund.m.coef = as.data.frame(fun.Mvabund.m$coefficients)

## mean-centering the contrasts
fun.Mvabund.m.coef.contrast = fun.Mvabund.m.coef - apply (fun.Mvabund.m.coef,2,mean)
### use the confidence intervas and coefficient for each species

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
###source*dust interactions plot
### Negative binomial glm for the species affected by the source* dust interaction

Microsphaeriopsis.model= glm.nb(MyAbund$Microsphaeriopsis_olivacea ~ locality + time + source*dust,
                               data= MetaData)
Micros.anova = anova (Microsphaeriopsis.model, test = "Chisq")                             
Micros.summary= summary(Microsphaeriopsis.model)
# ploting the effect
Micros.effec= effect("source:dust", Microsphaeriopsis.model)
Micros.effec.sum = summary (Micros.effec)
plot(Micros.effec, ylab = "Microsphaeriopsis_olivacea", xlab = "Dust Deposition (mg/cm2)", 
     ylim = c(-30,5))






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

# Ordination plot. Please make it similar as the NMDS.
dev.off()
par(mar=c(4,4,1,1))
ordibora= ordiplot(ModelOrd$lv.median, choices = c(1,2), type = "none", cex =0.5 
         ,display = "sites", xlim = c(-0.3,0.3))
points(ordibora,"sites", pch=20 ,col=as.numeric(MetaOrd$source))
#ordispider(ordibora,MetaOrd$source, col= "gray" )
mylegend = legend(0.7, 0.9, c("Leaf","Branch","Dust"), 
                  fill=c("gray", "red", "green"), border="white", bty="n")
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=.5, 
                          draw="polygon", col="gray",
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
### Because the NMDS plot likely has overdispersion effects.

