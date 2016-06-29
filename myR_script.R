
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
table(MetaData$locality, MetaData$time)
table(MetaData$locality, MetaData$source)
table(MetaData$time, MetaData$source)
table(MetaData$tree, MetaData$source)
table(MetaData$tree, MetaData$time)
table(MetaData$tree, MetaData$source)
table(MetaData$success, MetaData$source)

# Continuous variables
hist(MetaData$success, nclass=20)
summary(MetaData$success)
hist(MetaData$dust, xlab="Dust (mg/cm2)")

# Some simple plots
# Time?
boxplot(MetaData$success ~ MetaData$time)
# Locality and dust quantity
boxplot(MetaData$dust ~ MetaData$locality)
# Dust and sampling time
boxplot(MetaData$dust ~ MetaData$time)
# Dust quantity at localities at different times
boxplot(MetaData$dust[MetaData$time == 1] ~ 
          MetaData$locality[MetaData$time == 1])
boxplot(MetaData$dust[MetaData$time == 2] ~ 
          MetaData$locality[MetaData$time == 2])
boxplot(MetaData$dust[MetaData$time == 3] ~ 
          MetaData$locality[MetaData$time == 3])

# Q1: Dust deposition and isolation success
plot(MetaData$dust, MetaData$success, xlab=c("Dust (mg/cm2)"), 
         ylab=c("Isolation success"))

# Q2: Locality and isolation success
boxplot(MetaData$success ~ MetaData$locality, 
        ylab="Isolation success")

# Q4: Source of isolation
boxplot(MetaData$success ~ MetaData$source)

# different 
par(mar = c(12,4,3,1))
boxplot(MetaData$success ~ MetaData$locality*MetaData$source, 
        las=2)
dev.off()

# Models of culturing success

success.m1 = glm(success ~ time + locality + source + dust, data=MetaData, family = poisson(link = "log"))
summary(success.m1)
anova(success.m1, test="LRT")
# plot(success.m1)

### interaction model
success.m2=glm(success~locality+time+source*dust, data=MetaData, family = poisson(link = "log"))
success.sammary=summary(success.m2)
success.anova=anova(success.m2, test="LRT")
# plot(success.m2)

## plot the interaction
library(effects)
plot(effect("source:dust",success.m2 ,  multiline=TRUE))

AIC(success.m1)
AIC(success.m2)
## based on this the m2 model is better for describing our data

# Coefficients and 95% confidence intervals
Effs = cbind(Coefficients = coefficients(success.m2), confint(success.m2))

###########################
###########################
#### 2- Diversity analysis:
### Diversity Data input

fungal.abundance= read.csv(file = "morphotype_matrix_incubator.csv", header = T, sep=",", row.names = 1)
my.metadata=read.csv(file = "metadata.csv",header = T, sep = ",",row.names = 1)

# factors
#my.metadata$time=factor(my.metadata$time) 
locality=factor(my.metadata$locality)

# I changed source to Source, as source() seems to be an R function
Source = factor(my.metadata$source, levels = c("Leaf", "Branch", "Dust"))
my.metadata$time<-factor(my.metadata$time, levels = c("1","2","3"))
#dust=(my.metadata$dust)

#### Richness (Species number) model
# I removed the zero-observation samples.
Richness = apply(fungal.abundance,1,function(vec) sum(vec>0))
hist(log(Richness))

# We need to remove the zero-samples.
RichNotZero = Richness > 0
fungl.abun.not.zero=fungal.abundance[RichNotZero,]

# Richness with samples that have observations
RichnessNotZero = specnumber(fungl.abun.not.zero)
hist(RichnessNotZero)

# This removes the samples with one observation
MetaRich = my.metadata[RichNotZero,]
MetaRich$time<-factor(MetaRich$time, levels = c("1","2","3"))

## fitting the models for Richness

Richness.m1=lm(RichnessNotZero~time+locality+source*dust, data = MetaRich)
Richness.m2=glm(RichnessNotZero~time+locality+source+dust, data = MetaRich, 
                family=poisson(link = "log"))
Richness.m3= glm(RichnessNotZero~locality+time+source*dust, data = MetaRich, 
                 family=poisson(link = "log"))
 
# model diagnostic plots
# Change plotting parameters, so all diagnostic plots are on the same sheet
par(mfrow=c(2,2))
plot(Richness.m1)
plot(Richness.m2)
plot(Richness.m3)
# summary(Richness.m1)
# summary(Richness.m2)
# summary(Richness.m3)
# Look for the AICs: Akaike Information Criteria and model selection
AIC(Richness.m1)
AIC(Richness.m2)
AIC(Richness.m3)
# The AIC says that the simple linear model (m1) is the best fit 
# (although still crappy according to the diagnostic plot)

### Richness Model Evaluations
Richness.m1.anova=anova(Richness.m1, test="Chisq")
Richness.m1.summary=summary(Richness.m1)

# What does it mean? With predicted values
par(mfrow=c(1,1))
boxplot(fitted(Richness.m1) ~ MetaRich$source, xlab="Source", ylab="Richness")
# Richness in branch and dust similar, but richness stat. sign. lower in leaf

# Interaction plots:
Richness.effect = effect("source:dust",Richness.m1,multiline=TRUE,  ylim=c(-10,10))
Rich.eff.sum = summary(Richness.effect)

par(mfrow=c(1,3), mar = c(5,3,2,1))
for (i in levels(MetaOne$source)) {
  plot(c(min(MetaOne$dust), max(MetaOne$dust)),
       c(min(Rich.eff.sum$lower[i,]), max(Rich.eff.sum$upper[i,])), 
       type="n", xlab = paste(i), ylab = "Richness")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(Rich.eff.sum$effect[i,]), max(Rich.eff.sum$effect[i,])))
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(Rich.eff.sum$lower[i,]), max(Rich.eff.sum$lower[i,])), 
        lty="dashed")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(Rich.eff.sum$upper[i,]), max(Rich.eff.sum$upper[i,])),
        lty="dashed")}

#### Shannon and Simpson models.
# Keep only samples with at least two OTUs
RichNotOne = Richness > 1
fungl.abun.not.one=fungal.abundance[RichNotOne,]

# Richness with samples that have observations

# This keeps observations with at least two OTUs
MetaOne = my.metadata[RichNotOne,]
MetaOne$time<-factor(MetaOne$time, levels = c("1","2","3"))

# Calculate diversity indices
shannon = diversity(fungl.abun.not.one,index = "shannon")
simpson = diversity(fungl.abun.not.one,index = "simpson")
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
shannon.m1=lm(shannon ~ locality+time +source+dust, data=MetaOne)
shannon.m2=glm(shannon~time+locality+source+dust, data=MetaOne, family=Gamma(link="log"))
shannon.m3=glm(shannon~locality+time+source*dust, data=MetaOne, family=Gamma(link="log"))

# plot all diagnostics
par(mfrow=c(2,2))
plot(shannon.m1)
plot(shannon.m2)
plot(shannon.m3)

# AICs - the third is the best fit + it has the best diagnostic plots
AIC(shannon.m1)
AIC(shannon.m2)
AIC(shannon.m3)

# These are probably not necessary, as these models are worse fits
# summary(shannon.m1) 
# anova(shannon.m1, test="Chisq")
# summary(shannon.m2) 
# anova(shannon.m2, test = "Chisq")

# Best Shannon model statistics
shannon.summary=summary(shannon.m3) 
shannon.anova=anova(shannon.m3, test = "Chisq")

# Interaction plots:
shannon.effect = effect("source:dust",shannon.m3,multiline=TRUE,  ylim=c(-10,10))
shan.eff.sum = summary(shannon.effect)

par(mfrow=c(1,3), mar = c(5,3,2,1))
for (i in levels(MetaOne$source)) {
  plot(c(min(MetaOne$dust), max(MetaOne$dust)),
       c(min(shan.eff.sum$lower[i,]), max(shan.eff.sum$upper[i,])), 
       type="n", xlab = paste(i), ylab = "Shannon diversity")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(shan.eff.sum$effect[i,]), max(shan.eff.sum$effect[i,])))
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(shan.eff.sum$lower[i,]), max(shan.eff.sum$lower[i,])), 
        lty="dashed")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(shan.eff.sum$upper[i,]), max(shan.eff.sum$upper[i,])),
        lty="dashed")}
### I tried changing the order of the factors in the model and it doesn't really change anything. the model m3 is
###the best model describing our data

## Simpson models:
simpson.m1=lm(simpson ~ time + locality +source+dust, data=MetaOne)
simpson.m2=glm(simpson~time+locality+source+dust, data=MetaOne,family=Gamma(link = "log"))
simpson.m3=glm(simpson~locality+time+source*dust, data = MetaOne,family = Gamma(link="log"))

# Diagnostics
plot(simpson.m1)
plot(simpson.m2)
plot(simpson.m3)

# AICs - again everything says the m3 is the best.
AIC(simpson.m1)
AIC(simpson.m2)
AIC(simpson.m3)

# Model statistics
# summary(simpson.m1)
# anova(simpson.m1, test= "Chisq")
# summary(simpson.m2)
# anova(simpson.m2, test= "Chisq")

simpson.summary=summary(simpson.m3)
simpson.anova=anova(simpson.m3, test= "Chisq")

# Interaction plots:
simpson.effect = effect("source:dust",simpson.m3,multiline=TRUE,  ylim=c(-10,10))
Simp.eff.sum = summary(simpson.effect)

par(mfrow=c(1,3), mar = c(5,3,2,1))
for (i in levels(MetaOne$source)) {
  plot(c(min(MetaOne$dust), max(MetaOne$dust)),
       c(min(Simp.eff.sum$lower[i,]), max(Simp.eff.sum$upper[i,])), 
       type="n", xlab = paste(i), ylab = "Simpson diversity")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(Simp.eff.sum$effect[i,]), max(Simp.eff.sum$effect[i,])))
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(Simp.eff.sum$lower[i,]), max(Simp.eff.sum$lower[i,])), 
        lty="dashed")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(Simp.eff.sum$upper[i,]), max(Simp.eff.sum$upper[i,])),
        lty="dashed")}
# Actually, let's plot the model predictions.
#plot(MetaOne$dust, fitted(shannon.m3), xlab="Dust concentration")
par(mfrow =c(1,3))
boxplot(fitted(Richness.m1) ~ MetaRich$source, xlab="Source of isolation", ylab="Richness")
## why the ylab doesn't work here?
boxplot(fitted(shannon.m3) ~ MetaOne$source, xlab="Shannon")
boxplot(fitted(simpson.m3) ~ MetaOne$source, xlab="Simpson")

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
##### I didn't get this ERROR!!!
# Error message: MCMC fitting through JAGS failed. This is likely due to the
# prior on the dispersion (size) parameter of the negative binomial distribution
# been too uninformative (see below). Please consider a tougher prior or switch
# to a Poisson family for those response that don't appear to actually be
# overdispersed. [1] "Error : Error in node all.params[7,4]\nSlicer stuck at
# value with infinite density\n\n\n" attr(,"class") [1] "try-error" 
# attr(,"condition") <simpleError: Error in node all.params[7,4] Slicer stuck at
# value with infinite density

# I think this failed because the data is extremely overdispersed
# plot(0.1,0.1, type="n", xlim=c(0.1,5), ylim=c(0.1, max(apply(fung.abund.reduced,2,var))), xlab="Mean", ylab="Variance") # ann=FALSE, 
# points(apply(fung.abund.reduced,2,mean), apply(fung.abund.reduced,2,var), pch=20)
# title(main="Untransformed data")

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

### Why this plot is different from NMDS?
# Because the NMDS plot likely has overdispersion effects.
dev.off()
par(mar=c(4,4,1,1))
plot(ModelOrd$lv.median, col=as.numeric(MetaOrd$source), 
            pch=19, main="Latent variable model", las=1.5)

# The NMDS will not be necessary anymore. If you run the two lines below, 
# you will see how locality and dispersal are mixed up because of the overdispersion.
# # This is very strong for the green samples (leaves?).
# NMDS1<-metaMDS(fung.abun.reduced)
# plot(NMDS1$points, ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2", 
#      col=as.numeric(MetaOrd$source), pch=19)
# 
# 


### 3- Community composition models

png(file="mean-variance-plot.jpeg", units="mm", height=90, width=90, 
    pointsize=10, bg="white", res=1200)
par(mar = c(4,4,1,1))
plot(0.1,0.1, type="n", xlim=c(0.1,100), 
     ylim=c(0.1, max(apply(Corfun0,2,var))),
     xlab="Mean (log scale)", ylab="Variance (log scale)", log="xy", xaxt="n", yaxt="n")
for (i in 1:length(colnames(Corfun0))) {
  points(mean(Corfun0[,i]), var(Corfun0[,i]), pch=20, cex=0.7)}
axis(1, at=c(0.1, 1,5,10), labels=c(0,1,5,10))
ticks.2 = seq(1,9, by=2)
labels.2 <- sapply(ticks.2, function(i) as.expression(bquote(10^ .(i))))
axis(2, at=c(0.1, 10, 1000, 100000, 10000000), 
     labels=labels.2)


funcor.mvabund = mvabund(Corfun0)
plot(funcor.mvabund)

## Model selection
### if I remember correctly you said that I don't have to do this part and I should just go with the most complicated 
## model from diversity models

funcor.m1 = manyglm(funcor.mvabund ~ locality+time+source+dust, data= metacor0,
                 family="negative.binomial", show.residuals=T)
summary(funcor.m1)

funcor.m2 = manyglm(funcor.mvabund ~ locality+time+source*dust, data= metacor0,
                    family="negative.binomial", show.residuals=T)

funcor.m3 = manyglm(funcor.mvabund ~ locality*time+source*dust , data= metacor0, 
                    family="negative.binomial", show.residuals=T) 
summary(funcor.m3)
anova(funcor.m2, funcor.m3, nBoot=50)

funcor.m4= manyglm(funcor.mvabund ~ locality*dust+time+source*dust , data= metacor0, 
                   family="negative.binomial", show.residuals=T)
summary(funcor.m4)
anova(funcor.m3,funcor.m4, nBoot = 50)


### model 3 is the better option comparing to other models. but I can't find any explanation for the 
## interaction effect of time and locality! aparently somthing happend in those times in those localities but 
### we don't know what!!
### so I want to continue with model 2
## factors:
metacor0$time=factor(metacor0$time, levels = c("1","2","3"))
metacor0$source = factor(metacor0$source, levels = c("Leaf", "Branch", "Dust"))
metacor0$locality=factor(metacor0$locality)

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



