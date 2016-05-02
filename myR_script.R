
## Pakages we need for our project
library(mvabund)
library(vegan)
library(mgcv)

# Data input

fungal.abundance= read.csv(file = "morphotype_matrix_incubator.csv", header = T, sep=",", row.names = 1)
my.metadata=read.csv(file = "metadata.csv",header = T, sep = ",",row.names = 1)

# factors
#my.metadata$time=factor(my.metadata$time) 
#locality=factor(my.metadata$locality)
source = factor(my.metadata$source, levels = c("Leaf", "Branch", "Dust"))
#dust=(my.metadata$dust)

# Richness (Species number) model --  the zero samples were included in richness analysis
Richness = apply(fungal.abundance,1,function(vec) sum(vec>0))
hist(log(Richness))

Richness.m1=lm(Richness~time+locality+source*dust, data = my.metadata)
plot(Richness.m1) 
# This is a model diagnostic plot.
summary(Richness.m1)

Richness.m2=glm(Richness~time+locality+source+dust, data = my.metadata, family=poisson(link = "log"))
plot(Richness.m2)
summary(Richness.m2)

Richness.m3= glm(Richness~locality+ time +source*dust, data = my.metadata, family=poisson(link = "log"))
par(mfrow=c(2,2))
plot(Richness.m3)
summary(Richness.m3)
#when I moved the locality to the first the effect of time is less significant but the numbers are the same.
# The model diagnstic plots show that the third model has a much better fit - much less trend in the residuals

# Look for the AICs: Akaike Information Criteria and model selection
AIC(Richness.m1)
AIC(Richness.m2)
AIC(Richness.m3)
# The Poisson-GLM with the source*dust interactions looks like a reliable one, and it has a much better fit to the data according to AIC. 

# Evaluations
anova(Richness.m3, test="Chisq")
summary(Richness.m3)

# What does it mean?
boxplot(Richness ~ my.metadata$source, xlab="Source", ylab="Richness")
boxplot(Richness ~ my.metadata$time, xlab="Source", ylab="Richness")


#### Shannon and Simpson models. We need to remove the zero-samples.

speciesnum0 = Richness > 1
fungl.abun2=fungal.abundance[speciesnum0,]

#abun0=Richness>0
#fun.abun0=fungal.abundance[abun0,]
#Meta0=my.metadata[abun0,]
#row.names(fun.abun0)==row.names(Meta0)

# This removes the samples with one observation
MetaObs = my.metadata[speciesnum0,]


shannon= diversity(fungl.abun2,index = "shannon")
simpson=diversity(fungl.abun2,index = "simpson")
hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

## Fishers a log-series and evenness is also in my proposal...



## Shannon models:

shannon.m1=lm(shannon ~ locality+time +source+dust, data=MetaObs)
plot(shannon.m1)

summary(shannon.m1)
AIC(shannon.m1)
anova(shannon.m1, test="Chisq")


shannon.m2=glm(shannon~time+locality+source+dust, data=MetaObs, family=Gamma(link="log"))
plot(shannon.m2)

summary(shannon.m2) 
anova(shannon.m2, test = "Chisq")

shannon.m3=glm(shannon~locality+time+source*dust, data=MetaObs, family=Gamma(link="log"))
par(mfrow=c(2,2))
plot(shannon.m3)
summary(shannon.m3) 
anova(shannon.m3, test = "Chisq")


AIC(shannon.m1)
AIC(shannon.m2)
AIC(shannon.m3)
### I tried changing the order of the factors in the model and it doesn't really change anything. the model m3 is
###the best model describing our data


## Simpson models:

simpson.m1=lm(simpson ~ time + locality +source+dust, data=MetaObs)
plot(simpson.m1)
summary(simpson.m1)
anova(simpson.m1, test= "Chisq")



simpson.m2=glm(simpson~time+locality+source+dust, data=MetaObs,family=Gamma(link = "log"))
plot(simpson.m2)
summary(simpson.m2)
anova(simpson.m2, test= "Chisq")


simpson.m3=glm(simpson~locality+time+source*dust, data = MetaObs,family = Gamma(link="log"))
par(mfrow=c(2,2))
plot(simpson.m3)
summary(simpson.m3)
anova(simpson.m3, test= "Chisq")

AIC(simpson.m1)
AIC(simpson.m2)
AIC(simpson.m3)


#### Plot effects
plot(MetaObs$dust, shannon)
boxplot(shannon ~ MetaObs$source)

###Final models:
### I changed the order of the factors in the model and it did not really change the results and these three
### models are final
Richness.m3= glm(Richness~locality+ time +source*dust, data = my.metadata, family=poisson(link = "log"))
shannon.m3=glm(shannon~locality+time+source*dust, data=MetaObs, family=Gamma(link="log"))
simpson.m3=glm(simpson~time+locality+source*dust, data = MetaObs,family = Gamma(link="log"))

####source*dust interactions plot
library(effects)
plot(effect("source:dust",Richness.m3,multiline=TRUE))

plot(effect("source:dust",shannon.m3,multiline=TRUE))

plot(effect("source:dust",simpson.m3,multiline=TRUE))

#### Define core species:
## Summarize reads
funTotCount = apply(fun.abun0,2,sum)

## The average read number of OTUs
funMeanCount=apply(fun.abun0,2,function(vec) mean(vec[vec>0]))

## In how many samples is an OTU present?
funTotPresent = apply(fun.abun0,2,function(vec) sum(vec>0))

## The highest read number of an OTU in a sample
funMaxCount=apply(fun.abun0,2,max)

## Plotting incidence against abundance
plot(funTotPresent,funMaxCount, xlab="Incidence",
     ylab="Maximum Abundance", pch=20)

plot(funTotPresent, log(funMaxCount), xlab="Incidence",
     ylab="log(Maximum Abundance)", pch=20)

## Create a smoothed trendline
funGam1 = gam(log(funMaxCount)~s(funTotPresent))

plot(funGam1, residuals=T, shade=T, rug=F, cex=2.6,
     xlab="Incidence", ylab="logMean Abundance") # , xaxp=c(0,150,15)

## keep core OTUs

funFreq = funTotPresent > 15
Corfun= fungal.abundance[,funFreq]
length(Corfun)
Corname = colnames(Corfun)
corsamplename=row.names(Corfun)

## Core OTUs in leaf, branch and dust?
## Core OTUs in each locality?
## Core OTUs in each sampling time?

### 4. Visualize differences in community composition
## run NMDS
#MDS.all <- metaMDS(Corfun)
#MDS.all <- metaMDS(Corfun, previous = MDS.all)
#NMDS1=metaMDS(Corfun,k=2)

### Trying to fix the error!!!
#Corfun0=Corfun[abun0,]
#csum<-colSums(Corfun0)
#any(is.na(csum))
#which(is.na(csum))

### try to delete the rows with zero
rowCount=apply(Corfun, 1,FUN=sum)
nozero=rowCount>0
Corfun0=Corfun[nozero,]
colnames(Corfun0)
row.names(Corfun0)
metacor0=my.metadata[nozero,]
row.names(Corfun0)==row.names(metacor0)
### try to do the NMDS with the new matrix

NMDS1<-metaMDS(Corfun0)
NMDS1<-metaMDS(Corfun0, previous= NMDS1)

### yes its working :) .

### ploting the NMDS:
## plot NMDS
par(mar=c(4,4,1,1))
plot(NMDS1$points, type="n", ylim=c(-0.9,0.9), xlab="NMDS1", ylab="NMDS2")
plot(NMDS1)

ordispider(NMDS1, metacor0$source , col="grey")
points(NMDS1, pch=20, cex=exp(2*shannon)/100, col="black")
mylegend = legend(-1, 0.95, c("leaf","branch","dust"), 
                  fill=c("green","orange","gray"), border="white", bty="n")
with(metacor0,ordiellipse(NMDS1, metacor0$source,cex=.5, 
                          draw="polygon", col=c("green"),
                          alpha=100,kind="se",conf=0.95, 
                          show.groups=(c("Leaf"))))
with(metacor0,ordiellipse(NMDS1, metacor0$source,cex=.5, 
                             draw="polygon", col=c("orange"),
                             alpha=100,kind="se",conf=0.95, 
                             show.groups=(c("Branch"))))
with(my.metadata,ordiellipse(NMDS1, metacor0$source,cex=.5, 
                             draw="polygon", col=c("gray"),
                             alpha=100,kind="se",conf=0.95, 
                             show.groups=(c("Dust"))))
## can I do this NMDS plots for my localities and times?

## Three axes
NMDS.3 <- metaMDS(Corfun0, k=3, trymax=100)
NMDS.3 <- metaMDS(Corfun0, previous = NMDS.3, k=3, trymax=100)


## colors with experiment
## did not understand the meaning of this part in your Rscript?

### Axes 1 & 2
### why did we use shannon???
png(file="community_NMDS_1-2.jpeg", units="mm", height=60, width=180, 
    pointsize=10, bg="white", res=1200)
par(mfrow=c(1,3))
par(mar=c(4,4,1,1))
NMDS.1.2 = ordiplot(NMDS.3, choices=c(1,2), type="n", 
                       xlab="NMDS1", ylab="NMDS2")
ordispider(NMDS.1.2,metacor0$source, col="grey")
points(NMDS.3$points[,1], NMDS.3$points[,2], pch=20, 
       cex=exp(2*shannon)/100,col="gray")
ordiellipse(NMDS.1.2, metacor0$source,cex=.5, 
            draw="polygon", col=c("green"),
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("Leaf")), border="green")
ordiellipse(NMDS.1.2, metacor0$source,cex=.5, 
            draw="polygon", col=c("orange"),
            alpha=100,kind="se",conf=0.95,
            show.groups=(c("Branch")), border="orange")
ordiellipse(NMDS.1.2, metacor0$source,cex=.5, 
            draw="polygon", col=c("gray"),
            alpha=50,kind="se",conf=0.95,
            show.groups=(c("Dust")), border="black")

### Axes 1 & 3

png(file="community_NMDS_1-3.jpeg", units="mm", height=60, width=180, 
    pointsize=10, bg="white", res=1200)
par(mfrow=c(1,3))
par(mar=c(4,4,1,1))
NMDS.1.3 = ordiplot(NMDS.3, choices=c(1,3), type="n", 
                    xlab="NMDS1", ylab="NMDS3")
ordispider(NMDS.1.3,metacor0$source, col="grey")
points(NMDS.3$points[,1], NMDS.3$points[,3], pch=20, 
       cex=exp(2*shannon)/100,col="gray")
ordiellipse(NMDS.1.3, metacor0$source,cex=.5, 
            draw="polygon", col=c("green"),
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("Leaf")), border="green")
ordiellipse(NMDS.1.3, metacor0$source,cex=.5, 
            draw="polygon", col=c("orange"),
            alpha=100,kind="se",conf=0.95,
            show.groups=(c("Branch")), border="orange")
ordiellipse(NMDS.1.3, metacor0$source,cex=.5, 
            draw="polygon", col=c("gray"),
            alpha=50,kind="se",conf=0.95,
            show.groups=(c("Dust")), border="black")

### Axes 2 & 3

png(file="community_NMDS_2-3.jpeg", units="mm", height=60, width=180, 
    pointsize=10, bg="white", res=1200)
par(mfrow=c(1,3))
par(mar=c(4,4,1,1))
NMDS.2.3 = ordiplot(NMDS.3, choices=c(2,3), type="n", 
                    xlab="NMDS2", ylab="NMDS3")
ordispider(NMDS.2.3,metacor0$source, col="grey")
points(NMDS.3$points[,2], NMDS.3$points[,3], pch=20, 
       cex=exp(2*shannon)/100,col="gray")
ordiellipse(NMDS.2.3, metacor0$source,cex=.5, 
            draw="polygon", col=c("green"),
            alpha=100,kind="se",conf=0.95, 
            show.groups=(c("Leaf")), border="green")
ordiellipse(NMDS.2.3, metacor0$source,cex=.5, 
            draw="polygon", col=c("orange"),
            alpha=100,kind="se",conf=0.95,
            show.groups=(c("Branch")), border="orange")
ordiellipse(NMDS.2.3, metacor0$source,cex=.5, 
            draw="polygon", col=c("gray"),
            alpha=50,kind="se",conf=0.95,
            show.groups=(c("Dust")), border="black")

## did not understand this part
## PCA
fun.pca = rda(Corfun0)
fun.pca.scores = scores(fun.pca, choices=c(1,2,3))
fun.pca.eigenvals = eigenvals(fun.pca)

## explained by the first three axes: 52%
(fun.pca.eigenvals[1] + fun.pca.eigenvals[2] + fun.pca.eigenvals[3])/sum(fun.pca.eigenvals)

## axis distributions
fun.pca1 = fun.pca.scores$sites[,1]
fun.pca2 = fun.pca.scores$sites[,2]
fun.pca3 = fun.pca.scores$sites[,3]

hist(fun.pca1^2)
hist(fun.pca2^2)
hist(fun.pca3^2)


### why did you do this glm thing here:? 
## GLMs on scores
summary(glm(log(fun.pca1^2) ~ time+locality+source*dust, data = metacor0))
par(mfrow=c(2,2))
plot(glm(log(fun.pca1^2) ~ time+locality+source*dust, data = metacor0))
plot(glm(fun.pca1^2 ~ time+locality+source*dust, data = metacor0))

plot(glm(log(fun.pca2^2) ~ time+locality+source*dust, data = metacor0))
plot(glm(fun.pca2^2 ~ time+locality+source*dust, data = metacor0))

plot(glm(log(fun.pca3^2) ~ time+locality+source*dust, data = metacor0))
plot(glm(fun.pca3^2 ~ time+locality+source*dust, data = metacor0))

## dynamic 3D ordination plot

library(vegan3d)
ordirgl(NMDS.3)
orglspider(NMDS.3, metacor0$source)
orgltext(NMDS.3, rownames(Corfun0))
orgltext(NMDS.3, colnames(Corfun0))

### 5. Community composition models
## format the data

## Mean-abundance relationship
## why only this OTU??
#mean(fun.some$OTU_460)
#var(fun.some$OTU_460)
#boxplot(fun.some$OTU_460)

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

## Model selection
### if I remember correctly you said that I don't have to do this part and I should just go with the most complicated 
## model from diversity models

funcor.m1 = manyglm(funcor.mvabund ~ locality+time+source+dust, data= metacor0,
                 family="negative.binomial", show.residuals=T)
summary(funcor.m1)

funcor.m2 = manyglm(funcor.mvabund ~ locality+time+source*dust, data= metacor0,
                    family="negative.binomial", show.residuals=T)
summary(funcor.m2)
anova(funcor.m1,funcor.m2, nBoot = 50)

funcor.m3 = manyglm(funcor.mvabund ~ locality*time+source*dust , data= metacor0, 
                    family="negative.binomial", show.residuals=T) 
summary(funcor.m3)
anova(funcor.m2, funcor.m3, nBoot=50)

funcor.m4= manyglm(funcor.mvabund ~ locality*dust+time+source*dust , data= metacor0, 
                   family="negative.binomial", show.residuals=T)
summary(funcor.m4)
anova(funcor.m3,funcor.m4, nBoot = 50)

#funcor.m5=  manyglm(funcor.mvabund ~ locality*dust+locality*time+time*dust+source*dust , data= metacor0, 
                    #family="negative.binomial", show.residuals=T)
#anova(funcor.m4,funcor.m5, nBoot = 50)

#funcor.m6=manyglm(funcor.mvabund ~ locality*dust+time*dust+source*dust , data= metacor0, 
                 # family="negative.binomial", show.residuals=T)
#anova(funcor.m5,funcor.m6, nBoot = 50)

### model 3 is the better option comparing to other models. but I can't find any explanation for the 
## interaction effect of time and locality! aparently somthing happend in those times in those localities but 
### we don't know what!!
### so I want to continue with model 2


## Selected model: funcor.m2
funcor.m2 = manyglm(funcor.mvabund ~ locality+time+source*dust, data= metacor0,
                    family="negative.binomial", show.residuals=T)
m2.summary = summary(funcor.m2, nBoot=50, test="LR", p.uni="adjusted",
                      resamp="montecarlo")

## Analysis of variance explained by the predictors
funcor.anova.m2 = anova(funcor.m2, nBoot=300, test="LR", p.uni="adjusted", 
                     resamp="montecarlo")

# ## OTUs significantly explained by the affected by the source?? (at p<=0.001???)
m2.p.anova <- as.data.frame(funcor.anova.m2$uni.p)
fun.m2.exp3 = colnames(m2.p.anova)[m2.p.anova[3,]<=0.05]
fun.m2.exp2 = colnames(m2.p.anova)[m2.p.anova[2,]<=0.05]
fun.m2.exp1 = colnames(m2.p.anova)[m2.p.anova[1,]<=0.05]

## Visualization of source*dus interactions 
## Coefficients
funcor.m2.coef = as.data.frame(funcor.m2$coefficients)

## How do I know if i have an outlier or not?


## mean-centering the contrasts
coef.mean.contrast = funcor.m2.coef - apply(funcor.m2.coef,2,mean)

### 6. Post-hoc test of community model predictions
funcor.predict = fitted(funcor.m2)
funcor.predict = cbind(funcor.predict, experiment=metacor0$source)

Tukey.corcom = vector("list")
for (i in fun.m2.exp3){Tukey.corcom[[i]] = TukeyHSD(aov(log(funcor.predict[,colnames(funcor.predict) == i]) ~ 
                                                          metacor0$source))}

## extract the significances and test values for each OTU
Tukey.funcor.p = data.frame()
for (i in 1:length(Tukey.corcom)){
  Tukey.funcor.p = rbind(Tukey.funcor.p, Tukey.corcom[[i]]$source[,4]
}
colnames(Tukey.funcor.p) = c("Leaf","Branch","Dust")
                          
rownames(Tukey.funcore.p) = fun.exp.mc

Tukey.core.test = data.frame()
for (i in 1:length(Tukey.common)){
  Tukey.core.test = rbind(Tukey.core.test, Tukey.common[[i]]$experiment[,1])
}
colnames(Tukey.core.test) = c("South - North heated",
                              "South - North natural",
                              "North natural - North heated")
rownames(Tukey.core.test) = fun.exp.mc

Tukey.core.lwr = data.frame()
for (i in 1:length(Tukey.common)){
  Tukey.core.lwr = rbind(Tukey.core.lwr, Tukey.common[[i]]$experiment[,2])
}
colnames(Tukey.core.lwr) = c("South - North heated",
                             "South - North natural",
                             "North natural - North heated")
rownames(Tukey.core.lwr) = fun.exp.mc

Tukey.core.upr = data.frame()
for (i in 1:length(Tukey.common)){
  Tukey.core.upr = rbind(Tukey.core.upr, Tukey.common[[i]]$experiment[,3])
}
colnames(Tukey.core.upr) = c("South - North heated",
                             "South - North natural",
                             "North natural - North heated")
rownames(Tukey.core.upr) = fun.exp.mc

Tukey.full=cbind(Diff.S.NH=Tukey.core.test[,1],
                 lwr.S.NH=Tukey.core.lwr[,1],
                 upr.S.NH=Tukey.core.upr[,1],
                 p.S.NH = Tukey.core.p[,1],
                 Diff.S.NN=Tukey.core.test[,2],
                 lwr.S.NN=Tukey.core.lwr[,2],
                 upr.S.NN=Tukey.core.upr[,2],
                 p.S.NN = Tukey.core.p[,2],
                 Diff.NN.NH=Tukey.core.test[,3],
                 lwr.NN.NH=Tukey.core.lwr[,3],
                 upr.NN.NH=Tukey.core.upr[,3],
                 p.NN.NH = Tukey.core.p[,3])
rownames(Tukey.full) = rownames(Tukey.core.p)

write.csv(file="Tukey.csv", Tukey.full)
























