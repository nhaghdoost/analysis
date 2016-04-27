
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
my.metadata$source <- factor(my.metadata$source, levels = c("Leaf", "Branch", "Dust"))
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

abun0=Richness>0
fun.abun0=fungal.abundance[abun0,]
Meta0=my.metadata[abun0,]
row.names(fun.abun0)==row.names(Meta0)

# This removes the samples with one observation
MetaObs = my.metadata[speciesnum0,]


shannon= diversity(fungl.abun2,index = "shannon")
simpson=diversity(fungl.abun2,index = "simpson")
hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

## Fisher’s a log-series and evenness is also in my proposal...



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
Corfun= fun.abun0[,funFreq]
length(Corfun)
Corname = colnames(Corfun)


## Core OTUs in leaf, branch and dust?
## Core OTUs in each locality?
## Core OTUs in each sampling time?

### 4. Visualize differences in community composition
## run NMDS
MDS.all <- metaMDS(Corfun0)
MDS.all <- metaMDS(Corfun0, previous = MDS.all)

NMDS1=metaMDS(Corfun0,k=2)

### Trying to fix the error!!!
Corfun0=Corfun[abun0,]
csum<-colSums(Corfun0)
any(is.na(csum))
which(is.na(csum))
#### aparently I have NA in my data which i have no idea where did it come from!!!
corfun00=na.omit(Corfun0)
NMDS2=metaMDS(corfun00,k=2)
MDS.all <- metaMDS(corfun00)
MDS.all <- metaMDS(corfun00, previous = MDS.all)
### it didi not solve the problem!!!
