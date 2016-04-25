# This is a change

## Pakages we need for our project
library(mvabund)
library(vegan)
library(mgcv)

# Data input

fungal.abundance= read.csv(file = "morphotype_matrix_incubator.csv", header = T, sep=",", row.names = 1)
my.metadata=read.csv(file = "metadata.csv",header = T, sep = ",",row.names = 1)

# factors
time=factor(my.metadata$time) 
locality=factor(my.metadata$locality)
source= my.metadata$source <- factor(my.metadata$source, levels = c("Leaf", "Branch", "Dust"))
dust=factor(my.metadata$dust)

# Isolation success
Success=apply(fungal.abundance,1,sum)


# HILL Diversities and Shannon diversity index
Specieshill=renyi(fungal.abundance,scales=c(0,1,2),hill=T)
myhill.1=Specieshill$"0"
myhill.2=Specieshill$"1"
myhill.3=Specieshill$"2"
shannon= diversity(fungal.abundance,index = "shannon",MARGIN = 1,base = exp(1))

#First hill 
myhill.1.m1d=lm(myhill.1~dust)
myhill.1.m1=lm(myhill.1~source)
myhill.1.m2=lm(myhill.1~time+source)
myhill.1.m3=lm(myhill.1~time+locality+source)
myhill.1.m4=lm(myhill.1~time+locality+source+dust)
myhill.1.m5=lm(myhill.1~locality+time+source*dust)

anova.hill.1=anova(myhill.1.m1d,myhill.1.m1,myhill.1.m2,myhill.1.m3,myhill.1.m4,myhill.1.m5)
#model m2 is better?
summary(myhill.1.m2)

#Second Hill
myhill.2.m1d=lm(myhill.2~dust)
myhill.2.m1=lm(myhill.2~source)
myhill.2.m2=lm(myhill.2~time+source)
myhill.2.m3=lm(myhill.2~time+locality+source)
myhill.2.m4=lm(myhill.2~time+locality+source+dust)
myhill.2.m5=lm(myhill.2~locality+time+source*dust)

anova.hill.2=anova(myhill.2.m1d,myhill.2.m1,myhill.2.m2,myhill.2.m3,myhill.2.m4,myhill.2.m5)
summary(myhill.2.m2)

#Third hill
# an error is accuring in this part...What is it?
myhill.3.m1d=lm(myhill.3~dust)
myhill.3.m1=lm(myhill.3~source)
myhill.3.m2=lm(myhill.3~time+source)
myhill.3.m3=lm(myhill.3~time+locality+source)
myhill.3.m4=lm(myhill.3~time+locality+source+dust)
myhill.3.m5=lm(myhill.3~locality+time+source*dust)

##Shannon diversity
myshannon.m1d=lm(shannon~dust)
myshannon.m1=lm(shannon~source)
myshannon.m2=lm(shannon~time+source)
myshannon.m3=lm(shannon~time+locality+source)
myshannon.m4=lm(shannon~time+locality+source+dust)
myshannon.m5=lm(shannon~locality+time+source*dust)

Anova.shannon=anova(myshannon.m1d,myshannon.m1,myshannon.m2,myshannon.m3,myshannon.m4,myshannon.m5)


## Show the partial residuals of Hill's numbers 
## after accounting for sequence number differences (partial residuals)
# I need some explanation about this residuals. i am not sure that i am understanding it correctly...
##### why did you do the following only for model1 of every hill?

png(file="diversities.png", units="mm", height=70, width=90, pointsize=10, bg="white", res=1200)
par(mfrow=c(1,3))
par(mar = c(7,2.2,2,0.5))

boxplot(myhill.1.m1$residuals ~ source, boxwex=0.5, notch=T,xlim=c(0.5,3.5), ylab=NA, main="Hill's N0", xaxt="n", 
        lwd=0.7)
axis(side=1, at=c(1,2,3), lwd=0,labels=c("leaf", "branch", "dust"), las=2)

boxplot(myhill.2.m1$residuals ~ source, boxwex=0.5, notch=F, xlim=c(0.5,3.5), ylab=NA, main="Hill's N1", xaxt="n",
        lwd=0.7,col="lightgrey")
axis(side=1, at=c(1,2,3), lwd=0,labels=c("leaf", "branch", "dust"), las=2)

boxplot(myhill.3.m1$residuals ~ source, boxwex=0.5, notch=T, xlim=c(0.5,3.5), ylab=NA, main="Hill's N2", xaxt="n", 
        lwd=0.7,col="darkgrey")
axis(side=1, at=c(1,2,3), lwd=0,labels=c("leaf", "branch", "dust"), las=2)

## Post-hoc Tukey tests among the three experimental treatments with partial residuals, after accounting for 
## differential sequencing success

## Hill N0
TukeyHSD(aov(myhill.1.m1$residuals ~ source))

## Hill N1
TukeyHSD(aov(myhill.2.m1$residuals ~ source))

## Hill N2
TukeyHSD(aov(myhill.3.m1$residuals ~ source))

##### 3. Define core OTUs

## species abundance
abundance=apply(fungal.abundance,2,sum)

## The average read number of OTUs
Meanabundance=apply(fungal.abundance,2,function(vec) mean(vec[vec>0]))

## In how many samples is an OTU present? 
## I think in this part we should count the samples instead of sum
sample.present = apply(fungal.abundance,2,function(vec) sum(vec>0))

## The highest read number of an OTU in a sample
Maxobservation=apply(fungal.abundance,2,max)

## Plotting incidence against abundance
plot(sample.present, Maxobservation, xlab="Incidence",
     ylab="Maximum Abundance", pch=20)

plot(sample.present, log(Maxobservation), xlab="Incidence",
     ylab="log(Maximum Abundance)", pch=20)

## Create a smoothed trendline
mygam1 = gam(log(Maxobservation)~s(sample.present))

plot(mygam1, residuals=T, shade=T, rug=F, cex=2.6,
     xlab="Incidence", ylab="logMean Abundance") # , xaxp=c(0,150,15)

## keep core OTUs
##### Didn't understand this part...
## how do you know which one is a core OTU? where did this no 30 come from?
IsFreq = sample.present > 30
fun.some = fungal.abundance[,IsFreq]

## How many of these in the South, in the North Natural, and North Heated?
## can I do this for each source and each locality or time?
some.South = fun.some[experiment == "ControlWarm",]
length(colnames(some.South)[apply(some.South,2,sum) > 0])

some.NorthNatural = fun.some[experiment == "NotHeated",]
length(colnames(some.NorthNatural)[apply(some.NorthNatural,2,sum) > 0])

some.NorthHeated = fun.some[experiment == "Heated",]
length(colnames(some.NorthHeated)[apply(some.NorthHeated,2,sum) > 0])







