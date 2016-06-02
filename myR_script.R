
## Pakages we need for our project
library(mvabund)
library(vegan)
library(mgcv)

# Data input

fungal.abundance= read.csv(file = "morphotype_matrix_incubator.csv", header = T, sep=",", row.names = 1)
my.metadata=read.csv(file = "metadata.csv",header = T, sep = ",",row.names = 1)

# factors
#my.metadata$time=factor(my.metadata$time) 
locality=factor(my.metadata$locality)
source = factor(my.metadata$source, levels = c("Leaf", "Branch", "Dust"))
my.metadata$time<-factor(my.metadata$time, levels = c("1","2","3"))
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

# The model diagnstic plots show that the third model has a much better fit - much less trend in the residuals

# Look for the AICs: Akaike Information Criteria and model selection
AIC(Richness.m1)
AIC(Richness.m2)
AIC(Richness.m3)
# The Poisson-GLM with the source*dust interactions looks like a reliable one, and it has a much better fit to the data according to AIC. 

# Evaluations
Richness.m3.anova=anova(Richness.m3, test="Chisq")
Richness.m3.summary=summary(Richness.m3)

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
MetaObs$time<-factor(MetaObs$time, levels = c("1","2","3"))



shannon= diversity(fungl.abun2,index = "shannon")
simpson=diversity(fungl.abun2,index = "simpson")
hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

## Fishers alpha log-series
Fisheralpha= fisher.alpha(fungl.abun2)

## Evenness
Evenness= shannon/log(Richness)
plot(Evenness)
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
shannon.summary=summary(shannon.m3) 
shannon.anova=anova(shannon.m3, test = "Chisq")


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
simpson.summary=summary(simpson.m3)
simpson.anova=anova(simpson.m3, test= "Chisq")

AIC(simpson.m1)
AIC(simpson.m2)
AIC(simpson.m3)


#### Plot effects
plot(MetaObs$dust, shannon)
boxplot(shannon ~ MetaObs$source)

#### Fishers alpha log-series 
plot(Fisheralpha)
Fisheralpha.m=glm(Fisheralpha~time+source*dust+locality, data= MetaObs)
Fisher.summary=summary(Fisheralpha.m)
fisher.anova=anova(Fisheralpha.m)

###Final models:
### I changed the order of the factors in the model and it did not really change the results and these three
### models are final
Richness.m3= glm(Richness~locality+ time +source*dust, data = my.metadata, family=poisson(link = "log"))
shannon.m3=glm(shannon~locality+time+source*dust, data=MetaObs, family=Gamma(link="log"))
simpson.m3=glm(simpson~time+locality+source*dust, data = MetaObs,family = Gamma(link="log"))

####source*dust interactions plot
library(effects)
plot(effect("source:dust", Richness.m3, multiline=TRUE))

plot(effect("source:dust",shannon.m3,multiline=TRUE))

plot(effect("source:dust",simpson.m3,multiline=TRUE))

##### Species Accumulation Curves
acum=specaccum(fungal.abundance,method = "exact", permutations = 100,
          conditioned =TRUE, gamma = "jack1",  w = NULL)
plot(acum)
plot(acum, add = FALSE, random = FALSE, ci = 0, 
     ci.type = c("line"), col = "black",xlab = "Number of samples" , ylab = "Richness"
     , xvar = c("sites", "individuals", "effort"),ylim )

## Fit Lomolino model to the exact accumulation
acumodel=fitspecaccum(acum, "lomolino")
coef(acumodel)
fitted(acumodel)
plot(acum,random = FALSE, ci = 0,ci.type = c("line"),xlab = "Number of samples" , ylab = "Richness")

## Add Lomolino model using argument 'add'
plot(acumodel, add = TRUE, col=2, lwd=2)

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




### Visualize differences in community composition
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
plot(NMDS1$points, type="n", ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2")
plot(NMDS1)

# show overlapping samples
ordiplot(NMDS1, type="text")
ordispider(NMDS1, metacor0$source , col="grey")




ordispider(NMDS1, metacor0$source , col="grey")
points(NMDS1, pch=20, cex=exp(2*shannon)/100, col="grey")
mylegend = legend(2, 1.5, c("leaf","branch","dust"), 
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

## Three axes
NMDS.3 <- metaMDS(Corfun0, k=3, trymax=100)
NMDS.3 <- metaMDS(Corfun0, previous = NMDS.3, k=3, trymax=100)


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

## dynamic 3D ordination plot

library(vegan3d)
ordirgl(NMDS.3)
orglspider(NMDS.3, metacor0$source)
orgltext(NMDS.3, rownames(Corfun0))
orgltext(NMDS.3, colnames(Corfun0))

### NMDS plot for localities: 
par(mar=c(4,4,1,1))
plot(NMDS1$points, type="n", ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2")
ordispider(NMDS1, metacor0$locality , col="grey")
points(NMDS1, pch=20, cex=exp(2*shannon)/100, col="grey")
mylegend2= legend(2, 1.5, c("BISOTON","HASAN ABAD","KHOSRO ABAD", "KEREND", "SORKHE DIZE"), 
                  fill=c("green","orange","yellow","red","blue"), border="white", bty="n")
with(metacor0,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
                          draw="polygon", col=c("green"),
                          alpha=100,kind="se",conf=0.95, 
                          show.groups=(c("BISOTON"))))
with(metacor0,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
                          draw="polygon", col=c("orange"),
                          alpha=100,kind="se",conf=0.95, 
                          show.groups=(c("HASAN ABAD"))))
with(my.metadata,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
                             draw="polygon", col=c("yellow"),
                             alpha=100,kind="se",conf=0.95, 
                             show.groups=(c("KHOSRO ABAD"))))
with(my.metadata,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
                             draw="polygon", col=c("red"),
                             alpha=100,kind="se",conf=0.95, 
                             show.groups=(c("KEREND"))))
with(my.metadata,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
                             draw="polygon", col=c("blue"),
                             alpha=100,kind="se",conf=0.95, 
                             show.groups=(c("SORKHE DIZE"))))

### NMDS plot for time of sampling
par(mar=c(4,4,1,1))
plot(NMDS1$points, type="n", ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2")
ordispider(NMDS1, metacor0$time , col="grey")
points(NMDS1, pch=20, cex=exp(2*shannon)/100, col="grey")

mylegend3= legend(2, 1.5, c("First sampling","Second sampling","Third sampling"), 
                  fill=c("green","orange","blue"), border="white", bty="n")
with(metacor0,ordiellipse(NMDS1, metacor0$time,cex=.5, 
                          draw="polygon", col=c("green"),
                          alpha=100,kind="se",conf=0.95, 
                          show.groups=(c("1"))))
with(metacor0,ordiellipse(NMDS1, metacor0$time,cex=.5, 
                          draw="polygon", col=c("orange"),
                          alpha=100,kind="se",conf=0.95, 
                          show.groups=(c("2"))))
with(my.metadata,ordiellipse(NMDS1, metacor0$time,cex=.5, 
                             draw="polygon", col=c("blue"),
                             alpha=100,kind="se",conf=0.95, 
                             show.groups=(c("3"))))




### Community composition models

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
### Do a simple glm for the two species affected by the source* dust







### similarity analysis using anosim and adonis functions (corfun0=fungal abundance matrix for core OTUs and metacor0 
## is the metadata)
braysimi= anosim(Corfun0,metacor0$source,permutations = 999,distance = "bray")

jaccardsimi= anosim(Corfun0,metacor0$source, permutations = 999, distance = "jaccard")


corfunAdon= adonis(Corfun0~locality+time+source*dust,data= metacor0,permutations = 999, method = "bray")






