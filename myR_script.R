
## Pakages we need for our project
library(mvabund)
library(vegan)
library(boral)
# library(mgcv)
library(effects)
# Data input

fungal.abundance= read.csv(file = "morphotype_matrix_incubator.csv", header = T, sep=",", row.names = 1)
my.metadata=read.csv(file = "metadata.csv",header = T, sep = ",",row.names = 1)

# factors
#my.metadata$time=factor(my.metadata$time) 
locality=factor(my.metadata$locality)

# I changed source to Source, as source() seems to be an R function
Source = factor(my.metadata$source, levels = c("Leaf", "Branch", "Dust"))
my.metadata$time<-factor(my.metadata$time, levels = c("1","2","3"))
#dust=(my.metadata$dust)

# Richness (Species number) model
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

# Change plotting parameters, so all diagnostic plots are on the same sheet
par(mfrow=c(2,2))

# Please see what order do you want for the variables, and then modify all models
# accordingly.
Richness.m1=lm(RichnessNotZero~time+locality+source*dust, data = MetaRich)
plot(Richness.m1) 
# This is a model diagnostic plot.
# summary(Richness.m1)

Richness.m2=glm(RichnessNotZero~time+locality+source+dust, data = MetaRich, 
                family=poisson(link = "log"))
plot(Richness.m2)
# summary(Richness.m2)

Richness.m3= glm(RichnessNotZero~locality+time+source*dust, data = MetaRich, 
                 family=poisson(link = "log"))
plot(Richness.m3) # this model is still quite crappy, but no idea what else to include.

# Look for the AICs: Akaike Information Criteria and model selection
AIC(Richness.m1)
AIC(Richness.m2)
AIC(Richness.m3)
# The AIC says that the simple linear model is the best fit 
# (although still crappy according to the diagnostic plot)

# Evaluations
Richness.m1.anova=anova(Richness.m1, test="Chisq")
Richness.m1.summary=summary(Richness.m1)

# What does it mean? With predicted values
par(mfrow=c(1,1))
boxplot(fitted(Richness.m1) ~ MetaRich$source, xlab="Source", ylab="Richness")
# Richness in branch and dust similar, but richness stat. sign. lower in leaf

# Put the interaction plots here:

#### Shannon and Simpson models.
#abun0=Richness>0
#fun.abun0=fungal.abundance[abun0,]
#Meta0=my.metadata[abun0,]
#row.names(fun.abun0)==row.names(Meta0)

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
hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

## Fishers alpha log-series
# Check what Fisher's alpha means. 
# You get very high values with only a few species observations.
Fisheralpha= fisher.alpha(fungl.abun.not.one)
hist(Fisheralpha)

## Evenness
RichnessMinTwo = Richness[Richness > 1]
Evenness= shannon/log(RichnessMinTwo)
hist(Evenness)

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

#### Plot effects
par(mfrow=c(1,2))

# Actually, let's plot the model predictions.
# plot(MetaOne$dust, fitted(shannon.m3), xlab="Dust concentration")
boxplot(fitted(shannon.m3) ~ MetaOne$source, xlab="Source of isolation")

#### Fishers alpha log-series: do we really need this?
hist(Fisheralpha)
Fisheralpha.m=glm(Fisheralpha~time+source*dust+locality, data= MetaOne)

par(mfrow=c(2,2))
plot(Fisheralpha.m)

Fisher.summary=summary(Fisheralpha.m)
fisher.anova=anova(Fisheralpha.m, test="Chisq")

###Final models:
### I changed the order of the factors in the model and it did not really change the results and these three
### models are final
# Richness.m3= glm(Richness~locality+ time +source*dust, data = my.metadata, family=poisson(link = "log"))
# shannon.m3=glm(shannon~locality+time+source*dust, data=MetaObs, family=Gamma(link="log"))
# simpson.m3=glm(simpson~locality+time+source*dust, data = MetaObs,family = Gamma(link="log"))

####source*dust interactions plot
# library(effects)
# plot(effect("source:dust", Richness.m3, multiline=TRUE, ylim=c(0,1)))

shannon.effect = effect("source:dust",shannon.m3,multiline=TRUE,  ylim=c(-10,10))
shan.eff.sum = summary(shannon.effect)

# Adapt this for richness and simpson
par(mfrow=c(1,3), mar = c(5,3,2,1))
for (i in levels(MetaOne$source)) {
  plot(c(min(MetaOne$dust), max(MetaOne$dust)),
       c(min(shan.eff.sum$lower[i,]), max(shan.eff.sum$upper[i,])), 
       type="n", xlab = paste(i), ylab = "")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(shan.eff.sum$effect[i,]), max(shan.eff.sum$effect[i,])))
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(shan.eff.sum$lower[i,]), max(shan.eff.sum$lower[i,])), 
        lty="dashed")
  lines(c(min(MetaOne$dust), max(MetaOne$dust)),
        c(min(shan.eff.sum$upper[i,]), max(shan.eff.sum$upper[i,])),
        lty="dashed")
}

# ##### Species Accumulation Curves
# acum=specaccum(fungal.abundance,method = "exact", permutations = 100,
#           conditioned =TRUE, gamma = "jack1",  w = NULL)
# plot(acum)
# plot(acum, add = FALSE, random = FALSE, ci = 0, 
#      ci.type = c("line"), col = "black",xlab = "Number of samples" , ylab = "Richness"
#      , xvar = c("sites", "individuals", "effort"),ylim )
# 
# ## Fit Lomolino model to the exact accumulation
# acumodel=fitspecaccum(acum, "lomolino")
# coef(acumodel)
# fitted(acumodel)
# plot(acum,random = FALSE, ci = 0,ci.type = c("line"),xlab = "Number of samples" , ylab = "Richness")
# 
# ## Add Lomolino model using argument 'add'
# plot(acumodel, add = TRUE, col=2, lwd=2)

# #### Define core species:
# ## Summarize reads
# funTotCount = apply(fun.abun0,2,sum)
# 
# ## The average read number of OTUs
# funMeanCount=apply(fun.abun0,2,function(vec) mean(vec[vec>0]))
# 
# ## In how many samples is an OTU present?
# funTotPresent = apply(fun.abun0,2,function(vec) sum(vec>0))
# 
# ## The highest read number of an OTU in a sample
# funMaxCount=apply(fun.abun0,2,max)
# 
# ## Plotting incidence against abundance
# plot(funTotPresent,funMaxCount, xlab="Incidence",
#      ylab="Maximum Abundance", pch=20)
# 
# plot(funTotPresent, log(funMaxCount), xlab="Incidence",
#      ylab="log(Maximum Abundance)", pch=20)
# 
# ## Create a smoothed trendline
# funGam1 = gam(log(funMaxCount)~s(funTotPresent))
# 
# plot(funGam1, residuals=T, shade=T, rug=F, cex=2.6,
#      xlab="Incidence", ylab="logMean Abundance") # , xaxp=c(0,150,15)
# 
# ## keep core OTUs
# 
# funFreq = funTotPresent > 15
# Corfun= fungal.abundance[,funFreq]
# length(Corfun)
# Corname = colnames(Corfun)
# corsamplename=row.names(Corfun)
# 
# ## Core OTUs in leaf, branch and dust?
# 
# 
# ## Core OTUs in each locality?
# 
# 
# ## Core OTUs in each sampling time?
# 
# 
# 
# 

# Model of community composition with mvabund

funcor.mvabund = mvabund(fung.abun.reduced)
plot(funcor.mvabund, MetaOrd$source)
## Model selection

# funcor.m1 = manyglm(funcor.mvabund ~ locality+time+source+dust, data= metacor0,
#                  family="negative.binomial", show.residuals=T)
# summary(funcor.m1)
# 
# funcor.m2 = manyglm(funcor.mvabund ~ locality+time+source*dust, data= metacor0,
#                     family="negative.binomial", show.residuals=T)

# funcor.m3 = manyglm(funcor.mvabund ~ time+source*dust+locality, data= metacor0, 
#                     family="negative.binomial", show.residuals=T) 
# summary(funcor.m3)
# anova(funcor.m2, funcor.m3, nBoot=50)
# 
# funcor.m4= manyglm(funcor.mvabund ~ locality*dust+time+source*dust , data= metacor0, 
#                    family="negative.binomial", show.residuals=T)
# summary(funcor.m4)
# anova(funcor.m3,funcor.m4, nBoot = 50)
# 
# # I kept only the model from the diversity comparisons, for comparability.
# funcor.m3 = manyglm(funcor.mvabund ~ time+source*dust+locality, data= metacor0, 
#                     family="negative.binomial", show.residuals=T) 

# Finalize the predictor order!!!

### model 3 is the better option comparing to other models. but I can't find any explanation for the 
## interaction effect of time and locality! aparently somthing happend in those times in those localities but 
### we don't know what!!
### so I want to continue with model 2
## factors:

# metacor0$time=factor(metacor0$time, levels = c("1","2","3"))
# metacor0$source = factor(metacor0$source, levels = c("Leaf", "Branch", "Dust"))
# metacor0$locality=factor(metacor0$locality)

## Selected model: funcor.m2
funcor.m = manyglm(funcor.mvabund ~ locality+time+source*dust, data=MetaOrd,
                    family="negative.binomial", show.residuals=T)

# Looks good.
plot.manyglm(funcor.m, which=c(1:3))

### ask about the error

# Run these once with nBoot = 1000
Model.summary= summary.manyglm(funcor.m, nBoot=50, test="LR",p.uni="adjusted")

## Analysis of variance explained by the predictors
funcor.anova.m = anova.manyglm(funcor.m, nBoot=50, test="LR", p.uni="adjusted", 
                               resamp="montecarlo")

# ## OTUs significantly affected by the source?? (at p<=0.001???)
m2.p.anova <- as.data.frame(funcor.anova.m$uni.p)

# fun.m2.source = colnames(m2.p.anova)[m2.p.anova["source",]<=0.05]
colnames(m2.p.anova)[m2.p.anova["source",]<=0.05]
colnames(m2.p.anova)[m2.p.anova["locality",]<=0.05]
colnames(m2.p.anova)[m2.p.anova["time",]<=0.05]
colnames(m2.p.anova)[m2.p.anova["source:dust",]<=0.05]

# ### Visualize differences in community composition
# ## run NMDS
# #MDS.all <- metaMDS(Corfun)
# #MDS.all <- metaMDS(Corfun, previous = MDS.all)
# #NMDS1=metaMDS(Corfun,k=2)
# 
# ### Trying to fix the error!!!
# #Corfun0=Corfun[abun0,]
# #csum<-colSums(Corfun0)
# #any(is.na(csum))
# #which(is.na(csum))
# 
# ### try to delete the rows with zero
# rowCount=apply(Corfun, 1,FUN=sum)
# nozero=rowCount>0
# Corfun0=Corfun[nozero,]
# colnames(Corfun0)
# row.names(Corfun0)
# metacor0=my.metadata[nozero,]
# row.names(Corfun0)==row.names(metacor0)

# Model-based ordination

# Remove species present in only one sample
CountInSample = apply(fungl.abun.not.zero,2,sum)

# I needed to remove the species that were seen in a few samples.
fung.abun.reduced = fungl.abun.not.zero[,CountInSample > 10]
fung.abun.reduced = fung.abun.reduced[apply(fung.abun.reduced,1,sum) > 0,]

MetaOrd = MetaRich[rownames(MetaRich) %in% rownames(fung.abun.reduced),]

ModelOrd <- boral(fung.abun.reduced, family = "negative.binomial", num.lv = 2, 
                  n.burnin = 10, n.iteration = 100, n.thin = 1)

## model diagnostics
plot(ModelOrd, ask = FALSE, mfrow = c(2,2))

# Ordination plot. Please make it similar as the NMDS.
plot(ModelOrd$lv.median, col=as.numeric(MetaOrd$source), 
     pch=19, main="Latent variable model", las=1)
legend(0.8,-0.7, legend = levels(MetaOrd$source), cex = .8, pch=19,
       col=c(1,2,3), bty = "n")

# This should visualize the interaction effect between source*dust
# It seems that the abundance of at least some species abundance goes UP in dust 
# if there is more dust, while some species abundance goes DOWN in branches 
# if there is more dust.
# Probably it is a good idea to do an effect plot for each species,
# similarly to the shannon effect plots above.
ordisurf(ModelOrd$lv.median, MetaOrd$dust, add=T)

# The NMDS will not be necessary anymore. If you run the two lines below, 
# you will see how locality and dispersal are mixed up because of the overdispersion.
# This is very strong for the green samples (leaves?).
NMDS1<-metaMDS(fung.abun.reduced)
plot(NMDS1$points, ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2", 
     col=as.numeric(MetaOrd$source), pch=19)

# ### ploting the NMDS:
# ## plot NMDS
# par(mar=c(4,4,1,1), mfrow=c(1,1))
# # plot(NMDS1)
# plot(NMDS1$points, type="n", ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2", 
#      col=as.numeric(MetaOrd$source))
# ordispider(NMDS1,MetaRich$source, col="grey")
# points(NMDS1, pch=20, cex=exp(2*shannon)/100, col="grey")
# mylegend = legend(2, 1.5, c("leaf","branch","dust"), 
#                   fill=c("green","orange","blue"), border="white", bty="n")
# with(metacor0,ordiellipse(NMDS1, metacor0$source,cex=.5, 
#                           draw="polygon", col=c("green"),
#                           alpha=100,kind="se",conf=0.95, 
#                           show.groups=(c("Leaf"))))
# with(metacor0,ordiellipse(NMDS1, metacor0$source,cex=.5, 
#                              draw="polygon", col=c("orange"),
#                              alpha=100,kind="se",conf=0.95, 
#                              show.groups=(c("Branch"))))
# with(my.metadata,ordiellipse(NMDS1, metacor0$source,cex=.5, 
#                              draw="polygon", col=c("blue"),
#                              alpha=100,kind="se",conf=0.95, 
#                              show.groups=(c("Dust"))))
# # show overlapping samples
# ordiplot(NMDS1, type="text")
# ordispider(NMDS1, metacor0$source , col="grey")
# ## Three axes
# NMDS.3 <- metaMDS(Corfun0, k=3, trymax=100)
# NMDS.3 <- metaMDS(Corfun0, previous = NMDS.3, k=3, trymax=100)
# 
# 
# ### Axes 1 & 2
# ### why did we use shannon???
# png(file="community_NMDS_1-2.jpeg", units="mm", height=60, width=180, 
#     pointsize=10, bg="white", res=1200)
# par(mfrow=c(1,3))
# par(mar=c(4,4,1,1))
# NMDS.1.2 = ordiplot(NMDS.3, choices=c(1,2), type="n", 
#                        xlab="NMDS1", ylab="NMDS2")
# ordispider(NMDS.1.2,metacor0$source, col="grey")
# points(NMDS.3$points[,1], NMDS.3$points[,2], pch=20, 
#        cex=exp(2*shannon)/100,col="gray")
# ordiellipse(NMDS.1.2, metacor0$source,cex=.5, 
#             draw="polygon", col=c("green"),
#             alpha=100,kind="se",conf=0.95, 
#             show.groups=(c("Leaf")), border="green")
# ordiellipse(NMDS.1.2, metacor0$source,cex=.5, 
#             draw="polygon", col=c("orange"),
#             alpha=100,kind="se",conf=0.95,
#             show.groups=(c("Branch")), border="orange")
# ordiellipse(NMDS.1.2, metacor0$source,cex=.5, 
#             draw="polygon", col=c("gray"),
#             alpha=50,kind="se",conf=0.95,
#             show.groups=(c("Dust")), border="black")
# 
# ### Axes 1 & 3
# 
# png(file="community_NMDS_1-3.jpeg", units="mm", height=60, width=180, 
#     pointsize=10, bg="white", res=1200)
# par(mfrow=c(1,3))
# par(mar=c(4,4,1,1))
# NMDS.1.3 = ordiplot(NMDS.3, choices=c(1,3), type="n", 
#                     xlab="NMDS1", ylab="NMDS3")
# ordispider(NMDS.1.3,metacor0$source, col="grey")
# points(NMDS.3$points[,1], NMDS.3$points[,3], pch=20, 
#        cex=exp(2*shannon)/100,col="gray")
# ordiellipse(NMDS.1.3, metacor0$source,cex=.5, 
#             draw="polygon", col=c("green"),
#             alpha=100,kind="se",conf=0.95, 
#             show.groups=(c("Leaf")), border="green")
# ordiellipse(NMDS.1.3, metacor0$source,cex=.5, 
#             draw="polygon", col=c("orange"),
#             alpha=100,kind="se",conf=0.95,
#             show.groups=(c("Branch")), border="orange")
# ordiellipse(NMDS.1.3, metacor0$source,cex=.5, 
#             draw="polygon", col=c("gray"),
#             alpha=50,kind="se",conf=0.95,
#             show.groups=(c("Dust")), border="black")
# 
# ### Axes 2 & 3
# 
# png(file="community_NMDS_2-3.jpeg", units="mm", height=60, width=180, 
#     pointsize=10, bg="white", res=1200)
# par(mfrow=c(1,3))
# par(mar=c(4,4,1,1))
# NMDS.2.3 = ordiplot(NMDS.3, choices=c(2,3), type="n", 
#                     xlab="NMDS2", ylab="NMDS3")
# ordispider(NMDS.2.3,metacor0$source, col="grey")
# points(NMDS.3$points[,2], NMDS.3$points[,3], pch=20, 
#        cex=exp(2*shannon)/100,col="gray")
# ordiellipse(NMDS.2.3, metacor0$source,cex=.5, 
#             draw="polygon", col=c("green"),
#             alpha=100,kind="se",conf=0.95, 
#             show.groups=(c("Leaf")), border="green")
# ordiellipse(NMDS.2.3, metacor0$source,cex=.5, 
#             draw="polygon", col=c("orange"),
#             alpha=100,kind="se",conf=0.95,
#             show.groups=(c("Branch")), border="orange")
# ordiellipse(NMDS.2.3, metacor0$source,cex=.5, 
#             draw="polygon", col=c("gray"),
#             alpha=50,kind="se",conf=0.95,
#             show.groups=(c("Dust")), border="black")
# 
# ## PCA
# fun.pca = rda(Corfun0)
# fun.pca.scores = scores(fun.pca, choices=c(1,2,3))
# fun.pca.eigenvals = eigenvals(fun.pca)
# 
# ## explained by the first three axes: 52%
# (fun.pca.eigenvals[1] + fun.pca.eigenvals[2] + fun.pca.eigenvals[3])/sum(fun.pca.eigenvals)
# 
# ## axis distributions
# fun.pca1 = fun.pca.scores$sites[,1]
# fun.pca2 = fun.pca.scores$sites[,2]
# fun.pca3 = fun.pca.scores$sites[,3]
# 
# hist(fun.pca1^2)
# hist(fun.pca2^2)
# hist(fun.pca3^2)
# 
# ## dynamic 3D ordination plot
# 
# library(vegan3d)
# ordirgl(NMDS.3)
# orglspider(NMDS.3, metacor0$source)
# orgltext(NMDS.3, rownames(Corfun0))
# orgltext(NMDS.3, colnames(Corfun0))
# 
# ### NMDS plot for localities: 
# par(mar=c(4,4,1,1))
# plot(NMDS1$points, type="n", ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2")
# ordispider(NMDS1, metacor0$locality , col="grey")
# points(NMDS1, pch=20, cex=exp(2*shannon)/100, col="grey")
# mylegend2= legend(2, 1.5, c("BISOTON","HASAN ABAD","KHOSRO ABAD", "KEREND", "SORKHE DIZE"), 
#                   fill=c("green","orange","yellow","red","blue"), border="white", bty="n")
# with(metacor0,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
#                           draw="polygon", col=c("green"),
#                           alpha=100,kind="se",conf=0.95, 
#                           show.groups=(c("BISOTON"))))
# with(metacor0,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
#                           draw="polygon", col=c("orange"),
#                           alpha=100,kind="se",conf=0.95, 
#                           show.groups=(c("HASAN ABAD"))))
# with(my.metadata,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
#                              draw="polygon", col=c("yellow"),
#                              alpha=100,kind="se",conf=0.95, 
#                              show.groups=(c("KHOSRO ABAD"))))
# with(my.metadata,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
#                              draw="polygon", col=c("red"),
#                              alpha=100,kind="se",conf=0.95, 
#                              show.groups=(c("KEREND"))))
# with(my.metadata,ordiellipse(NMDS1, metacor0$locality,cex=.5, 
#                              draw="polygon", col=c("blue"),
#                              alpha=100,kind="se",conf=0.95, 
#                              show.groups=(c("SORKHE DIZE"))))
# 
# ### NMDS plot for time of sampling
# par(mar=c(4,4,1,1))
# plot(NMDS1$points, type="n", ylim=c(-1.5,1.5), xlab="NMDS1", ylab="NMDS2")
# ordispider(NMDS1, metacor0$time , col="grey")
# points(NMDS1, pch=20, cex=exp(2*shannon)/100, col="grey")
# 
# mylegend3= legend(2, 1.5, c("First sampling","Second sampling","Third sampling"), 
#                   fill=c("green","orange","blue"), border="white", bty="n")
# with(metacor0,ordiellipse(NMDS1, metacor0$time,cex=.5, 
#                           draw="polygon", col=c("green"),
#                           alpha=100,kind="se",conf=0.95, 
#                           show.groups=(c("1"))))
# with(metacor0,ordiellipse(NMDS1, metacor0$time,cex=.5, 
#                           draw="polygon", col=c("orange"),
#                           alpha=100,kind="se",conf=0.95, 
#                           show.groups=(c("2"))))
# with(my.metadata,ordiellipse(NMDS1, metacor0$time,cex=.5, 
#                              draw="polygon", col=c("blue"),
#                              alpha=100,kind="se",conf=0.95, 
#                              show.groups=(c("3"))))




### Community composition models
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



# ## Coefficients
# funcor.m2.coef = as.data.frame(funcor.m2$coefficients)
# 
# ## How do I know if i have an outlier or not?
# 
# 
# ## mean-centering the contrasts
# fun.coef.mean.contrast = funcor.m2.coef - apply(funcor.m2.coef,2,mean)
# 
# 
# ### plot effects
# 
# ?predict.glm
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
# 
# 
# ###source*dust interactions plot
# ### Do a simple glm for the two species affected by the source* dust
# library(MASS)
# Microsphaeriopsis.model= glm.nb(Corfun0$Microsphaeriopsis_olivacea ~ locality+time+source*dust, data= metacor0)
#                                
# plot(effect("source:dust",Microsphaeriopsis.model ,multiline=TRUE,confidence.level = 0.95))
# 
# anova(Microsphaeriopsis.model)
# 
# 
# Aureobasidium.model= glm.nb (Corfun0$Aureobasidium_sp_A30 ~ locality+time+source*dust, data= metacor0)
# 
# plot(effect("source:dust",Aureobasidium.model ,multiline=TRUE,confidence.level = 0.95))
# 
# anova(Aureobasidium.model)


### similarity analysis using anosim and adonis functions (corfun0=fungal abundance matrix for core OTUs and metacor0 
## is the metadata)
# braysimi= anosim(fung.abun.reduced,MetaOrd$source,permutations = 999,
#                  distance = "bray")
# 
# jaccardsimi= anosim(Corfun0,metacor0$source, permutations = 999, distance = "jaccard")
# 
# 
# corfunAdon= adonis(fung.abun.reduced~locality+time+source*dust,data=MetaOrd,
#                    permutations = 999, method = "bray")

### permanova for all of the species
# I don't trust this, just have a look on the NMDS with all the species.
# The permanova results are then only the statistics for that plot.
# permanova.total=adonis(fungl.abun2~locality+time+source*dust,data=MetaObs,permutations = 999,method = "bray")

# Compare distances among source classes (leaf, branch, dust)
# Use the sample placement according to the latent variable model
# Source for formatting code: http://stackoverflow.com/questions/17367277/how-to-extract-intragroup-and-intergroup-distances-from-a-distance-matrix-in-r


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
TukeyHSD(aov(value ~ groups, data = GroupedDist))

# It makes sense if you did a good surface sterilization: probably it is easier 
# for a fungus to go into a leaf than into a branch
boxplot(GroupedDist$value ~GroupedDist$groups)
