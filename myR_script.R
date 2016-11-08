####################################
####################################
## Community analysis of Persian oak fungal endophytes 
## under dust storm conditions shows context dependency##
## BY: N. Hagh Doust, M. Akbarinia, N. Safaie, H. Yousefzadeh and M. Bálint
###################################
###################################

## Pakages we need for the analysis
library(mvabund)
library(vegan)
library(coda)
library(rjags)
library(boral)
library(effects)
library(MASS)
library(reshape)
##################################
##################################

##### 1- Isolation success analysis:
## Isolation success data input

# Get the input data from here first:
# https://figshare.com/s/0f4283512f39cc4f3c0e
# https://figshare.com/s/fbf753647e6e86467707

MyAbund = read.csv(file="morphotype_matrix_incubator.csv", 
                   header = T, row.names = 1)
MetaData = read.csv(file="metadata.csv", header = T, row.names = 1)

##### 1- Isolation success calculation
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

### Isolation success model results
success.sammary = summary(success.m)
success.anova = anova(success.m, test="LRT")

# Visualize the culturing success model results
# Q1: Dust deposition and isolation success
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(MetaData$success ~ MetaData$dust, xlab=c("Dust Deposition (mg/cm2)"), 
      ylab=c("Isolation success"))

# Q2: Locality and isolation success
boxplot(fitted(success.m) ~ MetaData$locality, 
        ylab="Isolation success")

# Q4: Source of isolation
boxplot(fitted(success.m) ~ MetaData$source,
        ylab="Source of isolation")

## Creating the interaction object for ploting
success.effect = effect("source:dust",success.m , multiline=TRUE)

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

Source = factor(MetaData$source, levels = c("Leaf", "Branch", "Dust"))
MetaData$time <- factor(MetaData$time, levels = c("1","2","3"))

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
Richness.m= glm(Richness ~ locality + time + source*dust, data = MetaRich, 
                family=poisson(link = "log"))

### Richness Model Evaluations
Richness.anova = anova (Richness.m, test="LRT")
Richness.summary = summary (Richness.m)

# Visualize Richness changes
par(mfrow=c(1,1))
boxplot(fitted(Richness.m) ~ MetaRich$source, xlab="Source", ylab="Richness")
# Richness in branch and dust similar, but richness stat. sign. lower in leaf

# Creating the interaction object for ploting
Richness.effect = effect("source:dust",Richness.m, multiline=TRUE, confidence.level = 0.95)

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
# Shannon model
shannon.m = lm (shannon ~ locality + time + source*dust, data = MetaNotOne)

## Shannon model statistics
shannon.summary = summary (shannon.m) 
shannon.anova = anova (shannon.m, test = "LRT")

# Creating the interaction object for ploting
shannon.effect = effect("source:dust", shannon.m, multiline=TRUE)
shan.eff.sum = summary (shannon.effect)

## Simpson model:
simpson.m = lm (simpson ~ locality + time + source*dust, data=MetaNotOne)
simpson.summary = summary (simpson.m)
simpson.anova = anova (simpson.m, test= "Chisq")

## Creating the interaction object for ploting
simpson.effect = effect("source:dust", simpson.m, multiline=TRUE)

### PLOT all interaction objects together
success.effect = effect("source:dust",success.m , multiline=TRUE)
Richness.effect = effect("source:dust",Richness.m, multiline=TRUE, confidence.level = 0.95)
shannon.effect = effect("source:dust", shannon.m, multiline=TRUE)
simpson.effect = effect("source:dust", simpson.m, multiline=TRUE)

### USE the results from effect objects in the plots
### The plot
dev.off()
pdf(file = "final diversity interaction plot.pdf", paper = "a4", width = 7, height = 4)
par(mfrow=c(1,3),xpd=TRUE, oma= c(0,0,2,0))
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
reset <- function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}
reset()
legend("top",legend = c("Isolation success","Richness","Shannon diversity",
                        "Simpson diversity"),lwd=2, horiz=TRUE,
       col= c("black","red","blue","green"),  border="white", bty="n", cex =0.8)
dev.off()

#######################################################
#######################################################
### 3- Community Composition Analysis

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

### source*dust interactions 
### Negative binomial glm for the species affected by the source* dust interaction
Microsphaeriopsis.model= glm.nb(MyAbund$Microsphaeriopsis_olivacea ~ locality +
                                  time + source*dust,data= MetaData)
Micros.anova = anova (Microsphaeriopsis.model, test = "Chisq")                             
Micros.summary= summary(Microsphaeriopsis.model)

# ploting the interaction effect
Micros.effec= effect("source:dust", Microsphaeriopsis.model, multiline = FALSE)
Micros.effec.sum = summary (Micros.effec)

Mic.eff<-plot(Micros.effec, ylab = "M. olivacea abundance",
              xlab = "Dust Deposition (mg/cm2)",main= NULL, ylim = c(-30,5),
             se=FALSE,ci.style="none")
dev.off()
pdf(file = "M.oliv interaction plot.pdf", paper = "a4", width = 7, height = 4)

Mic.eff # plot original object
 
Mic.eff$x.scales$cex <- c(1.1, 1.1)
Mic.eff$xlab$cex <- 1.1
Mic.eff$ylab$cex <- 1.1
Mic.eff # plot edited object
dev.off()

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
### creating a negative binomial GLM model for each OTU
Alternaria_sp_A25.model= glm.nb(MyAbund$Alternaria_sp_A25 ~ locality +
                                  time + source*dust,data= MetaData)
Alternaria_sp_A9.model= glm.nb(MyAbund$Alternaria_sp_A9 ~ locality +
                            time + source*dust,data= MetaData)
Aureobasidium_sp_A30.model= glm.nb(MyAbund$Aureobasidium_sp_A30 ~ locality +
                               time + source*dust,data= MetaData)
Byssochlamys_spectabilis_.model= glm.nb(MyAbund$Byssochlamys_spectabilis_ ~ locality +
                           time + source*dust,data= MetaData)
Cladosporium_herbarum_A8.model= glm.nb(MyAbund$Cladosporium_herbarum_A8 ~ locality +
                                time + source*dust,data= MetaData)
Penicillium_sp_A21.model= glm.nb(MyAbund$Penicillium_sp_A21 ~ locality +
                          time + source*dust,data= MetaData)
Pleosporaceae_sp_A5.model= glm.nb(MyAbund$Pleosporaceae_sp_A5 ~ locality +
                            time + source*dust,data= MetaData)
## interaction effect objects
Alternaria_sp_A25.model.effec= effect("source:dust",Alternaria_sp_A25.model)
Alternaria_sp_A9.model.effec= effect("source:dust", Alternaria_sp_A9.model)
Aureobasidium_sp_A30.model.effec= effect("source:dust", Aureobasidium_sp_A30.model)
Byssochlamys.effec= effect("source:dust", Byssochlamys_spectabilis_.model)
Cladosporium_herbarum_A8.model.effec= effect("source:dust", Cladosporium_herbarum_A8.model)
Penicillium_sp_A21.model.effec= effect("source:dust", Penicillium_sp_A21.model)
Pleosporaceae_sp_A5.model.effec= effect("source:dust", Pleosporaceae_sp_A5.model)

### The reaction plot
dev.off()
pdf(file = "OTUs interaction plot.pdf", paper = "a4", width = 7, height = 4)
par(mfrow=c(1,3),xpd=TRUE, oma= c(0,0,2,0))
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
points(c(0.01,0.06), c(0.0009805355,5.550239e-09), type = "l",lwd=2,lty=3, col= "black")#A5
#Branch
plot(c(0.01,0.06), c(0.0,4.224453e-2), type = "n", xlab="Dust deposition (mg/cm2)",
     main = "Branch", ylab = "",cex.axis = 1.5,
     cex.lab = 1.5,cex.main= 1.5)
points(c(0.01,0.06), c(0.03120315,1.575250e-05), type = "l",lty=6, lwd=2, col="black")#25
points(c(0.01,0.06), c(0.039682698,5.224453e-04), type = "l",lty=5, lwd=2, col="blue")#A9
points(c(0.01,0.06), c(3.372862e-03,3.891702e-02), type = "l",lty=3, lwd=2, col= "red")#A30
points(c(0.01,0.06), c(0.0311248787,0.0273451531), type = "l",lty=1,lwd=2, col= "green")#Byssochlamys 
points(c(0.01,0.06), c(0.032479263,0.02415204), type = "l",lwd=2,lty=6, col="red")#A8
points(c(0.01,0.06), c(0.010226362,3.445531e-08), type = "l",lwd=2,lty=4, col= "blue")#A21
points(c(0.01,0.06), c(0.0144565757,7.301888e-04), type = "l",lwd=2,lty=3, col= "black")#A5
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
points(c(0.01,0.06), c(0.0238803385,8.354105e-01), type = "l",lwd=2,lty=3, col= "black")#A5
reset <- function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}
reset()
legend(0.04, 1.03, c("A25", "A9", "A30","Bys","A8","A21","A5"), col = c("black", "blue", 
                 "red","green","red","blue","black"),box.col="white",horiz=TRUE,
       lwd=2,lty = c(6,5,3,1,6,4,3), cex = 1)
dev.off()

#############################
#############################
## 4- Similarity analysis

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
pdf(file = "LVM.pdf", paper = "a4", width = 6.5, height = 6.5)
ordibora= ordiplot(ModelOrd$lv.median, choices = c(1,2),display = "sites",
                   type = "none", cex = 1 , xlim = c(-0.3,0.3), cex.axis = 1,
         cex.lab = 1.2)
points(ordibora,"sites", pch=20 ,col=as.numeric(MetaOrd$source))
#ordispider(ordibora,MetaOrd$source, col= "gray" )
mylegend = legend("topright",legend=c("Leaf","Branch","Dust"), 
                 fill =  c("black", "red", "green"), border="white", bty="n", cex = 1)
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=0.5, 
                          draw="polygon", col="black",
                          alpha=200,kind="se",conf=0.95, 
                          show.groups=(c("Leaf"))))
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=0.5, 
                          draw="polygon", col="red",
                          alpha=200,kind="se",conf=0.95,
                          show.groups=(c("Branch")))) 
with(MetaOrd,ordiellipse(ordibora, MetaOrd$source,cex=0.5, 
                          draw="polygon", col="green",
                          alpha=200,kind="se",conf=0.95,
                          show.groups=(c("Dust"))))
dev.off()
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

# It makes sense: probably it is easier for a fungus to go into a leaf than into a branch
boxplot(GroupedDist$value ~GroupedDist$groups)


