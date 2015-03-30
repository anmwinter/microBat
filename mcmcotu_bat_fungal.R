library(MCMC.OTU)

#Data needs to look like this:
#Rowname  Data  Data  sample  otu otu
#SampleID data  data  SampleID  otu otu

#Loads your own data, use functions read.table or read.csv or read.delim
setwd("~/Desktop/bat_microbiota/")

func_level1 <- read.delim("~/Desktop/bat_fungal/bat_fungal_mcmcotu.txt",sep="\t",row.names=1)
str(func_level1)
head(func_level1)
# removing low-count samples and OTUs (?purgeOutliers to see how to adjust cutoffs)
# also adjust otu.columns value for your dataset
func_level1=purgeOutliers(func_level1,count.columns=6:39,otu.cut=0.001)

# this is how the format of the data table should look
# note the 'sample' column giving a unique name for each sample - this one is required
# (row names are not required)
head(func_level1)

# what is the proportion of samples with data for these OTUs?
apply(func_level1[,5:length(func_level1[1,])],2,function(x){sum(x>0)/length(x)})

# what percentage of global total counts each OTU represents?
apply(func_level1[,6:length(func_level1[1,])],2,function(x){sum(x)/sum(func_level1[,6:length(func_level1[1,])])})

# stacking the data; adjust otu.columns and condition.columns values for your data
gss=otuStack(func_level1,count.columns=c(6:length(func_level1[1,])),condition.columns=c(1:5))
head(gss)

# fitting the model. Replace the formula specified in 'fixed' with yours, add random effects if present. 
# See ?mcmc.otu for these and other options. 
mm=mcmc.otu(
  fixed="REGION",
  data=gss,
  nitt=55000,thin=50,burnin=15000 # a long MCMC chain to improve modeling of rare OTUs
)

summary(mm)
# selecting the OTUs that were modeled reliably
# (OTUs that are too rare for confident parameter estimates are discarded) 
acpass=otuByAutocorr(mm,gss)
head(acpass)
# calculating differences and p-values between all pairs of factor combinations
smm0=OTUsummary(mm,gss,otus=acpass,summ.plot=TRUE) 

# adjusting p-values for multiple comparisons:
#?p.adjust for different methods
smmA=padjustOTU(smm0, method="BH")

# significant OTUs at FDR<0.05:
sigs=signifOTU(smmA,p.cutoff=0.05)
sigs

# plotting the significant ones
smm1=OTUsummary(mm,gss,otus=sigs)

# now plotting them by species
smm1=OTUsummary(mm,gss,otus=sigs,xgroup="color")

# table of log10-fold changes and p-values: this one goes into supplementary info in the paper
smmA$otuWise[sigs]
write.table(smmA$otuWise[sigs],file = "log10_fold_pvalues_by_place.csv",sep=",")

summary(mm)
#eff.samp values much smaller than the nominal sample size indicate poor mixing
#post.mean = posterier mean
#Intercept = intercept of the post

autocorr(mm$VCV)
#Ideally, all samples of the posterior distribution should be independent,
#and the autocorrelation for all lag values greater than zero should be near zero
#If not near zero you may want to run the model longer

posterior.mode(mm$VCV)
#Obtain estimates of the additive phylogenetic and residual variance
#by calculating the modes of the posterior distributions
#See proportion!

posterior.mode(mm$Sol)
#Most likely value is where the posterior distribution peaks

#############################################
# Principal coordinate analysis

library(vegan)
library(MCMC.OTU)

# purging under-sequenced samples and OTUs represented in less than 10% of all samples 
func_level1=purgeOutliers(func_level1,count.columns=3:56,zero.cut=0.1,otu.cut=0)

# creating a log-transfromed normalized dataset for PCoA:
func_level1.log=logLin(data=func_level1,count.columns=4:length(names(func_level1)))

# computing Bray-Curtis distances and performing PCoA:
func_level1.dist=vegdist(func_level1.log,method="bray")
func_level1.pcoa=pcoa(func_level1.dist)

# plotting by bank:
scores=func_level1.pcoa$vectors
conditions=func_level1[,1:3]
head(conditions)
margin=0.01
plot(scores[,1], scores[,2],type="n",
     xlim=c(min(scores[,1])-margin,max(scores[,1])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("Axis1 (",round(func_level1.pcoa$values$Relative_eig[1]*100,1),"%)",sep=""),
     ylab=paste("Axis2 (",round(func_level1.pcoa$values$Relative_eig[2]*100,1),"%)",sep=""),
     main="PCoA colored by bank")
points(scores[conditions$type=="Lavatubes",1],scores[conditions$bank=="Lavatubes",2])
points(scores[conditions$type=="Surface Soil",1],scores[conditions$bank=="Surface Soil",2],pch=19)
legend("bottomright", c("East","West"), pch=c(1, 19), cex=0.8)

# plotting by species:
margin=0.01
conditions=func_level1[,1:3]
plot(scores[,1], scores[,2],type="n",
     xlim=c(min(scores[,1])-margin,max(scores[,1])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("Axis1 (",round(func_level1.pcoa$values$Relative_eig[1]*100,1),"%)",sep=""),
     ylab=paste("Axis2 (",round(func_level1.pcoa$values$Relative_eig[2]*100,1),"%)",sep=""),
     main="PCoA colored by species")
points(scores[conditions$species=="franksi",1],scores[conditions$species=="franksi",2])
points(scores[conditions$species=="faveolata",1],scores[conditions$species=="faveolata",2],pch=19)
legend("bottomright", c("franksi","faveolata"), pch=c(1, 19), cex=0.8)


#############################################
# exploring correlations between OTUs 

library(MCMC.OTU)
#data(green.data)
#func_level1=purgeOutliers(green.data,count.columns=4:156,otu.cut=0.001,zero.cut=0.2)

# creating a log-transformed normalized dataset ignoring zero counts:
nl=startedLog(data=func_level1,count.columns=4:length(names(func_level1)),logstart=0)
names(nl)=c("I","II","III","IV","V","VI","VII")

# displaying a matrix of scatterplots and p-values of OTU correlations
# (onlu p-values better than 0.1 are displayed)
pairs(nl,lower.panel=panel.smooth,upper.panel=panel.cor.pval)

#----------------------------------------