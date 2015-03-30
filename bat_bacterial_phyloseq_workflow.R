library("ggplot2")
library("ape")
library("plyr")
library("phyloseq"); packageVersion("phyloseq")
library("vegan")
library("RColorBrewer")
library("reshape")
library("scales")
library("DESeq2")

# Set this to where your working directory is i.e where you biom, mapping file, and 
# rep.set.tree are at
setwd("~/Desktop/bat_biom/")

###These need to be the names of your files generated in QIIME 1.8
biom_file = "non_chimeric_bat.biom"    ###This is your .biom file
map_file = "mapping_8_18_14.csv"      ###This is your mapping file with all the metadata
#tree_file = "rep_set_tree.tre"  ### This is the tree built after assinging all taxonomy

# This reads in your newick tree file.
# Comment this out if you don't have a tree and all occuanaces of tree need to be removed
#tree <-read_tree(tree_file)

# Reads in a QIIME formatted mapping file
map <- import_qiime_sample_data(map_file)

# Creates the phyloseq object with all data in it
phylo <- import_biom(biom_file,parseFunction=parse_taxonomy_greengenes)
# for fungal parseFunction=parse_taxonomy_default
# phylo_fungal <- import_biom(biom_file,parseFunction=parse_taxonomy_default)

# Just for error checking
warnings(phylo)
# I get a lot of warnings all related to some taxonomy assignment issues. 
intersect(phylo)

# This merges everthing into one phyloseq object
phylo <- merge_phyloseq(phylo,map)

# What follows here are a series of checks to verfy that your data was read
# read in correctly
phylo
# You should see something like this
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5408 taxa and 26 samples ]
# sample_data() Sample Data:       [ 26 samples by 15 sample variables ]
# tax_table()   Taxonomy Table:    [ 5408 taxa by 8 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 5408 tips and 5406 internal nodes ]

ntaxa(phylo)
nsamples(phylo)
sample_names(phylo)
rank_names(phylo)
sample_variables(phylo)
otu_table(phylo)[1:10, 1:5]
tax_table(phylo)[1:10, 1:5]

#Show you all unique taxa
get_taxa_unique(phylo, "Rank2")

# Prune OTUs that are not present in any of the samples
# This should be the only filtering you do at this step! Resisit the urge to filter!!!
phylo <- prune_taxa(taxa_sums(phylo) > 0, phylo)

# Print out number of reads per sample for the first 10 samples
sample_sums(phylo)[1:10]

# Basic statistics on the number of reads in the sample collection
mean(sample_sums(phylo))

min(sample_sums(phylo))

max(sample_sums(phylo))

sd(sample_sums(phylo))


# This will plot thefrequency of occuance vs average relative abundance of the OTUs
# Calculate average abundance of OTUs
otu.abun = apply(otu_table(phylo),1,mean)
# Calculate the frequency of each OTU across all samples
# 
otu.freq = rowSums(otu_table(phylo) != 0)/nsamples(phylo)
# Reassign names of phyla so we only color by the top 5 phyla and mark all others as "other"
phyla = as.vector(data.frame(tax_table(phylo))$Rank2)
levels(phyla) = c(levels(phyla),"other")
keephyla = c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Tenericutes")
phyla[!(phyla %in% keephyla)] = "Other"
phyla = as.vector(phyla)
phyla=as.factor(phyla)

otuabun = data.frame(abundance=log(otu.abun),frequency=otu.freq,phyla)

# Use color brewer to pick a color scheme for the phyla
brew = brewer.pal(6, "Set1")

# Create a scatterplot of OTUs showing their average relative abundance and frequency 
# This plot shows how rare and abundant OTUs are distributed across the
# your study.

ggplot(otuabun, aes(x=abundance,y=frequency,color=phyla)) + geom_point(size=3) +
  xlab("Average relative abundance (log scale)") + 
  ylab("Frequency in bat hosts") + 
  scale_colour_brewer(palette="Set2") +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=6)
  )


# Relative abundance of the top five phyla in your data set
# This plot shows how rare and abundant OTUs are distributed across the
# your study.

phylum.sum = tapply(taxa_sums(phylo), tax_table(phylo)[, "Phylum"], sum, na.rm = TRUE)
topphyla = names(sort(phylum.sum, TRUE))[1:5]
top5phyla = prune_taxa((tax_table(phylo)[, "Phylum"] %in% topphyla), phylo)
top5phyla
#
# Calculate average abundance of OTUs
top_otu.abun = apply(otu_table(top5phyla),1,mean)
str(top_otu.abun)
# Calculate the frequency of each OTU across all samples
top_otu.freq = rowSums(otu_table(top5phyla) != 0)/nsamples(top5phyla)
str(top_otu.freq)
top_phyla = as.vector(data.frame(tax_table(top5phyla))$Phylum)
top_phyla=as.factor(top_phyla)

topotuabun = data.frame(abundance=log(top_otu.abun),frequency=top_otu.freq,top_phyla)
str(topotuabun)
head(topotuabun)

# Use color brewer to pick a color scheme for the phyla
brew = brewer.pal(6, "Set1")

# Create a scatterplot of OTUs showing their average relative abundance and frequency 
ggplot(topotuabun, aes(x=abundance,y=frequency,color=top_phyla)) + geom_point(size=3) +
  xlab("Average relative abundance (log scale)") + 
  ylab("Frequency in bat hosts") + 
  scale_colour_brewer(palette="Set2") +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=6)
  )


plot_bar(phylo, x = "REGION", fill = "Phylum")+
  guides(fill = guide_legend(ncol = 3))+ 
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=8)
  )

# This is a gnarly way to make a nice OTU barplot. You will need your game on if you want it
# SUPER FANCY RELATIVE ABUNDANCE BAR PLOT, it's painful

# Stolen wholesale from Michelle Berry
# An analysis of butterfly gut bacteria in r
# http://www.rpubs.com/michberr/25722 accessed 9/29/2014

## First, I have to average the OTU counts for each individual within a species

# Make each species its own physeq object
speciesList <- tapply(sample_names(phylo), get_variable(phylo, "REGION"), c)
speciesPhyseq <- lapply(speciesList, prune_samples, phylo)

## For each item in speciesPhyseq find the average of the taxa rows

# Make a list of OTU tables for each species
speciesOTUtable <- lapply(speciesPhyseq,otu_table)

# Make a list of average OTU table for each species
speciesAvg <- lapply(speciesOTUtable,rowMeans)

# Put every column of average otu counts back into a matrix where each column is a 
# species and each row is an OTU
pooledOTUtable = t(do.call(rbind,speciesAvg))
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable)

# Add in taxonomy info for OTUs
TT = tax_table(phylo)
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(phylo))
pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
str(pOTUtax)
# You will need to change the 6 to the number of columns of metadata you have
# Add in a column for total sequences
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:6]))

# Take only top X OTUs (user defined)
# Change 500 to whatever number of OTUs you want to look at
pOTU = pOTU[order(-pOTU$SeqTotal),]
pOTUtop = pOTU[1:500,]

# This calculaton will tell you what percentage of the data you are representing in the plot
sum(pOTUtop$SeqTotal)/sum(pOTU$SeqTotal)

# Plot bar chart of phylum level differences
# Change the 6 to the number of metadata you have AND
# 8 determines which level of OTUs you look at 
pOTU.phylum =pOTUtop[,c(2:6,8)]
melt.phylum = melt(pOTU.phylum,id.vars="Phylum")
colnames(melt.phylum)[2]="species"
agg.phylum=aggregate(.~Phylum+species,melt.phylum,sum)


ggplot(agg.phylum,aes(x=species,y=value,fill=Phylum)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(labels = percent_format())+ 
  xlab("By Region Captured") +
  ylab("Relative Abudance") +
  scale_fill_manual(values = c("grey26","chartreuse3","cyan",  "red", 
                               "darkorange","royalblue", "darkgreen", 
                               "blue4", "yellow1", "violetred", 
                               "deepskyblue", "mediumorchid3",
                               "#89C5DA", "#DA5724", "#74D944", "#C84248", 
                               "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                               "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", 
                               "#8A7C64", "#599861")) +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=6)
  )+
  theme(axis.title.x = element_text(face="bold",size=16),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        axis.text.y = element_text(colour = "black", size=14),
        axis.title.y = element_text(face="bold", size=16),
        plot.title = element_text(size = 18),
        legend.title = element_text(size=14),
        legend.text = element_text(size = 13),
        legend.position="right",
        #Manipulating the facet features
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="black")) + # Black rectangle around facet title
  guides(fill = guide_legend(reverse = TRUE))  +
  ggtitle("Relative abudance of phyla (top 500 OTUs)") 


# UNRARIFEID ALPHA DIVERSITY
# By rarefying you REMOVE OTUs that actaully occur in your samples. Your alpha indices will
# then indicted samples are less rich and diverse then they actualy are.

# LOOK AT THE FIRST LINE THAT STARTS WITH plot_richness
# Change x = "someword" to the metadata category you want to look at
# Change ggtitle to the name of your title

# Scatter plots of alpha diversity measures
plot_richness(phylo, measures = c("Observed","Chao1", "Shannon"), x = "REGION") +
  ggtitle("Alpha Diversity Indices for PROJECTNAME by METADATA") +
  guides(fill = guide_legend(ncol = 3))+
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=6)
  )

# Boxplots of alpha diversity measures
plot_richness(phylo, measures = c("Observed","Chao1", "Shannon"), x = "REGION") + geom_boxplot() +
  ggtitle("Alpha Diversity Indices for PROJECTNAME by METADATA") +
  guides(fill = guide_legend(ncol = 3))+
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
    ,axis.text.x  = element_text(angle=90, vjust=0.5, size=6)
  )

phylo_alpha_richness = (estimate_richness(phylo,measures="Shannon"))


# # You may not want to do these transformations. They are better then rarefying 
# # but you are still losing data
# 
# # Transforms to fractional counts so all sample abundances = 1, i.e. normalized
# phylo_frac <- transform_sample_counts(phylo, function(OTU) OTU/sum(OTU))
# 
# # Transforms to even sampling depth of 1,000,000
# phylo_even <- transform_sample_counts(phylo, function(x) 1e+06 * x/sum(x))


# NMDS Ordination
# Simplified NMDS

# Specify the number of dimensions m you want to use (into which you want to scale down the
# distribution of samples in multidimensional space - that's why it's scaling).
# 
# Construct initial configuration of all samples in m dimensions as a starting point of
# iterative process. The result of the whole iteration procedure may depend on this step, 
# so it's somehow crucial - the initial configuration could be generated by random, but better 
# way is to help it a bit, e.g. by using PCoA ordination as a starting position.
# 
# An iterative procedure tries to reshuffle the objects in given number of dimension in such a
# way that the real distances among objects reflects best their compositional dissimilarity. 
# Fit between these two parameters is expressed as so called stress value - the lower stress value
# the better. 
# 
# Algorithm stops when new iteration cannot lower the stress value - the solution has 
# been reached.
# 
# After the algorithm is finished, the final solution is rotated using PCA to ease its interpretation
# (that's why final ordination diagram has ordination axes, even if original algorithm doesn't produce
# any).

# NMDS uses non-linear mapping

# Why Brays-Curtis
# invariant to changes in units
# unaffected by additions/removals of species that are not present in two communities
# unaffected by the addition of a new community
# recognize differences in total abundances when relative abundances are the same

phylo.ord <- ordinate(phylo, "NMDS", "bray")

# stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great,
# >0.2 is good/ok, and stress > 0.3 provides a poor representation

# Shows the stress of the samples. If there is a lage variation around the line then
# your stress is to high. Try k=3. 
stressplot(phylo.ord)

# How well your samples fit. Larger circles are a worse fit. 
gof <- goodness(phylo.ord)
plot(phylo.ord, display = "sites", type = "n")
points(phylo.ord, display = "sites", cex = 2*gof/mean(gof))

# NMDS with type = samples. So a plot of the distances between samples.
# READ THE FIRST LINE THAT STARTS WITH p2 <-
# change color = "someword" to one of your metadata categories

p2 <- plot_ordination(phylo, phylo.ord, type = "samples", color = "REGION")+
  #label="X.SampleID" )+
  geom_point(size = 6, alpha = 0.75)+ 
  guides(fill = guide_legend(ncol = 3))+
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
  )
p2

# Density plot version of the sample NMDS
den_nmds_sample = ggplot(p2$data, p2$mapping) + geom_density2d() +
  ggtitle("NMDS Density Plot on Brays Distance by Susceptibility") +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
  )
den_nmds_sample

# PCoA Analysis
# PCoA on Unifrac Weighted
# Use the ordinate function to simultaneously perform weightd UniFrac and then perform
# a Principal Coordinate Analysis on that distance matrix (first line)

ordu = ordinate(phylo, "PCoA", "bray")

plot_ordination(phylo, ordu, color = "REGION") +
  geom_point(size = 7, alpha = 0.75)+
  ggtitle("MDS/PCoA on Weighted-UniFrac Distance")+
  geom_vline(xintercept=c(0,0), linetype="dotted")+
  geom_hline(xintercept=c(0,0), linetype="dotted")+
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.background = element_blank()
  )

# Network Analysis
# Creates a network by OTUs and colors by metadata category
ig <- make_network(phylo,type ="samples",
                   max.dist = 0.7, distance = "bray", keep.isolates=TRUE)

plot_network(ig, phylo, color="REGION")


#DESEQ2

diagdds = phyloseq_to_deseq2(phylo, ~ PLACE)

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.1
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(cont_wat_func)[rownames(sigtab), ], "matrix"))
head(sigtab)

# This looks at just which taxa are enriched between the control samples and the water samples
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "ontology1","ontology2")]
posigtab

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(ontology2))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$ontology2, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$ontology2 = factor(as.character(sigtabgen$ontology2), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$ontology2, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$ontology2 = factor(as.character(sigtabgen$ontology2), levels=names(x))
ggplot(sigtabgen, aes(x=ontology2, y=log2FoldChange, color=ontology1)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

