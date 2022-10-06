# getting all of the packages/programs needed
library(adegenet)
require(vegan)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(pairwiseAdonis)

#anything with 80% minInd ibs matrix is called OK.

# making matrix in R with values from 80% ibs matrix
inf=commandArgs(T)

inf="OK.ibsMat"
ma = as.matrix(read.table(inf))

# reading in file with list of samples (called 'bams'), cleaning up file names
bams=scan("bams_noclones",what="character") # list of bam files
bams=sub(".bam","",bams)
bams=sub("/.+/","",bams)
dimnames(ma)=list(bams,bams)

# reading sequencing depths (proportion of sites with coverage >5x, for each sample - output of plotQC.R)
# qc=read.table("quality.txt",header=F)
# names(qc)=c("srr","q")
# qc$srr=sub(".bam","",qc$srr)
# row.names(qc)=qc$srr
# qc=qc[bams,]

# plotting PCA with 80% ibs matrix
pp0=capscale(ma~1)
plot(pp0, scaling=1)

# doing this again with OK50 (50% minInd ibs matrix)

inf50="OK50.ibsMat"
ma50 = as.matrix(read.table(inf50))

dimnames(ma50)=list(bams,bams)

pp1=capscale(ma50~1)
plot(pp1, scaling=1)

# seeing if there is a significant difference between 50 and 80% ibs matrix - to see if there is any point in using 50% minInd matrix
plot(procrustes(pp0, pp1))

### MAKING CONDITIONS TABLE

# making vectors with species, spawning time and location
species <- substr(bams, 2, 2)
spawn <- substr(bams, 3,3)
location <- substr(bams,2,3)
location <- recode(location, MA = "NR", SA = "BI", SS = "BI", MS = "AR")

# making a table combining these vectors
conds=data.frame(cbind(species,spawn,location))
conds

# renaming the factors to make sense in the legends
conds$location <- gsub("NR", "Ningaloo Reef", conds$location)
conds$location <- gsub("AR", "Ashmore Reef", conds$location)
conds$location <- gsub("BI", "Barrow Island", conds$location)
conds$spawn <- gsub("A", "Autumn", conds$spawn)
conds$spawn <- gsub("S", "Spring", conds$spawn)
conds$species <- gsub("M", "Acropora millepora", conds$species)
conds$species <- gsub("S", "Acropora samoensis", conds$species)

# making a new factor, north and south
nthsth=location
nthsth <- recode(nthsth, south = "South",south = "South", north = "North")
conds=cbind(conds,nthsth)


### SAMOENSIS ONLY

# making a new matrix with only samoensis samples, then making pps and plotting it
sam=ma[25:49,25:49]
pps=capscale(sam~1)
plot(pps, scaling=1)

# making scoressam, which is a table of the values for MDS1 and MDS2 from samoensis PCA - for plotting
scoressam=data.frame(scores(pps,display="sites",choices=axes2plot))

# plotting samoensis samples coloured by spawning time - fix colours
plot(scoressam[,1:2],pch=16,col=myspawnColors[spawncolor],asp=1)


### ADONIS2

# running adonis2 with all of the samples, using different combinations of factors
adonis2(ma~spawn+species+location,conds)
adonis2(ma~species+spawn+location,conds)
adonis2(ma~location+species+spawn,conds)
adonis2(ma~location+spawn+species,conds)

adonis2(ma~location+spawn,conds)
adonis2(ma~species+location,conds)

adonis2(ma~species,conds)
adonis2(ma~location,conds)

pairwise.adonis(ma,conds$location)
pairwise.adonis2(ma~location,conds)

# adonis with north south
adonis2(ma~nthsth+spawn+species,conds)

# adonis with only samoensis
adonis2(sam~spawnsam)


# plot all samples with spawning time as the colours
ggplot(scores,aes(MDS1,MDS2,color=conds$spawn))+geom_point()+coord_equal()+theme_bw()+labs(colour="Spawning Time")

# plot all samples with species as the colours
ggplot(scores,aes(MDS1,MDS2,color=conds$species))+geom_point()+coord_equal()+theme_bw()+labs(colour="Species")

# plot all samples with location as the colours
ggplot(scores,aes(MDS1,MDS2,color=conds$location))+geom_point()+coord_equal()+theme_bw()+labs(colour="Location")+ theme(legend.title = element_text(face = "bold", size=10))+ theme(legend.text = element_text(size = 8, colour = "black"))

# plot all samples with spawning time as colour, location as symbol
ggplot(scores,aes(MDS1,MDS2,color=conds$spawn,))+geom_point((aes(shape=conds$location)))+coord_equal()+theme_bw()+labs(colour="Spawning Time",shape="Location")+ theme(legend.title = element_text(face = "bold", size=10))+ theme(legend.text = element_text(size = 8, colour = "black"))

# plot all samples with north/south as the colours

ggplot(scores,aes(MDS1,MDS2,color=conds$nthsth))+stat_ellipse()+geom_point()+coord_equal()+theme_bw()+labs(colour="Location")+theme(aspect.ratio=1,legend.title = element_blank(),legend.text = element_text(size = 9, colour = "black"),panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position=c(0.15,0.85),legend.background = element_blank(),legend.box.background = element_rect(colour = "grey"))
ordiellipse(scores[,1:2],group= conds$nthsth,draw="polygon",label=T)



### TESTING ASSUMPTIONS

# making histogram of number of sites 
hist(log(as.numeric(znum),10),xlab="Log(Number of Sites)",breaks=10, main=NULL, col="lightblue")

# testing for uneven coverag - plotting coverage
ggplot(scores,aes(MDS1,MDS2,color=logznum))+scale_color_continuous(type="viridis")+geom_point()+coord_equal()+theme_bw()+labs(colour="Log(Number of Sites)")+theme(aspect.ratio=1,legend.text = element_text(size = 9, colour = "black"),panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position="right",legend.background = element_blank(),legend.box.background = element_blank())


# testing for homogeneity of variance
distma=dist(ma, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

# permutation test of homogeneity of variance between locations
betaobj=betadisper(distma, as.factor(location), type = c("median","centroid"), bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
permutest(betaobj, pairwise = FALSE, permutations = 999,parallel = 1)

# boxplot of location (homogeneity of variance)
boxplot(betaobj, ylab = "Distance to centroid",xlab="Location", col="light grey", border="black", boxwex=0.5)
boxplot(betaobj,ylab= "Distance to centroid",xlab="Location",medcol = "black", medlty=1 ,boxlty = 0, whisklty = 1, staplelwd = 1, outpch = 2, outcex = 1)