
# makes an object called bams which is a list of the bam files
bams=read.table("clonebams")[,1] # list of bam files

# substitute "" with "" bams is a vector (string of text) - we want to remove .bam and just have the sample names instead of filenames
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("OKall.ibsMat"))

# changing the row and column names to the list of samples that is listed in "bams"
dimnames(ma)=list(bams,bams)

# looking at the first 6 rows and columns (looking at a section of the table) - need square brackets for this 
ma[1:6,1:6]

# as.dist - changing the matrix to distances (formatting)
# hclust - heirarchical clustering function - builds a phylogenetic tree - clusters them according to the distances (in this case this is genetic differences)
# second argument for hclust is the type of clustering we are asking for - ave = average 
# putting the result of this into an object called hc 
hc=hclust(as.dist(ma),"ave")

# look at this tree to identify clones - there were none