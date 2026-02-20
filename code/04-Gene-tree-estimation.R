# 04 Gene tree estimation

#Load libraries
library(ape)
library(phangorn)
library(phytools)
library(ggplot2)

#Check the working directory. we want to be in the results/RAxML folder
getwd() 

setwd("/Users/morven/Desktop/PLPATH875/results/RAxML") 

#List all .bestTree files. $ ensures the end of the name
tree_files <-list.files(pattern="\\.raxml.bestTree$")

# list with all the trees
gene_trees<- list() 

#make it a multiphylo object for ease of use with other 
class(gene_trees)<- "multiPhylo" 

i<-1
for(tree_file in tree_files){ ##go thru each file and read the tree
  gene_trees[[i]]<- read.tree(tree_file)
  i<-i+1
}

#gene_trees[[1]]
#plot(gene_trees[[1]])

#Find the appropriate root taxon
root_taxa <- c("H_vulgare_HVens23", "Er_bonaepartis_TB1", "S_vavilovii_Tr279", "Ta_caputMedusae_TB2")

# Extract the first matching species for each tree
gene_tree_outgroup<- rep(NA,length(gene_trees))
for(i in 1:length(gene_trees)){
  gene_tree <- gene_trees[[i]]
  found_taxa <- root_taxa[root_taxa %in% gene_tree$tip.label]
  
  # Return the first one if it exists
  if (length(found_taxa) > 0) {
    gene_tree_outgroup[i]<-found_taxa[1]
  }
}

#gene_tree_outgroup

#check if any gene trees don't have an outgroup taxon
any(is.na(gene_tree_outgroup)) ## should return FALSE

#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(gene_trees)){
  gene_trees[[i]]<- root(gene_trees[[i]],
                         outgroup = gene_tree_outgroup[i],
                         resolve.root=TRUE)
  gene_trees[[i]]<-chronos(gene_trees[[i]]) ## make ultrametric for nicer densitree
}

#plot(gene_trees) to plot all gene trees sequentially
plot(gene_trees[[2]]) # plots the second gene tree

# see which genes have which taxa
all_labels <- unlist(lapply(gene_trees, function(x) x$tip.label)) ##count how often each individual appears
df <- as.data.frame(table(all_labels) / length(gene_trees)) #organize as a nice data frame for plotting

ggplot(df, aes(x = all_labels, y = Freq)) +
  geom_col() +
  labs(x = "Individual", y = "Proportion")+
  theme(axis.text.x = element_text(angle = 90))

# Many genes nevertheless lack some taxa — not every taxon appears in every gene — which makes it harder to summarize and visualize species relationships.

# A simple maximum parsimony approach to create a supertree from all our genes to get a sense of the overall signal from the gene trees.
st<-superTree(gene_trees)
st<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st)

# Plot all of the genes as a densitree plot, a plot where all the gene tree topologies are overlaid over one another to give a high-level cloud of all the signal
densiTree(gene_trees,consensus=st,scaleX=T,type='cladogram', alpha=0.05)

common_tips <- Reduce(
  intersect,
  lapply(gene_trees, function(tr) tr$tip.label)
)

length(common_tips) ## 6

trees_pruned <- lapply(
  gene_trees[1:10],
  function(tr) drop.tip(tr, setdiff(tr$tip.label, common_tips))
)

## to remove NULLs:
## trees_pruned <- trees_pruned[!sapply(trees_pruned, is.null)]

densityTree(trees_pruned,use.edge.length=FALSE,type="cladogram",nodes="centered")

trees = read.tree(file="Ae_bicornis_Tr406_BIS2_Contig10132_simExt_macseNT_noFS_clean.aln.raxml.mlTrees")

rtrees <- lapply(trees, function(tr) {
  root(tr, outgroup = "H_vulgare_HVens23", resolve.root = TRUE)
})

# Plot the densiTree (before and after rooting)
densityTree(trees,type="cladogram",nodes="intermediate")
densityTree(rtrees,type="cladogram",nodes="intermediate")

# If we don’t care about differences in branch lengths and we want to focus on differences in topology
densityTree(trees,use.edge.length=FALSE,type="cladogram",nodes="centered")
densityTree(rtrees,use.edge.length=FALSE,type="cladogram",nodes="centered")






