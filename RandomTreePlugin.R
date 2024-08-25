library(microbiome)
library(ggplot2)
#library(phyloseq)
library(ape)
library(psadd)
library(picante)
library(httr)
library(foreach)
verbose()

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1]; 
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
   otu.path <<- paste(pfix, toString(parameters["otufile", 2]),sep="")
   tax.path <<- paste(pfix, toString(parameters["tax", 2]), sep="")
   map.path <<- paste(pfix, toString(parameters["mapping", 2]), sep="")
}

run <- function() {
   physeq <<- read_csv2phyloseq(otu.file=otu.path, taxonomy.file=tax.path, metadata.file=map.path)
   mytree <<- rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
   physeq <<- merge_phyloseq(physeq, mytree)
}
output <- function(outputfile) {
	ape::write.tree(phy_tree(physeq), outputfile)
}

