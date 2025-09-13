# install packages if needed
if(!requireNamespace("ape", quietly=TRUE)) install.packages("ape")
if(!requireNamespace("phangorn", quietly=TRUE)) install.packages("phangorn")
if(!requireNamespace("tidyverse", quietly=TRUE)) install.packages("tidyverse")

library(ape)
library(phangorn)
library(tidyverse)

counts_file <- "Orthogroups.GeneCount.tsv"     # OrthoFinder family counts
tree_file   <- "SpeciesTree_rooted.txt"        # OrthoFinder species tree (Newick)
out_prefix  <- "welwitschiae_subset"           # output prefix
# ----------------------------------------------------------------

#Read Orthogroups counts
og <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)

# Ensure first column is the family id (often "Orthogroup")| Here, I am keeping only my species of interest
fam_colname <- colnames(og)[1]
keep_species <- c(
  "Aspergillus_welwitschiae_cbs",
  "Aspergillus_welwitschiae_ccmb663",
  "Aspergillus_welwitschiae_ccmb674",
  "Aspergillus_welwitschiae_ihem",
  "Aspergillus_welwitschiae_imrd1",
  "Aspergillus_welwitschiae_item",
  "Aspergillus_welwitschiae_ocstreb1"
)

gene_counts_subset <- og %>%
  select(Orthogroup, all_of(keep_species))
head(gene_counts_subset)

#Drop families less than 4 counts across all selected species
gene_counts_subset <- gene_counts_subset[rowSums(gene_counts_subset[keep_species], na.rm = TRUE) > 4, ]

#Read and prune the tree
tr <- read.tree(tree_file)
message("Original tree tips: ", length(tr$tip.label))

#Root by outgroup if present
outgroup_name <- "Neurospora_crassa"
if(outgroup_name %in% tr$tip.label) {
  tr <- root(tr, outgroup = outgroup_name, resolve.root = TRUE)
  message("Rooted tree by outgroup: ", outgroup_name)
} else {
  message("Outgroup not found in tree tips - proceeding with input tree as-is.")
}

#prune tree to only the selected species (drop any extra species)
keep_tips <- intersect(tr$tip.label, keep_species)
message("Tips present in both tree and counts: ", length(keep_tips), " / ", length(keep_species))
missing_in_tree <- setdiff(keep_species, tr$tip.label)
if(length(missing_in_tree) > 0){
  message("WARNING: these species are in counts but not in tree -> ", paste(missing_in_tree, collapse=", "))
  # stop or decide to rename - we continue and will error if mismatch remains
}

tr_pruned <- drop.tip(tr, setdiff(tr$tip.label, keep_tips))

#check that names match exactly between pruned tree and counts columns
tree_tips <- tr_pruned$tip.label
count_cols <- keep_species

m1 <- setdiff(tree_tips, count_cols)
m2 <- setdiff(count_cols, tree_tips)

if(length(m1) || length(m2)){
  message("Name mismatches detected.")
  if(length(m1)) message("In tree but not in counts: ", paste(m1, collapse=", "))
  if(length(m2)) message("In counts but not in tree: ", paste(m2, collapse=", "))
  stop("Please fix naming differences (typos, underscores, suffixes) before proceeding.")
}

#make the pruned tree ultrametric (chronos). If chronos fails, fallback to Grafen
message("Attempting to make tree ultrametric using chronos (penalized likelihood).")
tr_ultra <- NULL
try({
  # chronos requires branch lengths > 0; ensure no zero-length branches
  if(any(tr_pruned$edge.length <= 0, na.rm = TRUE)){
    tr_pruned$edge.length[tr_pruned$edge.length <= 0] <- 1e-8
  }
  tr_ultra <- chronos(tr_pruned, lambda = 1)  # lambda can be tuned
  message("chronos completed.")
}, silent = TRUE)

if(is.null(tr_ultra)){
  message("chronos failed or unavailable — falling back to Grafen's method (approximate ultrametric).")
  # Grafen produces an ultrametric tree; we keep branch lengths
  tr_ultra <- compute.brlen(tr_pruned, method = "Grafen", power = 1)
  # ensure it's of class phylo
  class(tr_ultra) <- "phylo"
}

# ensure tree is rooted
if(!is.rooted(tr_ultra)){
  tr_ultra <- root(tr_ultra, outgroup = tr_ultra$tip.label[1], resolve.root = TRUE)
}

# reorder columns of count table to match tree tip order (very important)
ordered_species <- tr_ultra$tip.label
sub_ordered <- gene_counts_subset[, c(fam_colname, ordered_species)]

# final check: no NA and numeric
if(any(is.na(sub_ordered[ordered_species]))){
  message("NA detected in counts for selected species - converting NAs to 0")
  sub_ordered[ordered_species][is.na(sub_ordered[ordered_species])] <- 0
}
#Add the "Desc" column for CAFE which is required by CAFE 5.0
sub_ordered <- sub_ordered %>%
  mutate(Desc = "(null)") %>%
  select(Desc, everything())  # move Desc to first column


#write outputs for CAFE
cafe_counts_file <- paste0(out_prefix, "_cafe_counts_min4.tsv")
write.table(sub_ordered, file=cafe_counts_file, sep="\t", quote=FALSE, row.names=FALSE)

cafe_tree_file <- paste0(out_prefix, "_cafe_tree_ultra_min4.nwk")
write.tree(tr_ultra, file=cafe_tree_file)

# 8. also write a text file with the species list used (optional)
writeLines(ordered_species, paste0(out_prefix, "_species_list.txt"))

message("Wrote counts to: ", cafe_counts_file)
message("Wrote ultrametric tree to: ", cafe_tree_file)
message("Species list: ", paste0(out_prefix, "_species_list.txt"))
