install.packages("tidyverse")
install.packages("UpSetR")
library(tidyverse)
library(UpSetR)

# Read OrthoFinder orthogroups
orthogroups <- read.delim("Orthogroups.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
View(orthogroups)
# Use only columns 10 to 16 (i.e., columns 10:16 + the Orthogroup column)
welwitschiae_data <- orthogroups[, c(1, 10:16)]
head(welwitschiae_data)
# Optional: Rename strain columns for clarity
colnames(welwitschiae_data) <- c("Orthogroup", "CBS 139.54", 
                                 "CCMB 663", "CCMB 674", "IHEM 2864", 
                                 "IMRD1", "ITEM 11945", "AwOcstreb1")
colnames(welwitschiae_data)

# Create binary matrix
wel_matrix <- welwitschiae_data %>%
  mutate(across(-Orthogroup, ~ ifelse(.x == "" | is.na(.x), 0, 1))) %>%
  column_to_rownames("Orthogroup")
dim(wel_matrix)

# UpSet plot for presence/absence across strains
UpSetR::upset(as.data.frame(wel_matrix),
      sets = colnames(wel_matrix),
      nsets = 7,
      order.by = "freq",
      keep.order = TRUE,
      main.bar.color = "#377eb8",
      sets.bar.color = "#4daf4a")


# View how many unique OGs each strain has
sapply(strain_uniques, length)

# (Optional) Write the list for Endophyte strain
writeLines(as.character(strain_uniques$AwOcstreb1), "Endophyte_Unique_OGs.txt")

# Load list of endophyte-specific orthogroups
endophyte_ogs <- readLines("Endophyte_Unique_OGs.txt")
View(endophyte_ogs)

# Subset OrthoFinder table to only these orthogroups
endo_og_genes <- orthogroups %>%
  filter(Orthogroup %in% endophyte_ogs) %>%
  select(Orthogroup, Aspergillus_welwitschiae_ocstreb1) %>%
  separate_rows(Aspergillus_welwitschiae_ocstreb1, sep = ",\\s*")  # Split multiple genes in one cell
dim(endo_og_genes)
eggnog <- read.delim("eggnog_annotated_orthogroups.emapper.annotations", comment.char = "#", stringsAsFactors = FALSE)
# Merge with your endophyte-specific gene list
annotated_endo_genes <- endo_og_genes %>%
  rename(query = Aspergillus_welwitschiae_ocstreb1) %>%
  left_join(eggnog, by = c("query" = "query"))  # or "#query" if that's the column name in your file

# Check result
dim(annotated_endo_genes)
write.csv(annotated_endo_genes, "Annotated_Endophyte_Unique_OGs.csv", row.names = FALSE)
# Count top KEGG pathways
top_kegg <- annotated_endo_genes %>%
  filter(KEGG_ko != "" & !is.na(KEGG_ko)) %>%
  count(KEGG_ko, sort = TRUE)
top_kegg$KEGG_ko
# View top 10
head(top_kegg, 20)

# Barplot (optional)
library(ggplot2)
ggplot(top_kegg[1:10,], aes(x = reorder(KEGG_ko, n), y = n)) +
  geom_bar(stat = "identity", fill = "#377eb8") +
  coord_flip() +
  labs(title = "Top KEGG Pathways in Endophyte-unique Genes",
       x = "KEGG Pathway", y = "Gene Count")

## Add Presence_Count column before running the rest

wel_matrix$type <- case_when(
  wel_matrix$Presence_Count == 7 ~ "Core",
  wel_matrix$Presence_Count == 1 ~ "Unique",
  wel_matrix$Presence_Count < 7 & wel_matrix$Presence_Count > 1 ~ "Accessory"
)
wel_matrix <- wel_matrix[wel_matrix[, "Presence_Count"] != 0, ]
# Summary table
summary_table <- wel_matrix %>%
  count(type)

print(summary_table)

# Optional: pie chart
library(ggplot2)

ggplot(summary_table, aes(x = "", y = n, fill = type)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Pangenome Composition") +
  scale_fill_manual(values = c("Core" = "#1b9e77", "Accessory" = "#7570b3", "Unique" = "#d95f02"))



