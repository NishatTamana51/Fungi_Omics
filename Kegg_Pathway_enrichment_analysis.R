# ==== Load libraries ====
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
# 1. Install KEGGREST if not already installed
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("KEGGREST")
}

# 2. Load the KEGGREST library
library(KEGGREST)

#Load gene-to-KEGG data
# Assume: columns = Gene_ID, KEGG_ko (with ko/map IDs, comma-separated)
gene2kegg <- read.csv("kos_endophyte_all_og_genes.csv", stringsAsFactors = FALSE)

#Load expressed genes list
expressed_genes <- read.csv("kos_endophyte_unique_og_genes.csv", stringsAsFactors = FALSE)$GENE_ID
head(expressed_genes)

#Split KEGG annotations into rows
kegg_long <- gene2kegg %>%
  filter(!is.na(KEGG_ko) & KEGG_ko != "-") %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(
    KEGG_ko = trimws(KEGG_ko),
    Map_ID = str_extract(KEGG_ko, "map\\d{5}"),
    Status = if_else(GENE_ID %in% expressed_genes, "Expressed", "Unexpressed")
  ) %>%
  filter(!is.na(Map_ID))
View(kegg_long)
table(kegg_long$Status)

#Count gene usage per KEGG map ID
kegg_summary <- kegg_long %>%
  group_by(Map_ID, Status) %>%
  summarise(Gene_Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Status, values_from = Gene_Count, values_fill = 0) %>%
  mutate(Total = Expressed + Unexpressed)

#Get KEGG map descriptions
kegg_list_all <- keggList("pathway", "ko")
head(kegg_list_all)
kegg_df <- data.frame(
  Map_ID = str_extract(names(kegg_list_all), "\\d{5}"),
  Pathway_Name = as.character(kegg_list_all),
  stringsAsFactors = FALSE
)
View(kegg_df)
#Tag and remove broad/general categories|| As kegg always include some broad categories, I chose to skip them 
broad_terms <- c(
  "Global and overview", "Metabolic pathways",
  "Biosynthesis of secondary metabolites",
  "Biosynthesis of amino acids",
  "Microbial metabolism in diverse environments",
  "Overview", "General", "Miscellaneous"
)
# I also filtered disease terms which are not related to our organism (fungus) by choosing the terms from initial enriched list
disease_terms <- c("Huntington", "Parkinson", "Alzheimer", "Cancer", "Disease", "virus")
filtered_kegg <- kegg_df %>%
  filter(!str_detect(Pathway_Name, paste(disease_terms, collapse = "|")))

kegg_df <- filtered_kegg %>%
  mutate(IsBroad = str_detect(Pathway_Name, paste(broad_terms, collapse = "|")))

# Strip 'map' prefix from Map_ID before joining
kegg_summary <- kegg_summary %>%
  mutate(Map_ID = str_replace(Map_ID, "^map", ""))
View(kegg_df)

# Merge with pathway summary & remove broad pathways
kegg_cleaned <- left_join(kegg_summary, kegg_df, by = "Map_ID") %>%
  filter(!IsBroad)
View(kegg_cleaned)

# Extract filtered Map_IDs
valid_maps <- unique(kegg_cleaned$Map_ID)
head(valid_maps)

# Filter kegg_long using valid (cleaned) Map_IDs
head(kegg_long)
kegg_long <- kegg_long %>%
  mutate(Map_ID = str_replace(Map_ID, "^map", ""))

kegg_filtered_long <- kegg_long %>%
  filter(Map_ID %in% valid_maps)
head(kegg_filtered_long)

# Now build TERM2GENE and TERM2NAME
TERM2GENE <- kegg_filtered_long %>%
  select(Map_ID, Gene_ID) %>%
  distinct()

TERM2NAME <- kegg_cleaned %>%
  select(Map_ID, Pathway_Name) %>%
  distinct()

# Define gene universe
all_genes <- unique(kegg_filtered_long$Gene_ID)
expressed <- expressed_genes
unexpressed <- setdiff(all_genes, expressed)

# Run enrichment
kegg_enrich_expressed <- enricher(
  gene = expressed,
  universe = all_genes,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

write.csv(as.data.frame(kegg_enrich_expressed), "kegg_enrichment_expressed.csv", row.names = FALSE)

#Convert to dataframe
enrich_df <- as.data.frame(kegg_enrich_expressed)

#Add expressed and unexpressed gene counts for each pathway
pathway_counts <- kegg_filtered_long %>%
  group_by(Map_ID, Expression = if_else(Gene_ID %in% expressed, "Expressed", "Unexpressed")) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Expression, values_from = Count, values_fill = 0)

#Merge with enrichment result
plot_data <- enrich_df %>%
  rename(Map_ID = ID) %>%
  left_join(pathway_counts, by = "Map_ID") %>%
  mutate(Total = Expressed + Unexpressed) %>%
  slice_max(order_by = Count, n = 20)  # Optional: top 20

#Plot
ggplot(plot_data, aes(x = reorder(Description, Total), y = Total)) +
  geom_bar(aes(fill = "Total"), stat = "identity", alpha = 0.3) +
  geom_bar(aes(y = Expressed, fill = "Expressed"), stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Expressed" = "#330099", "Total" = "#660000")) +
  labs(
    x = "KEGG Pathway",
    y = "Gene Count",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10)
  )+
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.title.position = "plot"
  )
