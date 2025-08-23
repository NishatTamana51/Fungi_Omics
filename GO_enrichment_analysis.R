BiocManager::install("topGO")
library(topGO)
library(dplyr)
library(tidyr)
library(GO.db)
library(ggplot2)
library(stringr)
# GO file with just two columns GENE_ID and GO_IDs
gene2go_raw <- read.csv("gos_endophyte_all_og_genes.csv", stringsAsFactors = FALSE)
# Create vector of GO terms (GO_IDs column contains multiple GO separated by comma)
gene2go <- setNames(
  lapply(strsplit(gene2go_raw$GO_IDs, ","), trimws),
  gene2go_raw$GENE_ID
)
# Expressed genes of interest
genes_of_interest <- read.csv("gos_endophyte_unique_og_genes.csv", header = TRUE)$GENE_ID
View(genes_of_interest)
# All genes as background
all_genes <- unique(names(gene2go))
View(all_genes)
# Create factor indicating whether each gene is in the gene set
# Make sure 'geneList' is a factor with exactly two levels: 0 (not selected), 1 (selected)
geneList <- factor(as.integer(all_genes %in% genes_of_interest))
names(geneList) <- all_genes
table(geneList)
# Sanity check
stopifnot(all(names(geneList) %in% names(gene2go)))

# ===== Create topGO object =====
GOdata <- new("topGOdata",
              ontology = "CC",  # Change to "BP", "CC", or "MF" as needed
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.gene2GO,
              gene2GO = gene2go,
              nodeSize = 5)  # Filter small GO terms (<5 genes)
GOdata1 <- new("topGOdata",
               ontology = "BP",
               allGenes = geneList,
               geneSelectionFun = function(x)(x == 1),
               annot = annFUN.gene2GO,
               gene2GO = gene2go)

GOdata2 <- new("topGOdata",
               ontology = "MF",
               allGenes = geneList,
               geneSelectionFun = function(x)(x == 1),
               annot = annFUN.gene2GO,
               gene2GO = gene2go)
res_CC <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
res_BP <- runTest(GOdata1, algorithm = "weight01", statistic = "fisher")
res_MF <- runTest(GOdata2, algorithm = "weight01", statistic = "fisher")

# ===== Generate tables =====
top_n <- 20

tab_CC <- GenTable(GOdata, weightFisher = res_CC, topNodes = top_n)
tab_BP <- GenTable(GOdata1, weightFisher = res_BP, topNodes = top_n)
tab_MF <- GenTable(GOdata2, weightFisher = res_MF, topNodes = top_n)

# ===== Save CSVs =====
write.csv(tab_CC, "endo_unique_GO_CC_top20.csv", row.names = FALSE)
write.csv(tab_BP, "endo_unique_GO_BP_top20.csv", row.names = FALSE)
write.csv(tab_MF, "endo_unique_GO_MF_top20.csv", row.names = FALSE)

# Step 1: Load Data
bp <- read.csv("endo_unique_GO_BP_top20.csv") %>% mutate(Ontology = "Biological Process")
mf <- read.csv("endo_unique_GO_MF_top20.csv") %>% mutate(Ontology = "Molecular Function")
cc <- read.csv("endo_unique_GO_CC_top20.csv") %>% mutate(Ontology = "Cellular Component")

# Step 2: Combine and prepare data
all_go <- bind_rows(bp, mf, cc) %>%
  mutate(
    GO.Term = str_trunc(Term, 60),  # Truncate long GO terms
    logP = -log10(as.numeric(weightFisher)),
    Ontology = as.character(Ontology)  # Convert to character to match manual color names
  ) %>%
  filter(weightFisher < 0.05) %>%  # Filter only significant terms (p < 0.05)
  group_by(Ontology) %>%
  slice_max(order_by = logP, n = 15) %>%
  ungroup()

# Step 3: Order terms for plotting (by ontology and logP)
all_go <- all_go %>%
  arrange(Ontology, logP) %>%
  mutate(GO.Term = factor(GO.Term, levels = rev(unique(GO.Term))))  # y-axis order

# Step 4: Dot plot
ggplot(all_go, aes(x = logP, y = GO.Term, color = Ontology, size = Significant)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Biological Process" = "#1f77b4",
      "Molecular Function" = "#ff7f0e",
      "Cellular Component" = "#2ca02c"
    ),
    breaks = c(
      "Biological Process",
      "Molecular Function",
      "Cellular Component"
    )
  ) +
  labs(
    title = "(a)",
    x = "-log10(p-value)",
    y = NULL,
    size = "Gene Count",
    color = "Ontology"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )