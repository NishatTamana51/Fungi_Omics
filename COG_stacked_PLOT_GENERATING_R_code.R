# Load libraries
library(tidyverse)

df <- data.frame(
  COG_class = c("Function unknown","Carbohydrate transport and metabolism",
                "Secondary metabolites biosynthesis, transport and catabolism", "Amino acid transport and metabolism","Posttranslational modification, protein turnover, chaperones",
                "Transcription",
                "Energy production and conversion",
                "Intracellular trafficking, secretion, and vesicular transport",
                "Lipid transport and metabolism",
                "Translation, ribosomal structure and biogenesis",
                "Signal transduction mechanisms",
                "Inorganic ion transport and metabolism",
                "RNA processing and modification",
                "Coenzyme transport and metabolism",
                "Replication, recombination and repair",
                "Cell cycle control, cell division, chromosome partitioning",
                "Cell wall/membrane/envelope biogenesis",
                "Cytoskeleton",
                "Nucleotide transport and metabolism","Chromatin structure and dynamics","Defense mechanisms", "Nuclear structure", "Extracellular structures", "Cell motility"),
  Predicted = c(2989, 802, 724, 662, 610, 549, 529, 507, 464, 402, 399, 323, 298, 298, 289, 162, 143, 140, 136, 110, 65, 24, 6, 5),
  Active_Expression = c(67.82, 58.10, 45.03, 64.95, 82.46, 80.33, 65.97, 85.21, 65.73, 93.03, 82.96, 73.99, 92.95, 70.47, 78.55, 90.74, 58.74, 90.0, 76.47, 89.09, 50.77, 100.0, 66.67, 100.0),
  High = c(4.08, 6.23, 3.45, 8.91, 11.97, 4.01, 15.69, 7.30, 6.47, 28.61, 7.52, 6.19, 4.70, 9.06, 0.69, 5.56, 5.59, 10.71, 9.56, 4.55, 1.54, 0.0, 0.0, 20.0),
  Moderate = c(27.80, 21.95, 15.88, 30.06, 49.34, 39.16, 25.90, 56.02, 30.60, 38.56, 52.38, 37.77, 43.96, 29.53, 26.30, 56.17, 24.48, 55.00, 40.44, 41.82, 21.54, 91.67, 50.0, 40.0),
  Low = c(35.93, 29.93, 25.69, 25.98, 21.15, 37.16, 24.39, 21.89, 28.66, 25.87, 23.06, 30.03, 44.30, 31.88, 51.56, 29.01, 28.67, 24.29, 26.47, 42.73, 27.69, 8.33, 16.67, 40.0)
)

# Pivot to long format
df_long <- df %>%
  pivot_longer(cols = c(High, Moderate, Low), names_to = "Expression_Level", values_to = "Percentage")

# Set factor order for coloring
df_long$Expression_Level <- factor(df_long$Expression_Level, levels = c("Low", "Moderate", "High"))

# Merge predicted and active values for annotations
df_annotate <- df %>%
  mutate(BarPosition = 103, LabelPosition = 96)

# Plot
ggplot(df_long, aes(x = COG_class, y = Percentage, fill = Expression_Level)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(data = df_annotate,
            aes(x = COG_class, y = BarPosition, label = Predicted),
            size = 3, angle = 45, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = df_annotate,
            aes(x = COG_class, y = LabelPosition, label = paste0(round(Active_Expression, 1), "%")),
            color = "red", size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Low" = "#440154", "Moderate" = "#31688e", "High" = "#35b779")) +
  labs(
    title = "Expression Level Distribution with Predicted Gene Count and Active Expression",
    x = "COG Category",
    y = "Percentage of Expressed Genes",
    fill = "Expression Level"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylim(0, 110)

