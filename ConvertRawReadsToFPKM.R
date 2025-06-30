# Load the data
data <- read.table("ocstreb1_19.fcount.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
View(data)
head(data)
total_mapped_reads <- sum(data$ocstreb1_against_19.bam)
data$Length_kb <- data$Length / 1000
# Calculate FPKM
data$FPKM <- (data$ocstreb1_against_19.bam / data$Length_kb) / (total_mapped_reads / 1e6)
head(data)
# Preview result
head(data[, c("Geneid", "FPKM")])

# Optionally, write to a new file
write.table(data, file = "featurecounts_with_FPKM_ocstreb1_19.txt", sep = "\t", row.names = FALSE, quote = FALSE)
dim(data)

data18 <- read.table("featurecounts_with_FPKM_ocstreb1_18.txt", header = TRUE, sep = "\t")
data19 <- read.table("featurecounts_with_FPKM_ocstreb1_19.txt", header = TRUE, sep = "\t")

head(data18)
head(data19)
fpkm19 <- data19[,c("Geneid", "FPKM")]
colnames(fpkm19)[2] <- "fpkm_19"
head(fpkm19)
# Merge with the full first file (data18)
merged_data <- merge(data18, fpkm19, by = "Geneid")
head(merged_data)

# Calculate average FPKM
merged_data$FPKM_avg <- rowMeans(merged_data[, c("FPKM", "fpkm_19")], na.rm = TRUE)

##split them in three files based on expression level

# Define groups based on FPKM ranges
low_expr <- subset(merged_data, FPKM_avg >= 1 & FPKM_avg < 10)
moderate_expr <- subset(merged_data, FPKM_avg >= 10 & FPKM_avg<= 100)
high_expr <- subset(merged_data, FPKM_avg > 100)
dim(low_expr)
dim(moderate_expr)
dim(high_expr)

# Write to separate files
write.table(low_expr, file = "low_expressed_genes_18_19.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(moderate_expr, file = "moderate_expressed_genes_18_19.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(high_expr, file = "high_expressed_genes_18_19.txt", sep = "\t", row.names = FALSE, quote = FALSE)