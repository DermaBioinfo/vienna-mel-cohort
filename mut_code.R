library("readxl")
library("tidyverse")

mut <- c("BRAF", "NRAS", "NF1", "PTEN", "CDKN2A", "KIT", "TP53", "ATM", "CCND1", "CDK4", "CTNNB1", "FGFR2", "RB1", "APC", "MET", "SMO", "HNF1A")

s3 <- read.table("NIHMS362881-supplement-2.txt", header = T, sep = "\t")
s3 <- s3[, c(5, 6, 7, 8)]
s3 <- unique(subset(s3, s3[, 4] %in% mut))
s3 <- s3[, c(4, 3, 1, 2)]

s4 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S4", skip = 5, col_names = T))
s4 <- s4[, c(2, 7, 9, 10, 20)]
s4 <- unique(subset(s4, s4[, 1] %in% mut))

s6 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S6", skip = 5, col_names = T))
s6 <- s6[, c(2, 8, 10, 11, 19)]
s6 <- unique(subset(s6, s6[, 1] %in% mut))

s8 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S8", skip = 6, col_names = T))
s8 <- unique(s8[, c(2, 14, 16, 18)])
s8 <- subset(s8, s8[, 1] %in% mut)
s8$`Variant_Classification` <- "Fusion"
names(s8)[1] <- "Hugo_Symbol"
s8 <- s8[, c(1, 5, 4, 2, 3)]

s4A <- read.table("NIHMS393895-supplement-03.txt", header = T, sep = "\t", fill=T)
s4A <- s4A[, c(1, 7, 9, 10, 11, 24)]
s4A <- unique(subset(s4A, s4A[, 1] %in% mut))
s4A[, 5] <- NULL
names(s4A)[4] <- "Tumor_Seq_Allele"

s4B <- as.data.frame(read_excel("NIHMS393895-supplement-02.xlsx", sheet = "TABLE S4B", skip = 5, col_names = T))
s4B <- s4B[, c(1, 7, 9, 10, 24)]
s4B <- unique(subset(s4B, s4B[, 1] %in% mut))
names(s4B)[4] <- "Tumor_Seq_Allele"

s6_1 <- as.data.frame(read_excel("NIHMS393895-supplement-02.xlsx", sheet = "TABLE S6", skip = 5, col_names = T))
s6_1 <- s6_1[, c(2, 6, 9)]
s6_1 <- unique(subset(s6_1, s6_1[, 1] %in% mut))
names(s6_1)[2] <- "Variant_Classification"
names(s6_1)[3] <- "Protein_Change"
names(s6_1)[1] <- "Hugo_Symbol"

result1 <- rbind.fill(s4, s6, s4A, s4B, s3, s6_1)
result1 <- result1[!duplicated(result1[1:4]),]
result1 <- result1[!duplicated(result1[1, 2, 5]),]
result2 <- rbind.fill(result1, s8)
result2 <- arrange(result2, Hugo_Symbol, Variant_Classification)
write.xlsx(result2, "MUT_res.xlsx", sheetName = "mut", col.names = TRUE, row.names = TRUE, append = FALSE)