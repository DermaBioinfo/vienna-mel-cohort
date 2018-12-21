library("readxl")
library("plyr")
library("xlsx")

mut <- c("BRAF", "NRAS", "NF1")

s6 <- as.data.frame(read_excel("NIHMS393895-supplement-02.xlsx", sheet = "TABLE S6", skip = 5, col_names = T))
s6 <- s6[, c(1, 2)]
s6 <- subset(s6, grepl("ME0", s6[, 1]) & (s6[, 2] %in% mut))
s6[, c("BRAF", "NRAS", "NF1")] <- c(grepl("BRAF", s6[, 2]), grepl("NRAS", s6[, 2]), grepl("NF1", s6[, 2]))
s6[, 2] <- NULL
s6 <- aggregate(s6, by=list(s6[, 1]), FUN = any)
s6[, 2] <- NULL
names(s6)[1] <- "Sample"

s4 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S4", skip = 5, col_names = T))
s4 <- s4[, c(1, 2)]
s4 <- subset(s4, grepl("ME0", s4[, 1]) & (s4[, 2] %in% mut))
s4[, c("BRAF", "NRAS", "NF1")] <- c(grepl("BRAF", s4[, 2]), grepl("NRAS", s4[, 2]), grepl("NF1", s4[, 2]))
s4[, 2] <- NULL
s4 <- aggregate(s4, by=list(s4[, 1]), FUN = any)
s4[, 2] <- NULL
names(s4)[1] <- "Sample"

s10 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S10", skip = 6, col_names = T))
s10 <- s10[, c(1, 6, 7)]
s10 <- subset(s10, grepl("ME0", s10[, 1]))
s10[, 1] <- gsub("T", "", s10[,1])
s10[, c(2, 3)] <- !is.na(s10[, c(2, 3)])
names(s10)[c(2, 3)] <- c("BRAF", "NRAS")
s10$NF1 <- F

s1 <- as.data.frame(read_excel("NIHMS362881-supplement-3.xlsx", sheet = "TABLE S1", skip = 6, col_names = T))
s1 <- s1[, c(1, 16, 17)]
s1[, 1] <- gsub("T", "", s1[,1])
s1[, c(2, 3)] <- c(!grepl("wild type", s1[, 2]), !grepl("wild type", s1[, 3]))
names(s1)[c(2, 3)] <- c("BRAF", "NRAS")
s1$NF1 <- F

result <- aggregate(. ~ Sample, rbind.fill(s1, s4, s6, s10), FUN = any)
result$twt <- !(result$BRAF | result$NRAS | result$NF1)
write.xlsx(result, "TWT_res.xlsx", sheetName = "twt", col.names = TRUE, row.names = TRUE, append = FALSE)