#!/bin/bash
data <- read.table("combined_info_scores.txt", header = FALSE, sep = " ")
filtered_data <- data[data$V7 < 0.3, ]
result <- filtered_data$V2

write.table(result, file = "eli_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
