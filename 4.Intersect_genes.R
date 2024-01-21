library(VennDiagram)
##### FOR ALL
# Finding the common genes
B73_noP_all <- read.csv("result/all_analysis/B73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])
B73_PvsB73_noP <- read.csv("result/all_analysis/B73_PvsB73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])
MICH_noPvsB73_noP <- read.csv("result/all_analysis/MICH_noPvsB73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])
MICH_PvsB73_noP <- read.csv("result/all_analysis/MICH_PvsB73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])

head(B73_noP_all)
head(B73_PvsB73_noP)
head(MICH_noPvsB73_noP)
head(MICH_PvsB73_noP)



# Extract 'gene' columns
genes1 <- B73_noP_all$gene
genes2 <- B73_PvsB73_noP$gene
genes3 <- MICH_noPvsB73_noP$gene
genes4 <- MICH_PvsB73_noP$gene

# Compute intersections and unique parts manually
venn_list <- list(
  "B73_noP_all_ONLY" = setdiff(genes1, union(genes2, union(genes3, genes4))),
  "B73_PvsB73_noP_ONLY" = setdiff(genes2, union(genes1, union(genes3, genes4))),
  "MICH_noPvsB73_noP_ONLY" = setdiff(genes3, union(genes1, union(genes2, genes4))),
  "MICH_PvsB73_noP_ONLY" = setdiff(genes4, union(genes1, union(genes2, genes3))),
  "B73_noP_all_B73_PvsB73_noP" = intersect(genes1, genes2),
  "B73_noP_all_MICH_noPvsB73_noP" = intersect(genes1, genes3),
  "B73_noP_all_MICH_PvsB73_noP" = intersect(genes1, genes4),
  "B73_PvsB73_noP_MICH_noPvsB73_noP" = intersect(genes2, genes3),
  "B73_PvsB73_noP_MICH_PvsB73_noP" = intersect(genes2, genes4),
  "MICH_noPvsB73_noP_MICH_PvsB73_noP" = intersect(genes3, genes4),
  "ALL" = Reduce(intersect, list(genes1, genes2, genes3, genes4))
)

# Print the results
venn_list


# Compute the maximum length among all lists
max_length <- max(sapply(venn_list, length))

# Function to pad lists with NA to make their lengths equal to max_length
pad_with_na <- function(x) {
  c(x, rep(NA, max_length - length(x)))
}

# Apply the function to all elements of the list and combine into a dataframe
venn_df <- as.data.frame(lapply(venn_list, pad_with_na), stringsAsFactors = FALSE)

# Print the resulting dataframe
print(venn_df)

# Saving
write.csv(venn_df,"result/all_analysis/Venn_results.csv",row.names=F)

