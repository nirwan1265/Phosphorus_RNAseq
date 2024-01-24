library(VennDiagram)
##### FOR ALL
# Finding the common genes
B73_noP_all <- read.csv("result/all_analysis/B73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])
B73_PvsB73_noP <- read.csv("result/all_analysis/B73_PvsB73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])
MICH_noPvsB73_noP <- read.csv("result/all_analysis/MICH_noPvsB73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])
MICH_PvsB73_noP <- read.csv("result/all_analysis/MICH_PvsB73_noP.csv") %>% dplyr::select(c(1,2)) %>% dplyr::arrange(.[,2])

# Extract unique values from each vector
B73_noP_all <- unique(B73_noP_all$gene)
B73_PvsB73_noP <- unique(B73_PvsB73_noP$gene)
MICH_noPvsB73_noP <- unique(MICH_noPvsB73_noP$gene)
MICH_PvsB73_noP <- unique(MICH_PvsB73_noP$gene)

# 1. Common to all four
common_all <- Reduce(intersect, list(B73_noP_all, B73_PvsB73_noP, MICH_noPvsB73_noP, MICH_PvsB73_noP))

# 2-5. Triple combinations
common_A_B_C <- setdiff(intersect(intersect(B73_noP_all, B73_PvsB73_noP), MICH_noPvsB73_noP), MICH_PvsB73_noP)
common_A_B_D <- setdiff(intersect(intersect(B73_noP_all, B73_PvsB73_noP), MICH_PvsB73_noP), MICH_noPvsB73_noP)
common_A_C_D <- setdiff(intersect(intersect(B73_noP_all, MICH_noPvsB73_noP), MICH_PvsB73_noP), B73_PvsB73_noP)
common_B_C_D <- setdiff(intersect(intersect(B73_PvsB73_noP, MICH_noPvsB73_noP), MICH_PvsB73_noP), B73_noP_all)

# 6-11. Pair combinations
common_A_B <- setdiff(intersect(B73_noP_all, B73_PvsB73_noP), union(MICH_noPvsB73_noP, MICH_PvsB73_noP))
common_A_C <- setdiff(intersect(B73_noP_all, MICH_noPvsB73_noP), union(B73_PvsB73_noP, MICH_PvsB73_noP))
common_A_D <- setdiff(intersect(B73_noP_all, MICH_PvsB73_noP), union(B73_PvsB73_noP, MICH_noPvsB73_noP))
common_B_C <- setdiff(intersect(B73_PvsB73_noP, MICH_noPvsB73_noP), union(B73_noP_all, MICH_PvsB73_noP))
common_B_D <- setdiff(intersect(B73_PvsB73_noP, MICH_PvsB73_noP), union(B73_noP_all, MICH_noPvsB73_noP))
common_C_D <- setdiff(intersect(MICH_noPvsB73_noP, MICH_PvsB73_noP), union(B73_noP_all, B73_PvsB73_noP))

# 12-15. Unique elements
unique_A <- setdiff(B73_noP_all, union(B73_PvsB73_noP, union(MICH_noPvsB73_noP, MICH_PvsB73_noP)))
unique_B <- setdiff(B73_PvsB73_noP, union(B73_noP_all, union(MICH_noPvsB73_noP, MICH_PvsB73_noP)))
unique_C <- setdiff(MICH_noPvsB73_noP, union(B73_noP_all, union(B73_PvsB73_noP, MICH_PvsB73_noP)))
unique_D <- setdiff(MICH_PvsB73_noP, union(B73_noP_all, union(B73_PvsB73_noP, MICH_noPvsB73_noP)))


# Just the MICH null and P
MICH <- c(unique_C,unique_D,common_C_D)
MICH <- data.frame(MICH)

write.csv(MICH,"result/all_analysis/MICH_only.csv", row.names=F)
