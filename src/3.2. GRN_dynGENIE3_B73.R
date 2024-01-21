# Loading the packages
## Note: These files have to be downloaded from their github
## https://github.com/vahuynh/dynGENIE3/tree/master/dynGENIE3_R_C_wrapper
system("ls")
system("R CMD SHLIB dynGENIE3.c")
source("dynGENIE3.R")

# Preparing the data
TS1 <- B73 %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,B73_P_null_T1_R1:B73_P_null_T1_R1)),
                Sample_2 = rowMeans(dplyr::select(.,B73_P_null_T2_R1:B73_P_null_T2_R3)),
                Sample_3 = rowMeans(dplyr::select(.,B73_P_null_T3_R1:B73_P_null_T3_R2)),
                Sample_4 = rowMeans(dplyr::select(.,B73_P_null_T4_R1:B73_P_null_T4_R2)),
                ) %>%
  dplyr::select(21:24)

# Selecting genes
controlvsP_B73 <- read.csv("result/ControlvsP_significant_genes_B73.csv")
controlvsP_B73 <- read.csv("result/onlyP_significant_genes_B73.csv")
TS1 <- TS1[controlvsP_B73$x, ]
# Removing zeroes and replacing with lowest/2 value
TS1 <- TS1 %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1 <- log(TS1)

time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS1 <- rbind(time_points_row, TS1)
rownames(TS1)[1] <- "time_points"
TS1 <- as.matrix(TS1)


TS2 <- B73 %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,B73_P_add_T1_R1:B73_P_add_T1_R3)),
                Sample_2 = rowMeans(dplyr::select(.,B73_P_add_T2_R1:B73_P_add_T2_R3)),
                Sample_3 = rowMeans(dplyr::select(.,B73_P_add_T3_R1:B73_P_add_T3_R3)),
                Sample_4 = rowMeans(dplyr::select(.,B73_P_add_T4_R1:B73_P_add_T4_R3)),
  ) %>%
  dplyr::select(21:24)

# Selecting genes
controlvsP_B73 <- read.csv("result/ControlvsP_significant_genes_B73.csv")
controlvsP_B73 <- read.csv("result/onlyP_significant_genes_B73.csv")
TS2 <- TS2[controlvsP_B73$x, ]
# Removing zeroes and replacing with lowest/2 value
TS2 <- TS2 %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS2 <- log(TS2)

time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS2 <- rbind(time_points_row, TS2)
rownames(TS2)[1] <- "time_points"
TS2 <- as.matrix(TS2)


time.points <- list(TS1[1,],TS2[1,])
TS.data <- list(TS1[2:nrow(TS1), ],TS2[2:nrow(TS2), ])

res <- dynGENIE3(TS.data, time.points)

link.list <- get.link.list(res$weight.matrix)
head(link.list)
res$weight.matrix

# Gene degradation rates
decay_rates <- 0.02
# Or alternatively
gene.names <- rownames(TS.data[[1]])
decay_rates <- rep(0.02, length(gene.names))
decay_rates <- setNames(decay_rates, gene.names)

# Run dynGENIE3
res2 <- dynGENIE3(TS.data, time.points, alpha=decay_rates)
link.list <- get.link.list(res2$weight.matrix)
head(link.list)
