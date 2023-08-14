# Preparing the data
names(MICH)
TS1 <- MICH %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,MICH_P_null_T1_R1:MICH_P_null_T1_R4)),
                Sample_2 = rowMeans(dplyr::select(.,MICH_P_null_T2_R1:MICH_P_null_T2_R2)),
                Sample_3 = rowMeans(dplyr::select(.,MICH_P_null_T3_R1:MICH_P_null_T3_R3)),
                Sample_4 = rowMeans(dplyr::select(.,MICH_P_null_T4_R1:MICH_P_null_T4_R2)),
  ) %>%
  dplyr::select(24:27)

# Selecting genes
controlvsP_MICH <- read.csv("result/ControlvsP_significant_genes_MICH.csv")
controlvsP_MICH <- read.csv("result/onlyP_significant_genes_MICH.csv")

TS1 <- TS1[controlvsP_MICH$x, ]
# Removing zeroes and replacing with lowest/2 value
TS1 <- TS1 %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1 <- log(TS1)

time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS1 <- rbind(time_points_row, TS1)
rownames(TS1)[1] <- "time_points"
TS1 <- as.matrix(TS1)


TS2 <- MICH %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,MICH_P_add_T1_R1:MICH_P_add_T1_R3)),
                Sample_2 = rowMeans(dplyr::select(.,MICH_P_add_T2_R1:MICH_P_add_T2_R3)),
                Sample_3 = rowMeans(dplyr::select(.,MICH_P_add_T3_R1:MICH_P_add_T3_R3)),
                Sample_4 = rowMeans(dplyr::select(.,MICH_P_add_T4_R1:MICH_P_add_T4_R3)),
  ) %>%
  dplyr::select(24:27)

# Selecting genes
controlvsP_MICH <- read.csv("result/ControlvsP_significant_genes_MICH.csv")
controlvsP_MICH <- read.csv("result/onlyP_significant_genes_MICH.csv")

TS2 <- TS2[controlvsP_MICH$x, ]
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

adjacency_matrix <- res$weight.matrix
library(igraph)
network <- graph.adjacency(adjacency_matrix, mode = "undirected", weighted = TRUE)
quartz()
plot(network, vertex.label = V(network)$name, vertex.size = 10, vertex.color = "lightblue")
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
