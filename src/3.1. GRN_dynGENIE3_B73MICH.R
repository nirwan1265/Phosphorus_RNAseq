# Loading the packages
## Note: These files have to be downloaded from their github
## https://github.com/vahuynh/dynGENIE3/tree/master/dynGENIE3_R_C_wrapper
system("ls")
system("R CMD SHLIB dynGENIE3.c")
source("dynGENIE3.R")



#Getting the data
B73MICH <- read.csv("data/MICHvsB73_noP.csv")
rownames(B73MICH) <- B73MICH[,1]
B73MCH <- B73MICH[,-1]
names(B73MICH)
# Preparing the data
TS1 <- B73MICH %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,B73_P_null_T1_R1:B73_P_null_T1_R1)),
                Sample_2 = rowMeans(dplyr::select(.,B73_P_null_T2_R1:B73_P_null_T2_R3)),
                Sample_3 = rowMeans(dplyr::select(.,B73_P_null_T3_R1:B73_P_null_T3_R2)),
                Sample_4 = rowMeans(dplyr::select(.,B73_P_null_T4_R1:B73_P_null_T4_R2)),
  ) %>%
  dplyr::select(21:24)



# Selecting genes
onlyP <- read.csv("result/onlyP_significant_genes_B73MCH.csv")
CvP <- read.csv("result/ControlvsP_significant_genes_B73MCH.csv")
all_P_CvP <- read.csv("result/all_significant_genes_B73MCH.csv")

onlyP <- read.csv("onlyP_significant_genes_B73MCH.csv")
CvP <- read.csv("ControlvsP_significant_genes_B73MCH.csv")
all_P_CvP <- read.csv("all_significant_genes_B73MCH.csv")

TS1_onlyP <- TS1[onlyP$gene, ]
TS1_CvP <- TS1[CvP$gene, ]
TS1_P_CvP <- TS1[all_P_CvP$gene, ]


#####
# Only P
#####

# Removing zeroes and replacing with lowest/2 value
TS1_onlyP <- TS1_onlyP %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1_onlyP <- log(TS1_onlyP)

time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS1_onlyP <- rbind(time_points_row, TS1_onlyP)
rownames(TS1_onlyP)[1] <- "time_points"
TS1_onlyP <- as.matrix(TS1_onlyP)

#####
# C vs P
#####

# Removing zeroes and replacing with lowest/2 value
TS1_CvP <- TS1_CvP %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1_CvP <- log(TS1_CvP)

time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS1_CvP <- rbind(time_points_row, TS1_CvP)
rownames(TS1_CvP)[1] <- "time_points"
TS1_CvP <- as.matrix(TS1_CvP)


#####
# All genes
#####

# Removing zeroes and replacing with lowest/2 value
TS1_P_CvP <- TS1_P_CvP %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1_P_CvP <- log(TS1_P_CvP)

time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS1_P_CvP <- rbind(time_points_row, TS1_P_CvP)
rownames(TS1_P_CvP)[1] <- "time_points"
TS1_P_CvP <- as.matrix(TS1_P_CvP)

# Second condition
TS2 <- B73MICH %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,MICH_P_null_T1_R1:MICH_P_null_T1_R4)),
                Sample_2 = rowMeans(dplyr::select(.,MICH_P_null_T2_R1:MICH_P_null_T2_R2)),
                Sample_3 = rowMeans(dplyr::select(.,MICH_P_null_T3_R1:MICH_P_null_T3_R3)),
                Sample_4 = rowMeans(dplyr::select(.,MICH_P_null_T4_R1:MICH_P_null_T4_R2)),
  ) %>%
  dplyr::select(21:24)

TS2_onlyP <- TS2[onlyP$gene, ]
TS2_CvP <- TS2[CvP$gene, ]
TS2_P_CvP <- TS2[all_P_CvP$gene, ]

#####
# Only P
#####
# Removing zeroes and replacing with lowest/2 value
TS2_onlyP <- TS2_onlyP %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS2_onlyP <- log(TS2_onlyP)

time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS2_onlyP <- rbind(time_points_row, TS2_onlyP)
rownames(TS2_onlyP)[1] <- "time_points"
TS2_onlyP <- as.matrix(TS2_onlyP)


#####
# Only P
#####

# Removing zeroes and replacing with lowest/2 value
TS2_CvP <- TS2_CvP %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS2_CvP <- log(TS2_CvP)


time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS2_CvP <- rbind(time_points_row, TS2_CvP)
rownames(TS2_CvP)[1] <- "time_points"
TS2_CvP <- as.matrix(TS2_CvP)

#####
# All sig genes
#####

# Removing zeroes and replacing with lowest/2 value
TS2_P_CvP <- TS2_P_CvP %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS2_P_CvP <- log(TS2_P_CvP)


time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)
TS2_P_CvP <- rbind(time_points_row, TS2_P_CvP)
rownames(TS2_P_CvP)[1] <- "time_points"
TS2_P_CvP <- as.matrix(TS2_P_CvP)


##### Running the GRN
#### Only P
time.points_onlyP <- list(TS1_onlyP[1,],TS2_onlyP[1,])
TS.data_onlyP <- list(TS1_onlyP[2:nrow(TS1_onlyP), ],TS2_onlyP[2:nrow(TS2_onlyP), ])

res_onlyP <- dynGENIE3(TS.data_onlyP, time.points_onlyP)
link.list_onlyP <- get.link.list(res_onlyP$weight.matrix)
res_onlyP$weight.matrix

# Gene degradation rates
decay_rates_onlyP <- 0.5
# Or alternatively
gene.names_onlyP <- rownames(TS.data_onlyP[[1]])
decay_rates_onlyP <- rep(0.02, length(gene.names_onlyP))
decay_rates_onlyP <- setNames(decay_rates_onlyP, gene.names_onlyP)

# Run dynGENIE3
res_onlyP2 <- dynGENIE3(TS.data_onlyP, time.points_onlyP, alpha=decay_rates_onlyP)
res_onlyP2 <- get.link.list(res_onlyP2$weight.matrix)
head(res_onlyP2)



#### CvP
time.points_CvP <- list(TS1_CvP[1,],TS2_CvP[1,])
TS.data_CvP <- list(TS1_CvP[2:nrow(TS1_CvP), ],TS2_CvP[2:nrow(TS2_CvP), ])

res_CvP <- dynGENIE3(TS.data_CvP, time.points_CvP)
link.list_CvP <- get.link.list(res_CvP$weight.matrix)
res_CvP$weight.matrix

# Gene degradation rates
decay_rates_CvP <- 0.02
# Or alternatively
gene.names_CvP <- rownames(TS.data_CvP[[1]])
decay_rates_CvP <- rep(0.02, length(gene.names_CvP))
decay_rates_CvP <- setNames(decay_rates_CvP, gene.names_CvP)

# Run dynGENIE3
res_CvP2 <- dynGENIE3(TS.data_CvP, time.points_CvP, alpha=decay_rates_CvP)
res_CvP2 <- get.link.list(res_CvP2$weight.matrix)
head(res_CvP2)


### All sig genes

time.points_P_CvP <- list(TS1_P_CvP[1,],TS2_P_CvP[1,])
TS.data_P_CvP <- list(TS1_P_CvP[2:nrow(TS1_P_CvP), ],TS2_P_CvP[2:nrow(TS2_P_CvP), ])

res_P_CvP <- dynGENIE3(TS.data_P_CvP, time.points_P_CvP)
link.list_P_CvP <- get.link.list(res_P_CvP$weight.matrix)
res_P_CvP$weight.matrix

# Gene degradation rates
decay_rates_P_CvP <- 0.02
# Or alternatively
gene.names_P_CvP <- rownames(TS.data_P_CvP[[1]])
decay_rates_P_CvP <- rep(0.02, length(gene.names_P_CvP))
decay_rates_P_CvP <- setNames(decay_rates_P_CvP, gene.names_P_CvP)

# Run dynGENIE3
res_P_CvP2 <- dynGENIE3(TS.data_P_CvP, time.points_P_CvP, alpha=decay_rates_P_CvP)
res_P_CvP2 <- get.link.list(res_P_CvP2$weight.matrix)
head(res_P_CvP2)


# Set a threshold for strong links
threshold <- 0.001  # Adjust this threshold value as needed


# Convert values below the threshold to 0, keeping only strong links
adjacency_matrix <- as.matrix(res_onlyP$weight.matrix)
adjacency_matrix[adjacency_matrix <= threshold] <- 0

# Create a graph from the adjacency matrix
graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")

# Layout
LO <- layout_with_fr(graph)

# Removing empty nodes
isolated <- which(igraph::degree(graph) == 0)
graph <- delete.vertices(graph, isolated)
LO <- LO[-isolated, ]

quartz()
plot(graph)
