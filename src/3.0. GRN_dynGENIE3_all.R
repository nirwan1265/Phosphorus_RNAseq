# Loading the packages
## Note: These files have to be downloaded from their github
## https://github.com/vahuynh/dynGENIE3/tree/master/dynGENIE3_R_C_wrapper
library(igraph)
system("ls")
system("R CMD SHLIB dynGENIE3.c")
source("dynGENIE3.R")


#Getting the data
data <- read.csv("data/B73&MICH.csv")
head(data)
rownames(data) <- data[,1]
data <- data[,-1]
data[1:4,1:5]
names(data)

# Preparing the data
TS1_B73_null <- data %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,B73_P_null_T1_R1:B73_P_null_T1_R1)),
                Sample_2 = rowMeans(dplyr::select(.,B73_P_null_T2_R1:B73_P_null_T2_R3)),
                Sample_3 = rowMeans(dplyr::select(.,B73_P_null_T3_R1:B73_P_null_T3_R2)),
                Sample_4 = rowMeans(dplyr::select(.,B73_P_null_T4_R1:B73_P_null_T4_R2)),
  ) %>%
  dplyr::select(44:47)

TS1_B73_P <- data %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,B73_P_add_T1_R1:B73_P_add_T1_R3)),
                Sample_2 = rowMeans(dplyr::select(.,B73_P_add_T2_R1:B73_P_add_T2_R3)),
                Sample_3 = rowMeans(dplyr::select(.,B73_P_add_T3_R1:B73_P_add_T3_R3)),
                Sample_4 = rowMeans(dplyr::select(.,B73_P_add_T4_R1:B73_P_add_T4_R3)),
  ) %>%
  dplyr::select(44:47)

TS1_MICH_null <- data %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,MICH_P_null_T1_R1:MICH_P_null_T1_R4)),
                Sample_2 = rowMeans(dplyr::select(.,MICH_P_null_T2_R1:MICH_P_null_T2_R2)),
                Sample_3 = rowMeans(dplyr::select(.,MICH_P_null_T3_R1:MICH_P_null_T3_R3)),
                Sample_4 = rowMeans(dplyr::select(.,MICH_P_null_T4_R1:MICH_P_null_T4_R2)),
  ) %>%
  dplyr::select(44:47)

TS1_MICH_P <- data %>%
  dplyr::mutate(Sample_1 = rowMeans(dplyr::select(.,MICH_P_add_T1_R1:MICH_P_add_T1_R3)),
                Sample_2 = rowMeans(dplyr::select(.,MICH_P_add_T2_R1:MICH_P_add_T2_R3)),
                Sample_3 = rowMeans(dplyr::select(.,MICH_P_add_T3_R1:MICH_P_add_T3_R3)),
                Sample_4 = rowMeans(dplyr::select(.,MICH_P_add_T4_R1:MICH_P_add_T4_R3)),
  ) %>%
  dplyr::select(44:47)



# Selecting differentiated genes
B73_null <- read.csv("result/all_analysis/B73_noP.csv")
B73_P <- read.csv("result/all_analysis/B73_PvsB73_noP.csv")
MICH_null <- read.csv("result/all_analysis/MICH_noPvsB73_noP.csv")
MICH_P <- read.csv("result/all_analysis/MICH_PvsB73_noP.csv")


# Subsetting the important genes only
# For MICH_null
#TS1_B73_null <- TS1_B73_null[MICH_null$gene, ]
#TS1_B73_P <- TS1_B73_P[MICH_null$gene, ]
#TS1_MICH_null <- TS1_MICH_null[MICH_null$gene, ]
#TS1_MICH_P <- TS1_MICH_null[MICH_null$gene, ]

# For B73_null
TS1_B73_null <- TS1_B73_null[B73_null$gene, ]
TS1_B73_P <- TS1_B73_P[B73_null$gene, ]
TS1_MICH_null <- TS1_MICH_null[B73_null$gene, ]
TS1_MICH_P <- TS1_MICH_null[B73_null$gene, ]


unique_MICH_null <- c("Zm00001eb084030",
"Zm00001eb342100",
"Zm00001eb170400",
"Zm00001eb197050",
"Zm00001eb326430",
"Zm00001eb011330",
"Zm00001eb168200",
"Zm00001eb108340",
"Zm00001eb196430",
"Zm00001eb293200",
"Zm00001eb291010",
"Zm00001eb303340")

MICH_nullP_highR <- c("Zm00001eb190240",
                      "Zm00001eb191790",
                      "Zm00001eb191820",
                      "Zm00001eb188750",
                      "Zm00001eb188070",
                      "Zm00001eb189350",
                      "Zm00001eb189910",
                      "Zm00001eb196780",
                      "Zm00001eb196760",
                      "Zm00001eb193260",
                      "Zm00001eb193620",
                      "Zm00001eb192850",
                      "Zm00001eb334460",
                      "Zm00001eb241410",
                      "Zm00001eb196600",
                      "Zm00001eb190350",
                      "Zm00001eb264460"
                      )

MICH_only <- c("Zm00001eb197050","Zm00001eb170400","Zm00001eb108340","Zm00001eb084030",
               "Zm00001eb342100","Zm00001eb196430","Zm00001eb291010","Zm00001eb011330",
               "Zm00001eb168200","Zm00001eb326430","Zm00001eb361560","Zm00001eb293200",
               "Zm00001eb303340","Zm00001eb298230","Zm00001eb351110","Zm00001eb140380",
               "Zm00001eb337150","Zm00001eb090640","Zm00001eb304570","Zm00001eb154770",
               "Zm00001eb250820","Zm00001eb298910","Zm00001eb377030","Zm00001eb407220",
               "Zm00001eb264460","Zm00001eb334460","Zm00001eb186400","Zm00001eb241410",
               "Zm00001eb187890","Zm00001eb187280","Zm00001eb249670","Zm00001eb189710",
               "Zm00001eb189580","Zm00001eb186390","Zm00001eb189570","Zm00001eb189820",
               "Zm00001eb060540","Zm00001eb213070","Zm00001eb190750","Zm00001eb308810",
               "Zm00001eb196880",
               "Zm00001eb191990",
               "Zm00001eb196100",
               "Zm00001eb192330",
               "Zm00001eb194380",
               "Zm00001eb112470",
               "Zm00001eb189640",
               "Zm00001eb276970",
               "Zm00001eb189700",
               "Zm00001eb195790",
               "Zm00001eb186860",
               "Zm00001eb192580",
               "Zm00001eb000540",
               "Zm00001eb186610",
               "Zm00001eb196900",
               "Zm00001eb048840",
               "Zm00001eb186980",
               "Zm00001eb194200",
               "Zm00001eb186790",
               "Zm00001eb193120",
               "Zm00001eb187870",
               "Zm00001eb249540",
               "Zm00001eb191150",
               "Zm00001eb195180",
               "Zm00001eb186410",
               "Zm00001eb195900",
               "Zm00001eb188030",
               "Zm00001eb013800",
               "Zm00001eb191620",
               "Zm00001eb025970",
               "Zm00001eb190090",
               "Zm00001eb195860",
               "Zm00001eb195370",
               "Zm00001eb189720",
               "Zm00001eb195340",
               "Zm00001eb187610",
               "Zm00001eb187460",
               "Zm00001eb194180",
               "Zm00001eb157030",
               "Zm00001eb187290",
               "Zm00001eb007720",
               "Zm00001eb188350",
               "Zm00001eb187400",
               "Zm00001eb191330",
               "Zm00001eb190230",
               "Zm00001eb190880",
               "Zm00001eb225320",
               "Zm00001eb183030")

# TS1_B73_null <- TS1_B73_null[unique_MICH_null, ]
# TS1_B73_P <- TS1_B73_P[unique_MICH_null, ]
# TS1_MICH_null <- TS1_MICH_null[unique_MICH_null, ]
# TS1_MICH_P <- TS1_MICH_P[unique_MICH_null, ]

#TS1_B73_null <- TS1_B73_null[MICH_nullP_highR, ]
#TS1_B73_P <- TS1_B73_P[MICH_nullP_highR, ]
#TS1_MICH_null <- TS1_MICH_null[MICH_nullP_highR, ]
#TS1_MICH_P <- TS1_MICH_P[MICH_nullP_highR, ]
# TS1_B73_null <- TS1_B73_null[MICH_only, ]
# TS1_B73_P <- TS1_B73_P[MICH_only, ]
# TS1_MICH_null <- TS1_MICH_null[MICH_only, ]
# TS1_MICH_P <- TS1_MICH_P[MICH_only, ]

# Removing zeroes and replacing with lowest/2 value and taking log value
TS1_B73_null <- TS1_B73_null %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1_B73_null <- log(TS1_B73_null)

TS1_B73_P <- TS1_B73_P %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1_B73_P <- log(TS1_B73_P)

TS1_MICH_null <- TS1_MICH_null %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1_MICH_null <- log(TS1_MICH_null)

TS1_MICH_P <- TS1_MICH_P %>% 
  mutate(across(everything(), ~ ifelse(. == 0, min(.[. > 0], na.rm = TRUE) / 2, .)))
TS1_MICH_P <- log(TS1_MICH_P)



# Adding the time points
time_points_row <- data.frame(Sample_1 = 1, Sample_2 = 2, Sample_3 = 3, Sample_4 = 4)

TS1_B73_null <- rbind(time_points_row, TS1_B73_null)
rownames(TS1_B73_null)[1] <- "time_points"
TS1_B73_null <- as.matrix(TS1_B73_null)

TS1_B73_P <- rbind(time_points_row, TS1_B73_P)
rownames(TS1_B73_P)[1] <- "time_points"
TS1_B73_P <- as.matrix(TS1_B73_P)

TS1_MICH_null <- rbind(time_points_row, TS1_MICH_null)
rownames(TS1_MICH_null)[1] <- "time_points"
TS1_MICH_null <- as.matrix(TS1_MICH_null)

TS1_MICH_P <- rbind(time_points_row, TS1_MICH_P)
rownames(TS1_MICH_P)[1] <- "time_points"
TS1_MICH_P <- as.matrix(TS1_MICH_P)


##### Running the GRN
#### Only P
## All 4
time.points_all4 <- list(TS1_B73_null[1,],TS1_B73_P[1,],TS1_MICH_null[1,],TS1_MICH_P[1,])
TS.data_all4 <- list(TS1_B73_null[2:nrow(TS1_B73_null), ],TS1_B73_P[2:nrow(TS1_B73_P), ],TS1_MICH_null[2:nrow(TS1_MICH_null), ],TS1_MICH_P[2:nrow(TS1_MICH_P), ])

#MICH_null MICH_P and B73_null
time.points_MICH_NULL_P_B73_NULL <- list(TS1_B73_null[1,],TS1_MICH_null[1,],TS1_MICH_P[1,])
TS.data_B73_MICH_NULL_P_B73_NULL <- list(TS1_B73_null[2:nrow(TS1_B73_null), ],TS1_MICH_null[2:nrow(TS1_MICH_null), ],TS1_MICH_P[2:nrow(TS1_MICH_P), ])

## MICH_null and MICH_P
time.points_MICH_null_MICH_P <- list(TS1_MICH_null[1,],TS1_MICH_P[1,])
TS.data_B73_MICH_null_MICH_P <- list(TS1_MICH_null[2:nrow(TS1_MICH_null), ],TS1_MICH_P[2:nrow(TS1_MICH_P), ])

## MICH_null and B73_null
time.points_B73_null_MICH_null <- list(TS1_MICH_null[1,],TS1_B73_null[1,])
TS.data_B73_null_MICH_null <- list(TS1_MICH_null[2:nrow(TS1_MICH_null), ],TS1_B73_null[2:nrow(TS1_B73_null), ])

## B73_P and B73_null
time.points_B73_null_B73_P <- list(TS1_B73_null[1,],TS1_B73_P[1,])
TS.data_B73_null_B73_P <- list(TS1_B73_null[2:nrow(TS1_B73_null), ],TS1_B73_P[2:nrow(TS1_B73_P), ])

## Just MICH_null
time.points_MICH_null <- list(TS1_MICH_null[1,])
TS.data_MICH_null <- list(TS1_MICH_null[2:nrow(TS1_MICH_null), ])

## Just B73_null
time.points_B73_null <- list(TS1_B73_null[1,])
TS.data_B73_null <- list(TS1_B73_null[2:nrow(TS1_B73_null), ])


# Running DynGENIE3
res_all4 <- dynGENIE3(TS.data_all4, time.points_all4)
link.list_all4 <- get.link.list(res_all4$weight.matrix)

res_MICH_NULL_P_B73_NULL <- dynGENIE3(TS.data_B73_MICH_NULL_P_B73_NULL, time.points_MICH_NULL_P_B73_NULL)
link.list_MICH_NULL_P_B73_NULL <- get.link.list(res_MICH_NULL_P_B73_NULL$weight.matrix)

res_MICH_null_MICH_P <- dynGENIE3(TS.data_B73_MICH_null_MICH_P, time.points_MICH_null_MICH_P)
link.list_MICH_null_MICH_P <- get.link.list(res_MICH_null_MICH_P$weight.matrix)

res_B73_null_MICH_null <- dynGENIE3(TS.data_B73_null_MICH_null, time.points_B73_null_MICH_null)
link.list_B73_null_MICH_null <- get.link.list(res_B73_null_MICH_null$weight.matrix)

res_MICH_null <- dynGENIE3(TS.data_MICH_null, time.points_MICH_null)
link.list_MICH_null <- get.link.list(res_MICH_null$weight.matrix)

res_B73_null <- dynGENIE3(TS.data_B73_null, time.points_B73_null)
link.list_B73_null <- get.link.list(res_MICH_null$weight.matrix)

res_B73_null_B73_P <- dynGENIE3(TS.data_B73_null_B73_P, time.points_B73_null_B73_P)
link.list_B73_null_B73_P <- get.link.list(res_B73_null_B73_P$weight.matrix)


### Plotting:
# Open a new Quartz device
quartz()

# Set up a 2x2 plotting layout
par(mfrow=c(2,2))

network_all4 <- graph_from_data_frame(link.list_all4[1:20,], directed = TRUE)
isolated_nodes <- which(degree(network_all4) == 0)
network_all4 <- delete_vertices(network_all4, isolated_nodes)

network_MICH_NULL_P_B73_NULL <- graph_from_data_frame(link.list_MICH_NULL_P_B73_NULL[1:10,], directed = TRUE)
isolated_nodes <- which(igraph::degree(network_MICH_NULL_P_B73_NULL) == 0)
network_MICH_NULL_P_B73_NULL <- delete_vertices(network_MICH_NULL_P_B73_NULL, isolated_nodes)

network_MICH_null_MICH_P <- graph_from_data_frame(link.list_MICH_null_MICH_P[1:20,], directed = TRUE)
isolated_nodes <- which(igraph::degree(network_MICH_null_MICH_P) == 0)
network_MICH_null_MICH_P <- delete_vertices(network_MICH_null_MICH_P, isolated_nodes)

network_B73_null_MICH_null <- graph_from_data_frame(link.list_B73_null_MICH_null[1:10,], directed = TRUE)
isolated_nodes <- which(igraph::degree(network_B73_null_MICH_null) == 0)
network_B73_null_MICH_null <- delete_vertices(network_B73_null_MICH_null, isolated_nodes)

network_MICH_null <- graph_from_data_frame(link.list_MICH_null[1:20,], directed = TRUE)
isolated_nodes <- which(igraph::degree(network_MICH_null) == 0)
network_MICH_null <- delete_vertices(network_MICH_null, isolated_nodes)

network_B73_null <- graph_from_data_frame(link.list_B73_null[1:20,], directed = TRUE)
isolated_nodes <- which(igraph::degree(network_B73_null) == 0)
network_B73_null <- delete_vertices(network_B73_null, isolated_nodes)
B73_null_unique <- c(unique(link.list_B73_null[1:20,1]),unique(link.list_B73_null[1:20,2]))

plot(network_all4, 
     layout=layout_with_fr, 
     edge.arrow.size=.5, 
     vertex.color="skyblue", 
     vertex.size=igraph::degree(network)*5,  # Using a constant size for all vertices
     edge.width=1)   # Using a constant width for all edges


plot(network_MICH_NULL_P_B73_NULL, 
     layout=layout_with_fr, 
     edge.arrow.size=.5, 
     vertex.color="skyblue", 
     vertex.size=4, 
     edge.width=2)

plot(network_MICH_null_MICH_P, 
     layout=layout_with_fr, 
     edge.arrow.size=.5, 
     vertex.color="skyblue", 
     #vertex.size=igraph::degree(network)*2, 
     vertex.size=5, 
     edge.width=3)

plot(network_B73_null_MICH_null, 
     layout=layout_with_fr, 
     edge.arrow.size=.5, 
     vertex.color="skyblue", 
     vertex.size=igraph::degree(network)*5, 
     edge.width=1)

plot(network_MICH_null, 
     layout=layout_with_fr, 
     edge.arrow.size=.5, 
     vertex.color="skyblue", 
     vertex.size=5, 
     edge.width=3)

plot(network_B73_null, 
     layout=layout_with_fr, 
     edge.arrow.size=.5, 
     vertex.color="skyblue", 
     vertex.size=5, 
     edge.width=3)

write_graph(network_MICH_null, file="network_MICH_null.gml", format="gml")
write_graph(network_MICH_null_MICH_P, file="network_MICH_null_MICH_P.gml", format="gml")
write_graph(network_B73_null, file="network_B73_null.gml", format="gml")


dev.off()
