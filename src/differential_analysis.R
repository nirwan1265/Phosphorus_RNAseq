# Data directory
data_dir <- "data/"
system("ls data/")


# Load exp data
B73 <- read.csv(paste0(data_dir,"B73.csv"))
rownames(B73) <- B73[,1]
B73 <- B73[,-1]
colnames(B73)

MICH <- read.csv(paste0(data_dir,"MICH.csv"))
rownames(MICH) <- MICH[,1]
MICH <- MICH[,-1]
colnames(MICH)


# Design_matrix
# B73
design_matrix_B73 <- data.frame(
  colnames = c("B73_P_null_T1_R1", "B73_P_null_T2_R1", "B73_P_null_T2_R2", "B73_P_null_T2_R3",
               "B73_P_null_T3_R1", "B73_P_null_T3_R2", "B73_P_null_T4_R1", "B73_P_null_T4_R2",
               "B73_P_add_T1_R1", "B73_P_add_T1_R2", "B73_P_add_T1_R3", "B73_P_add_T2_R1",
               "B73_P_add_T2_R2", "B73_P_add_T2_R3", "B73_P_add_T3_R1", "B73_P_add_T3_R2",
               "B73_P_add_T3_R3", "B73_P_add_T4_R1", "B73_P_add_T4_R2", "B73_P_add_T4_R3"),
  Time = c(1, 2, 2, 2, 3, 3, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
  Replicate = c(1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8),
  Phosphorus = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Control = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

rownames(design_matrix_B73) <- design_matrix_B73[,1]
design_matrix_B73 <- design_matrix_B73[,-1]

#MICH
design_matrix_MICH <- data.frame(
  colnames = c("MICH_P_null_T1_R1", "MICH_P_null_T1_R2", "MICH_P_null_T1_R3", "MICH_P_null_T1_R4",
               "MICH_P_null_T2_R1", "MICH_P_null_T2_R2", "MICH_P_null_T3_R1", "MICH_P_null_T3_R2", "MICH_P_null_T3_R3",
               "MICH_P_null_T4_R1", "MICH_P_null_T4_R2",
               "MICH_P_add_T1_R1", "MICH_P_add_T1_R2", "MICH_P_add_T1_R3",
               "MICH_P_add_T2_R1", "MICH_P_add_T2_R2", "MICH_P_add_T2_R3",
               "MICH_P_add_T3_R1", "MICH_P_add_T3_R2", "MICH_P_add_T3_R3",
               "MICH_P_add_T4_R1", "MICH_P_add_T4_R2", "MICH_P_add_T4_R3"),
  Time = c(1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4,
        1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
  Replicate = c(1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4,
        5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8),
  Phosphorus = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Control = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

rownames(design_matrix_MICH) <- design_matrix_MICH[,1]
design_matrix_MICH <- design_matrix_MICH[,-1]





# Creating a regression matrix for the full regression model, degrees = 3 for 4 time points
design_B73 <- make.design.matrix(design_matrix_B73, degree = 3)
design_MICH <- make.design.matrix(design_matrix_MICH, degree = 3)

# Significant genes
fit_B73 <- p.vector(B73, design_B73, Q = 0.05, MT.adjust = "BH", min.obs = 20)
fit_MICH <- p.vector(MICH, design_MICH, Q = 0.05, MT.adjust = "BH", min.obs = 20)


# Significant differences
tstep_B73 <- T.fit(fit_B73, step.method = "backward", alfa = 0.05)
tstep_MICH <- T.fit(fit_MICH, step.method = "backward", alfa = 0.05)



# List of Significant genes
sig_B73 <- get.siggenes(tstep_B73, rsq = 0.6, vars = "groups")
sig_MICH <- get.siggenes(tstep_MICH, rsq = 0.6, vars = "groups")

# Venn
suma2Venn(sig_B73$summary[, c(1:2)])
suma2Venn(sig_MICH$summary[, c(1:2)])


list_all_sig_B73 <- tstep_B73[["sol"]]
list_all_sig_MICH <- tstep_MICH[["sol"]]


list_phosphorus_sig_B73 <- sig_B73[["sig.genes"]][["Phosphorus"]][["sig.pvalues"]]
list_phosphorus_sig_B73 <- list_phosphorus_sig_B73 %>% rownames_to_column("gene")
list_phosphorus_sig_MICH <- sig_MICH[["sig.genes"]][["Phosphorus"]][["sig.pvalues"]]
list_phosphorus_sig_MICH <- list_phosphorus_sig_MICH %>% rownames_to_column("gene")




list_controlvsphosphorus_sig_B73 <- sig_B73[["sig.genes"]][["ControlvsPhosphorus"]][["sig.pvalues"]]
list_controlvsphosphorus_sig_B73 <- list_controlvsphosphorus_sig_B73 %>% rownames_to_column("gene")
list_controlvsphosphorus_sig_MICH <- sig_MICH[["sig.genes"]][["ControlvsPhosphorus"]][["sig.pvalues"]]
list_controlvsphosphorus_sig_MICH <- list_controlvsphosphorus_sig_MICH %>% rownames_to_column("gene")

list_onlyP_B73 <-list_phosphorus_sig_B73$gene[!(list_phosphorus_sig_B73$gene %in% list_controlvsphosphorus_sig_B73$gene )]
list_onlyP_B73 <- list_phosphorus_sig_B73[list_phosphorus_sig_B73$gene %in% list_onlyP_B73, ]
list_onlyP_MICH <-list_phosphorus_sig_MICH$gene[!(list_phosphorus_sig_MICH$gene %in% list_controlvsphosphorus_sig_MICH$gene )]
list_onlyP_MICH <- list_phosphorus_sig_MICH[list_phosphorus_sig_MICH$gene %in% list_onlyP_MICH, ]


list_onlyControlvsPhosphorus_B73 <- list_controlvsphosphorus_sig_B73$gene[!(list_controlvsphosphorus_sig_B73$gene %in% list_phosphorus_sig_B73$gene)]
list_onlyControlvsPhosphorus_B73 <- list_controlvsphosphorus_sig_B73[list_controlvsphosphorus_sig_B73$gene %in% list_onlyControlvsPhosphorus_B73, ]
list_onlyControlvsPhosphorus_MICH <- list_controlvsphosphorus_sig_MICH$gene[!(list_controlvsphosphorus_sig_MICH$gene %in% list_phosphorus_sig_MICH$gene)]
list_onlyControlvsPhosphorus_MICH <- list_controlvsphosphorus_sig_MICH[list_controlvsphosphorus_sig_MICH$gene %in% list_onlyControlvsPhosphorus_MICH, ]


# Venn
suma2Venn(sig_B73$summary[, c(1:2)])
suma2Venn(sig_B73$summary[, c(1:2)])


# Gene list
sig_B73$sig.genes$ControlvsPhosphorus$g

# Cluster
quartz()
see.genes(sig_B73$sig.genes$ControlvsPhosphorus, show.fit = T, dis =design_B73$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 4)

see.genes(sig_MICH$sig.genes$ControlvsPhosphorus, show.fit = T, dis =design_MICH$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 4)
