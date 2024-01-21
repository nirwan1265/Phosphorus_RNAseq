# Data directory
data_dir <- "data/"
system("ls data/")

# Time series RNA seq data
data <- read.csv("data/B73&MICH.csv")
head(data)
rownames(data) <- data[,1]
data <- data[,-1]
data[1:4,1:5]


# Design_matrix
design_matrix <- data.frame(
  colnames = c("B73_P_null_T1_R1", "B73_P_null_T2_R1", "B73_P_null_T2_R2", "B73_P_null_T2_R3",
               "B73_P_null_T3_R1", "B73_P_null_T3_R2", "B73_P_null_T4_R1", "B73_P_null_T4_R2",
               "B73_P_add_T1_R1", "B73_P_add_T1_R2", "B73_P_add_T1_R3", "B73_P_add_T2_R1",
               "B73_P_add_T2_R2", "B73_P_add_T2_R3", "B73_P_add_T3_R1", "B73_P_add_T3_R2",
               "B73_P_add_T3_R3", "B73_P_add_T4_R1", "B73_P_add_T4_R2", "B73_P_add_T4_R3",
               "MICH_P_null_T1_R1", "MICH_P_null_T1_R2", "MICH_P_null_T1_R3", "MICH_P_null_T1_R4",
               "MICH_P_null_T2_R1", "MICH_P_null_T2_R2", "MICH_P_null_T3_R1", "MICH_P_null_T3_R2", "MICH_P_null_T3_R3",
               "MICH_P_null_T4_R1", "MICH_P_null_T4_R2",
               "MICH_P_add_T1_R1", "MICH_P_add_T1_R2", "MICH_P_add_T1_R3",
               "MICH_P_add_T2_R1", "MICH_P_add_T2_R2", "MICH_P_add_T2_R3",
               "MICH_P_add_T3_R1", "MICH_P_add_T3_R2", "MICH_P_add_T3_R3",
               "MICH_P_add_T4_R1", "MICH_P_add_T4_R2", "MICH_P_add_T4_R3"),
  Time = c(1, 2, 2, 2, 3, 3, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4,
           1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4,1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
  Replicate = c(1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8,
                1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
  B73_noP = c(1, 1, 1, 1, 1, 1, 1, 1, 
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  B73_P = c(0, 0, 0, 0, 0, 0, 0, 0, 
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  MICH_noP = c(0, 0, 0, 0, 0, 0, 0, 0, 
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  MICH_P = c(0, 0, 0, 0, 0, 0, 0, 0, 
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)
rownames(design_matrix) <- design_matrix[,1]
design_matrix <- design_matrix[,-1]




# Creating a regression matrix for the full regression model, degrees = 3 for 4 time points
design <- make.design.matrix(design_matrix, degree = 3)

# Fit the first model
fit <- p.vector(data, design, Q = 0.05, MT.adjust = "BH", min.obs = 10)

# Significant differences using second step wise regression
tstep <- T.fit(fit, step.method = "forward", alfa = 0.05)

# Get significant genes
sig_each <- get.siggenes(tstep, rsq = 0.6, vars = "each")
sig_group <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
sig_all <- get.siggenes(tstep, rsq = 0.6, vars = "all")


# Significant genes
B73_noP <- sig_group$sig.genes$B73_noP$sig.pvalues
B73_PvsB73_noP <- sig_group$sig.genes$B73_PvsB73_noP$sig.pvalues
MICH_noPvsB73_noP <- sig_group$sig.genes$MICH_noPvsB73_noP$sig.pvalues
MICH_PvsB73_noP <- sig_group$sig.genes$MICH_PvsB73_noP$sig.pvalues

# saving:
write.csv(B73_noP,"result/all_analysis/B73_noP.csv")
write.csv(B73_PvsB73_noP,"result/all_analysis/B73_PvsB73_noP.csv")
write.csv(MICH_noPvsB73_noP,"result/all_analysis/MICH_noPvsB73_noP.csv")
write.csv(MICH_PvsB73_noP,"result/all_analysis/MICH_PvsB73_noP.csv")


# Venn Diagram
suma2Venn(sig_each$summary[, c(5,9,13)])
suma2Venn(sig_group$summary[, c(2:4)])


# Cluster
quartz()
see.genes(sig_group$sig.genes$B73_noP, show.fit = T, dis =design$dis, groups.vector = design$groups.vector,
          cluster.method="hclust" ,cluster.data = 1, k.mclust=TRUE,distance="euclidean")
?see.genes()


see.genes(data,design_matrix,how.fit = T, dis =design$dis, groups.vector = design$groups.vector,
          cluster.method="hclust" ,cluster.data = 1, k.mclust=TRUE,distance="euclidean")
?see.genes())
# Get sig genes



x

PlotGroups(design, edesign = design_matrix)
data("data.abiotic")
data("edesign.abiotic")

typeof(edesign.abiotic)
typeof(design_matrix)
head(edesign.abiotic)
head(design_matrix)
class(edesign.abiotic)
class(design_matrix)

design_matrix <- apply(design_matrix, 2, as.integer)
design <- data.frame(design)

typeof(design)
typeof(data.abiotic)
class(design)
class(data.abiotic)


design_matrix
STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]


# For individual gene. replace data with individual gene
quartz()
PlotGroups (data, edesign = design_matrix, show.fit = T,
             dis = design$dis, groups.vector = design$groups.vector)


