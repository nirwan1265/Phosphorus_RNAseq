# Load
B73MCH <- read.csv("data/MICHvsB73_noP.csv")
rownames(B73MCH) <- B73MCH[,1]
B73MCH <- B73MCH[,-1]
names(B73MCH)
# Design Matrix
design_matrix_B73MCH <- data.frame(
  colnames = c("B73_P_null_T1_R1", "B73_P_null_T2_R1", "B73_P_null_T2_R2", "B73_P_null_T2_R3",
               "B73_P_null_T3_R1", "B73_P_null_T3_R2", "B73_P_null_T4_R1", "B73_P_null_T4_R2",
               "MICH_P_null_T1_R1", "MICH_P_null_T1_R2","MICH_P_null_T1_R3","MICH_P_null_T1_R4",
               "MICH_P_null_T2_R1", "MICH_P_null_T2_R2", "MICH_P_null_T3_R1", "MICH_P_null_T3_R2",
               "MICH_P_null_T3_R3", "MICH_P_null_T4_R1","MICH_P_null_T4_R2"),
  Time = c(1, 2, 2, 2, 3, 3, 4, 4, 1,1,1,1,2,2,3,3,3,4,4),
  Replicate = c(1, 2, 2, 2, 3, 3, 4, 4, 5,5,5,5,6,6,7,7,7,8,8),
  B73_noP = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  MICH_noP = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)
rownames(design_matrix_B73MCH) <- design_matrix_B73MCH[,1]
design_matrix_B73MCH <- design_matrix_B73MCH[,-1]

design_B73MCH <- make.design.matrix(design_matrix_B73MCH, degree = 3)


# Regression fit
fit_B73MCH <- p.vector(B73MCH, design_B73MCH, Q = 0.05, MT.adjust = "BH", min.obs = 10, counts = TRUE)

# Stepwise regression fit
tstep_B73MCH <- T.fit(fit_B73MCH, step.method = "backward", alfa = 0.05,family = negative.binomial() )


# Significant genes
sig_B73MCH_groups <- get.siggenes(tstep_B73MCH, rsq = 0.6, vars = "groups")
sig_B73MCH_all <- get.siggenes(tstep_B73MCH, rsq = 0.6, vars = "all")
sig_B73MCH_each <- get.siggenes(tstep_B73MCH, rsq = 0.6, vars = "each")

suma2Venn(sig_B73MCH_groups$summary[, c(1:2)])
suma2Venn(sig_B73MCH_each$summary[, c(1:8)])
suma2Venn(sig_B73MCH_each$summary[, c(1:8)])


# Saving
B73_noP <- sig_B73MCH_groups$sig.genes$B73_noP$sig.pvalues
MICH_noPvsB73_noP <- sig_B73MCH_groups$sig.genes$MICH_noPvsB73_noP$sig.pvalues

write.csv(B73_noP,"result/B73_noP_VS_MICH_noP/B73_noP.csv")
write.csv(MICH_noPvsB73_noP,"result/B73_noP_VS_MICH_noP/MICH_noPvsB73_noP.csv")


list_all_sig_B73MCH <- tstep_B73MCH[["sol"]]

list_phosphorus_sig_B73MCH <- sig_B73MCH[["sig.genes"]][["Phosphorus"]][["sig.pvalues"]]

list_controlvsphosphorus_sig_B73MCH <- sig_B73MCH[["sig.genes"]][["ControlvsPhosphorus"]][["sig.pvalues"]]

list_onlyP_B73MCH <-as.data.frame(rownames(list_phosphorus_sig_B73MCH)[!(rownames(list_phosphorus_sig_B73MCH) %in% rownames(list_controlvsphosphorus_sig_B73MCH))])

list_onlyCvsP_B73MCH <-as.data.frame(rownames(list_controlvsphosphorus_sig_B73MCH)[!(rownames(list_controlvsphosphorus_sig_B73MCH) %in% rownames(list_phosphorus_sig_B73MCH))])


write.csv(list_all_sig_B73MCH,"result/B73MICH_noP/all_significant_genes_B73MCH_noP.csv", row.names = T)

write.csv(list_onlyP_B73MCH,"result/B73MICH_noP/noP_significant_genes_B73MCH.csv", row.names = F)

write.csv(list_onlyCvsP_B73MCH,"result/B73MICH_noP/MICHnoPVSB73noP_significant_genes_B73MCH.csv", row.names = F)




