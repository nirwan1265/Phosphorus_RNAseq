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
  Phosphorus = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Control = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)
rownames(design_matrix_B73MCH) <- design_matrix_B73MCH[,1]
design_matrix_B73MCH <- design_matrix_B73MCH[,-1]

design_B73MCH <- make.design.matrix(design_matrix_B73MCH, degree = 3)

fit_B73MCH <- p.vector(B73MCH, design_B73MCH, Q = 0.05, MT.adjust = "BH", min.obs = 10, counts = TRUE)

tstep_B73MCH <- T.fit(fit_B73MCH, step.method = "backward", alfa = 0.05)

sig_B73MCH <- get.siggenes(tstep_B73MCH, rsq = 0.6, vars = "groups")
suma2Venn(sig_B73MCH$summary[, c(1:2)])


list_all_sig_B73MCH <- tstep_B73MCH[["sol"]]

list_phosphorus_sig_B73MCH <- sig_B73MCH[["sig.genes"]][["Phosphorus"]][["sig.pvalues"]]

list_controlvsphosphorus_sig_B73MCH <- sig_B73MCH[["sig.genes"]][["ControlvsPhosphorus"]][["sig.pvalues"]]

list_onlyP_B73MCH <-as.data.frame(rownames(list_phosphorus_sig_B73MCH)[!(rownames(list_phosphorus_sig_B73MCH) %in% rownames(list_controlvsphosphorus_sig_B73MCH))])

list_onlyCvsP_B73MCH <-as.data.frame(rownames(list_controlvsphosphorus_sig_B73MCH)[!(rownames(list_controlvsphosphorus_sig_B73MCH) %in% rownames(list_phosphorus_sig_B73MCH))])


write.csv(list_all_sig_B73MCH,"all_significant_genes_B73MCH.csv", row.names = T)

write.csv(list_onlyP_B73MCH,"onlyP_significant_genes_B73MCH.csv", row.names = F)

write.csv(list_onlyCvsP_B73MCH,"ControlvsP_significant_genes_B73MCH.csv", row.names = F)




