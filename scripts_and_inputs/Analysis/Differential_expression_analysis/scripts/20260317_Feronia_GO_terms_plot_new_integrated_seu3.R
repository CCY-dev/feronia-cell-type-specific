rm(list=ls())
library(tidyverse)
library(ggrepel)
library(readxl)
library(ComplexHeatmap)
library(circlize)
# Inputs are read from this repository; set FERONIA_REPO to your clone (default: current directory).
repo <- Sys.getenv("FERONIA_REPO", ".")

### Up-regulated
# Read in the GO terms Michelle selected
up_tbl_terms_tbl <- read_xlsx(path = file.path(repo, "scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/GO_term_selection/Upregulated_genes_GO_table-15 GO terms.xlsx"))
up_tbl_terms <- up_tbl_terms_tbl$Go_term[up_tbl_terms_tbl$`Selected?` %in% c("Yes", "yes")]
length(up_tbl_terms)

# Read in the GO table from the results of the newly integrated Seurat 3 obj
up_tbl <- read.table(file.path(repo, "scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/GO/Upregulated_genes_GO_table.txt"),
                     sep = "\t", header = T, check.names = F)


up_tbl_sel <- up_tbl %>%
  filter(Go_term %in% up_tbl_terms)

# The GO terms Michelle selected but didn't appear here
setdiff(up_tbl_terms, up_tbl_sel$Go_term)

# Convert to data.frame
mat <- as.data.frame(up_tbl_sel)

# Set GO terms as rownames
rownames(mat) <- mat$Go_term
mat$Go_term <- NULL

# Convert to numeric matrix
mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

col_fun <- colorRamp2(
  breaks = c(
    0, 10
  ),
  colors = c("beige", "red2")
)

#pdf("20260317_Feronia_up_selected_GO_fromNewIntegrated_seu3obj.pdf", width = 9, height = 5)
ht <- Heatmap(
  mat,
  name = "Enrichment",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "-log10_pval",
    legend_height = unit(2, "cm")
  )
)

# BLTR
draw(
  ht,
  heatmap_legend_side = "left",
  padding = unit(c(8, 10, 2, 2), "mm") 
)

dev.off()


### Down-regulated
# Read in the GO terms Michelle selected
down_tbl_terms_tbl <- read_xlsx(path = file.path(repo, "scripts_and_inputs/Analysis/Differential_expression_analysis/inputs/GO_term_selection/Downregulated_genes_GO_table-15 GO term.xlsx"))
down_tbl_terms <- down_tbl_terms_tbl$Go_term[down_tbl_terms_tbl$`Selected?` %in% c("Yes", "yes")]
length(down_tbl_terms)

# Read in the GO table from the results of the newly integrated Seurat 3 obj
down_tbl <- read.table(file.path(repo, "scripts_and_inputs/Analysis/Differential_expression_analysis/outputs/GO/Downregulated_genes_GO_table.txt"),
                     sep = "\t", header = T, check.names = F)

down_tbl_sel <- down_tbl %>%
  filter(Go_term %in% down_tbl_terms)

# The GO terms Michelle selected but didn't appear here
setdiff(down_tbl_terms, down_tbl_sel$Go_term)

# Convert to data.frame
mat <- as.data.frame(down_tbl_sel)

# Set GO terms as rownames
rownames(mat) <- mat$Go_term
mat$Go_term <- NULL

# Convert to numeric matrix
mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

col_fun <- colorRamp2(
  breaks = c(
    0, 10
  ),
  colors = c("beige", "blue3")
)

#pdf("20260317_Feronia_down_selected_GO_fromNewIntegrated_seu3obj.pdf", width = 9, height = 5)
ht <- Heatmap(
  mat,
  name = "Enrichment",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  heatmap_legend_param = list(
    title = "-log10_pval",
    legend_height = unit(2, "cm")
  )
)

# BLTR
draw(
  ht,
  heatmap_legend_side = "left",
  padding = unit(c(8, 10, 2, 2), "mm") 
)

dev.off()

