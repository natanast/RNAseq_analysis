


rm(list = ls())
gc()


# libraries -------

library(data.table)
library(stringr)
library(tidyverse)

library(ComplexHeatmap)
library(colorRamp2)

library(ggplotify)
library(ggplot2)

library(pathfindR)
library(genekitr)


# load data ------------

sample_metadata = "sample_metadata.xlsx" |> readxl::read_xlsx() |> setDT()


# heatmap  (overall) ------------------

x = list()

x[[1]] = "DESeq2/DESeq2_Trisomies_vs_mutated_no_subset.xlsx" |>
    readxl::read_xlsx(sheet = 1 ) |>
    setDT()

x[[2]] = "DESeq2/DESeq2_Trisomies_vs_subset_4_IgG.xlsx" |>
    readxl::read_xlsx(sheet = 1) |>
    setDT()

x[[3]] = "DESeq2/DESeq2_Trisomies_vs_subset_8.xlsx" |>
    readxl::read_xlsx(sheet = 1) |>
    setDT()


x = x |> rbindlist(use.names = TRUE, fill = TRUE)




df = fread("DESeq2/gene_counts_deseq2_normalized.txt")



x = x[which(abs(log2FoldChange) > 2 & padj <= 0.01)]



df = df[which(gene_id %in% x$GeneID)]


ann = sample_metadata[which(Sample %in% colnames(df))]

ann = ann[Group1 != "subset_6" & Group1 != "subset_4_IgM"]


# Define order 
ann$Group1 <- factor(ann$Group1, levels = c("Trisomies", "subset_4_IgG", "subset_8", "mutated_no_subset"))

# Reorder the samples 
ann <- ann[order(ann$Group1), ]



# z-score ------------
mm = df[, ann$Sample, with = FALSE] |> setDF(rownames = df$gene_id) |> as.matrix()


zmat = mm |> t() |> scale(scale = TRUE, center = TRUE) |> t()




# Heatmap -------

col = colorRamp2(breaks = c(-4, -2, 0, 2, 4), colors =  c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a'))

ha = HeatmapAnnotation("Group" = ann$Group1, col = list("Group" = c("Trisomies" = "#00429d", "subset_4_IgG" =  '#73a2c6', "subset_8" = "#93003a", "mutated_no_subset" = "#f4777f")))


ht = Heatmap(
    zmat, name = "Z-score",

    col = col,
    border = TRUE,
    

    clustering_distance_columns = "pearson",
    clustering_distance_rows = "pearson",

    clustering_method_columns = "ward.D2",
    clustering_method_rows = "ward.D2",

    # cluster_columns = FALSE,
    column_split = ann$Group1,
    # column_split = 4,
    row_split = 5,


    top_annotation = ha,
    show_row_names = FALSE,
    # row_names_gp = gpar(fontsize = 3),
    column_names_gp = gpar(fontsize = 10)
)


gr = ht |>
    draw() |>
    grid.grabExpr() |>
    as.ggplot()

gr


ggsave(
    plot = gr, filename = "heatmap_overall_supervised.jpeg",
    width = 12, height = 12, units = "in", dpi = 600
)
