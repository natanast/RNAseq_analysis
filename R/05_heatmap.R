

rm(list = ls())
gc()


# load libraries --------

library(data.table)
library(stringr)
library(ComplexHeatmap)
library(colorRamp2)
library(ggplotify)
library(ggplot2)


# load data ----------

d = "DESeq2/DESeq2_Group_M_vs_G.xlsx" |> 
    readxl::read_xlsx() |> 
    setDT()

mat <- "DESeq2/gene_counts_deseq2_normalized.txt" |> fread()

sm <- "sample_metadata.xlsx" |> readxl::read_xlsx() |> setDT()


# filtering ----------

d <- d[which(padj <= .05)]

index <- d$GeneID |> unique()

mat <- mat[which(Geneid %in% index)]


# z-score ------------

m <- mat[, -1] |> setDF(rownames = mat$Geneid) |> as.matrix()

zmat <- m |> t() |> scale(center = TRUE, scale = TRUE) |> t()



# Heatmap ------------

zmat <- zmat[, sm$Sample]



col = colorRamp2(
    breaks = c(-4, -2, 0, 2, 4), 
    colors =  c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a')
)


ha = HeatmapAnnotation(
    "Group" = sm$Group, 
    col = list(
        "Group" = c("M" = "#00429d", 
                    "G" = "#93003a")
    )
)



ht <- Heatmap(
    zmat, use_raster = TRUE, name = "Z-score",
    
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    
    col = col,
    
    border = TRUE,
    

    column_split = 2,
    row_split = 2,
    
    # column_split = sm$Group,
    top_annotation = ha,
    # right_annotation = ha_row,

    
    # top_annotation = ha,
    # right_annotation = ha_row,
    
    show_row_names = FALSE
)



gr <- ht |>
    draw(merge_legends = TRUE) |>
    grid.grabExpr() |>
    as.ggplot()


gr


ggsave(
    plot = gr, filename = paste0("heatmap_unsupervised_M_vs_G_without.png"),
    width = 10, height = 10, units = "in", dpi = 600
)    

