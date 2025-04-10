

rm(list = ls())
gc()


# libraries -------

library(data.table)
library(stringr)

library(ggplot2)
library(ggsci)
library(ggforce)
library(colorspace)



# load data -------

sample_metadata = "sample_metadata.xlsx" |> readxl::read_xlsx() |> setDT()


df1 = "DESeq2/DESeq2_Trisomies_vs_subset_8.xlsx" |> 
    readxl::read_xlsx() |>
    setDT()


df2 = paste0("DESeq2/gene_counts_deseq2_normalized.txt") |>
    fread()



# clean data ----------------

sample_metadata = sample_metadata[Group1 %in% c("Trisomies", "subset_8"),]


index = which( (df2 |> colnames()) %in% c("gene_id", sample_metadata$Sample) )

df2 = df2[, index, with = FALSE]



significant_results <- df1[which(padj <= 0.01 & abs(log2FoldChange) >= 1)]
significant_genes   <- significant_results$GeneID

filtered_counts <- df2[which(gene_id %in% significant_genes)]


filtered_counts_df = filtered_counts |> transpose(make.names = "gene_id", keep.names = "Sample")

filtered_counts_df = filtered_counts_df[, 2:ncol(filtered_counts_df)] |> 
  setDF(rownames = filtered_counts_df$Sample)



# PCA analysis --------

res.pca <- prcomp(filtered_counts_df, scale = TRUE, center = TRUE)


df = res.pca$x |> as.data.frame() |> setDT(keep.rownames = "Sample")

sample_metadata = sample_metadata[which(Sample %in% df$Sample)]

df = df |> merge(sample_metadata, by = "Sample")



# plot -----------

summary(res.pca)


prop_var <- summary(res.pca)$importance[2, ]  # Extract Proportion of Variance
pc1_label <- paste0("PC1 (", round(prop_var[1] * 100, 2), "%)")
pc2_label <- paste0("PC2 (", round(prop_var[2] * 100, 2), "%)")



gr = ggplot(df, aes(PC1, PC2)) +
    
    geom_mark_ellipse(aes(fill = Group1, label = Group1), alpha = .1, expand = unit(2, "mm")) +
    
    geom_point(aes(fill = Group1), shape = 21, size = 3, stroke = .25, color = "white") +  
    
    scale_fill_manual(
        values = c(
            "subset_4_IgG"          = "#1170AA",
            "Trisomies"             = "#F28E2B",
            "mutated_no_subset"     = "#E15759",
            "subset_8"              = "#76B7B2"
        )
    ) +
    
    geom_text(aes(label = Sample)) +
    
    scale_x_continuous(limits = c(-150, 150)) +
    scale_y_continuous(limits = c(-150, 100)) +
    
    
    theme_minimal() +
    
    theme(
        legend.position = "none",
        plot.margin = margin(20, 20, 20, 20)
    ) +
    
    labs(x = pc1_label, y = pc2_label)
  


gr

ggsave(
  plot = gr, filename = paste0("pca_trisomies_vs_subset_8_labels.jpeg"), 
  width = 8, height = 8, units = "in",dpi = 600
)






# # PCA 2 (all)--------------

df2 = paste0("DESeq2/gene_counts_deseq2_normalized.txt") |>
    fread()

# data
sample_metadata = sample_metadata[Group1 != "subset_6" & Group1 != "subset_4_IgM",]


index = which( (df2 |> colnames()) %in% c("gene_id", sample_metadata$Sample) )

df2 = df2[, index, with = FALSE]


df1 = list()

df1[["Trisomies_vs_mutated_no_subset"]] = "DESeq2/DESeq2_Trisomies_vs_mutated_no_subset.xlsx" |>
  readxl::read_xlsx(sheet = 1) |>
  setDT()

df1[["Trisomies_vs_subset_4_IgG"]] = "DESeq2/DESeq2_Trisomies_vs_subset_4_IgG.xlsx" |>
  readxl::read_xlsx(sheet = 1) |>
  setDT()

df1[["Trisomies_vs_subset_8"]] = "DESeq2/DESeq2_Trisomies_vs_subset_8.xlsx" |>
  readxl::read_xlsx(sheet = 1) |>
  setDT()


df1 = df1 |> rbindlist(use.names = TRUE, fill = TRUE)




significant_results <- df1[which(padj <= 0.05), ]

significant_genes <- significant_results$GeneID |> unique()



filtered_counts <- df2[which(gene_id %in% significant_genes), ]

filtered_counts_df = transpose(filtered_counts, make.names = "gene_id", keep.names = "Sample")

filtered_counts_df = filtered_counts_df[, 2:ncol(filtered_counts_df)] |>
  setDF(rownames = filtered_counts_df$Sample)




# pca
res.pca <- prcomp(filtered_counts_df, scale = TRUE, center = TRUE)


df = res.pca$x |> as.data.frame() |> setDT(keep.rownames = "Sample")

df = merge(df, sample_metadata, by = "Sample")


# plot ---------

prop_var <- summary(res.pca)$importance[2, ]  # Extract Proportion of Variance
pc1_label <- paste0("PC1 (", round(prop_var[1] * 100, 2), "%)")
pc2_label <- paste0("PC2 (", round(prop_var[2] * 100, 2), "%)")



gr <- ggplot(df, aes(PC1, PC2)) +

  geom_mark_ellipse(aes(fill = Group1, label = Group1), alpha = .1, expand = unit(3, "mm")) +

  geom_point(aes(fill = Group1), shape = 21, size = 3, stroke = .25, color = "white") +

  scale_fill_manual(
    values = c(
      "subset_4_IgG"          = "#1170AA",
      "Trisomies"             = "#F28E2B",
      "mutated_no_subset"     = "#E15759",
      "subset_8"              = "#76B7B2"
    )
  ) +

  scale_x_continuous(limits = c(-150, 120)) +
  scale_y_continuous(limits = c(-200, 100)) +

  theme_minimal() +

  theme(
      
    #legend.position = "none",
      
    plot.margin = margin(20, 20, 20, 20)
  ) +

    labs(x = pc1_label, y = pc2_label)



gr


ggsave(
  plot = gr, filename = paste0("pca_total_PC1_PC2.jpeg"),
  width = 8, height = 8, units = "in",dpi = 600
)




# PCA viz total-------- 

df2 <- df |> melt(
    id.vars = c("Sample", "Group1", "Patient", "Batch"),
    value.factor = FALSE,
    variable.factor = FALSE
) |> setDT()


library(ggforce)

gr1 = df2[which(variable %in% paste0("PC", 1:9))]|>
    ggplot(aes(variable, value)) +
    geom_point(aes(color = Group1, fill = Group1), 
               shape = 21, stroke = .15, size = 3,
               position = position_jitternormal(sd_y = 0, sd_x = .05)) +
    
    scale_fill_manual(
        values = c(
            "subset_4_IgG" = "#1170AA",
            "Trisomies"    = "#F28E2B",
            "mutated_no_subset"     = "#E15759",
            "subset_8"     = "#76B7B2"
        ) |> lighten(.25)
    ) +
    
    scale_color_manual(
        values = c(
            "subset_4_IgG" = "#1170AA",
            "Trisomies"    = "#F28E2B",
            "mutated_no_subset"     = "#E15759",
            "subset_8"     = "#76B7B2"
        ) |> darken(.25)
    ) +
    
    theme_minimal() +
    
    theme(
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        
        plot.margin = margin(20, 20, 20, 20)
    ) +
    
    labs(x = "", y = "Value")


gr1



ggsave(
    plot = gr1, filename = paste0("pca_total.jpeg"),
    width = 8, height = 8, units = "in",dpi = 600
)
