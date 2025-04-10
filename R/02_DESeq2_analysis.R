

rm(list = ls())
gc()


# libraries ---------

library(data.table)
library(stringr)
library(DESeq2)


# load data ---------

sample_metadata = "sample_metadata.xlsx" |> readxl::read_xlsx() |> setDT()

df = fread("gene_counts_merged_clean.txt")


# clean data --------

index = which( (df |> colnames()) %in% c("gene_name", sample_metadata$Sample) )
df = df[, index, with = FALSE]


index = which( sample_metadata$Sample %in% colnames(df) )
sample_metadata = sample_metadata[index]


sample_metadata = sample_metadata[order(Group1)]


sample_metadata$Group1  = sample_metadata$Group1 |> 
    factor(levels = c("Trisomies", "mutated_no_subset", "subset_4_IgG", "subset_4_IgM", "subset_6", "subset_8"))



# differential expression analysis 1 ----------

index = sample_metadata$Sample

mm <- df[, index, with = FALSE] |> setDF(rownames = df$gene_name) |> as.matrix()


dds <- DESeqDataSetFromMatrix(
  countData = mm,
  colData   = sample_metadata,
  design    = ~ 0 + Group1
)

dds = estimateSizeFactors(dds)


normalized_counts = counts(dds, normalized = TRUE)

normalized_counts = normalized_counts |>
  as.data.frame() |>
  setDT(keep.rownames = "gene_id")


fwrite(
  normalized_counts, row.names = FALSE, quote = FALSE, sep = "\t",
  paste0("gene_counts_deseq2_normalized.txt")
)



dds = DESeq(dds, test = "Wald")


resultsNames(dds)


res = dds |>
  results(contrast = c("Group1", "Trisomies", "subset_8")) |>
  as.data.frame() |>
  setDT(keep.rownames = "GeneID")

res = res[order(padj)]


writexl::write_xlsx(res, paste0("DESeq2_Trisomies_vs_subset_8.xlsx"))



# differential expression analysis 2 ----------


index = sample_metadata$Sample

mm <- df[, index, with = FALSE] |> setDF(rownames = df$gene_name) |> as.matrix()


# ~ Group
state <- list(
    c("Trisomies", "subset_8"),
    c("Trisomies", "mutated_no_subset"),
    c("Trisomies", "subset_4_IgG")
)


for(i in seq_along(state)) {
    
    q = state[[i]]
    
    
    s0 = sample_metadata[Group1 %in% c(q[1], q[2])]
    
    s0$Group1 <- s0$Group1 |> factor(levels = c(q[2], q[1]))
    
    index = s0$Sample
    
    mm = df[, index, with = FALSE] |> setDF(rownames = df$gene_name)
    
    index2 = mm |> rowSums()
    
    mm = mm[which( index2 >= 10 ), ]
    
    
    dds = DESeqDataSetFromMatrix(
        countData = mm,
        colData   = s0,
        design    = ~ Group1
    )
    
    dds               = estimateSizeFactors(dds)
    normalized_counts = counts(dds, normalized = TRUE)
    
    normalized_counts = normalized_counts |> 
        as.data.frame() |> 
        setDT(keep.rownames = "GeneID")
    
    
    fwrite(
        normalized_counts, row.names = FALSE, quote = FALSE, sep = "\t",
        paste0("gene_counts_clean_DESEq2_normalized_", q[1], "_vs_", q[2], ".txt")
    )
    
    
    dds = DESeq(dds, test = "Wald")
    
    
    # resultsNames(dds)
    
    res = dds |> 
        results(pAdjustMethod = "fdr", name = paste0("Group1_", q[1],"_vs_", q[2])) |>
        as.data.frame() |>
        setDT(keep.rownames = "GeneID")
    
    
    res = res[order(padj), ]
    
    res$comparison = paste0(q[1], "_vs_", q[2])
    
    writexl::write_xlsx(res, paste0("DESeq2_", q[1], "_vs_", q[2], ".xlsx"))
    
}
