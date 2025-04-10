

rm(list = ls())
gc()


# libraries ---------

library(data.table)
library(stringr)
library(DESeq2)


# load data ---------

sm = "sample_metadata.xlsx" |> readxl::read_xlsx(sheet = 1) |> setDT()

df = fread("gene_count_clean/gene_counts_clean_all.txt")

sm = sm[order(Group)]


# clean data ---------

index = which( (df |> colnames()) %in% c("gene_name", sm$Sample) )

df = df[, index, with = FALSE]


# ~ Group
state <- list(
    c("Subsets_1_99", "Other_U_CLL"),
    c("Subsets_1_99", "M_CLL"),
    c("Other_U_CLL", "M_CLL")
)


# differential expression analysis -----------

for(i in seq_along(state)) {
    
    q = state[[i]]
    
    
    s0 = sm[Group %in% c(q[1], q[2])]
    
    s0$Group <- s0$Group |> factor(levels = c(q[2], q[1]))
    
    index = s0$Sample
    
    mm = df[, index, with = FALSE] |> setDF(rownames = df$gene_name)
    
    index2 = mm |> rowSums()
    
    mm = mm[which( index2 >= 10 ), ]
    
    
    dds = DESeqDataSetFromMatrix(
        countData = mm,
        colData   = s0,
        design    = ~ Group
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
        results(pAdjustMethod = "fdr", name = paste0("Group_", q[1],"_vs_", q[2])) |>
        as.data.frame() |>
        setDT(keep.rownames = "GeneID")
    
    
    res = res[order(padj), ]
    
    res$comparison = paste0(q[1], "_vs_", q[2])
    
    writexl::write_xlsx(res, paste0("DESeq2_", q[1], "_vs_", q[2], ".xlsx"))
    
}


# part 2 -------

# load data ---------

sm = "sample_metadata.xlsx" |> readxl::read_xlsx(sheet = 2) |> setDT()

df = fread("gene_count_clean/gene_counts_clean_without_last.txt")

index = which( (df |> colnames()) %in% c("gene_name", sm$Sample) )

df = df[, index, with = FALSE]



sm = sample_metadata[order(Group)]

sm$Group  = sm$Group |> factor(levels = c("M_CLL", "Other_U_CLL", "Subsets_1_99"))


# differential expression analysis ----------

index = sm$Sample

mm <- df[, index, with = FALSE] |> setDF(rownames = df$gene_name) |> as.matrix()


dds <- DESeqDataSetFromMatrix(
    countData = mm,
    colData   = sm,
    design    = ~ 0 + Group
)

dds = estimateSizeFactors(dds)


normalized_counts = counts(dds, normalized = TRUE)

normalized_counts = normalized_counts |> 
    as.data.frame() |> 
    setDT(keep.rownames = "gene_id")


fwrite(
    normalized_counts, row.names = FALSE, quote = FALSE, sep = "\t",
    paste0("gene_counts_deseq2_normalized_without_last.txt")
)

