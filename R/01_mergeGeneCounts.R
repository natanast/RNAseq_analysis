


rm(list = ls())
gc()


# libraries -------

library(data.table)
library(stringr)


# 1 ---------------------


d1 <- "../analysis_20241014_trisomy/gene-counts-clean-featureCounts.txt" |> fread()
d2 <- "../analysis_20241014_trisomy/salmon.merged.gene_counts_1.tsv" |> fread()
d3 <- "../analysis_20241014_trisomy/salmon.merged.gene_counts_2.tsv" |> fread()



d1 <- d1[, !c("12557IgG", "13240IgG", "13941IgG", "23556IgGhigh", "23556IgMlow", "36924IgG", "46184IgG", "24133_IgG", "23556IgG"), with = FALSE]


colnames(d1)[1] <- "gene_name"

d2$gene_id <- NULL
d3$gene_id <- NULL

d2 <- d2[, by = gene_name, lapply(.SD, sum)]
d3 <- d3[, by = gene_name, lapply(.SD, sum)]

for(i in 2:ncol(d2)) d2[[i]] <- d2[[i]] |> round()
for(i in 2:ncol(d3)) d3[[i]] <- d3[[i]] |> round()


# 2 --------------------------------------------------


d4 <- "../analysis_20241206_trisomy_puplic/counts.txt" |> fread()

mapping <- d4$geneID |>
    genekitr::transId(transTo = "symbol", unique = TRUE, hgVersion = "v38") |>
    setDT()

index <- match(d4$geneID, mapping$input_id)

d4$gene_name <- mapping[index]$symbol

d4$gene_name <- ifelse(is.na(d4$gene_name), d4$geneID, d4$gene_name)

d4$geneID <- NULL

d4 <- d4[, by = gene_name, lapply(.SD, sum)]

# 3 ----------------------------------

d5 <- "../analysis_20250408/gene-counts_subset_8.txt" |> fread()


d5 <- d5[, c(1, 7:ncol(d5)), with = FALSE]

colnames(d5) <- colnames(d5) |> str_split_i("\\.", 1)
colnames(d5)[1] <- "gene_name"


# clean environment -------------

rm(i, index, mapping)


# merge gene counts --------------------------

df <- d1 |>
    merge(d2, by = "gene_name", all = TRUE) |>
    merge(d3, by = "gene_name", all = TRUE) |>
    merge(d4, by = "gene_name", all = TRUE) |>
    merge(d5, by = "gene_name", all = TRUE)


for(i in 2:ncol(df)) df[[i]] <- ifelse( is.na(df[[i]]), 0, df[[i]] )

rm(i)
gc()


# find chromosome ----------------


mapping <- df$gene_name |> 
    genekitr::genInfo(unique = TRUE, hgVersion = "v38") |>
    setDT()



# filter -----------

keep <- which(mapping$chr %in% c( as.character(1:22), "X", "Y" ))


mapping <- mapping[keep]


df <- df[which(gene_name %in% mapping$input_id)]

chr12 <- mapping[chr == 12, ]

chr12$symbol 


df_filtered <- df[!(gene_name %in% chr12$symbol)]



fwrite(
    df_filtered, "gene_counts_merged.txt",
    row.names = FALSE, quote = FALSE, sep = "\t"
)












