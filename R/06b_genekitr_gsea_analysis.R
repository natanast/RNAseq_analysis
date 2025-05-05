

rm(list = ls())
gc()


# libraries --------

library(data.table)
library(stringr)
library(genekitr)
library(geneset)


# load data --------

fls = "DESeq2/" |> 
  list.files(full.names = TRUE, pattern = "DESeq2") |>
  str_subset("xlsx")



# gene sets ---------

go_ontology <- getGO(org = "human", ont = "mf")
kegg        <- getKEGG(org = "hsa", category = "pathway")
reactome    <- getReactome(org = "human")


gene_sets = list(
  "GO" = go_ontology, 
  "KEGG" = kegg, 
  "reactome" = reactome
)


# gsea -----------


for(i in seq_len(length(fls))) {
    
    df = fls[i] |> readxl::read_xlsx() |> setDT()
    
    comparison = fls[i] |>
        basename() |>
        str_sub(8, -1) |>
        str_sub(1, -6)
    
    
    df = df[, c("GeneID", "log2FoldChange", "pvalue", "padj"), with = FALSE]
    
    gene_map = df$GeneID |> transId(transTo = c("symbol", "entrez", "ensembl", "uniprot"), unique = TRUE) |> setDT()

    gene_map = gene_map[which(
        !is.na(symbol) & !is.na(entrezid) & !is.na(ensembl) & !is.na(uniprot)
    )]
    
    df = df |> merge(gene_map, by.x = "GeneID", by.y = "input_id")
    
    
    df = df[order(symbol, -abs(log2FoldChange))]
    
    df = df[, by = symbol, head(.SD, 1)]
    
    df = df[, by = uniprot, head(.SD, 1)]
    df = df[, by = ensembl, head(.SD, 1)]
    df = df[, by = entrezid, head(.SD, 1)]
    
    # ranking
    df = df[order(-log2FoldChange)]
    
    # df = df[which(!is.na(pvalue) & pvalue != 0)]
    # 
    # df[, rank_score := sign(log2FoldChange) * (-log10(pvalue))]
    # 
    # df = df[order(-rank_score)]  # Sort genes from highest to lowest rank
    
    
    out = list()
    
    for (j in names(gene_sets)) {
        
        # Create index
        # index <- df$rank_score
        index <- df$log2FoldChange
        
        names(index) <- df$entrezid
        
        # run genGSEA and handle errors 
        tmp <- tryCatch(
            {
                # Perform GSEA analysis
                result <- genGSEA(index, geneset = gene_sets[[j]])
                
                # Save the result as an RDS file
                saveRDS(result, file = paste0(comparison, "-", j, ".rds"))
                
                # Extract GSEA results as a data.table
                result$gsea_df |> setDT()
            },
            error = function(e) {
                message("Error in genGSEA for ", j, ": ", e$message)
                NULL # Return NULL in case of error
            }
        )
        
        # Skip if tmp is NULL
        if (is.null(tmp)) next
        
        # Store results if no error
        out[[j]] <- tmp
    }
    
    
    out = out |> rbindlist(use.names = TRUE, fill = TRUE, idcol = "db")
    
    
    writexl::write_xlsx(out, paste0("gsea_genekitr_", comparison, ".xlsx"))
    
    
}



