

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


# gene sets --------

go_ontology <- getGO(org = "human", ont = "mf")
kegg        <- getKEGG(org = "hsa", category = "pathway")
reactome    <- getReactome(org = "human")


gene_sets = list(
    "GO" = go_ontology, 
    "KEGG" = kegg, 
    "reactome" = reactome
)


# ORA --------

for(i in seq_len(length(fls))) {
    
      df = fls[i] |> readxl::read_xlsx() |> setDT()
      
      df = df[which( pvalue <= .05 )]
      
      if( nrow(df) == 0 ) next
      
      df = df[, c("GeneID", "log2FoldChange", "padj", "pvalue"), with = FALSE]
      
      comparison = fls[i] |> basename() |> str_sub(8, -1) |> str_sub(1, -6)
      
      
      out = list()
      

      for (j in names(gene_sets)) {
          
          
          tmp <- tryCatch(
              {
                  genORA(df$GeneID, geneset = gene_sets[[j]]) |> setDT()
              },
              error = function(e) {
                  message("No terms enriched for ", j, ". Skipping...")
                  NULL # Return NULL in case of error
              }
          )
          
          # Skip if tmp is NULL
          if (is.null(tmp)) next
          
          # Store results if no error
          out[[j]] <- tmp
      }
      
      
      out = out |> rbindlist(use.names = TRUE, fill = TRUE, idcol = "db")
      
      
      writexl::write_xlsx(out, paste0("ora_genekitr_", comparison, "_pvalue.xlsx"))
      

    
}
