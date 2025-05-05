


rm(list = ls())
gc()


# libraries -------

library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)
library(colorspace)
library(shadowtext)


# load data -------

df = "DESeq2/DESeq2_Group_M_vs_G.xlsx" |>
    readxl::read_xlsx(sheet = 1) |>
    setDT()


# data cleaning ------

df = df[which( !is.na(df$padj) )]

df$y = -log10(df$padj)


df$ann = ifelse(
    df$padj > .05, "Not significant",
    ifelse(
        df$log2FoldChange > 0, "Up regulated", "Down regulated"
    )
)


df$ann = ifelse(
  df$padj <= 0.05 & df$log2FoldChange > -1 & df$log2FoldChange < 1, 
  paste0(df$ann, " (low)"),
  df$ann
)



df2 = df[which(padj <= .05 & abs(log2FoldChange) > 1)]

df2 = df2[order( abs(log2FoldChange), decreasing = TRUE )]

df2 = df2[, by = ann, head(.SD, 10) ]



# plot -----------

gr = ggplot(data = df) +
  
  geom_point(aes(x = log2FoldChange, y = -log10(padj), fill = ann),
             shape = 21, stroke = NA, size = 2, alpha = .5) +
  
  geom_vline(xintercept = c(-1, 1), linewidth = .3, linetype = "dashed", lineend = "round") +
  geom_hline(yintercept = -log10(.05), linewidth = .3, linetype = "dashed", lineend = "round") +
  
  geom_point(data = df2, aes(x = log2FoldChange, y = -log10(padj), fill = ann), 
             shape = 21, stroke = .15, size = 2, color = "white") +
  
  geom_text_repel(
    data = df2, aes(x = log2FoldChange, y = -log10(padj), label = GeneID),
    max.overlaps = Inf, # label.size = NA, fill = "transparent",
    fontface = "bold", size = 2.5, bg.color = "white", bg.r = .05
  ) +
  
  scale_fill_manual(
    values = c(
      "Up regulated" = "#990000",
      "Up regulated (low)" = lighten("#990000", .5),
      
      "Down regulated" = "#004d99",
      "Down regulated (low)" = lighten("#004d99", .5),
      
      "Not significant" = "grey"
    ),
    
    breaks = c("Up regulated", "Not significant", "Down regulated"),
    
    guide = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    ),
    
  ) +
  
  #scale_x_continuous(trans = scales::pseudo_log_trans()) +
  
  scale_x_continuous(
    #limits = c(-10, 5),
    breaks = c(-5, -2.5, -1, 0, 1, 2.5, 5),
    trans = scales::pseudo_log_trans()
    ) +
  scale_y_continuous(expand = c(0, 0), breaks = c(2, 5, 10, 20, 30, 40),
                     trans = scales::pseudo_log_trans()
                     ) +
  
  coord_cartesian(clip = "off") +
  
  theme_minimal() +
  
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    
    # axis.title = element_text(size = 14),
    # axis.text = element_text(size = 14),
    
    axis.line = element_line(linewidth = .3, color = "black"),
    axis.ticks = element_line(linewidth = .3, color = "black"),
    
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey85"),
    
    plot.margin = margin(20, 20, 20, 20),
    
    
  ) +
  
  labs(y = "-log10(padj)", x = "log2(Fold Change)")


gr


ggsave(
  plot = gr, filename = paste0("volcano_M_vs_G.jpeg"), 
  width = 10, height = 10, units = "in", dpi = 600
)


