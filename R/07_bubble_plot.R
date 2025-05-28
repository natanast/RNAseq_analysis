

rm(list = ls())
gc()



library(data.table)
library(stringr)

library(ggplot2)

library(tidyverse)


x = "gsea_genekitr_Trisomies_vs_subset_4_IgG.xlsx" |>
  readxl::read_xlsx() |>
  setDT()

x$log_padj = -log10(x$p.adjust)

x = x[which(x$p.adjust <= 0.05)]


# bubble plot-----
gr = ggplot(data = x) +
  
  geom_point(
    aes(x = NES, y = Description,
        size = `Count`, 
        fill = log_padj
        ),
    
    shape = 21, stroke = .25
  ) +
  
  scale_size_continuous(range = c(5, 12)) +
  
  scale_fill_gradient(
    low = "grey", high = "#f30000",
    guide = guide_colorbar(
      title = "-log10(p.adj)",
      barheight = unit(10, "lines"),
      barwidth = unit(.75, "lines")
    )
  ) +
  
  theme_minimal() +
  
  theme(
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x =  element_text(size = 12),
    axis.text.y =  element_text(size = 12),
    
    panel.grid.major = element_line(linewidth = .4),
    
    panel.border = element_rect(fill = NA, linewidth = .4),
    axis.ticks = element_line(linewidth = .4),
    
    plot.margin = margin(20, 20, 20, 20)
  )
  
  #facet_wrap(vars(Group), nrow = 4, scales = "free_y")
  # facet_grid(
  #   #rows = vars(Group1),
  #   scales = "free_y",
  #   space = "free_y")

gr


ggsave(
  plot = gr, filename = paste0("bubble_plot_Trisomies_vs_subset_4_IgG.jpeg"), 
  width = 12, height = 12, units = "in", dpi = 600
)



# box plot------

gr = x |>
  ggplot(aes(y = NES, x = Description, fill = log_padj)) +

  geom_bar(stat = "identity", alpha = .8, width = .45) +
  
  scale_fill_gradient(
    low = "grey", high = "#f30000",
    guide = guide_colorbar(
      title = "-log10(p.adj)",
      barheight = unit(10, "lines"),
      barwidth = unit(.75, "lines")
    )
  ) +
  
  coord_flip() +
  
  theme_minimal() +
  
  theme(
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x =  element_text(size = 10),
    axis.text.y =  element_text(size = 10),
    
    panel.grid.major = element_line(linewidth = .4),
    
    panel.border = element_rect(fill = NA, linewidth = .4),
    axis.ticks = element_line(linewidth = .4),
    
    plot.margin = margin(20, 20, 20, 20)
  )


ggsave(
  plot = gr, filename = paste0("bar_plot_plot_Trisomies_vs_subset_4_IgG.jpeg"), 
  width = 12, height = 12, units = "in", dpi = 600
)

