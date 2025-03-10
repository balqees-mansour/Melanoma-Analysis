#============================================================================
# Load required libraries
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

# Convert gene symbols to ENTREZ IDs for upregulated genes
entrez_up <- mapIds(org.Hs.eg.db,
                    keys = rownames(res_upregulated),
                    keytype = "SYMBOL",
                    column = "ENTREZID")

# Convert gene symbols to ENTREZ IDs for downregulated genes
entrez_down <- mapIds(org.Hs.eg.db,
                      keys = rownames(res_downregulated),
                      keytype = "SYMBOL",
                      column = "ENTREZID")

# Remove any NA values
entrez_up <- na.omit(entrez_up)
entrez_down <- na.omit(entrez_down)

# Perform enrichment analysis
## Biological Process
ego_up <- enrichGO(gene = entrez_up,
                   OrgDb = org.Hs.eg.db,
                   ont = "ALL",
                   keyType = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)

ego_down <- enrichGO(gene = entrez_down,
                     OrgDb = org.Hs.eg.db,
                     ont = "ALL",
                     keyType = "ENTREZID",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)

# KEGG pathway analysis
ekegg_up <- enrichKEGG(gene = entrez_up,
                       organism = "hsa",
                       pvalueCutoff = 0.05)

ekegg_down <- enrichKEGG(gene = entrez_down,
                         organism = "hsa",
                         pvalueCutoff = 0.05)

# Create visualizations
## GO terms
p1 <- dotplot(ego_up, showCategory = 10, title = "Upregulated GO Terms") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 8))

p2 <- dotplot(ego_down, showCategory = 10, title = "Downregulated GO Terms") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 8))

## KEGG pathways
p3 <- dotplot(ekegg_up, showCategory = 10, title = "Upregulated KEGG Pathways") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 8))

p4 <- dotplot(ekegg_down, showCategory = 5, title = "Downregulated KEGG Pathways") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 8))
# Create improved dotplot
dotplot(ekegg_down, 
        showCategory = 10,
        title = "KEGG Pathways Enrichment Analysis") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  scale_color_gradient(low = "red", high = "blue")
# Arrange all plots
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
















# Convert gene symbols to ENTREZ IDs and remove NA values/duplicates
entrez_up <- mapIds(org.Hs.eg.db,
                    keys = rownames(res_upregulated),
                    keytype = "SYMBOL",
                    column = "ENTREZID",
                    multiVals = "first")
entrez_up <- unique(na.omit(entrez_up))

entrez_down <- mapIds(org.Hs.eg.db,
                      keys = rownames(res_downregulated),
                      keytype = "SYMBOL",
                      column = "ENTREZID",
                      multiVals = "first")
entrez_down <- unique(na.omit(entrez_down))

# Perform GO enrichment analysis
ego_up <- enrichGO(gene = entrez_up,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   keyType = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

ego_upmf <- enrichGO(gene = entrez_up,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   keyType = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

ego_upcc <- enrichGO(gene = entrez_up,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   keyType = "ENTREZID",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

ego_down <- enrichGO(gene = entrez_down,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",
                     keyType = "ENTREZID",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)

# Perform KEGG pathway analysis
ekegg_up <- enrichKEGG(gene = entrez_up,
                       organism = "hsa",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       use_internal_data = FALSE)

ekegg_down <- enrichKEGG(gene = entrez_down,
                         organism = "hsa",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         use_internal_data = FALSE)



# Enhanced Barplot for GO enrichment analysis (Upregulated)
p1 <- barplot(ego_up,
              showCategory = 10,  # Show top 10 categories
              title = "Upregulated BP GO Terms") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),  # Larger and bold text for categories
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_gradient(low = "#FF4B4B", high = "#4B7BFF")

p1mf <- barplot(ego_upmf,
              showCategory = 10,  # Show top 10 categories
              title = "Upregulated MF GO Terms") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),  # Larger and bold text for categories
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_gradient(low = "#FF4B4B", high = "#4B7BFF")


p1cc <- barplot(ego_upmf,
                showCategory = 10,  # Show top 10 categories
                title = "Upregulated CC GO Terms") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),  # Larger and bold text for categories
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_gradient(low = "#FF4B4B", high = "#4B7BFF")
# Enhanced Barplot for KEGG pathway analysis (Upregulated)
p3 <- barplot(ekegg_up,
              showCategory = 10,  # Show top 10 pathways
              title = "Upregulated KEGG Pathways") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),  # Larger and bold text for categories
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_gradient(low = "#FF4B4B", high = "#4B7BFF")

gridExtra::grid.arrange(p1, p1mf, p1cc, p3, ncol=2)

# Combine all enhanced barplots into a single figure
p_combined <- gridExtra::grid.arrange(p1, p1mf, p1cc, p3, ncol=2)

# Save the combined enhanced barplots as an image file
ggsave("enhanced_enrichment_analysis_barplots.png",
       plot=p_combined,
       width=14,
       height=12,
       dpi=300,
       units="in")




# Enhanced Barplot for GO enrichment analysis (Downregulated)
p2 <- barplot(ego_down,
              showCategory = 10,  # Show top 10 categories
              title = "Downregulated BP GO Terms") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),  # Larger and bold text for categories
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_gradient(low = "#FF4B4B", high = "#4B7BFF")



# Enhanced Barplot for KEGG pathway analysis (Downregulated)
p4 <- barplot(ekegg_down,
              showCategory = 10,  # Show top 10 pathways
              title = "Downregulated KEGG Pathways") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),   # Larger and bold text for categories
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_fill_gradient(low="#FF4B4B", high="#4B7BFF")

gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

# Combine all enhanced barplots into a single figure
p_combined <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

# Save the combined enhanced barplots as an image file
ggsave("enhanced_enrichment_analysis_barplots.png",
       plot=p_combined,
       width=14,
       height=12,
       dpi=300,
       units="in")





