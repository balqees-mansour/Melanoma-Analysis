# Load required libraries
install.packages("svglite")

library(svglite)
library(ggplot2)

# Assuming 'res' is your DESeq2 results object
# Convert to a data frame
res_df <- as.data.frame(res)

# Add -log10(p-value) for the y-axis
res_df$negLog10Pvalue <- -log10(res_df$pvalue)

# Define new significance thresholds
log2FC_threshold <- 3  # Log2 fold change threshold
pvalue_threshold <- 0.05 # P-value threshold

# Add a column for significance classification
res_df$Significance <- ifelse(
  res_df$pvalue < pvalue_threshold & abs(res_df$log2FoldChange) > log2FC_threshold,
  ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"),
  "Not Significant"
)

# Define colors for points
color_palette <- c("Downregulated" = "green", "Not Significant" = "black", "Upregulated" = "red")

# Create the volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10Pvalue, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = color_palette) +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plot",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~italic(P)~value),
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(volcano_plot)
# Save the plot as an SVG file with DPI set to 300
# Save the plot as a PNG file with 300 DPI
ggsave(
  filename = "volcano_plot.png",  # File name
  plot = volcano_plot,           # Plot object
  device = "png",                # Save as PNG
  dpi = 300,                     # Resolution (DPI)
  width = 8,                     # Width in inches
  height = 6                     # Height in inches
)















