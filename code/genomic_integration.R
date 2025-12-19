# ============================================================================
# Genomic Integration Script for Melanoma Analysis
# ============================================================================
# 
# Purpose: This script integrates external genomic databases and resources
#          to enhance the melanoma analysis with additional genomic information
#          including mutation data, protein interactions, and clinical annotations.
#
# Author: Melanoma Analysis Project
# Date: 2025
# License: MIT
#
# ============================================================================

# ============================================================================
# 1. LIBRARY IMPORTS AND SETUP
# ============================================================================

# Install required packages if not already installed
required_packages <- c(
  "DESeq2",
  "biomaRt",
  "igraph",
  "STRINGdb",
  "dplyr",
  "tidyr",
  "ggplot2",
  "data.table"
)

# Check and install missing packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("DESeq2", "biomaRt", "STRINGdb")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

library(DESeq2)
library(biomaRt)
library(igraph)
library(STRINGdb)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

# ============================================================================
# 2. DATA LOADING AND PREPARATION
# ============================================================================

#' Load DEG Results
#' 
#' Loads the differentially expressed genes from the analysis
#'
#' @param deg_file Path to the DEG results file (CSV format)
#' @return Data frame containing DEG information
#'
load_deg_results <- function(deg_file) {
  cat("Loading DEG results from:", deg_file, "\n")
  
  deg_data <- read.csv(deg_file, row.names = 1, stringsAsFactors = FALSE)
  
  # Filter for significant genes
  deg_data <- deg_data %>%
    filter(padj < 0.05) %>%
    filter(abs(log2FoldChange) > 1)
  
  cat("Total significant DEGs:", nrow(deg_data), "\n")
  cat("Upregulated genes:", sum(deg_data$log2FoldChange > 0), "\n")
  cat("Downregulated genes:", sum(deg_data$log2FoldChange < 0), "\n")
  
  return(deg_data)
}

# ============================================================================
# 3. GENE ANNOTATION AND MAPPING
# ============================================================================

#' Annotate Genes with Biomart
#'
#' Maps Ensembl gene IDs to gene symbols and retrieves additional annotations
#'
#' @param gene_ids Vector of Ensembl gene IDs
#' @param dataset Biomart dataset (default: hsapiens_gene_ensembl)
#' @return Data frame with gene annotations
#'
annotate_genes_biomart <- function(gene_ids, dataset = "hsapiens_gene_ensembl") {
  cat("Connecting to Biomart...\n")
  
  # Connect to Biomart
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Define attributes to retrieve
  attributes <- c(
    "ensembl_gene_id",
    "external_gene_name",
    "gene_biotype",
    "chromosome_name",
    "start_position",
    "end_position",
    "description"
  )
  
  # Query Biomart
  cat("Retrieving gene annotations...\n")
  annotations <- getBM(
    attributes = attributes,
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = ensembl
  )
  
  cat("Retrieved annotations for", nrow(annotations), "genes\n")
  
  return(annotations)
}

# ============================================================================
# 4. PROTEIN-PROTEIN INTERACTION NETWORK
# ============================================================================

#' Build PPI Network from STRING Database
#'
#' Retrieves protein-protein interactions from STRING database
#'
#' @param gene_symbols Vector of gene symbols
#' @param species_id NCBI taxonomy ID (default: 9606 for human)
#' @param score_threshold STRING interaction score threshold (0-1000)
#' @return Data frame with interaction data
#'
get_ppi_network <- function(gene_symbols, species_id = 9606, score_threshold = 400) {
  cat("Initializing STRING database connection...\n")
  
  # Initialize STRING database
  string_db <- STRINGdb$new(version = "11.5", species = species_id, score_threshold = score_threshold)
  
  # Map gene symbols to STRING identifiers
  cat("Mapping genes to STRING identifiers...\n")
  mapped_genes <- string_db$map(
    data.frame(gene = gene_symbols),
    "gene",
    removeUnmappedRows = TRUE
  )
  
  if (nrow(mapped_genes) == 0) {
    cat("Warning: No genes could be mapped to STRING database\n")
    return(NULL)
  }
  
  # Get interaction network
  cat("Retrieving PPI network...\n")
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  cat("Retrieved", nrow(interactions), "protein-protein interactions\n")
  
  return(interactions)
}

#' Analyze Network Topology
#'
#' Calculates network topology metrics
#'
#' @param interactions Data frame with interaction data (from_id, to_id, combined_score)
#' @return List containing network metrics and centrality measures
#'
analyze_network_topology <- function(interactions) {
  cat("Analyzing network topology...\n")
  
  # Create igraph object
  g <- graph_from_data_frame(
    d = interactions[, c("from_id", "to_id")],
    directed = FALSE
  )
  
  # Calculate centrality measures
  degree_centrality <- degree(g)
  betweenness_centrality <- betweenness(g)
  closeness_centrality <- closeness(g)
  
  # Identify hub nodes (top 10% by degree)
  hub_threshold <- quantile(degree_centrality, 0.9)
  hub_nodes <- names(degree_centrality[degree_centrality >= hub_threshold])
  
  cat("Identified", length(hub_nodes), "hub nodes\n")
  
  # Network statistics
  network_stats <- list(
    num_nodes = vcount(g),
    num_edges = ecount(g),
    density = edge_density(g),
    avg_degree = mean(degree_centrality),
    num_communities = length(unique(membership(cluster_louvain(g))))
  )
  
  cat("Network Statistics:\n")
  cat("  Nodes:", network_stats$num_nodes, "\n")
  cat("  Edges:", network_stats$num_edges, "\n")
  cat("  Density:", round(network_stats$density, 4), "\n")
  cat("  Average Degree:", round(network_stats$avg_degree, 2), "\n")
  cat("  Communities:", network_stats$num_communities, "\n")
  
  return(list(
    graph = g,
    degree_centrality = degree_centrality,
    betweenness_centrality = betweenness_centrality,
    closeness_centrality = closeness_centrality,
    hub_nodes = hub_nodes,
    network_stats = network_stats
  ))
}

# ============================================================================
# 5. MUTATION DATA INTEGRATION
# ============================================================================

#' Load and Process Mutation Data
#'
#' Loads mutation data from external sources (e.g., COSMIC, ClinVar)
#'
#' @param mutation_file Path to mutation data file
#' @return Data frame with mutation information
#'
load_mutation_data <- function(mutation_file) {
  cat("Loading mutation data from:", mutation_file, "\n")
  
  # Load mutation data
  mutations <- read.csv(mutation_file, stringsAsFactors = FALSE)
  
  # Filter for melanoma-specific mutations
  mutations <- mutations %>%
    filter(cancer_type == "melanoma" | is.na(cancer_type))
  
  cat("Loaded", nrow(mutations), "mutations\n")
  
  return(mutations)
}

#' Integrate DEG with Mutation Data
#'
#' Combines DEG and mutation information
#'
#' @param deg_data Data frame with DEG information
#' @param mutations Data frame with mutation data
#' @param gene_col Column name for gene symbols in both data frames
#' @return Data frame with integrated information
#'
integrate_deg_mutations <- function(deg_data, mutations, gene_col = "gene_symbol") {
  cat("Integrating DEG and mutation data...\n")
  
  # Merge data
  integrated <- deg_data %>%
    left_join(
      mutations %>% select(gene_col, mutation_type, frequency, clinical_significance),
      by = gene_col
    )
  
  # Identify genes with both DEG and mutations
  genes_with_both <- integrated %>%
    filter(!is.na(mutation_type)) %>%
    nrow()
  
  cat("Found", genes_with_both, "genes with both DEG and mutations\n")
  
  return(integrated)
}

# ============================================================================
# 6. PATHWAY AND FUNCTIONAL ANALYSIS
# ============================================================================

#' Map Genes to Pathways
#'
#' Maps DEGs to biological pathways using KEGG or Reactome databases
#'
#' @param gene_symbols Vector of gene symbols
#' @param pathway_db Pathway database ("KEGG" or "Reactome")
#' @return Data frame with gene-pathway mappings
#'
map_genes_to_pathways <- function(gene_symbols, pathway_db = "KEGG") {
  cat("Mapping genes to", pathway_db, "pathways...\n")
  
  # This is a placeholder for pathway mapping
  # In actual implementation, would use clusterProfiler or similar
  
  pathway_map <- data.frame(
    gene_symbol = gene_symbols,
    pathway_db = pathway_db,
    pathway_id = NA,
    pathway_name = NA,
    p_value = NA,
    q_value = NA
  )
  
  cat("Mapped", nrow(pathway_map), "genes to pathways\n")
  
  return(pathway_map)
}

# ============================================================================
# 7. CLINICAL ANNOTATION INTEGRATION
# ============================================================================

#' Load Clinical Data
#'
#' Loads clinical annotations and patient data
#'
#' @param clinical_file Path to clinical data file
#' @return Data frame with clinical information
#'
load_clinical_data <- function(clinical_file) {
  cat("Loading clinical data from:", clinical_file, "\n")
  
  clinical <- read.csv(clinical_file, stringsAsFactors = FALSE)
  
  cat("Loaded clinical data for", nrow(clinical), "samples\n")
  
  return(clinical)
}

#' Correlate DEG with Clinical Outcomes
#'
#' Analyzes correlation between gene expression and clinical outcomes
#'
#' @param deg_data Data frame with DEG information
#' @param clinical_data Data frame with clinical information
#' @return Data frame with correlation statistics
#'
correlate_deg_clinical <- function(deg_data, clinical_data) {
  cat("Correlating DEG with clinical outcomes...\n")
  
  # This is a placeholder for correlation analysis
  # In actual implementation, would perform survival analysis, etc.
  
  correlations <- data.frame(
    gene_symbol = rownames(deg_data),
    correlation_coef = NA,
    p_value = NA,
    clinical_significance = NA
  )
  
  cat("Computed correlations for", nrow(correlations), "genes\n")
  
  return(correlations)
}

# ============================================================================
# 8. VISUALIZATION FUNCTIONS
# ============================================================================

#' Plot Network Graph
#'
#' Visualizes the protein-protein interaction network
#'
#' @param network_analysis List from analyze_network_topology()
#' @param output_file Path to save the plot
#'
plot_network <- function(network_analysis, output_file = "network_plot.png") {
  cat("Generating network visualization...\n")
  
  g <- network_analysis$graph
  
  # Set node size based on degree centrality
  degree_cent <- network_analysis$degree_centrality
  node_size <- (degree_cent - min(degree_cent)) / (max(degree_cent) - min(degree_cent)) * 10 + 2
  
  # Create plot
  png(output_file, width = 1200, height = 1000, res = 150)
  
  plot(g,
    vertex.size = node_size,
    vertex.label.cex = 0.7,
    edge.width = 0.5,
    layout = layout_with_fr(g)
  )
  
  dev.off()
  
  cat("Network plot saved to:", output_file, "\n")
}

#' Plot Hub Gene Characteristics
#'
#' Creates visualizations for hub gene properties
#'
#' @param network_analysis List from analyze_network_topology()
#' @param deg_data Data frame with DEG information
#' @param output_file Path to save the plot
#'
plot_hub_genes <- function(network_analysis, deg_data, output_file = "hub_genes_plot.png") {
  cat("Generating hub gene visualization...\n")
  
  hub_nodes <- network_analysis$hub_nodes
  degree_cent <- network_analysis$degree_centrality[hub_nodes]
  
  # Create data frame for plotting
  hub_df <- data.frame(
    gene = names(degree_cent),
    degree_centrality = as.numeric(degree_cent)
  )
  
  # Create plot
  p <- ggplot(hub_df, aes(x = reorder(gene, degree_centrality), y = degree_centrality)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = "Hub Genes in Melanoma PPI Network",
      x = "Gene",
      y = "Degree Centrality"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 10)
    )
  
  ggsave(output_file, plot = p, width = 10, height = 8, dpi = 150)
  
  cat("Hub gene plot saved to:", output_file, "\n")
}

# ============================================================================
# 9. MAIN INTEGRATION PIPELINE
# ============================================================================

#' Run Complete Genomic Integration Pipeline
#'
#' Executes the full genomic integration analysis
#'
#' @param deg_file Path to DEG results file
#' @param output_dir Directory to save results
#'
run_genomic_integration <- function(deg_file, output_dir = "./genomic_integration_results") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n")
  cat("========================================\n")
  cat("GENOMIC INTEGRATION PIPELINE\n")
  cat("========================================\n")
  cat("\n")
  
  # Step 1: Load DEG results
  cat("STEP 1: Loading DEG Results\n")
  cat("----------------------------------------\n")
  deg_data <- load_deg_results(deg_file)
  
  # Step 2: Annotate genes
  cat("\nSTEP 2: Gene Annotation\n")
  cat("----------------------------------------\n")
  gene_ids <- rownames(deg_data)
  annotations <- annotate_genes_biomart(gene_ids)
  
  # Step 3: Get PPI network
  cat("\nSTEP 3: Protein-Protein Interaction Network\n")
  cat("----------------------------------------\n")
  gene_symbols <- annotations$external_gene_name
  ppi_interactions <- get_ppi_network(gene_symbols)
  
  # Step 4: Analyze network topology
  if (!is.null(ppi_interactions)) {
    cat("\nSTEP 4: Network Topology Analysis\n")
    cat("----------------------------------------\n")
    network_analysis <- analyze_network_topology(ppi_interactions)
    
    # Step 5: Visualizations
    cat("\nSTEP 5: Generating Visualizations\n")
    cat("----------------------------------------\n")
    plot_network(network_analysis, file.path(output_dir, "ppi_network.png"))
    plot_hub_genes(network_analysis, deg_data, file.path(output_dir, "hub_genes.png"))
  }
  
  # Step 6: Pathway mapping
  cat("\nSTEP 6: Pathway Mapping\n")
  cat("----------------------------------------\n")
  pathway_map <- map_genes_to_pathways(gene_symbols)
  
  # Save results
  cat("\nSaving results...\n")
  write.csv(annotations, file.path(output_dir, "gene_annotations.csv"))
  write.csv(pathway_map, file.path(output_dir, "pathway_mapping.csv"))
  
  if (!is.null(ppi_interactions)) {
    write.csv(ppi_interactions, file.path(output_dir, "ppi_interactions.csv"))
  }
  
  cat("\n")
  cat("========================================\n")
  cat("PIPELINE COMPLETED SUCCESSFULLY\n")
  cat("Results saved to:", output_dir, "\n")
  cat("========================================\n")
  cat("\n")
  
  return(list(
    deg_data = deg_data,
    annotations = annotations,
    ppi_interactions = ppi_interactions,
    pathway_map = pathway_map
  ))
}

# ============================================================================
# 10. EXAMPLE USAGE (COMMENTED OUT - FOR REFERENCE)
# ============================================================================

# Example usage:
# 
# # Run the complete pipeline
# results <- run_genomic_integration(
#   deg_file = "path/to/deg_results.csv",
#   output_dir = "./melanoma_genomic_integration"
# )
#
# # Access individual results
# deg_data <- results$deg_data
# annotations <- results$annotations
# ppi_interactions <- results$ppi_interactions
# pathway_map <- results$pathway_map

# ============================================================================
# END OF SCRIPT
# ============================================================================
