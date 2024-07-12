# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

necessary_packages <- c("clusterProfiler", "ggplot2", "topGO", "DOSE", "ggsci", "forcats")
for (pkg in necessary_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
from_type <- args[2]
org_db <- args[3]
pvalue_cutoff <- as.numeric(args[4])
qvalue_cutoff <- as.numeric(args[5])
p_adjust_method <- args[6]

# Set gene ID type and database based on parameters
if (from_type == "E") {
  gene_column <- 1
  from_type_r <- "ENSEMBL"
} else if (from_type == "S") {
  gene_column <- 2
  from_type_r <- "SYMBOL"
} else {
  stop("Invalid fromType. Use 'E' for ENSEMBL or 'S' for SYMBOL.")
}

# Dynamically load organism annotation database
load_org_db <- function(org_db) {
  org_db_mapping <- list(
    mmu = "org.Mm.eg.db",
    hsa = "org.Hs.eg.db",
    bta = "org.Bt.eg.db",
    cfa = "org.Cf.eg.db",
    cel = "org.Ce.eg.db",
    gga = "org.Gg.eg.db",
    rno = "org.Rn.eg.db",
    tae = "org.Ta.eg.db",
    ptr = "org.Pt.eg.db",
    dme = "org.Dm.eg.db",
    dre = "org.Dr.eg.db",
    sce = "org.Sc.sgd.db",
    gma = "org.Gmax.eg.db",
    zma = "org.Zm.eg.db",
    ssc = "org.Ss.eg.db",
    ath = "org.At.tair.db",
    ecb = "org.Ec.eg.db",
    oar = "org.Oa.eg.db",
    osa = "org.Os.eg.db"
  )
  
  if (!org_db %in% names(org_db_mapping)) {
    stop("Invalid OrgDb.")
  }
  
  pkg <- org_db_mapping[[org_db]]
  
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
  
  return(pkg)
}

org_db_r <- load_org_db(org_db)

# Extract file prefix and output directory
output_prefix <- tools::file_path_sans_ext(basename(input_file))
output_dir <- output_prefix

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Read CSV file, assuming data is in the input_file
data <- read.csv(input_file)
gene_list <- data[, gene_column]

# Gene ID conversion
gene_symbol <- bitr(geneID = gene_list,
                    fromType = from_type_r,
                    toType = c("ENTREZID"),
                    OrgDb = get(org_db_r))

# Extract converted ENTREZID gene list
gene <- gene_symbol[, 2]

# Check if any gene IDs were mapped
if (length(gene) == 0) {
  stop("No gene IDs were mapped to ENTREZID. Please check your input data and fromType parameter.")
}

# Define a function to perform enrichment analysis and save results and plots
perform_enrichment <- function(gene, ont, orgDb, output_prefix, pvalue_cutoff, qvalue_cutoff, p_adjust_method) {
  # GO enrichment analysis
  enrichment_result <- enrichGO(gene = gene,
                                keyType = "ENTREZID",
                                OrgDb = orgDb,
                                ont = ont,
                                pvalueCutoff = pvalue_cutoff,
                                pAdjustMethod = p_adjust_method,
                                minGSSize = 1,
                                maxGSSize = 500,
                                qvalueCutoff = qvalue_cutoff,
                                readable = TRUE)
  
  # Check if there are results, if not create NA entry
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    result_df <- data.frame(
      ID = NA, Description = NA, GeneRatio = NA, BgRatio = NA, pvalue = NA, p.adjust = NA, qvalue = NA,
      geneID = NA, Count = NA, Category = ont
    )
  } else {
    result_df <- enrichment_result@result
    result_df$Category <- ont
  }
  
  # Save GO enrichment results to CSV file
  write.csv(result_df, file.path(output_dir, paste0(output_prefix, "_", ont, ".csv")), row.names = FALSE)
  
  # Plot GO enrichment results as dotplot
  if (!is.na(result_df$ID[1])) {
    dotplot_plot <- dotplot(enrichment_result,
                            x = "GeneRatio",
                            color = "p.adjust",
                            showCategory = 10,
                            size = NULL,
                            title = paste0(ont, "_dotplot"))
    ggsave(file.path(output_dir, paste0(output_prefix, "_dotplot_", ont, ".pdf")), plot = dotplot_plot, width = 8, height = 8)
    ggsave(file.path(output_dir, paste0(output_prefix, "_dotplot_", ont, ".png")), plot = dotplot_plot, width = 8, height = 8)
  } else {
    empty_plot <- ggplot() + 
      geom_text(aes(0.5, 0.5, label = "No significant enrichment")) + 
      ggtitle(paste0(ont, "_dotplot")) +
      theme_void() +
      theme(plot.background = element_rect(color = "black", size = 1))
    ggsave(file.path(output_dir, paste0(output_prefix, "_dotplot_", ont, ".pdf")), plot = empty_plot, width = 8, height = 8)
    ggsave(file.path(output_dir, paste0(output_prefix, "_dotplot_", ont, ".png")), plot = empty_plot, width = 8, height = 8)
  }
  
  # Plot GO enrichment results as barplot
  if (!is.na(result_df$ID[1])) {
    barplot_plot <- barplot(enrichment_result,
                            x = "Count",
                            color = "p.adjust",
                            showCategory = 10,
                            size = NULL,
                            title = paste0(ont, "_barplot"))
    ggsave(file.path(output_dir, paste0(output_prefix, "_barplot_", ont, ".pdf")), plot = barplot_plot, width = 8, height = 8)
    ggsave(file.path(output_dir, paste0(output_prefix, "_barplot_", ont, ".png")), plot = barplot_plot, width = 8, height = 8)
  } else {
    empty_plot <- ggplot() + 
      geom_text(aes(0.5, 0.5, label = "No significant enrichment")) + 
      ggtitle(paste0(ont, "_barplot")) +
      theme_void() +
      theme(plot.background = element_rect(color = "black", size = 1))
    ggsave(file.path(output_dir, paste0(output_prefix, "_barplot_", ont, ".pdf")), plot = empty_plot, width = 8, height = 8)
    ggsave(file.path(output_dir, paste0(output_prefix, "_barplot_", ont, ".png")), plot = empty_plot, width = 8, height = 8)
  }
  
  return(result_df)
}

# Perform GO enrichment analysis for CC, MF, and BP
CC_result <- perform_enrichment(gene, "CC", get(org_db_r), output_prefix, pvalue_cutoff, qvalue_cutoff, p_adjust_method)
MF_result <- perform_enrichment(gene, "MF", get(org_db_r), output_prefix, pvalue_cutoff, qvalue_cutoff, p_adjust_method)
BP_result <- perform_enrichment(gene, "BP", get(org_db_r), output_prefix, pvalue_cutoff, qvalue_cutoff, p_adjust_method)

# Combine all GO results for overall enrichment analysis and save results and plots
all_go_results <- rbind(CC_result, MF_result, BP_result)
write.csv(all_go_results, file.path(output_dir, paste0(output_prefix, "_BP_MF_CC_combined.csv")), row.names = FALSE)

# Convert GeneRatio from fraction to decimal and round to 2 decimal places
all_go_results$GeneRatio <- sapply(all_go_results$GeneRatio, function(x) round(eval(parse(text = x)), 2))

# Save the modified CSV file
write.csv(all_go_results, file.path(output_dir, paste0(output_prefix, "_BP_MF_CC_combined_decimal.csv")), row.names = FALSE)

# Select top 10 significant terms from each category
top_CC <- head(CC_result[order(CC_result$pvalue), ], 10)
top_MF <- head(MF_result[order(MF_result$pvalue), ], 10)
top_BP <- head(BP_result[order(BP_result$pvalue), ], 10)
top_results <- rbind(top_CC, top_MF, top_BP)

# Create combined GO bar plot
ggplot(top_results) + 
  geom_bar(aes(Count, fct_reorder(Description, Category), fill = Category), stat = "identity") + 
  theme_bw() + 
  theme(text = element_text(size = 8), 
        axis.text = element_text(size = 8, colour = 'black')) + 
  geom_text(data = subset(top_results, p.adjust < 0.05, select = c('Count', 'Description')), 
            aes(x = Count + 3, y = as.factor(Description), label = '*')) + 
  labs(x = 'Number of genes', y = 'GO term') + 
  scale_fill_manual(values = c(BP = "#79B494", CC = "#D67E56", MF = "#848CBD")) +
  theme(axis.text = element_text(size = 8))  # Adjust font size for long descriptions
ggsave(file.path(output_dir, paste0(output_prefix, "_top_BP_MF_CC_barplot.pdf")), width = 8, height = 8)
ggsave(file.path(output_dir, paste0(output_prefix, "_top_BP_MF_CC_barplot.png")), width = 8, height = 8)

# Create combined GO bubble plot
goinput <- read.csv(file.path(output_dir, paste0(output_prefix, "_BP_MF_CC_combined_decimal.csv")))

# Select top 10 from each category
goinput <- rbind(
  head(goinput[goinput$Category == "CC", ], 10),
  head(goinput[goinput$Category == "MF", ], 10),
  head(goinput[goinput$Category == "BP", ], 10)
)

x <- goinput$GeneRatio
y <- factor(goinput$Description, levels = goinput$Description)
p <- ggplot(goinput, aes(x, y))
p1 <- p + geom_point(aes(size = Count, color = -0.5 * log(p.adjust), shape = Category)) +
  scale_color_gradient(low = "SpringGreen", high = "DeepPink") + theme_minimal()
p2 <- p1 + labs(color = expression(-log[10](p.adjust)),
                size = "Count",
                x = "GeneRatio",
                y = "GO term",
                title = "GO enrichment") + 
                theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave(file.path(output_dir, paste0(output_prefix, "_top_BP_MF_CC_bubbleplot.pdf")), plot = p2, width = 10, height = 8)
ggsave(file.path(output_dir, paste0(output_prefix, "_top_BP_MF_CC_bubbleplot.png")), plot = p2, width = 10, height = 8)

# Perform KEGG enrichment analysis
KEGG <- tryCatch({
  enrichKEGG(gene = gene,
             organism = org_db,
             keyType = "ncbi-geneid",
             minGSSize = 1,
             maxGSSize = 500,
             pvalueCutoff = pvalue_cutoff,
             pAdjustMethod = p_adjust_method,
             qvalueCutoff = qvalue_cutoff)
}, error = function(e) {
  NULL
})

# Check if KEGG enrichment results are available
if (is.null(KEGG) || nrow(KEGG@result) == 0) {
  KEGG_result <- data.frame(
    ID = NA, Description = NA, GeneRatio = NA, BgRatio = NA, pvalue = NA, p.adjust = NA, qvalue = NA,
    geneID = NA, Count = NA
  )
  warning("No significant KEGG enrichment results found.")
  
  # Plot empty KEGG enrichment results as dotplot
  empty_plot <- ggplot() + 
    geom_text(aes(0.5, 0.5, label = "No significant enrichment")) + 
    ggtitle("KEGG_dotplot") +
    theme_void() +
    theme(plot.background = element_rect(color = "black", size = 1))
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_dotplot.pdf")), plot = empty_plot, width = 8, height = 8)
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_dotplot.png")), plot = empty_plot, width = 8, height = 8)
  
  # Plot empty KEGG enrichment results as barplot
  empty_plot <- ggplot() + 
    geom_text(aes(0.5, 0.5, label = "No significant enrichment")) + 
    ggtitle("KEGG_barplot") +
    theme_void() +
    theme(plot.background = element_rect(color = "black", size = 1))
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_barplot.pdf")), plot = empty_plot, width = 8, height = 8)
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_barplot.png")), plot = empty_plot, width = 8, height = 8)
} else {
  # Display KEGG enrichment results
  KEGG_result <- KEGG@result
  print(KEGG_result)
  
  # Save KEGG enrichment results to CSV file
  write.csv(KEGG_result, file.path(output_dir, paste0(output_prefix, "_KEGG.csv")), row.names = FALSE)
  
  # Plot KEGG enrichment results as dotplot
  kegg_dotplot <- dotplot(KEGG,
                          x = "GeneRatio",
                          color = "p.adjust",
                          showCategory = 10,
                          size = NULL,
                          title = "KEGG_dotplot")
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_dotplot.pdf")), plot = kegg_dotplot, width = 8, height = 8)
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_dotplot.png")), plot = kegg_dotplot, width = 8, height = 8)
  
  # Plot KEGG enrichment results as barplot
  kegg_barplot <- barplot(KEGG,
                          x = "Count",
                          color = "p.adjust",
                          showCategory = 10,
                          size = NULL,
                          title = "KEGG_barplot")
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_barplot.pdf")), plot = kegg_barplot, width = 8, height = 8)
  ggsave(file.path(output_dir, paste0(output_prefix, "_KEGG_barplot.png")), plot = kegg_barplot, width = 8, height = 8)
}

# Save KEGG enrichment results to CSV file with NA if no significant enrichment
write.csv(KEGG_result, file.path(output_dir, paste0(output_prefix, "_KEGG.csv")), row.names = FALSE)
