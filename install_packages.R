if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("topGO")
BiocManager::install("DOSE")
BiocManager::install("org.Mm.eg.db")

