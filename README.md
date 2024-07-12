# Usage
Use the command line to run the clusterProfiler R package for rapid or batch GO and KEGG enrichment analysis, while also generating various scientific images.

Before running the analysis, various R packages that clusterProfiler depends on need to be installed, including clusterProfiler, topGO, DOSE, and annotation packages (depending on the species being analyzed), such as mmu (Mouse), hsa (Human), bta (Cow), cfa (Dog), cel (Worm), gga (Chicken), rno (Rat), tae (Wheat), ptr (Chimpanzee), dme (Fruit Fly), dre (Zebrafish), sce (Yeast), gma (Soybean), zma (Maize), ssc (Pig), ath (Arabidopsis), ecb (Horse), oar (Sheep), osa (Rice). If R has been added to the environment variables, you can use BiocManager to install directly them with the following command:
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("clusterProfiler")
    BiocManager::install("topGO")
    BiocManager::install("DOSE")
    BiocManager::install("org.Mm.eg.db") #depend the species analyzed




