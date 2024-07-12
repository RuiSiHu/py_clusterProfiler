# Usage
Use the command line to run the clusterProfiler R package for rapid or batch GO and KEGG enrichment analysis, while also generating various scientific images.

Before running the analysis, various R packages that clusterProfiler depends on need to be installed, including clusterProfiler, topGO, DOSE, and annotation packages (depending on the species being analyzed), such as mmu (Mouse), hsa (Human), bta (Cow), cfa (Dog), cel (Worm), gga (Chicken), rno (Rat), tae (Wheat), ptr (Chimpanzee), dme (Fruit Fly), dre (Zebrafish), sce (Yeast), gma (Soybean), zma (Maize), ssc (Pig), ath (Arabidopsis), ecb (Horse), oar (Sheep), osa (Rice). If R has been added to the environment variables, you can use BiocManager to install directly them with the following command:
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install("clusterProfiler")
    BiocManager::install("topGO")
    BiocManager::install("DOSE")
    BiocManager::install("org.Mm.eg.db") #depend the species analyzed

Or:
Rscript install_packages.R

All dependencies have been installed, the following command can be used to run:
    
    python3 py_clusterProfiler.py -i test.csv -fromType E -OrgDb mmu  -pvalueCutoff 0.01 -qvalueCutoff 0.01 -pAdjustMethod BH
    
Where:

    python3 py_clusterProfiler.py -h 
    
The format of the CSV file is as follows: the first column should be the ENSEMBL ID, and the second column should be the SYMBOL ID, represented as E and S respectively in the command line.

    ENSMUSG00000002257	Def6
    ENSMUSG00000036918	Ttc7
    ENSMUSG00000024066	Xdh
    ENSMUSG00000050107	Gsg2
    ENSMUSG00000032501	Trib1
    ENSMUSG00000017652	Cd40
    ENSMUSG00000031490	Eif4ebp1
    ENSMUSG00000022454	Nell2
    ENSMUSG00000031304	Il2rg
    ENSMUSG00000022218	Tgm1

py_clusterProfiler.py [-h] -i INPUT -fromType {E,S} -OrgDb
                             {mmu,hsa,bta,cfa,cel,gga,rno,tae,ptr,dme,dre,sce,gma,zma,ssc,ath,ecb,oar,osa}
                             [-pvalueCutoff PVALUECUTOFF] [-qvalueCutoff QVALUECUTOFF]
                             [-pAdjustMethod {BH,fdr,BY,holm,hochberg,hommel,bonferroni}]

Run GO and KEGG enrichment analysis.

Options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV file
  -fromType {E,S}       Type of gene IDs in input file (E for ENSEMBL, S for SYMBOL)
  -OrgDb {mmu,hsa,bta,cfa,cel,gga,rno,tae,ptr,dme,dre,sce,gma,zma,ssc,ath,ecb,oar,osa}
                        Organism database: mmu (Mouse), hsa (Human), bta (Cow), cfa (Dog), cel (Worm), gga (Chicken),
                        rno (Rat), tae (Wheat), ptr (Chimpanzee), dme (Fruit Fly), dre (Zebrafish), sce (Yeast), gma
                        (Soybean), zma (Maize), ssc (Pig), ath (Arabidopsis), ecb (Horse), oar (Sheep), osa (Rice)
  -pvalueCutoff PVALUECUTOFF
                        P-value cutoff for enrichment analysis (**default: 0.05**)
  -qvalueCutoff QVALUECUTOFF
                        Q-value cutoff for enrichment analysis (**default: 0.05**)
  -pAdjustMethod {BH,fdr,BY,holm,hochberg,hommel,bonferroni}
                        P-value adjustment method (**default: BH**): BH (Benjamini-Hochberg), fdr (False Discovery Rate,
                        equivalent to BH), BY (Benjamini-Yekutieli), holm (Holm-Bonferroni), hochberg (Hochberg),
                        hommel (Hommel), bonferroni (Bonferroni)

# Results
After run, the enrichment results for GO (CC, MF, and BP) and KEGG will be generated and presented in CSV file format. The image files will include dot plots and bubble plots for GO, and bar charts and bubble plots for KEGG.

Such as:
bubbleplot for GO enrichment (CC, MF, and BP)

![image](https://raw.githubusercontent.com/RuiSiHu/py_clusterProfiler/main/test/test_top_BP_MF_CC_bubbleplot.png)

barplot for GO enrichment (CC, MF, and BP)

![image](https://raw.githubusercontent.com/RuiSiHu/py_clusterProfiler/tree/main/test/test_top_BP_MF_CC_barplot.png)

dotplot for KEGG enrichment

![image](https://raw.githubusercontent.com/RuiSiHu/py_clusterProfiler/tree/main/test/test_KEGG_dotplot.png)

barplot for KEGG enrichment

![image](https://raw.githubusercontent.com/RuiSiHu/py_clusterProfiler/tree/main/test/test_KEGG_barplot.png)


# Reference
Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012, 16(5): 284-287. doi: 10.1089/omi.2011.0118.
Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu X, Liu S, Bo X, Yu G. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation (Camb). 2021, 2(3): 100141. doi: 10.1016/j.xinn.2021.100141.





