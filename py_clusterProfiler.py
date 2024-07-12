
import subprocess
import argparse
import os

def run_r_script(input_file, from_type, org_db, pvalue_cutoff, qvalue_cutoff, p_adjust_method):

    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = base_name

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    r_script_path = os.path.join("scripts", "clusterProfiler_analysis.R")
    rscript_path = r"C:\Program Files\R\R-4.4.0\bin\Rscript.exe"
    command = [
        rscript_path, r_script_path,
        input_file, from_type, org_db, str(pvalue_cutoff), str(qvalue_cutoff), p_adjust_method
    ]

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error in R script execution:")
        print(result.stderr)
    else:
        print(result.stdout)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run GO and KEGG enrichment analysis.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-fromType", required=True, choices=["E", "S"], help="Type of gene IDs in input file (E for ENSEMBL, S for SYMBOL)")
    parser.add_argument(
        "-OrgDb",
        required=True,
        choices=[
            "mmu", "hsa", "bta", "cfa", "cel", "gga", "rno", "tae", "ptr",
            "dme", "dre", "sce", "gma", "zma", "ssc", "ath", "ecb", "oar", "osa"
        ],
        help=(
            "Organism database: "
            "mmu (Mouse), hsa (Human), bta (Cow), cfa (Dog), cel (Worm), "
            "gga (Chicken), rno (Rat), tae (Wheat), ptr (Chimpanzee), "
            "dme (Fruit Fly), dre (Zebrafish), sce (Yeast), gma (Soybean), "
            "zma (Maize), ssc (Pig), ath (Arabidopsis), ecb (Horse), "
            "oar (Sheep), osa (Rice)"
        )
    )
    parser.add_argument("-pvalueCutoff", required=False, type=float, default=0.05, help="P-value cutoff for enrichment analysis (default: %(default)s)")
    parser.add_argument("-qvalueCutoff", required=False, type=float, default=0.05, help="Q-value cutoff for enrichment analysis (default: %(default)s)")
    parser.add_argument(
        "-pAdjustMethod",
        required=False,
        choices=["BH", "fdr", "BY", "holm", "hochberg", "hommel", "bonferroni"],
        default="BH",
        help=(
            "P-value adjustment method (default: %(default)s): "
            "BH (Benjamini-Hochberg), fdr (False Discovery Rate, equivalent to BH), "
            "BY (Benjamini-Yekutieli), holm (Holm-Bonferroni), hochberg (Hochberg), "
            "hommel (Hommel), bonferroni (Bonferroni)"
        )
    )

    args = parser.parse_args()
    run_r_script(args.input, args.fromType, args.OrgDb, args.pvalueCutoff, args.qvalueCutoff, args.pAdjustMethod)
