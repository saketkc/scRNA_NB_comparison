from itertools import chain, combinations
from os.path import abspath, join
import glob
import re
import operator
import pandas as pd

RAWDATA_DIR_RDS = abspath("../data/rds_filtered")
OUTPUT_DIR = abspath("../output/snakemake_output")
SCRIPTS_DIR = abspath("../code")


SAMPLE_NAMES = ["PBMC__Smart-seq3"]



workdir: OUTPUT_DIR


os.makedirs(join(OUTPUT_DIR, "slurm-logs"), exist_ok=True)


KNNS = [200, 300, 400, 500]


rule all:
    input:
        expand(
            "metacell_knn_seurat_output/{sample_name}/{knn}/seurat_sct_object.rds",
            sample_name=SAMPLE_NAMES,
            knn=KNNS,
        ),


rule run_metacell:
    input:
        rds=RAWDATA_DIR_RDS + "/" + "{sample_name}.rds",
        scriptx=SCRIPTS_DIR + "/" + "05_run_metacell.R",
    output:
        raw="metacell_knn_output/{sample_name}/{knn}/seurat_metacell_aggregatedcounts_raw.rds",
        filtered="metacell_knn_output/{sample_name}/{knn}/seurat_metacell_aggregatedcounts_filtered.rds",
        subset="metacell_knn_output/{sample_name}/{knn}/seurat_metacell_filtered_subset.rds",
    params:
        prefix="metacell_knn_output/{sample_name}/{knn}",
        knn="{knn}",
    shell:
        r"""Rscript --vanilla {input.scriptx} {input.rds} {params.knn} {params.prefix}"""


rule run_seurat_metacell_filtered:
    input:
        rds="metacell_knn_output/{sample_name}/{knn}/seurat_metacell_aggregatedcounts_filtered.rds",
        scriptx=SCRIPTS_DIR + "/" + "06_run_sct.R",
    output:
        objects="metacell_knn_seurat_output/{sample_name}/{knn}/seurat_sct_object.rds",
    params:
        outputprefix="metacell_knn_seurat_output/{sample_name}/{knn}",
        method="vst2",
    shell:
        r"""Rscript --vanilla {input.scriptx} {input.rds} {params.method} {params.outputprefix}"""
