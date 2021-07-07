from itertools import chain, combinations
from os.path import abspath, join
import glob
import re
import operator
import pandas as pd

sample_name = "PBMC__Smart-seq3_sampled_counts"
RAWDATA_DIR_RDS = abspath("../data/sampled_counts/{}".format(sample_name))
OUTPUT_DIR = abspath("../output/snakemake_output/{}".format(sample_name))

SCRIPTS_DIR = abspath("../code")

NCELLS = [1000]
METHODS = ["vst2"]

ALL_DATASETS = glob.glob("{}/*.rds".format(RAWDATA_DIR_RDS))
SAMPLE_NAMES = [x.split("/")[-1].replace(".rds", "") for x in ALL_DATASETS]


workdir: OUTPUT_DIR


os.makedirs(join(OUTPUT_DIR, "slurm-logs"), exist_ok=True)


rule all:
    input:
        expand(
            "poisson_glm_output/{ncells}/{sample_name}/model_fit.csv",
            sample_name=SAMPLE_NAMES,
            ncells=NCELLS,
        ),
        expand(
            "seurat_output/{sample_name}/{method}/gene_attr.csv",
            sample_name=SAMPLE_NAMES,
            method=METHODS,
        ),


rule run_poisson_glm:
    input:
        rds=RAWDATA_DIR_RDS + "/" + "{sample_name}.rds",
        scriptx=SCRIPTS_DIR + "/" + "01_run_glm.R",
    output:
        "poisson_glm_output/{ncells}/{sample_name}/model_fit.csv",
    params:
        ncells="{ncells}",
        glmtype="poisson",
    shell:
        r"""Rscript --vanilla {input.scriptx} {input.rds} {params.ncells} {params.glmtype} {output}"""

rule run_seurat:
    input:
        rds=RAWDATA_DIR_RDS + "/" + "{sample_name}.rds",
        scriptx=SCRIPTS_DIR + "/" + "02_run_seurat.R",
    output:
        objects="seurat_output/{sample_name}/{method}/gene_attr.csv",
    params:
        outputprefix="seurat_output/{sample_name}/{method}",
        method="{method}",
    shell:
        r"""Rscript --vanilla {input.scriptx} {input.rds} {params.method} {params.outputprefix}"""
