from itertools import chain, combinations
from os.path import abspath, join
import glob
import re
import operator
import pandas as pd

RAWDATA_DIR_RDS = abspath("../data/rds_filtered")
OUTPUT_DIR = abspath("../output/snakemake_output/")
SCRIPTS_DIR = abspath("../code")

NCELLS = [
    "1000",
    "2000",
    "5000",
    "10000",
    "25000",
    "50000",
    "100000",
    "200000",
    "300000",
    "Full",
]

ALL_DATASETS = glob.glob("{}/*.rds".format(RAWDATA_DIR_RDS))
SAMPLE_NAMES = ["Fetal__sci-RNA-seq3"]
# [x.split("/")[-1].replace(".rds", "") for x in ALL_DATASETS]


workdir: OUTPUT_DIR


os.makedirs(join(OUTPUT_DIR, "slurm-logs"), exist_ok=True)


rule all:
    input:
        expand(
            "sct_ncells_benchmarks/{sample_name}/{ncells}/times.csv",
            ncells=NCELLS,
            sample_name=SAMPLE_NAMES,
        ),


rule run_vst_ncells:
    input:
        rds=RAWDATA_DIR_RDS + "/" + "{sample_name}.rds",
        scriptx=SCRIPTS_DIR + "/" + "04_run_vst_ncells.R",
    output:
        "sct_ncells_benchmarks/{sample_name}/{ncells}/times.csv",
    params:
        ncells="{ncells}",
        prefix="sct_ncells_benchmarks/{sample_name}/{ncells}",
    shell:
        r"""Rscript --vanilla {input.scriptx} {input.rds} {params.ncells} {params.prefix}"""
