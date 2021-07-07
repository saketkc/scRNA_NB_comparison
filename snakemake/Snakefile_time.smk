from itertools import chain, combinations
from os.path import abspath, join
import glob
import re
import operator
import pandas as pd

RAWDATA_DIR_RDS = abspath("../data/rds_filtered")
OUTPUT_DIR = abspath("../output/snakemake_output/")
SCRIPTS_DIR = abspath("../code")

SEEDS = [20141015, 20160823, 20161003, 20170523, 20170726]

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
]

ALL_DATASETS = glob.glob("{}/*.rds".format(RAWDATA_DIR_RDS))
SAMPLE_NAMES = ["Fetal__sci-RNA-seq3"]


workdir: OUTPUT_DIR


os.makedirs(join(OUTPUT_DIR, "slurm-logs"), exist_ok=True)


rule all:
    input:
        expand(
            "vst2_time_benchmarks/{sample_name}/{ncells}/{seed}/times.csv",
            ncells=NCELLS,
            seed=SEEDS,
            sample_name=SAMPLE_NAMES,
        ),


rule run_vst:
    input:
        rds=RAWDATA_DIR_RDS + "/" + "{sample_name}.rds",
        scriptx=SCRIPTS_DIR + "/" + "03_run_vst2_downsample.R",
    output:
        "vst2_time_benchmarks/{sample_name}/{ncells}/{seed}/times.csv",
    params:
        ncells="{ncells}",
        seed="{seed}",
        prefix="vst2_time_benchmarks/{sample_name}/{ncells}/{seed}",
    shell:
        r"""Rscript --vanilla {input.scriptx} {input.rds} {params.seed} {params.ncells} {params.prefix}"""
