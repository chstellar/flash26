"""
metaSPLASH pipeline (for bacterial genomes)
kmer-based embeddings + predictions

Daniel Cotter
10/2024

This pipeline is designed to take SPLASH results and use them to predict on phenotypic metadata.
The pipeline is designed to be run on Sherlock, and is written in snakemake.
"""

## Importing necessary modules
from pathlib import Path
import pandas as pd
import csv

## Define the config files for the pipeline
configfile: "config.yaml"
metadata_table_path = Path(config["data_table"])
TEMP_DIR = Path(config["temp_dir"])

## read the metadata table with the short names as indices
metadata_table = pd.read_csv(metadata_table_path, index_col = "dataset_short_name")

## Define the wildcards on which the pipeline will be run
# TODO: Dynamically generate {dataset} based on the metadata table
# TODO: Define the other wildcards based on the config file
DATASETS = list(metadata_table.index)
#DATASETS = ["eFaecium-CollEtAl"]
SELECT_TYPES = ["filter1", "filter2", "filter3"]
#SELECT_TYPES = ["filter1"]
CLUSTER_TYPES = ["shiftDist-keepTopES", "shiftDist-keepMostAbundant", "shiftDist-levFilter", "mmseqs-levFilter"]
#CLUSTER_TYPES = ["shiftDist-keepTopES"]
NUM_CLUSTERS = [4000]
KMER_WIDTH = [54]
KMER_STEP = [54]
MODELS = ["esm"]
NORMALIZE = ["normalized", "unnormalized"]

## constrain the wildcards of the pipeline
wildcard_constraints:
    dataset=r"[A-Za-z\d-]+",
    model=r"[A-Za-z\d-]+",
    select_type=r"[A-Za-z\d-]+",
    cluster_type=r"[A-Za-z\d-]+",
    num_clusters=r"\d+",
    kmer_width=r"\d+",
    kmer_step=r"\d+",
    normalize=r"[A-Za-z]+"


## Define the rules for the pipeline
rule all:
    input:
        # all nonzero coefficients files
        expand(Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", 
                    "{dataset}_{model}_glmnet_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_{FILE}"),
               dataset=DATASETS,
               select_type=SELECT_TYPES,
               cluster_type=CLUSTER_TYPES,
               model=MODELS,
               num_clusters=NUM_CLUSTERS,
               kmer_width=KMER_WIDTH,
               kmer_step=KMER_STEP,
               normalize=NORMALIZE,
               FILE = ["nonzero_coefficients.tsv", "confusion_matrices.pdf"]) ,
        # all kmer mapping to clusters files
        expand(Path("results", "{dataset}", "{select_type}", "{cluster_type}", 
                    "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_k{kmer_width}_s{kmer_step}.tsv"),
               dataset=DATASETS,
               select_type=SELECT_TYPES,
               cluster_type=CLUSTER_TYPES,
               num_clusters=NUM_CLUSTERS,
               kmer_width=KMER_WIDTH,
               kmer_step=KMER_STEP),
        # all ohe glmnet results
        expand(Path("results", "{dataset}", "{select_type}", "{cluster_type}", "ohe", 
                    "{dataset}_ohe_glmnet_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_{FILE}"),
               dataset=DATASETS,
               select_type=SELECT_TYPES,
               cluster_type=CLUSTER_TYPES,
               num_clusters=NUM_CLUSTERS,
               kmer_width=KMER_WIDTH,
               kmer_step=KMER_STEP,
               FILE = ["nonzero_coefficients.tsv", "confusion_matrices.pdf"])


rule choose_anchors:
    """
    This rule selects the top anchors from the SPLASH results based on various criteria (select_type)
    New scripts that select different anchors can be added to the config file under the "anchor_select_script" key
    """
    input:
        lambda wildcards: Path(metadata_table.loc[wildcards.dataset, "SPLASH_results"],
                               "result.after_correction.scores.tsv") # this is the path to the SPLASH results
    params:
        script = lambda wildcards: Path(config["scripts"]["anchor_select_script"][wildcards.select_type]),
        lookup_table = lambda wildcards: Path(metadata_table.loc[wildcards.dataset, "lookup_table"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type)
    threads: 2
    resources:
        # 64 GB of memory
        mem_mb = 64000
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}.txt")
    shell:"""
        ml R/4.3.2
        Rscript --vanilla {params.script} --input {input} --output {output} \
        --lookup_table {params.lookup_table} --temp_dir {params.tmp_dir}
    """


rule cluster_anchors:
    """
    This rule clusters the selected anchors based on the selected clustering method (cluster_type)
    New scripts that cluster the anchors can be added to the config file under the "cluster_script" key
    This rule has a mix of R and Python code, so it is necessary to load all appropriate modules
    The call to the python or Rscript is stored in the config file under the "cluster_script" key
    """
    input:
        Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}.txt")
    params:
        script = lambda wildcards: Path(config["scripts"]["cluster_script"][wildcards.cluster_type]),
        python_env = Path(config["envs"]["default_python"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type)
    threads: 4
    resources:
        # dynamically allocate memory based on the attempt
        mem_mb = lambda _, attempt: 32000 + ((attempt - 1) * 32000),
        time = "3:00:00"
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_clustered_anchors_{select_type}_{cluster_type}.txt")
    shell:"""
        ml R/4.3.2
        ml python/3.9.0
        source {params.python_env}
        {params.script} --input {input} --output {output} --temp_dir {params.tmp_dir}
    """


rule reorder_clusters:
    """
    This rule reorders the clusters based on the SPLASH results and the selected anchors and clustering method
    One reordering example is to sort the clusters and only grab 1 anchor per cluster.
    New scripts that reorder the clusters can be added to the config file under the "reorder_script" key
    """
    input:
        clusters = Path(TEMP_DIR, "{dataset}", "{dataset}_clustered_anchors_{select_type}_{cluster_type}.txt"),
        splash_results = lambda wildcards: Path(metadata_table.loc[wildcards.dataset, "SPLASH_results"],
                                                "result.after_correction.scores.tsv")
    params:
        script = lambda wildcards: Path(config["scripts"]["reorder_script"][wildcards.cluster_type]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type)
    threads: 4
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_reordered_clusters_{select_type}_{cluster_type}.txt")
    shell:"""
        ml R/4.3.2
        Rscript --vanilla {params.script} --input_anchor_clusters {input.clusters} \
        --splash_stats {input.splash_results} --output {output} --temp_dir {params.tmp_dir} \
        --num_cores {threads}
    """

rule select_N_clusters:
    """
    After reordering the clusters, this rule selects the top N clusters based on their id (since they are already reordered)
    """
    input:
        Path(TEMP_DIR, "{dataset}", "{dataset}_reordered_clusters_{select_type}_{cluster_type}.txt")
    output:
        clusters = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_clusters_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt"),
        anchors = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt")
    shell:"""
        awk '$1 <= {wildcards.num_clusters}' {input} > {output.clusters}
        cut -f2 {output.clusters} > {output.anchors}
    """


rule prepare_sequences:
    """
    This rule prepares the sequences for the selected clusters and anchors for the dataset.
    It formats the sequences into concatmers for each sample where all anchor-target pairs are concatenated
    together. The output is a fasta file and a tsv file with the sequences.
    """
    input:
        cluster_file = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_clusters_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt"),
        anchor_file = Path(TEMP_DIR, "{dataset}", "{dataset}_selected_anchors_{select_type}_{cluster_type}_top{num_clusters}-clusters.txt"),
        id_mapping = lambda wildcards: Path(metadata_table.loc[wildcards.dataset, "SPLASH_results"], "sample_name_to_id.mapping.txt")
    params:
        script = Path(config["scripts"]["prepare_sequences"]),
        satc_dir = lambda wildcards: Path(metadata_table.loc[wildcards.dataset, "SPLASH_results"], "result_satc"),
        output_prefix = lambda wildcards: Path(TEMP_DIR, f"{wildcards.dataset}", f"{wildcards.dataset}_prepared_sequences_{wildcards.select_type}_{wildcards.cluster_type}_top{wildcards.num_clusters}"),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters"),
    threads: 16
    resources:
        # 128 GB of memory
        mem_mb = 128000
    output:
        fasta = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_sample_sequences.fasta"),
        tsv = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_sample_sequences.tsv")
    shell:"""
        ml R/4.3.2
        Rscript --vanilla {params.script} --anchor_file {input.anchor_file} \
        --cluster_file {input.cluster_file} --id_mapping {input.id_mapping} \
        --satc_files {params.satc_dir} --output_prefix {params.output_prefix} \
        --temp_dir {params.tmp_dir} --num_cores {threads}
    """


rule decompose_kmers:
    """
    Process the sample sequences to decompose them into kmers of width kmer_width and step kmer_step
    The outputs are 1) a fasta file with the unique kmers and 2) a tsv file with the ordering of the kmers
    for each sample
    """
    input:
        Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_sample_sequences.fasta")
    params:
        script = Path(config["scripts"]["decompose_kmers"]),
        output_prefix = lambda wildcards: Path(TEMP_DIR, f"{wildcards.dataset}", f"{wildcards.dataset}_decomposed_kmers_{wildcards.select_type}_{wildcards.cluster_type}_top{wildcards.num_clusters}"),
        kmer_width = lambda wildcards: wildcards.kmer_width,
        kmer_step = lambda wildcards: wildcards.kmer_step,
        python_env = Path(config["envs"]["default_python"])
    output:
        unique_kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta"),
        order = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv")
    resources:
        # dynamically allocate memory based on the attempt
        mem_mb = lambda _, attempt: 16000 + ((attempt - 1) * 16000),
        time = "3:00:00"
    shell:"""
        ml python/3.9.0
        source {params.python_env}
        python {params.script} -k {wildcards.kmer_width} -s {wildcards.kmer_step} \
        {input} {params.output_prefix}
    """


rule match_kmers_to_clusters:
    """
    Process the ordering file to produce a tsv file mapping clusters to
    their component kmers
    """
    input:
        order = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv"),
        unique_kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta")
    params:
        script = Path(config["scripts"]["match_kmers_to_clusters"])
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{dataset}_sequences_per_cluster_top{num_clusters}-clusters_k{kmer_width}_s{kmer_step}.tsv")
    shell:"""
    ml R/4.3.2
    Rscript --vanilla {params.script} --ordering {input.order} --kmers {input.unique_kmers} --output {output}
    """

rule translate_kmers_ESM:
    """
    This rule translates the kmers using a provided genetic code so that they can be fed into the
    ESM2 model (or any other protein-based language model).
    """
    input:
        unique_kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta")
    params:
        script = Path(config["scripts"]["translate_script"]),
        translation_table = lambda wildcards: metadata_table.loc[wildcards.dataset, "translation_table"],
        python_env = Path(config["envs"]["default_python"])
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_translated_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta")
    shell:"""
        ml python/3.9.0
        source {params.python_env}
        python {params.script} -t {params.translation_table} {input} {output}
    """


rule embed_kmers_ESM:
    """
    This rule embeds the TRANSLATED kmers into a pre-trained language model to get averaged embeddings 
    for each kmer. Downstream, these kmers are recombined into their order and used as features
    to predict on the metadata.
    """
    input:
        translated_kmers = Path(TEMP_DIR, "{dataset}", "{dataset}_translated_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_unique_kmers.fasta"),
    params:
        torch_dir = Path(TEMP_DIR, "torch_cache"),
        extract_embeddings = Path(config["scripts"]["extract_embeddings"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters", 
                                         "k" + wildcards.kmer_width + "_s" + wildcards.kmer_step, "esm_embeddings", "raw_embeddings"),
        python_env = Path(config["envs"]["esm_env"])
    threads: 8
    resources:
        # 64 GB of memory
        time = "3:00:00",
        mem_mb = 32000,
        partition = "gpu,owners",
        slurm_extra = "-G 1 -C 'GPU_GEN:AMP|GPU_GEN:VLT|GPU_GEN:TUR'"
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_esm-embeddings_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}.tsv")
    shell:"""
        ml python/3.9.0
        source {params.python_env}
        export TORCH_HOME={params.torch_dir}
        esm-extract esm2_t33_650M_UR50D {input} {params.tmp_dir} --include mean per_tok
        python {params.extract_embeddings} {params.tmp_dir} {output}
    """


rule prepare_data_for_glmnet_top_variance:
    """
    This rule processes the embeddings to fit into a glmnet model by grabbing the top variance embeddings per cluster
    and then saves the resulting data frame as a feather object to be used in the glmnet model.
    """
    input:
        embeddings = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_{model}-embeddings_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}.tsv"),
        ordering = lambda wildcards: Path(TEMP_DIR, "{dataset}", "{dataset}_decomposed_kmers_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_kmer_ordering.tsv")
    params:
        script = Path(config["scripts"]["format_embeddings_variance"]),
        tmp_dir = lambda wildcards: Path(TEMP_DIR, wildcards.dataset, wildcards.select_type, wildcards.cluster_type, wildcards.num_clusters + "-clusters", 
                                         "k" + wildcards.kmer_width + "_s" + wildcards.kmer_step, wildcards.model + "_embeddings", wildcards.normalize),
        normalized_flag = lambda wildcards: "--normalized_embeddings" if wildcards.normalize =="normalized" else ""
    threads: 32
    resources:
        # dynamically allocate memory based on the attempt
        mem_mb = lambda _, attempt: 64000 + ((attempt - 1) * 64000),
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_{normalize}.feather")
    shell:"""
        ml R/4.3.2
        Rscript --vanilla {params.script} --embeddings {input.embeddings} --ordering {input.ordering} \
        --output {output} --temp_dir {params.tmp_dir} --num_threads {threads} --num_to_keep 100 {params.normalized_flag}
    """

rule run_glmnet:
    """
    This rule uses preprocessed embeddings for running the glmnet model to predict on the metadata.
    """
    input:
        embeddings = Path(TEMP_DIR, "{dataset}", "{dataset}_{model}_top_variance_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}_{normalize}.feather"),
        metadata = lambda wildcards: metadata_table.loc[wildcards.dataset, "metadata_file"]
    params:
        script = Path(config["scripts"]["glmnet_script"]),
        output_prefix = lambda wildcards: Path("results", f"{wildcards.dataset}", f"{wildcards.select_type}", f"{wildcards.cluster_type}", f"{wildcards.model}", f"{wildcards.normalize}", f"{wildcards.dataset}_{wildcards.model}_glmnet_results_top{wildcards.num_clusters}_k{wildcards.kmer_width}_s{wildcards.kmer_step}")
    threads: 16
    resources:
        # dynamically allocate memory based on the attempt
        mem_mb = lambda _, attempt: 256000 + ((attempt - 1) * 64000),
        time = "4:00:00"
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_glmnet_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_nonzero_coefficients.tsv"),
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "{model}", "{normalize}", "{dataset}_{model}_glmnet_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_confusion_matrices.pdf"),
    shell:"""
        ml R/4.3.2
        Rscript --vanilla {params.script} --embeddings {input.embeddings} \
        --metadata {input.metadata} --output_prefix {params.output_prefix} \
        --even_classes
    """


rule prepare_data_for_glmnet_ohe:
    input:
        sample_sequences = Path(TEMP_DIR, "{dataset}", "{dataset}_prepared_sequences_{select_type}_{cluster_type}_top{num_clusters}_sample_sequences.tsv")
    params:
        script = Path(config["scripts"]["format_sequences_ohe"]),
        kmer_width = lambda wildcards: wildcards.kmer_width,
        kmer_step = lambda wildcards: wildcards.kmer_step
    output:
        Path(TEMP_DIR, "{dataset}", "{dataset}_ohe_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}.feather")
    threads: 4
    shell:"""
        ml R/4.3.2
        Rscript --vanilla {params.script} --input {input} --output {output} --kmer_width {params.kmer_width}
    """

rule run_glmnet_ohe:
    input:
        features = Path(TEMP_DIR, "{dataset}", "{dataset}_ohe_features_for_glmnet_{select_type}_{cluster_type}_top{num_clusters}_k{kmer_width}_s{kmer_step}.feather"),
        metadata = lambda wildcards: metadata_table.loc[wildcards.dataset, "metadata_file"]
    params:
        script = Path(config["scripts"]["glmnet_script"]),
        output_prefix = lambda wildcards: Path("results", f"{wildcards.dataset}", f"{wildcards.select_type}", f"{wildcards.cluster_type}", f"ohe", f"{wildcards.dataset}_ohe_glmnet_results_top{wildcards.num_clusters}_k{wildcards.kmer_width}_s{wildcards.kmer_step}")
    threads: 16
    resources:
        # dynamically allocate memory based on the attempt
        mem_mb = lambda _, attempt: 256000 + ((attempt - 1) * 64000),
        time = "4:00:00"
    output:
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_glmnet_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_nonzero_coefficients.tsv"),
        Path("results", "{dataset}", "{select_type}", "{cluster_type}", "ohe", "{dataset}_ohe_glmnet_results_top{num_clusters}_k{kmer_width}_s{kmer_step}_confusion_matrices.pdf"),
    shell:"""
        ml R/4.3.2
        Rscript --vanilla {params.script} --embeddings {input.features} \
        --metadata {input.metadata} --output_prefix {params.output_prefix} \
        --even_classes
    """