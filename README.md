# *metaSPLASH* Pipeline

Daniel Cotter -- 10/07/2024

## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction and attribution on SPLASH results.

### Different filters used for anchor selection3

These filters are used to select the set of all anchors from the SPLASH results that will be sent for clustering and reordering. See `config.yaml` for the script used for the filtering.

#### Filter 1

1. effect size >= 0.6
2. number of nonzero samples > 10th percentile
3. no lookup table hits to artifacts
4. select top 150,000 anchors by number of nonzero samples

#### Filter 2

1. effect size >= 0.6
2. number of nonzero samples > 10th percentile
3. no lookup table hits to artifacts
4. select top 150,000 anchors by effect size

#### Filter 3

1. effect size >= 0.9
2. no lookup table hits to artifacts
3. further filter out the bottom 60% of anchors by number nonzero samples
4. filter 3: further filter out the bottom 60% of anchors by effect size
5. filter 3: select top 100,000 anchors by effect size * number of nonzero samples

### Different clustering methods

These are different methods for clustering the anchors from step 1 and then informed reordering/pruning of the clusters. See `config.yaml` for the scripts used for these steps.

#### shiftDist-keepTopES

1. Cluster using Tavor's shift distance code. shiftDist = 5
2. Arrange each cluster by effect size bin and keep the **one** highest effect size anchor per cluster.

#### shiftDist-keepMostAbundant

1. Cluster using Tavor's shift distance code. shiftDist = 5
2. Arrange each cluster by number nonzero samples and keep the **one** most abundant anchor per cluster.

## Running the workflow

### 1. Install `mamba` and `snakemake`

You must have conda installed with the `mamba` package manager. See [here](https://github.com/conda-forge/miniforge) for info on how to install mamba. After that, create a snakemake environment using the following command:

```{python}
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### 2. Get R configured to have the right packages

In lieu of having a virutal environment for R, one must have the required modules installed in their environment for `R/4.3.2` on Sherlock. First you must load the required modules on the command line:

```{bash}
# ON A DEV / INTERACTIVE NODE
ml R/4.3.2
ml system
ml fribidi
ml freetype
ml harfbuzz
ml fontconfig
```

Then inside an `R` session, install the modules listed below:

```{R}
# now inside R
install.packages("data.table") # you'll need to select a mirror
install.packages("tidyverse")
install.packages("glmnet")
install.packages("optparse")
install.packages("furrr")
install.packages("caret")
install.packages("resample")
install.packages("mltools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

q()
```

Now `R` is configured to have all of the required packages moving forward.

TODO: Create an environment or Docker images with all of these packages and just point at that to avoid this initial step.

### 3. Ensure the correct datasets are in the `dataset_table.csv` file

This file contains all of the datasets on which we want to run Snakemake. Replace the current datasets with ones you would like to run. To run on a new dataset, add a new row with a **short name**, a path to the **SPLASH run folder**, a path to the **metadata file**, a path to a **lookup table for artifact filtering** (you don't need to change this from the other rows), and a **translation table** (this is an integer corresponding to the correct genetic code to use for translation).

#### Inside the SPLASH run folder, the script is looking for the following files

1. `result.after_correction.scores.tsv`
2. `sample_name_to_id.mapping.txt`
3. `result_satc/` containing all of the satc files for the run

#### The metadata file needs to be formatted as follows

1. Column 1 contains the same sample ids as were used to run `SPLASH` and is named `sample_name` (the script will attempt to assign the first column as `sample_name` otherwise)
2. The other columns contain a simple column name describing the metadata and the metadata.

### 4. Ensure that the wildcards you want to use are specified at the top of the Snakefile

The current wildcards are defined as follows:

```{python}
## Define the wildcards on which the pipeline will be run
# TODO: Dynamically generate {dataset} based on the metadata table
# TODO: Define the other wildcards based on the config file
DATASETS = list(metadata_table.index)
SELECT_TYPES = ["filter1", "filter2", "filter3"]
CLUSTER_TYPES = ["shiftDist-keepTopES", "shiftDist-keepMostAbundant"]
NUM_CLUSTERS = [4000]
KMER_WIDTH = [54]
KMER_STEP = [54]
MODELS = ["esm"]
NORMALIZE = ["normalized", "unnormalized"]
```

You can change these by specifying new paramaters.In order to change the `SELECT_TYPES` or the `CLUSTER_TYPES`, you must also provide new scripts performing the desired steps.

### 4. Change out any necessary file paths in the `config.yaml` file

In the config file, you should change out the temp directory and point it at your personal `$SCRATCH` folder as the pipeline will create all necessary intermediate files in this folder.

### 5. Run `snakemake`

You can run the workflow by typing the following in an interactive node:

```{bash}
snakemake --profile slurm_profile/
```

This will then submit and watch all the jobs along the way until the pipeline is complete. You can also use the provied `run_snakemake.sbatch` command.

## Input Files 

Input files and paths are details in `dataset_table.csv` and the columns are descriptive. Metadata paths as well as SPLASH results paths are stored there.

## Output Files

Output files are available in the `results/` folder with the dataset name as a subfolder.
The files that the script generates are as follows: 

- Annotated nonzero coefficients files: results/dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_glmnet_results_top4000_k54_s54_nonzero_coefficients_annotated.tsv
- Pdfs: results/dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_glmnet_results_top4000_k54_s54_confusion_matrices.pdf
- Random Forests PDFS: results/dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_randomForests_results_top4000_k54_s54_confusion_matrics.pdf
- Random forests important features: results/dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_randomForests_results_top4000_k54_s54_important_features.tsv

## Additional details

### Adding new filters and new clustering scripts

You can add new filters or clustering scripts by adding a line to the config.yaml file:

```{yaml}
scripts:
  anchor_select_script: 
    filter4: src/anchor_selection/NEW_SCRIPT_HERE.R
```

The same can be said for the clustering script though in that case you need to add a subentry both under `cluster_script` and under `reorder_script` with the same name. Use the currently existing code as an example for what input and outputs should look like and how to process them on the command line. See the following rules: `choose_anchors`, `cluster_anchors`, `reorder_clusters`

### Sometimes jobs fail

Use the slurm log for a given job to find out why. All of the logs that slurm produces will be organized in the `logs/` folder with a subfolder for the date and then a subfolder for each rule.

### TODO LIST

- Add a rule to rclone everything automatically to Google Drive
- Add a rule to grab everything for a given dataset and produce a summary plot of the accuracies.
- Add in the attribution code
- Add in GPN/Hyena
