# *metaSPLASH* Pipeline

Daniel Cotter -- 10/07/2024

## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction and attribution on SPLASH results.

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

Similarly, if there are any issues with python environments

### 5. Run `snakemake`

You can run the workflow by typing the following in an interactive node:

```{bash}
snakemake --profile slurm_profile/
```

This will then submit and watch all the jobs along the way until the pipeline is complete. You can also use the provied `run_snakemake.sbatch` command.

## Input Files 

Input files and paths are detailed in `dataset_table.csv` and the columns are descriptive. Metadata paths as well as SPLASH results paths are stored there.

## Output Files

Output files are available in the `results/` folder with the dataset name as a subfolder.
The files that the script generates are as follows: 

- Annotated nonzero coefficients files: `results/{dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_adelie_results_top4000_k54_s54_nonzero_coefficients_annotated.tsv`
- Pdfs: `results/{dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_adelie_results_top4000_k54_s54_confusion_matrices.pdf`
- Random Forests PDFS: `results/{dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_randomForests_results_top4000_k54_s54_confusion_matrics.pdf`
- Random forests important features: `results/{dataset}/{select_type}/{cluster_type}/{model}/{normalize}/{dataset}_{model}_randomForests_results_top4000_k54_s54_important_features.tsv`

## Details on Pipeline steps
1. Filter anchors from the SPLASH stats file to decide which anchors to include downstream. See below for specific information on what each Filter is doing. This step takes in SPLASH stats and outputs a list of anchors to a text file that matches the filter criteria (and number). 
2. Cluster and report anchor clusters. These steps take in a set of anchors (From step 1) and output a two column file with cluster id and anchor after filtering and reordering the clusters. 
  - Cluster anchors using a specified algorithm. Some options we have explored (see cluster for more details):
    
    - shiftDist clustering
    - mmseqs for clusterings
    - AA-based clustering
Once clusters are determined, we reorder them by taking the average effect size across all anchors in a cluster and place a cluster first if it has the highest average effect size. 
Within each cluster, we reorder the anchors so that the most abundant anchor across all samples (as defined by SPLASH is first, and so on). 
We then take the top (N=20000) clusters to simplify downstream analysis since the tasks can get unwieldy
We prepare the sequences from all of the samples for downstream analysis by finding each sample’s highest count target(s) for each anchor. This step takes in a set of anchor clusters (From step 2) as well as a set of SATC files from the SPLASH run and dumps only those anchors (and their top count targets). It then formats these into anchor + target concatenated format for output to further steps. The final output is the kmer decomposition of this set of sample sequences that slides across these concatenated sequences and identifies all unique anchors+targets and puts them in a file to embed. 
The prepare sequences step uses the anchor clusters file and the SATC files to take the first anchor that each sample has (in each cluster) and it’s highest count target and then concatenates these anchor-target pairs into one long fasta entry per sample. If a given sample has NO anchors from a cluster present, then it will instead receive a 54mer composed of the first anchor in the cluster ( the “representative” ) and 27 Ns. 
The decompose kmers step slides along this sequence and decomposes it back into it’s constituent kmers. It will create a mapping file to identify where each 54mer appears in a given sequence and it will create a fasta file with ALL of the unique 54mers that will need to be embedded
For ESM only, this step translates the 54mers into an amino acid alphabet and returns a translated fasta prior to embedding.
Each unique 54mer receives a set of embeddings by processing it through a pre-trained language model. For example, in ESM (the Evolutionary State Model []), we embed each translated 54mer and then average the embeddings to get a vector that is 1xM (where M is the number of dimensions in the model).
In this way we receive a set of embeddings for each unique 54mer that can then be reassembled to represent each sample’s embeddings for a different anchor-target pair.
To prepare these embeddings for downstream use, we merge them back with the ordering file (generated from the decompose kmers step) and pivot this matrix to get one long vector of embeddings per sample. The dimension of this matrix is num_samples X num_clusters*num_dimensions. Because num_clusters * num_dimensions is a very large number, we also want to drop some dimensions here to reduce the scale of any downstream modeling. To do this, we keep only the highest K variance columns per cluster so that each cluster has even representation in the downstream models. 
Finally, we perform prediction on this data structure using glmnet and splitting the embedding matrix from 4 into two sets, train and test. To train, we use cross-validation and then to test we pick the model lambda that resulted in the minimum mean cross-validated error. We test by predicting the metadata directly from the test embeddings and compute accuracy by comparing the predictions to the known results. As a result, for every piece of metadata we have two things: 1) an accuracy for the prediction under this model and parameters, and 2) a set of nonzero coefficients that were determined to drive this prediction. Because these coefficients are for embeddings for specific anchor clusters, we can attribute a model to specific clusters and identify which constituent sequences are therein.



## Details on Pipeline paramaters

### Different filters used for anchor selection

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

#### shiftDist-levFilter

1. Cluster using Tavor's shift distance code. shiftDist = 5
2. Pick the highest ES anchor in each clusters as the Rep and remove any other anchors that are more than lev > 5 away from it.

#### shiftDist-hamFilter

1. Cluster using Tavor's shift distance code. shiftDist = 5
2. Pick the highest ES anchor in each clusters as the Rep and remove any other anchors that are more than hamming > 5 away from it.

## Additional details

### Notes on paths and environments

Currently there are a few python virtual environment paths that are hardcoded into the `config.yaml` file. These are hosted on $OAK and should be usable by anyone on our lab's partition. Additionally I have provided `requirements.txt` files in the `envs/` folder for each of these. Note that each of these environments is built on top of the `python/3.9.0` module provided by Sherlock, but if creating them directly this should not matter.

One problem with embedding in the hyena models is that the Docker image, checkpoint and config file are all located on our partition's OAK storage. The `src/hyena_utils/embed_fasta*` scripts can be modified to accomadate this and I can provide a copy of the Docker image if neccessary. Otherwise, the pipeline can still be run with the model set as `esm` or `ohe`.

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
