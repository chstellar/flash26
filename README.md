# *FLASH*

## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction and attribution on SPLASH results.

### IMPORTANT NOTE

This pipeline has been modified to run entirely removed from Sherlock (unless using the `--profile` command). The only thing that needs to be fixed are some of the paths to SLASH repos and some of the paths to the hyena embedding objects. Need to figure out how to wrap this into a more packageable format.

Also the large artifact lookup table is too big to add to github. May need to provide code to build it locally for the user.

Also the SPLASH output directory is currently hard-coded to the location on sherlock. This will need to be fixed.

## Running the workflow

### 1. Install `mamba` and `snakemake`

You must have conda installed with the `mamba` package manager. See [here](https://github.com/conda-forge/miniforge) for info on how to install mamba. After that, create a snakemake environment using the following command:

```{bash}
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

To use the executor profile for slurm you also must run the following while the snakemake environment is active:

```{bash}
mamba install snakemake-executor-plugin-cluster-generic
```

### 2. Ensure the correct datasets are in the `dataset_table.csv` file

This file contains all of the datasets on which we want to run the FLASH Snakemake pipeline. Replace the current datasets with ones you would like to run. To run on a new dataset, add a new row with a **short name**, a path to the **SPLASH run folder** (where you have run SPLASH with the option to dump SATC files), a path to the **metadata file**, a path to a **lookup table for artifact filtering** (you don't need to change this from the other rows), and a **translation table** (this is an integer corresponding to the correct genetic code to use for translation).

#### Inside the SPLASH run folder, the script is looking for the following files

1. `result.after_correction.scores.tsv`
2. `sample_name_to_id.mapping.txt`
3. `result_satc/` containing all of the satc files for the run

#### The metadata file needs to be formatted as follows

1. Column 1 contains the same sample ids as were used to run `SPLASH` and is named `sample_name` (the script will attempt to assign the first column as `sample_name` otherwise)
2. The other columns contain a simple column name describing the metadata as well as the observations of the metadata.

### 3. Ensure that the wildcards you want to use are specified at the top of the Snakefile

The current wildcards are defined as follows:

```{python}
## Define the wildcards on which the pipeline will be run
DATASETS = list(dataset_table.index)
SELECT_TYPES = ["filter1"]
CLUSTER_TYPES = ["shiftDist-levFilter", "masked-aa-clustered"]
NUM_CLUSTERS = [5000]
ANCHOR_LENGTH = 27 # this is the length of the anchor in nucleotides
TARGET_LENGTH = 27 # this is the length of the target in nucleotides
KMER_WIDTH = [ANCHOR_LENGTH + TARGET_LENGTH] # this is Anchor Length + Target Length
KMER_STEP = [ANCHOR_LENGTH + TARGET_LENGTH] # this can be used to let the steps be variable, but for now we will use a single value
MODELS = ["hyena"]
NORMALIZE=["normalized"]
TRAIN_PROPORTION = 0.5 # this is the proportion of the data to use for training, the rest will be used for testing
```

You can change these by specifying new paramaters.In order to change the `SELECT_TYPES` or the `CLUSTER_TYPES`, you must also provide new scripts on which to perform the desired steps.

### 4. Change out any necessary file paths in the `config.yaml` file

In the config file, you should change out the temp directory and point it at a temporary folder where the pipeline will create all necessary intermediate files in this folder.

### 5. Run `snakemake`

You can run the workflow by typing the following in an interactive node:

```{bash}
export NUM_CORES=1
snakemake --sdm conda -j $NUM_CORES all_ohe
```
`--sdm conda` is important because it will tell snakemake to build the required conda environments

If you want to run FLASH using the embedding mode or genomes mode, you will need a GPU in the environment in which you submit snakemake or you will need to use a profile and submission script that submits the embedding job to a node with a GPU resource.

For running on a cluster you can use a profile using --profile slurm_profile/config.v8+.yaml. **NOTE THIS IS UNTESTED IN THE CURRENT ITERATION**

## Input Files

Input files and paths are detailed in `dataset_table.csv` and the columns are descriptive. Metadata paths as well as SPLASH results paths are stored there.

## Output Files

Output files are available in the `results/` folder with the dataset name as a subfolder and subsequent subfolders defining different branching points along the way.