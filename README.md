# *FLASH*

## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction and attribution on SPLASH results.

## Running the workflow

### 1. Install `mamba` and `snakemake`

You must have conda installed with the `mamba` package manager. See [here](https://github.com/conda-forge/miniforge) for info on how to install mamba. After that, create a snakemake environment using the following command:

```{python}
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### 3. Ensure the correct datasets are in the `dataset_table.csv` file

This file contains all of the datasets on which we want to run the FLASH Snakemake pipeline. Replace the current datasets with ones you would like to run. To run on a new dataset, add a new row with a **short name**, a path to the **SPLASH run folder** (where you have run SPLASH with the option to dump SATC files), a path to the **metadata file**, a path to a **lookup table for artifact filtering** (you don't need to change this from the other rows), and a **translation table** (this is an integer corresponding to the correct genetic code to use for translation).

#### Inside the SPLASH run folder, the script is looking for the following files

1. `result.after_correction.scores.tsv`
2. `sample_name_to_id.mapping.txt`
3. `result_satc/` containing all of the satc files for the run

#### The metadata file needs to be formatted as follows

1. Column 1 contains the same sample ids as were used to run `SPLASH` and is named `sample_name` (the script will attempt to assign the first column as `sample_name` otherwise)
2. The other columns contain a simple column name describing the metadata as well as the observations of the metadata.

### 4. Ensure that the wildcards you want to use are specified at the top of the Snakefile

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
snakemake -j $NUM_CORES all
```

This will then submit and watch all the jobs along the way until the pipeline is complete. For running on a cluster you can use a profile using --profile slurm_profile/config.v8+.yaml

## Input Files

Input files and paths are detailed in `dataset_table.csv` and the columns are descriptive. Metadata paths as well as SPLASH results paths are stored there.

## Output Files

Output files are available in the `results/` folder with the dataset name as a subfolder and subsequent subfolders defining different branching points along the way.