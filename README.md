# *FLASH*: FunctionaL Assigning Sequence Homing

## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction and attribution on SPLASH results. Once SPLASH has been run on data and Snakemake has been installed, it should take *5-10 minutes* to configure the metadata table and check on input files, to prepare FLASH to run on a new dataset.

### IMPORTANT NOTE

This pipeline has been modified to run locally but it will require a large amount of resources (unless using the `--profile` command). 

#### To do

The only thing that needs to be fixed are some of the paths to SLASH repos and some of the paths to the hyena embedding objects. Need to figure out how to wrap this into a more packageable format.

Also the large artifact lookup table is too big to add to github. May need to provide code to build it locally for the user.

Also the SPLASH output directory is currently hard-coded to the location on sherlock. This will need to be fixed.

## Running the workflow

### 1. Install `mamba`, `snakemake`, and `SPLASH`

You must have conda installed with the `mamba` package manager. See [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge) for info on how to install mamba. After that, create a snakemake environment using the following command:

```{bash}
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

To use the executor profile for cluster job submission (with sbatch) you also must run the following while the new conda environment is active (to activate: `conda activate snakemake`):

```{bash}
mamba install snakemake-executor-plugin-cluster-generic
```

To install `SPLASH`, download the relase here: [github.com/refresh-bio/SPLASH/releases/tag/v2.11.6](https://github.com/refresh-bio/SPLASH/releases/tag/v2.11.6). Unzip it into a folder in the project root and change the `splash_bin` paramater in `config.yaml` to match. You can also run the script `get_splash.sh` inside the `splash` folder if you are using x64 linux.

### 2. Ensure `dataset_table.csv` is filled out correctly

This file contains all of the datasets on which we want to run the FLASH Snakemake pipeline. Replace the current datasets with ones you would like to run. To run on a new dataset, add a new row with a **short name**, a path to the **SPLASH run folder** (where you have run SPLASH with the option to dump SATC files), a path to the **metadata file**, a path to a **lookup table for artifact filtering** (you don't need to change this from the other rows), and a **translation table** (this is an integer corresponding to the correct genetic code to use for translation).

#### Inside the SPLASH run folder, the script is looking for the following files

1. `result.after_correction.scores.tsv`
2. `sample_name_to_id.mapping.txt`
3. `result_satc/` containing all of the satc files for the run.

There is an example script for running SPLASH in the `resources/utility_scripts` folder (run_splash.bash).

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
TRAIN_PROPORTION = 0.5 # this is the proportion of the data to use for training, the rest will be used for testing.
```

You can change these by specifying new values for the paramaters.In order to change the `SELECT_TYPES` or the `CLUSTER_TYPES`, you must also provide new scripts on which to perform the desired steps.

### 4. Change out any necessary file paths in the `config.yaml` file

In the config file, you should change the path to the temp directory and point it at a temporary folder where the pipeline will create all necessary intermediate files. Depending on the size of the input data, these intermediate files can be quite large.

### 5. Run `snakemake`

You can run the workflow by typing the following in an interactive node:

```{bash}
export NUM_CORES=1
snakemake --sdm conda -j $NUM_CORES all_ohe
```

This will execute `snakemake` to run jobs using up to `$NUM_CORES` maximum in parralel and will produce the files specifice by the rule `all_ohe`..

The flag `--sdm conda` is important because it will tell snakemake to build and use the required conda environments specified in the `envs/` directory. *Note: This can be quite slow whenever the environment is installed for the first time or is changed.*

If you want to run FLASH using the embedding mode or genomes mode, you will need a GPU in the environment where you run `snakemake` or you will need to use a profile.

The included example submission script (`resources/utility_scripts/run_snakemake.sbatch`) and profile (`slurm_profile/config.v8+.yaml`) submit the pipeline and request the required resources including GPU resources. For running on a cluster you can use a profile using --profile slurm_profile/config.v8+.yaml. **NOTE: This is untested outside of some specific cases. You will probably need to change the general resources specified in the profile as well as the rule-specific resources.**

## Input Files

Input files and paths are detailed in `dataset_table.csv` and the columns are descriptive. Metadata paths as well as SPLASH results paths are stored there.

- `dataset_short_name`: A short name for the dataset (Don't include underscores).
- `SPLASH_results`: Path to the SPLASH run folder containing results. Should contain `result.after_correction.scores.tsv`, `sample_name_to_id.mapping.txt`, and the folder `result_satc`. (An example script for running SPLASH to produce these outputs is in the `resources/helper_scripts` folder).
- `metadata_file`: Path to the metadata file associated with the dataset.
- `lookup_table`: Path to the lookup table for annotation.
- `translation_table`: Integer corresponding to the correct genetic code for translation.
- `genomes_folder`: Path to the folder containing the genome fasta files.
- `genome_list`: List of "testing data" genomes for extracting. Names must map to the filenames in `genomes folder` (e.g. name:`SAMPLE1` -> file: `SAMPLE1.fasta`)
- `genome_metadata_file`: Path to the genome metadata file. Takes the same form as the metadata file above. Column names and labels must match the original metadata file.
- `taxid`: Taxonomic ID associated with the dataset. Used to restrict the BLAST-based annotation pipeline to a specific taxon.

## Target rules and output files

Output files are available in the `results/` folder with the dataset name as a subfolder and subsequent subfolders defining different branching points along the way. The flag `GENERATE_PLOTS` can be changed to `True` to force the script to perform the annotation and plotting steps.

- The target rule **`rule_all_embeddings`** genarates all embedding-based predictions and summary fules using the nucleotide language model.
- The target rule **`rule_all_genomes`** generates models based on embeddings data and then predicts on these data using test data supplied in the genome directory as specified in the `dataset_table.csv` file.
- The target rule **`rule_all_ohe`** formats the data using One Hot Encoding and predicts on the data this way. This does not require setting up the gpu-based embedding pipeline.
- The target rule **`rule_all_umap`** generates the unsupervised umap clusters using the top variance embedding per cluster in the formatted data.
