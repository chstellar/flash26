# *FLASH*: FunctionaL Assigning Sequence Homing

## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction and attribution on SPLASH results. Once SPLASH has been run on data and Snakemake has been installed, it should take *5-10 minutes* to configure the metadata table and check on input files, to prepare FLASH to run on a new dataset.

### IMPORTANT NOTE

This pipeline has been modified to run locally but it requires extensive resources. It can be submitted to a scheduler by providing a `--profile` that contains a job submission command. See `slurm_profile/` for an example.

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

#### Installing SPLASH

To install `SPLASH`, download the relase here: [github.com/refresh-bio/SPLASH/releases/tag/v2.11.6](https://github.com/refresh-bio/SPLASH/releases/tag/v2.11.6). Unzip it into a folder in the project root and change the `splash_bin` parameter in `config.yaml` to match. You can also run the script `get_splash.sh` if you are using x64 linux.

### 2. Ensure `dataset_table.csv` is filled out correctly

This file contains all of the datasets on which we want to run the FLASH Snakemake pipeline. Replace the current datasets with ones you would like to run. To run on a new dataset, add a new row with a **short name**, a path to the **SPLASH run folder** (where you have run SPLASH with the option to dump SATC files), a path to the **metadata file**, a path to a **lookup table for artifact filtering** (you don't need to change this from the other rows), and a **translation table** (this is an integer corresponding to the correct genetic code to use for translation).

#### Inside the SPLASH run folder, the script is looking for the following files

1. `result.after_correction.scores.tsv`
2. `sample_name_to_id.mapping.txt`
3. `result_satc/` containing all of the satc files for the run.

There is an example script for running SPLASH in the `resources/utility_scripts` folder (run_splash.sh).

#### The metadata file needs to be formatted as follows

1. Column 1 contains the same sample ids as were used to run `SPLASH` and is named `sample_name` (the script will attempt to assign the first column as `sample_name` otherwise)
2. The other columns contain a simple column name describing the metadata as well as the observations of the metadata.

### 3. Ensure that the wildcards you want to use are specified at the top of the Snakefile

The current wildcards are defined as follows:

```{python}
## Define the wildcards on which the pipeline will be run, for example:
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

You can change these by specifying new values for the parameters.
In order to use your own `SELECT_TYPES` or `CLUSTER_TYPES`, you must also provide new scripts on which to perform the desired steps.

#### Generating annotations and plots

To generate annotations and plots you can set the flag `GENERATE_PLOTS = True` in the `Snakefile`. This will run BLAST on the output sequences for each analysis and generate a series of visualizations. This step can be very slow as the default behavior is to run blast in `remote` mode. By pre-downloading the blast databases for `core_nt` and `refseq_protein` you can speed this step up.

To point at the folder containing the local blast databases, modify the config file `config.yaml` and change the field `blast_db_path` to point to the folder containing the local blast databases.

These can be downloaded with `blast+` installed using `update_blastdb.pl --decompress core_nt refseq_protein`. For further instructions on downloading local copies of these databases, see <https://www.ncbi.nlm.nih.gov/books/NBK569850/>.

To enable taxonomic filtering of the blast features, using the `taxid` field in the file `dataset_table.csv`, you must be using local blast databases as remote blast does not support taxonomic filtering. You also must additionally download and decompress the file `taxdb.tar.gz` from <https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz> in the local blast database folder. Follow the instructions in the *Taxonomic filtering for BLAST databases* section of <https://www.ncbi.nlm.nih.gov/books/NBK569839/> for further details.

### 4. Change out any necessary file paths in the `config.yaml` file

The `config.yml` file is used by the Snakemake workflow to specify the correct locations of all paths and some paramaters.
It is necessary to change some of the fields in this config file to match your setup.

1. Change the `entrez_email` entry if you intend to run the pipeline with blast plots enabled.
2. Change the `temp_dir` entry to point to a location that has a large amount of storage and fast I/O speeds. I set this to the temporary SCRATCH storage on our file systems. Depending on the size of the input data, these intermediate files can be quite large.
3. Change `splash_bin` to point to the installed location of SPLASH. See note above on how to install SPLASH.

### 5. Run `snakemake`

#### One Hot Encoding

You can run the workflow for a given target rule (for example using One Hot Encoding) by typing the following in an interactive node:

```{bash}
export NUM_CORES=1
snakemake --sdm conda -j $NUM_CORES all_ohe
```

This will execute `snakemake` to run jobs using up to `$NUM_CORES` maximum in parallel and will produce the files specified by the rule `all_ohe`.

The flag `--sdm conda` is important because it will tell snakemake to build and use the required conda environments specified in the `envs/` directory. *Note: This can be quite slow whenever the environment is installed for the first time or is changed.*

#### Embedding mode or genome predictions

**IMPORTANT** Before you run in embeddings mode you will need to download the singularity image containing the model and embedding code into the folder `containers/`. See the `README` file in that folder for information.

If you want to run FLASH using the embedding mode or genome prediction mode, you will need a GPU in the environment where you run `snakemake` or you will need to use a profile that can allocate GPUs and other resources (like the one provided for slurm schedulers in the directory `slurm_profile/`).

The command for running the code in an interactive environment (where a GPU is present) is similar to above:

```{bash}
export NUM_CORES=1
snakemake --sdm conda -j $NUM_CORES all_embeddings
```

The code can also be run using an automatic scheduler. The included example submission script (`resources/utility_scripts/run_snakemake.sbatch`) and profile (`slurm_profile/config.v8+.yaml`) can submit the pipeline and request the required resources including GPU resources. For running on a cluster you can use a profile using --profile slurm_profile/config.v8+.yaml. **NOTE: You will also need to modify the profile to match the correct partitions, resources, and constraints necessary for your cluster. The provided is an example for our local clusters resources. This can be adapted to different schedulers by modifying the profile but is currently only set to work for `slurm`.**

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

- The target rule **`rule_all_embeddings`** genarates all embedding-based predictions and summary files using the nucleotide language model.
- The target rule **`rule_all_genomes`** generates models based on embeddings data and then predicts on these data using test data supplied in the genome directory as specified in the `dataset_table.csv` file.
- The target rule **`rule_all_ohe`** formats the data using One Hot Encoding and predicts on the data this way. This does not require setting up the gpu-based embedding pipeline.
- The target rule **`rule_all_umap`** generates the unsupervised umap clusters using the top variance embedding per cluster in the formatted data.

## Guide: Running SPLASH and preparing inputs for FLASH

### SPLASH

To run the pipeline, you must first run SPLASH on your raw sequencing data to generate the necessary input files. An example script for running SPLASH is provided in the `resources/utility_scripts` folder (run_splash.sh). This script requires a file called `sample_sheet.txt` that contains two columns: `sample_name` and `file_path`, where `sample_name` is a unique identifier for each sample and `file_path` is the path to the raw sequencing data file (in FASTQ format). Typically, when using paired end data, we run SPLASH only on the forward reads. The options in the example script are set to generate the necessary SATC files for FLASH as they are not generated by default. You can modify the anchor and target lengths as needed, but these must match the values used in the FLASH pipeline. By default, we typically use 27 for both anchor and target lengths, which corresponds to ~9 amino acids each when translated.

### Required Outputs

After running SPLASH, ensure that the output folder contains the following files:

1. `result.after_correction.scores.tsv`
2. `sample_name_to_id.mapping.txt`
3. `result_satc/` containing all of the satc files for the run.

These files are the necessary inputs for the FLASH pipeline. You can then proceed to fill out the `dataset_table.csv` file with the path to the SPLASH output folder containing these files and to the metadata file associated with your samples. The metadata file should have a column named `sample_name` that matches the sample names used in the SPLASH run.

### Running FLASH

You should not need to modify any of the paramaters in the `Snakefile` unless you change the anchor or target lengths when running SPLASH as mentioned above. If you do change these lengths, ensure that the `ANCHOR_LENGTH` and `TARGET_LENGTH` variables at the top of the `Snakefile` are updated accordingly. After confirming that all paths and parameters are correctly set, you can run the FLASH pipeline using Snakemake as described in the previous sections.

Note that if you use a very short anchor length, you should also modify the `CLUSTER_TYPES` variable in the `Snakefile` to avoid using clustering methods that rely on longer anchors, such as `shiftDist-levFilter`. This should instead be set to use `noCluster`:

```{python}
# If using default anchor length of 27, use:
ANCHOR_LENGTH = 27
TARGET_LENGTH = 27
CLUSTER_TYPES = ["shiftDist-levFilter"]

# If using a shorter anchor length, e.g., 9, use:
ANCHOR_LENGTH = 9
TARGET_LENGTH = 27
CLUSTER_TYPES = ["noCluster"]
```

## Example Run of FLASH on H5N1 Sample Data

An example dataset containing SPLASH results and metadata for H5N1 samples is provided in the `resources/metadata/` folder. There are several utility scripts in the project root that will download the example data and run SPLASH. FLASH can then be run on this example data by following the instructions above and modifying the `dataset_table.csv` file to point to the example data paths and the included metadata file.

The script `generate_example_data.sh` will download the example data into a folder `example_data/` and generate a file `sample_sheet.txt` for running SPLASH. The script `run_splash_example.sh` will run SPLASH on the example data (provided SPLASH has been installed and the path to the binaries has been set correctly in the script). After running SPLASH, you can modify the `dataset_table.csv` file to point to the SPLASH output folder and the metadata file `resources/metadata/H5N1_example_metadata.csv`. You can then run FLASH using Snakemake as described above.
