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

### 3. Ensure that the parameters you want to use are specified in the `config.yaml` file. Snakemake will use these to fill out the wildcards in the `Snakefile`

The current wildcards are defined as follows:

```yaml
options:
  seqs_from_raw_data: true # set to true if starting from raw FASTQ files, false if starting from pre-processed SPLASH results
  generate_plots: false # set to true to generate plots after model training
  num_clusters: 20000 # number of clusters to use for clustering anchors
  filters: ["filter1"] # define which anchor selection filters to use
  cluster_types: ["shiftDist-levFilter", "masked-nucleotide-clustered"] # define which clustering methods to use
  models: ["hyena"] # currently only 'hyena' is supported
  normalize_embeddings: ['normalized', 'unnormalized'] # define whether to use normalized or unnormalized embeddings, this will do both
  train_proportion: 0.8 # proportion of samples to use for training GLMnet models
  anchor_length: 27 # this is the length of the anchor in nucleotides, change as needed for different SPLASH runs
  target_length: 27 # this is the length of the target in nucleotides, change as needed for different SPLASH runs
  target_rank: 1 # this is the rank of the target to use when assembling the anchor-target pairs, change as needed

# you should not need to change anything below this line but you can alter these if needed
extended_options:
  num_anchors_to_select: 3000000 # number of anchors to select from each dataset for clustering
  effect_size_cutoff: 0.6 # this is the cutoff for selecting an anchor from the SPLASH results
  distance_threshold: 5 # distance threshold for filtering anchors after clustering
  num_embedding_features_to_keep:
    glmnet: 100 # number of top variance embedding features to keep for GLMnet
    umap: 1 # number of top variance embedding features to keep for UMAP
  num_PCs_umap: 10 # number of principal components to use for UMAP plotting
  min_samples_adelie: 28 # minimum number of samples in the smallest class for running adelie (requires at least N samples per class)

```

You can change these by specifying new values for the parameters.
In order to use your own filter clustering type, you must also provide new scripts (in the `scripts:` section) containing modified code for anchor selection and clustering. See the existing scripts for examples of how to do this.

#### Generating annotations and plots

To generate annotations and plots you can set the flag `generate_plots:` in the config file. This will run BLAST on the output sequences for each analysis and generate a series of visualizations. This step can be very slow as the default behavior is to run blast in `remote` mode. By pre-downloading the blast databases for `core_nt` and `refseq_protein` you can speed this step up.

To point at the folder containing the local blast databases, modify the config file, `config.yaml`, and change the field `blast_db_path` to point to the folder containing the local blast databases.

These can be downloaded with `blast+`; they are installed using `update_blastdb.pl --decompress core_nt refseq_protein`. For further instructions on downloading local copies of these databases, see <https://www.ncbi.nlm.nih.gov/books/NBK569850/>.

To enable taxonomic filtering of the blast features, using the `taxid` field in the file `dataset_table.csv`, you must A) be using local blast databases as remote blast does not support taxonomic filtering and B) You also must additionally download and decompress the file `taxdb.tar.gz` from <https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz> in the local blast database folder. Follow the instructions in the *Taxonomic filtering for BLAST databases* section of <https://www.ncbi.nlm.nih.gov/books/NBK569839/> for further details.

### 4. Update all paths in the config file, `config.yaml`

The `config.yaml` file is also used by the Snakemake workflow to specify the correct locations of all paths and some parameters.

It is necessary to change the fields below to match your setup and computing environment:

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

If you want to run FLASH using the embedding mode or genome prediction mode, you will need a GPU in the environment where you run `snakemake` or you will need to use a profile that can allocate GPUs and other resources (like the one provided for slurm schedulers in the directory `slurm_profile/`).

---

**IMPORTANT:** *Before you run in embeddings mode you will need to download the singularity image containing the model state and necessary code into the folder `containers/`. See the `containers/README.md` file in that folder for further information.*

---

The command for running the code in an interactive environment (where a GPU is present) is similar to above:

```{bash}
export NUM_CORES=1
snakemake --sdm conda -j $NUM_CORES all_embeddings
```

The code can also be run using an automatic scheduler. The included example submission script (`resources/utility_scripts/run_snakemake.sbatch`) and profile (`slurm_profile/config.v8+.yaml`) can submit the pipeline and request the required resources including GPU resources. For running on a cluster using `slurm` you can use this config using the --profile slurm_profile/config.v8+.yaml.

**NOTE: You must modify the profile to match the correct partitions, resources, and constraints necessary for your cluster. The provided config is an example for our local clusters resources. This can be adapted to different schedulers by modifying the profile but is currently only set to work for `slurm`.**

## Input Files

Input files and paths are detailed in `dataset_table.csv` and the columns are descriptive. Metadata paths as well as SPLASH results paths are stored there.

- `dataset_short_name`: A short name for the dataset (Don't include underscores).

- `SPLASH_results`: Path to the SPLASH run folder containing results. Should contain `result.after_correction.scores.tsv`, `sample_name_to_id.mapping.txt`, and the folder `result_satc`. (An example script for running SPLASH to produce these outputs is in the `resources/helper_scripts` folder).
  - If using the flag `seqs_from_raw_data: true` in the config file, this step will additionaly allow the script to find the `sample_sheet.tsv` file. This file is the default one used to run SPLASH and contains two columns (`sample_name` and `file_path`) that map sample names to raw FASTQ files. This will allow the script to extract sequences directly from the raw data instead of from the SPLASH outputs which circumvents the limitations of SPLASH in terms of outputting sequences that appear less frequently.

- `metadata_file`: Path to the metadata file associated with the dataset.
  - The metadata file must contain a column named `sample_name` that matches the sample names used in the SPLASH run.

- `lookup_table`: Path to the lookup table for annotation.

- `translation_table`: Integer corresponding to the correct genetic code for translation.

- `taxid`: Taxonomic ID associated with the dataset. Used to restrict the BLAST-based annotation pipeline to a specific taxon.

### Additional inputs for genome predictions

For running genome-based predictions, additional fields must be filled out in the `dataset_table.csv` file. If you are not running genome-based predictions, you can ignore set these fields to NA.

- `genomes_folder`: Path to the folder containing the genome fasta files.

- `genome_list`: List of "testing data" genomes for extracting. Names must map to the filenames in `genomes folder` (e.g. name:`SAMPLE1` -> file: `SAMPLE1.fasta`)

- `genome_metadata_file`: Path to the genome metadata file. Takes the same form as the metadata file above. Column names and labels must match the original metadata file.


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
4. `sample_sheet.txt` (if using the `seqs_from_raw_data: true` option in the config file)

These files are the necessary inputs for the FLASH pipeline. You can then proceed to fill out the `dataset_table.csv` file with the path to the SPLASH output folder containing these files and to the metadata file associated with your samples. The metadata file should have a column named `sample_name` that matches the sample names used in the SPLASH run.

### Running FLASH

You should not need to modify any of the parameters in the `config.yaml` unless you change the anchor or target lengths when running SPLASH as mentioned above. If you do change these lengths, ensure that the `anchor_length` and `target_length` variables in the `config.yaml` file are updated accordingly. After confirming that all paths and parameters are correctly set, you can run the FLASH pipeline using Snakemake as described in the previous sections.

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
