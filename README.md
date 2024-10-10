# *metaSPLASH* Pipeline

Daniel Cotter -- 10/07/2024

## Description

The following pipeline uses `Snakemake` to run metadata-driven prediction
and attribution on SPLASH results.

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

### Install `mamba` and `snakemake`

You must have conda installed with the `mamba` package manager. See [here](https://github.com/conda-forge/miniforge)for info on how to install mamba. After that, create a snakemake environment using the following command:

```{python}
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```


