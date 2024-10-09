# *metaSPLASH* Pipeline

Daniel Cotter -- 10/07/2024

## Description 

The following pipeline uses `Snakemake` to run metadata-driven prediction
and attribution on SPLASH results.

### Different filters used for anchor selection3

These filters are used to select the set of all anchors from the SPLASH results that will be sent for clustering and reordering. See `config.yaml` for the script used for the filtering.

 | **Filter** | **Description** |
 |:-----------|:---------------:|
 | filter1    | NA              |
 | filter2    | NA              |
 | filter3    | NA              |

### Different clustering methods

These are different methods for clustering the anchors from step 1 and then informed reordering/pruning of the clusters. See `config.yaml` for the scripts used for these steps.

 | **Cluster Method** | **Description** |
 |:-------------------|:---------------:|
 | shiftDist-keep1    | NA              |
 | filter2            | NA              |
 | filter3            | NA              |
