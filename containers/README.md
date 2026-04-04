# Singularity Container for Hyena Embedding Methods

In order to run this workflow with the embedding methods you must download the pre-compiled singularity image containing the model and modified code base.

The singularity image can be downloaded using the provided `containers/setup.sh script`. It requires the `singularity` is available and in your path.

It must be named `hyena_embedder.sif`.

This image was generated using modified code from the repository `https://github.com/HazyResearch/hyena-dna` and the weights for a trained model on SPLASH anchor and targets from thousands of microbial samples. 
