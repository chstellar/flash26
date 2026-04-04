# Get the directory where the script is located (should be /containers)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Configuration ---
# Ensure the URL is lowercase and uses the oras:// protocol
IMAGE_URL="oras://ghcr.io/djcotter/flash-smk/hyena_embedder:latest"
IMAGE_FILE="$SCRIPT_DIR/hyena_embedder.sif"

# 1. Ensure we are in the correct directory (optional but safe)
cd "$SCRIPT_DIR" || exit

# 2. Check if the image already exists
if [ -f "hyena_embedder.sif" ]; then
    echo "Image already exists at $IMAGE_FILE. Skipping download."
else
    echo "Downloading Singularity image from GHCR to $SCRIPT_DIR..."
    
    # Pull the image. The first argument is the output name, second is the source.
    singularity pull "hyena_embedder.sif" "$IMAGE_URL"
    
    # 3. Final verification
    if [ $? -eq 0 ]; then
        echo "Success: Image ready at $IMAGE_FILE"
    else
        echo "Error: Failed to pull image. Ensure image is Public or you are logged in."
        exit 1
    fi
fi