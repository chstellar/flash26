#!/bin/bash

SPLASH_REPO="https://github.com/refresh-bio/SPLASH/releases/download/v2.11.6/splash-2.11.6.linux.x64.tar.gz"

if [ ! -d "splash-2.11.6" ]; then
    mkdir splash-2.11.6
fi

if [ ! -f "splash-2.11.6/splash-2.11.6.linux.x64.tar.gz" ]; then
    echo "Downloading SPLASH..."
    wget $SPLASH_REPO -O splash-2.11.6/splash-2.11.6.linux.x64.tar.gz
else
    echo "SPLASH already downloaded."
fi

if [ ! -d "splash-2.11.6/splash-2.11.6" ]; then
    echo "Extracting SPLASH..."
    tar -xzf splash-2.11.6/splash-2.11.6.linux.x64.tar.gz -C splash-2.11.6
else
    echo "SPLASH already extracted."
fi