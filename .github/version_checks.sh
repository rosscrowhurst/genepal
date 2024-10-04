#!/usr/bin/env bash

config_version=$(grep -A 10 'manifest {' nextflow.config | sed -n "s/.*version.*= '\(.*\)'.*/\1/p")

# Check CHANGELOG version

head -n 10 CHANGELOG.md | grep "## $config_version - " >/dev/null \
    || (echo 'Failed to match CHANGELOG version'; exit 1)

# Check multiqc_config.yml version
config_version_plain=$(echo $config_version | sed 's/.*\+dev$/dev/; t; s/\+.*//')
echo "Looking for '$config_version_plain' in multiqc_config.yml"

head -n 5 assets/multiqc_config.yml | grep "https://github.com/plant-food-research-open/genepal/blob/$config_version_plain/docs/usage.md" >/dev/null \
    || (echo 'Failed to match multiqc_config.yml version'; exit 1)
