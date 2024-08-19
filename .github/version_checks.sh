#!/usr/bin/env bash

config_version=$(sed -n "s/.*version.*= '\(.*\)'.*/\1/p" nextflow.config)

# Check CHANGELOG version

head -n 10 CHANGELOG.md | grep "## $config_version - " >/dev/null \
    || (echo 'Failed to match CHANGELOG version'; exit 1)
