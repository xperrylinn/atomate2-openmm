# Underconstruction...

# atomate2-openmm
An add-on package to atomate2 for OpenMM

# Installation

This package will eventually be distributed using conda-forge. For now follow see the development 
environmentvsection.

# Development Environment Setup

## Python Environment Setup

A Python3.8 environment setup using conda according to the steps below:

1. `conda env create -f environment.yaml`
2. `conda activate atomate2-openmm`
3. `pip install -e .` from repository root

## Environment variables:

- ATLAS_USERNAME: <username for MongoDB Atlas>
- ATLAS_PASSWORD: <database password. See database access page.>

## Configuration Files

For AWS S3 IO, make sure you have an `~/.aws/credentials` INI file present with the following secion and key values pairs:

```
[atomate2-openmm-dev]
aws_access_key_id=<your AWS access key>
aws_secret_access_key=<your AWS secret key>

```

## Running Test

From root of repo, run the following command:

`python -m pytest tests`

# Contributing
- Open A PR or file a GitHub issue
- Numpy style documentation
- Formatting with black
