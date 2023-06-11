# Underconstruction...

# atomate2-openmm
An add-on package to atomate2 for OpenMM

# Manifest

## Scripts
- additional_stores_test_with_memory_store.py: small example of using additional stores in JobFlow API
- atlas_clear_collection.py: deletes objects in collection
- atlas_connection.py: test connection to MongoDB Atlas
- boto3_s3.py - test connection to S3 bucket
- s3_store.py: test Maggma S3 store
- generate_erd_diagram.py: script for generating ERD diagrams of PyDantic data models
- nvt_flow_simple_using_external_storage.py: simple NVT simulation Flow using MongoDB Atlas and S3 storage
- nvt_flow_simple_using_memory_store.py: simple NVT simulation Flow using local memory storage
- production_flow_using_external_storage.py: production Flow using MongoDB Atlas and S3 storage
- production_flow_using_memory_storage.py: production Flow using local memory storage

# Development Environment Setup

## Python Environment Setup

A Python3.8 environment setup using conda according to the steps below:

1. `conda env create -f environment.yaml`
2. `conda activate atomate2-openmm`
3. `pip install -e .` from repository root

## Environment variables:

- ATLAS_USERNAME: \<username for MongoDB Atlas\>
- ATLAS_PASSWORD: \<database password. See database access page\>

## Configuration Files

For AWS S3 IO, make sure you have an `~/.aws/credentials` INI file present with the following secion and key values pairs:

```
[atomate2-openmm-dev]
aws_access_key_id=<your AWS access key>
aws_secret_access_key=<your AWS secret key>

```



# Contributing

- Open A PR or file a GitHub issue
- Numpy style documentation
- Formatting with black

## Running Test

From root of repo, run the following command:

`pytest .`

All PRs must maintain 100% testing coverage. To confirm full test coverage run the following command from the root:

`pytest --cov=tests --cov-report term-missing`

