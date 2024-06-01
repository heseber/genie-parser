## Preparations for using the genie package

Before using the package, the GENIE data needs to be dowloaded and some
additional derived data files need to be created. For example, the mutations in
GENIE are not normalized (left-shifted). Therefore, the same mutation can be
specified differently for different panels.

### Step 1: Setting up the infrastructure

[Set up the environment](python_env.md), which mostly means creating a mamba
environment and installing some Python packages and bioinformatics tools.

### Step 2: Preparing GENIE data

This is required whenever a new GENIE release is published. Required steps:

1. [Download the data](download.md)
2. [Annotate mutations](annotate.md)


