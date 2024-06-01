# Downloading data

## Downloading GENIE data

The GENIE data need to be downloaded from
[Synapse](https://www.synapse.org/#!Synapse:syn7222066/wiki/410924).

It is recommended to use the command line tool `synapse` for downloading the
data. If you haven't done yet, install the `synapseclient` Python package.

!!! example "install_synapseclient.sh"
    ```sh
    pip install synapseclient
    ```

And then download the GENIE data into your current directory, here for version
15.0 of GENIE. Make sure you have [registered for Synapse GENIE
downloads](https://help.synapse.org/docs/Managing-Your-Account.2055405596.html)
and [created a download token](https://www.synapse.org/#!PersonalAccessTokens:)
for Synapse, you will need it for authentication.

!!! example "download_genie.sh"
    ```sh
    synapse get -r syn53210170 # the number is different for each GENIE release
    ```

!!! warning "WARNING: Fixing the PROV-TRISEQ-V2 panel required."

    There is an open issue with the PROV-TRISEQ-V2 panel, please see [the
    discussion on GENIE Synapse](
    https://www.synapse.org/#!Synapse:syn7222066/discussion/threadId=10045 ). To
    fix this issue, the script [`fix_gene_panels.py`](
    https://github.com/Bayer-Group/genie.parser/blob/main/fix_gene_panels.py )
    needs to be run once before using the `genie` package. This script will add
    a missing file for this panel in the `gene_panels` directory of the GENIE
    data.

    You have to edit the line `genie_dir` before calling `fix_gene_panels.py`.

## Downloading MANE transcript definitions

Download the latest release of *MANE_\human* from NCBI to the same directory
where the genie data is stored. This is required for reannotating variants with
MANE transcripts and for filtering for MANE transcripts.

!!! example "download_mane.sh"

    ```sh
    FILE="MANE.GRCh38.v1.3.summary.txt.gz"
    URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/$FILE"
    curl $URL > $FILE
    ```

!!! warning "Don't forget to update `default_config.json`"

    The file name of the MANE definition file is configured in
    `default_config.json` in the `genie` GitHub repository. If you download a
    new release of MANE, the filename in the config file needs to be updated.

