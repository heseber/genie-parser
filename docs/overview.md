# Package overview

## Main class `Genie`

The class `Genie` is the main class of the `genie` package. A `Genie` object contains all data of a Genie release and the results of a re-annotation of mutations based on *Illumina Connected Annotations* (previously known as *Nirvana*).

!!! example "How to create a Genie object"

    ```python
    # Loading all data takes about 2 minutes.
    import genie.genie as gd
    genie_dir = "/my/path/to/genie-<version>"
    g = gd.Genie(genie_dir)
    print(g.summary())
    ```


## Other data classes 

Other classes defined by the `genie` package store the different data types
provided by the Genie consortium. Objects of the main class `Genie` have
instances of these classes as members (attributes).

* `PatientInfo`: clinical data of patients.
* `SampleInfo`: sample annotation.
* `Mutations`: mutation data (SNVs and InDels) and mutation annotations.
* `CNA`: copy number alteration data.
* `Panel`: information about a single panel (annotation and content).
* `PanelSet`: all `Panel`s included in Genie.
* `TestedPositions`: all positions tested by any `Panel` from the `PanelSet`.
* `Meta`: meta information on data files.

## Configuration

The `Configuration` class stores the names of original data files provided by
Genie and of auxiliary files derived from these original data files. The names
of these files are defined in a JSON configuration file that is part of the
`genie` module. Users can provide a custom configuration file to override some
of the predefined file names if a new Genie version changes names.

The `Configuration` class also provides some functions for loading data files
and auxiliary files as Pandas DataFrames. The data is read from cache files in
Parquet format. If a cache file does not exist yet, it is automatically created
from the original TSV file on the fly.

## Details of data analysis

Please see the [Analysis Details](analysis_details.md) for a description how
individual mutations are aggregated to gene level and how co-mutations are
detected.





