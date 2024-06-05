import hashlib
import itertools
import json
import os
import re
from typing import Type, Union

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from Bio.SeqUtils import seq1
from importlib_resources import as_file, files
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from statsmodels.stats.proportion import proportion_confint as ci  # noqa: E402

import genie.data


class Merger:
    """Deep merging of dictionary hierarchies.

    This class can be used to merge two different dictionary trees.
    This is useful, for example, to load a default configuration from
    a JSON file and merge it with a user configuration file where the defaults
    are kept for all values not overwritten by the user configuration file.
    """

    @staticmethod
    def deep_merge(dict1: dict, dict2: dict) -> dict:
        """Deep merging of two hierarchical dictionary trees.

        The `dict1` dictionary is updated with any information from `dict2`.
        The value of a dictionary item for a particular key can be another
        dictionary, and a *flat* merge would simply replace the dictionary item
        for that key in `dict1` with the dictionary item from `dict2`.
        Such a *flat* merge would therefore lose the information for that
        dictionary item from `dict1`. A deep merge does not drop the value of a
        dictionary item in `dict1` and replace it with the dictionary item from
        `dict2` with the same key, but instead updates the dictionary item in
        `dict1` with the information from `dict2`, changing or adding only those
        keys that are part of `dict2`.

        Args:
            dict1: the dictionary to be updated.
            dict2: the dictionary with the update information.

        Returns:
            The updated dictionary (`dict1` updated with `dict2`).

        Examples:
            >>> genie.Merger.deep_merge(default_config, user_config)
        """
        for key, value in dict2.items():
            if (
                key in dict1
                and isinstance(dict1[key], dict)
                and isinstance(value, dict)
            ):
                Merger.deep_merge(dict1[key], value)
            else:
                dict1[key] = value
        return dict1


class Configuration:
    """Genie source data configuration.

    The default configuration is read from a default JSON config file
    which is part of the package. If a user config file is specified,
    it is used to update the information from the default config file.
    """

    def __init__(self, genie_dir: str, config_file: str = None):
        """Creates a new Configuration object.

        Args:
            genie_dir: the base directory name of the GENIE data files.
            config_file: optional user config file (JSON).

        Examples:
            >>> import os
            >>> import genie
            >>> home = os.getenv("HOME")
            >>> genie_dir = os.path.join(home, "genie_data")
            >>> user_config_file = os.path.join(home, ".genie_config.json")

            >>> # Without a user config file (use only defaults)
            >>> config = genie.Configuration(genie_dir)

            >>> # With a user config file, overwriting some defaults
            >>> config = genie.Configuration(genie_dir, user_config_file)
        """
        source = files(genie.data).joinpath("default_config.json")
        with as_file(source) as f:
            with open(f) as fp:
                default_config = json.load(fp)
        if config_file is not None:
            with open(config_file) as f:
                user_config = json.load(f)
                self.config = Merger.deep_merge(default_config, user_config)
        else:
            self.config = default_config

        self.config["dir"]["genie"] = genie_dir

    def _file_names(self, dir: str, type: str, key: str) -> Union[str, dict]:
        """Internal auxiliary function for getting file names

        This function returns the names of files as specified in the
        configuration file. The configuration file has key-value pairs with file
        names as values and stable identifiers as keys. The value can be a
        regular expression by prepending "re:". If it is a regular expression,
        files with matching names in `dir` are identified and returned. The
        `type` specifies whether this is for a "dir" (directory), "data" (data
        file), "meta" (meta information file), or "aux" (auxiliary file).
        Full pathnames are returned. If a regular expression was specified and
        there is more than one match, a list of file names is returned,
        otherwise a scalar string is returned.

        The advantage of specifying file names in a config file is that they
        do not need to be hard-coded into the Python package. If the GENIE
        consortium changes file names, a user configuration file can be used
        to adapt to such changes.

        Args:
            dir:  the directory where to search for files.
            type: the type to search for (these are the keys of the
                  highest level of the config JSON file, one of "dir", "data",
                  "meta", "doc").
            key:  the identifier for the file name.

        Returns:
            a file name or a list of file names or None.
        """
        file_names = self.config[type][key]
        if file_names.startswith("re:"):
            pattern = re.compile(file_names[3:])
            file_names = [x for x in os.listdir(dir) if re.fullmatch(pattern, x)]
        else:
            file_names = [file_names]
        file_names = [os.path.join(dir, x) for x in file_names]
        if len(file_names) == 0:
            return None
        elif len(file_names) == 1:
            return file_names[0]
        else:
            return file_names

    def _list_to_dict_by_pattern(self, pattern: str, names: list) -> dict:
        """Convert a list of file names to a dict based on regular expression.

        File names can contain an identifier (panel id, disease). This function
        takes a list of file names and a regular expression to convert the list
        of file names to a dictionary with the extracted identifier as key and
        the file name as value. For this to work, the pattern must include
        exactly one group, that is, a part of the regular expression surrounded
        by parentheses.

        Args:
            pattern: a regular expression with one group
            names: list of file names

        Returns:
            a dictionary with the extracted idenifiers as keys and file names
                as values
        """
        result = {}
        compiled_pattern = re.compile(pattern)
        for x in names:
            m = re.fullmatch(compiled_pattern, os.path.basename(x))
            if m:
                key = m.group(1)
                result[key] = x
        return result

    def get_config(self) -> dict:
        """Get the complete configuration data.

        Returns:
            complete configuration as dictionary.
        """
        return self.config

    def get_genie_dir(self) -> str:
        """Get the base directory of GENIE data.

        Returns:
            base directory name for GENIE data.
        """
        return self.config["dir"]["genie"]

    def get_dir(self, key: str) -> str:
        """Get a directory name from the configuration.

        Args:
            key: directory identfier used in the config file.

        Returns:
            directory name.
        """
        return self._file_names(self.get_genie_dir(), "dir", key)

    def get_meta_file_name(self, key: str) -> str:
        """Get the file name of a GENIE meta data file.

        Args:
            key: meta data file identfier used in the config file.

        Returns:
            the file name for this meta file.
        """
        return self._file_names(self.get_genie_dir(), "meta", key)

    def get_data_file_name(self, key: str) -> str:
        """Get the file name of a GENIE  data file.

        Args:
            key: data file identfier used in the config file.

        Returns:
            the file name for this data file.
        """
        return self._file_names(self.get_genie_dir(), "data", key)

    def get_aux_file_name(self, key: str) -> str:
        """Get the file name of an auxiliary data file.

        Args:
            key: auxiliary data file identfier used in the config file.

        Returns:
            the file name for this auxiliary data file.
        """
        return self._file_names(self.get_genie_dir(), "aux", key)

    def load_data_file(self, key: str) -> pd.DataFrame:
        """Load GENIE data from a file or file cache.

        Loads a dataframe from a file. The dataframe is read from the cache
        Parquet cache file if it exists. If the cache file does not exist yet,
        it will be created so that next time loading the data will be faster.

        Args:
            key: data file identfier used in the config file.

        Returns:
            the data loaded from the file or its cache.
        """
        file_name = self.get_data_file_name(key)
        return self.load_file(file_name)

    def load_aux_file(self, key: str) -> pd.DataFrame:
        """Load auxiliary data from a file or file cache.

        Loads a dataframe from a file. The dataframe is read from the cache
        Parquet cache file if it exists. If the cache file does not exist yet,
        it will be created so that next time loading the data will be faster.

        Args:
            key: auxiliary file identfier used in the config file.

        Returns:
            the data loaded from the file or its cache.
        """
        file_name = self.get_aux_file_name(key)
        return self.load_file(file_name)

    def get_case_list_file_names(self) -> dict:
        """Get all files with cast lists.

        Returns:
            dictionary with disease names as keys and file names as values.
        """
        pattern = self.config["data"]["case_lists"][3:]
        d = self.get_dir("case_lists")
        names = self._file_names(d, "data", "case_lists")
        return self._list_to_dict_by_pattern(pattern, names)

    def get_gene_panel_file_names(self) -> dict:
        """Get all files with gene panel descriptions.

        Returns:
            dictionary with panel IDs as keys and file names as values.
        """
        pattern = self.config["data"]["gene_panels"][3:]
        d = self.get_dir("gene_panels")
        names = self._file_names(d, "data", "gene_panels")
        return self._list_to_dict_by_pattern(pattern, names)

    def get_genie_version(self) -> str:
        """Get the version of the GENIE release.

        The version is obtained from the study meta file.

        Returns:
            release version of GENIE data.
        """
        m = Meta(self, "study")
        version = m["description"]
        return version

    def load_file(self, file_name: str) -> pd.DataFrame:
        """Load a dataframe from a file or file cache.

        Loads a dataframe from a file. The dataframe is read from the cache
        Parquet cache file if it exists. If the cache file does not exist yet,
        it will be created so that next time loading the data will be faster.

        Args:
            file_name: the name of the data file.

        Returns:
            the data loaded from the file or its cache.
        """
        # If this is already a parquet file, just read it and return
        if file_name.endswith(".parquet"):
            return pd.read_parquet(file_name)
        # Continue if file_name is not a parquet file.
        parquet_name = self.get_cache_file_name(file_name)
        # If the cache file exists, read it and return. Otherwise, create the
        # cache file.
        if os.path.isfile(parquet_name):
            df = pd.read_parquet(parquet_name)
        else:
            df = pd.read_table(
                file_name, low_memory=False, comment="#"
            ).drop_duplicates()
            df.to_parquet(parquet_name)
        return df

    def get_cache_file_name(self, file_name: str) -> str:
        """Get the name of the Parquet cache file.

        The name of the cache file is created from the original file name by
        stripping the extensions and appending the suffix '.parquet'.

        Args:
            file_name: the name of the data file

        Returns:
            the name of the cache file
        """
        parquet_name = re.sub("\\.(tsv|txt)(\\.gz)?", "", file_name)
        parquet_name += ".parquet"
        return parquet_name


class Meta(dict):
    """This class is a dictionary initialized with data from a meta file.

    GENIE comes with several files with meta information about the actual data
    files. This class is a dictionary that is initialized with key-value pairs
    from such a meta information file.
    """

    def __init__(self, config: Type["Configuration"], key: str):
        """Create a new Meta class.

        Args:
            config: the GENIE source data configuration.
            key: the file name key from the config file.
        """
        super().__init__()
        file_name = config.get_meta_file_name(key)
        with open(file_name, "rt") as fh:
            for line in fh.readlines():
                line = line.rstrip()
                k, v = line.split(": ")
                self[k] = v


class TestedPositions:
    """This class holds all tested positions for all panels."""

    def __init__(self, config: Type["Configuration"]):
        """Creates a new TestedPositions object.

        Args:
            config: the GENIE source data configuration.

        Attributes:
            data (pd.DataFrame): tested positions for all panels.
        """
        self.data = config.load_data_file("genomic_information").set_index(
            "SEQ_ASSAY_ID", drop=True
        )

    def positions_for_panel(self, panel_id: str) -> pd.DataFrame:
        """Get all tested positions for a particular panel.

        Args:
            panel_id: the panel identifier.

        Returns:
            tested positions for specified panel.
        """
        pos = self.data.loc[panel_id, :].reset_index()
        pos = pos[pos.includeInPanel]
        return pos


class TestedMutations:
    """This class knows for each panel-mutation pair if the mutation was tested.

    Checking if a particular mutation was tested by a panel could be done with
    the `TestedPositions` class. However, for thousands of mutations and
    hundreds of panels this takes quite long (about a day for all mutations
    detected in GENIE samples and for all panels). Therefore, in a data
    preparation step done once for each new GENIE release, all mutations found
    in any sample are checked against all panels and a panel-mutation matrix is
    created where matrix elements are _True_ if a mutation is tested by a panel,
    _False_ otherwise.

    The `TestedMutations` class uses this precomputed cache file, loading it
    takes only about 0.5 seconds compared to the about 1 day for computing this
    on the fly.
    """

    def __init__(self, config: Type["Configuration"]):
        """Creates a new TestedMutations object.

        Args:
            config: the GENIE source data configuration.

        Attributes:
            parquet_file (str): name of parquet cache file.
        """
        self.parquet_file = config.get_aux_file_name("panel_mut_tested")

    def is_tested_matrix(
        self, hgvsgs: list = None, panels: list = None
    ) -> pd.DataFrame:
        """Get table of mutations versus panels telling if mutation is tested.

        Get a table with mutations (HGVSGs) as rows and panels as columns where
        table cells are True or False depending on whether a particular mutation
        was tested by a particular panel.

        Args:
            hgvsgs: return table for subset of these mutations (or all mutations
                if None).
            panels: return table for subset of these panels (or all
                panels if None). This can be provided as a list of panel
                identifiers or a list of Panel objects.

        Returns:
            table with information which mutation was tested by which panel.
        """
        if panels is not None and isinstance(panels[0], Panel):
            panels = [p.id for p in panels]
        df = pq.read_table(
            self.parquet_file, columns=panels, use_pandas_metadata=True
        ).to_pandas()
        if hgvsgs is not None:
            df = df.loc[hgvsgs]
        return df


class Panel:
    """This class represents a panel (assay) included in GENIE."""

    def __init__(
        self,
        config: Type["Configuration"],
        global_tested_positions: Type["TestedPositions"],
        global_panel_info: pd.DataFrame,
        panel_id: str,
    ):
        """Create a new Panel object.

        Args:
            config: the GENIE source data configuration.
            global_tested_positions: tested positions for all panels.
            global_panel_info: assay information for all panels.
            panel_id: the panel identifier.

        Attributes:
            id (str): panel identifier.
            description (str): panel description.
            genes (list): genes on panel.
            tested_positions (TestedPositions): tested position for this panel.
            panel_info (dict): assay information for this panel.
        """

        # Read the panel content
        file_name = config.get_gene_panel_file_names()[panel_id]
        with open(file_name, "rt") as fh:
            id, description, genes = fh.readlines()
            self.id = id.rstrip().replace("stable_id: ", "", 1)
            self.description = description.rstrip().replace("description: ", "", 1)
            self.genes = genes.rstrip().split("\t")[1:]
        pos = global_tested_positions.positions_for_panel(panel_id)
        self.tested_positions = pos[pos.includeInPanel]
        self.panel_info = global_panel_info.loc[self.id].to_dict()

    def tested_positions(self) -> pd.DataFrame:
        """Get all tested positions for this panel.

        Returns:
            tested positions for this panel.
        """

        pos = self.global_tested_positions.positions_for_panel(self.id)
        pos = pos[pos.includeInPanel]
        return pos

    def gene_is_on_panel(self, gene_symbol: str) -> bool:
        """Check if a gene is included in this panel.

        Args:
            gene_symbol: HGNC gene symbol.

        Returns:
            True if gene is on panel, else False.
        """
        return gene_symbol in self.genes

    def range_is_tested(self, chr: str, start_pos: int, end_pos: int) -> bool:
        """Check if a genomic range is tested by this panel.

        Args:
            chr: chromosome.
            start_pos: start position of range on chromosome.
            end_pos: end position of range on chromosome.

        Returns:
            True if entire range is tested by this panel, else False.
        """
        chr = chr.replace("chr", "")
        if chr == "M":
            chr = "MT"
        elif chr == "22":
            chr = "X"
        elif chr == "23":
            chr == "Y"
        tested_ranges_including_range = self.tested_positions[
            (self.tested_positions.Chromosome == chr)
            & (self.tested_positions.Start_Position <= start_pos)
            & (self.tested_positions.End_Position >= end_pos)
        ]
        return not tested_ranges_including_range.empty

    def position_is_tested(self, chr: str, pos: int) -> bool:
        """Check if a genomic position is tested by this panel.

        Args:
            chr: chromosome.
            pos: position on chromosome.

        Returns:
            True if position is tested by this panel, else False.
        """
        return self.range_is_tested(chr, pos, pos)

    def get_total_range(self) -> int:
        """Get the overall size of genomic ranges tested by this panel.

        The overall size of tested genomic regions can be used to estimate TMB.

        Returns:
            Sum of lengths of all ranges tested by this panel.
        """
        total_length = sum(
            self.tested_positions.End_Position
            - self.tested_positions.Start_Position
            + 1
        )
        return total_length

    def __str__(self) -> str:
        s = f"{type(self).__name__}: {self.description}"
        return s


class PanelSet:
    """Set of all panels included in GENIE."""

    def __init__(self, config: Type["Configuration"]):
        """Create a new PanelSet.

        Args:
            config: the GENIE source data configuration.

        Attributes:
            panels (dict): all panels used by GENIE, with the panel identifier
                as key and `Panel` as value.
            tested_positions (TestedPositions): tested positions for all panels.
            tested_mutations (TestedMutations): tested mutations for all panels.
            panel_info (pandas.DataFrame): assay information for all panels.

        """
        self.panels = {}
        self.tested_positions = TestedPositions(config)
        self.tested_mutations = TestedMutations(config)
        self.panel_info = pd.read_table(
            config.get_data_file_name("assay_information"),
            index_col="SEQ_ASSAY_ID",
            low_memory=False,
        )
        panel_ids = set(self.tested_positions.data.index)
        for panel_id in panel_ids:
            panel = Panel(config, self.tested_positions, self.panel_info, panel_id)
            self.panels[panel.id] = panel

    def panels_for_position(self, chr: str, pos: int) -> list:
        """Get a list of all panels testing a genomic locaton.

        Args:
            chr: chromosome.
            pos: position on chromosome.

        Returns:
            List of panels probing this location.
        """
        return [p for p in self.panels.values() if p.position_is_tested(chr, pos)]

    def panels_for_range(self, chr: str, start_pos: int, end_pos: int) -> list:
        """Get a list of all panels testing a genomic range.

        Args:
            chr: chromosome.
            start_pos: start position of range.
            end_pos: end position of range.

        Returns:
            List of panels probing this range.
        """
        return [
            p
            for p in self.panels.values()
            if p.range_is_tested(chr, start_pos, end_pos)
        ]

    def panels_for_mutation(self, hgvsg: str) -> dict:
        """Get a list of all panels testing a particular mutation.

        This is done using the precomputed

        Args:
            hgvsg: mutation to be checked.

        Returns:
            dictionary with panel ids as key and `Panel` objects as values
                containing all panels testing this mutation.
        """
        df = self.tested_mutations.is_tested_matrix(hgvsgs=[hgvsg])
        panel_ids = df.columns[df.loc[hgvsg] is True]
        panels = {pid: self.panels[pid] for pid in panel_ids}
        return panels


class PatientInfo:
    """Patient information for all subjects included in GENIE.

    Patient information includes *sex*, *primary race*, *ethnicity*,
    *clinical center*, *contact id*, *dod(?) id*, *year of contact*,
    *dead or alive*, *year of death*.
    """

    def __init__(self, config: Type["Configuration"]):
        """Create a new Patients object.

        Args:
            config: the GENIE source data configuration.

        Attributes:
            data (pd.DataFrame): patient information for all subjects.
        """
        self.data = config.load_data_file("clinical_patient").set_index(
            "PATIENT_ID", drop=True
        )

    def get_info_for_patient(self, patient_id: str) -> dict:
        """Get patient information for a single patient.

        Args:
            patient_id: patient identifier.

        Returns:
            patient information.
        """
        return self.data.loc[patient_id].to_dict()

    def append_info(self, more_info: pd.DataFrame) -> None:
        """Add more columns to the patient information data.

        Args:
            more_info: table with additional columns to be added to the patient
                information. This dataframe must have the `PATIENT_ID` as index.

        Returns:
            nothing
        """
        self.data = self.data.join(more_info)


class SampleInfo:
    """Sample information for all samples included in GENIE.

    Sample information includes *patient id*, *age at sequencing*,
    *Oncotree code*, *sample type*, *sequencing assay id*, *cancer type*,
    *cancer type detailed*, *sample type detailed*.
    """

    def __init__(self, config: Type["Configuration"]):
        """Create a new Samples object.

        Args:
            config: the GENIE source data configuration.

        Attributes:
            data (pd.DataFrame): sample information for all samples.
        """
        self.data = config.load_data_file("clinical_sample").set_index(
            "SAMPLE_ID", drop=True
        )

    def get_info_for_sample(self, sample_id: str) -> dict:
        """Get sample information for for a single sample.

        Args:
            sample_id: sample identifier.

        Returns:
            sample information.
        """
        return self.data.loc[sample_id].to_dict()

    def get_sample_ids(self) -> list:
        """Get all sample sample identifiers.

        Returns:
            sample identifiers
        """
        return self.data.index.to_list()

    def append_info(self, more_info: pd.DataFrame) -> None:
        """Add more columns to the sample information data.

        Args:
            more_info: table with additional columns to be added to the sample
                information. This dataframe must have the `SAMPLE_ID` as index.

        Returns:
            nothing
        """
        self.data = self.data.join(more_info)


class CNA:
    """Copy number alteration data for all samples included in GENIE."""

    def __init__(self, config: Type["Configuration"]):
        """Create a new CNA object.

        Args:
            config: the GENIE source data configuration.

        Attributes:
            data (pd.DataFrame): copy number data for all samples.
                HGNC gene symbol is row index, sample id is column index. If
                CNA information is missing for a gene-sample pair, the value
                is NaN.
        """
        self.data = config.load_data_file("CNA").set_index("Hugo_Symbol", drop=True)

    def get_cna(self, sample_id: str, gene_symbol: str) -> float:
        """Get copy number data for a single sample.

        Args:
            sample_id: sample identifier.

        Returns:
            copy number (NaN if not measured).
        """
        return self.data.loc[gene_symbol, sample_id]


class Mutations:
    """Mutation data for all samples included in GENIE.

    Objects of this class contain all SNVs and indels found in any sample in
    the GENIE data base.
    """

    def __init__(self, config: Type["Configuration"]):
        """Create a new Mutations object.

        Args:
            config: the GENIE source data configuration.

        Attributes:
            data (pd.DataFrame): mutation data for all samples.
            annot (dict): annotations of mutations derived from re-annotation
                by ICA. The dictionary has two keys - "ensembl" and "mane".
                annot["ensembl"] contains annotations for all Ensembl
                transcripts in GENIE (GENIE is based on Ensembl transcripts).
                annot["mane"] contains annotations for all MANE transcripts
                mapping to the mutations in GENIE.
        """

        # Read mutation data
        self.data = config.load_data_file("mutations_extended")

        # Remove entries where REF and ALT are identical (there are four such
        # entries in GENIE 15.0)
        self.data = self.data[self.data.Reference_Allele != self.data.Tumor_Seq_Allele2]

        # Read and merge VCF ID to ICA `vid` mapping table and to `hgvsg`
        id_mapping = config.load_aux_file("vcf_id_to_hgvsg").loc[
            :, ["ID", "vid", "hgvsg"]
        ]
        self.data["ID"] = self.data.apply(self._get_sha256, axis=1)
        self.data = self.data.merge(id_mapping, on="ID")
        self.data = self.data.drop(columns="ID")

        # Read the ICA annotations
        self.annot = {}
        ens_file = config.get_aux_file_name("annot_ensembl")
        mane_file = config.get_aux_file_name("annot_mane")
        if os.path.isfile(ens_file) and os.path.isfile(mane_file):
            self.annot["ensembl"] = config.load_aux_file("annot_ensembl")
            self.annot["mane"] = config.load_aux_file("annot_mane")
        else:
            df = config.load_aux_file("annot")
            df.drop(
                columns=["sample", "genotype", "variantFrequency", "mutationStatus"],
                inplace=True,
            )
            df["transcriptIdNoVersion"] = df["transcriptId"].str.split(".").str.get(0)
            df["hgvsp_short"] = df["hgvsp"].astype(str).map(self.hgvsp3_to_hgvsp1)

            # Get annotations for ENST transcripts originally contained in GENIE
            enst = {x for x in self.data["Transcript_ID"] if isinstance(x, str)}
            self.annot["ensembl"] = (
                df[df["transcriptIdNoVersion"].isin(enst)]
                .drop(columns="transcriptIdNoVersion")
                .drop_duplicates()
            )
            self.annot["ensembl"].to_parquet(ens_file)

            # Get MANE definitions
            mane = pd.read_table(config.get_aux_file_name("mane"), low_memory=False)
            mane["transcriptIdNoVersion"] = mane["RefSeq_nuc"].str.split(".").str.get(0)
            self.annot["mane"] = (
                df.merge(
                    mane[["transcriptIdNoVersion", "MANE_status"]],
                    how="inner",
                    on="transcriptIdNoVersion",
                )
                .drop(columns="transcriptIdNoVersion")
                .drop_duplicates()
            )
            self.annot["mane"].to_parquet(mane_file)

    def _get_sha256(self, row):
        vid = ":".join(
            [
                row.Chromosome,
                str(row.Start_Position),
                row.Reference_Allele,
                row.Tumor_Seq_Allele2,
            ]
        )
        return hashlib.sha256(vid.encode("UTF-8")).hexdigest()

    def _check_universe(self, universe: str) -> None:
        """Check if universe is one of *ensembl* and *mane*.

        This function checks if the specified universe is one of *ensembl* and
        *mane*. If it is not, it raises a `ValueError` exception. Otherwise, it
        does nothing.
        """
        if universe not in {"ensembl", "mane"}:
            raise ValueError("universe must be one of 'ensembl' and 'mane'.")

    def get_detected_mutations(
        self, gene_symbols: list, universe: str = "ensembl"
    ) -> pd.DataFrame:
        """Get detected mutations for a list of genes.

        This function returns mutations that are found for a the specified genes
        in all samples included in GENIE. If a mutation is not included in the
        returned values for a particular sample, this does not mean that the
        gene is of wild type for this sample because the mutation may not be on
        the panel used for that sample.

        Args:
            gene_symbols: list of HGNC gene symbols.
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.

        Returns:
            all mutations of these genes found in GENIE samples.
        """
        self._check_universe(universe)
        hgvsgs = self.get_unique_mutations(gene_symbols, universe)
        # Need to filter for gene symbols again because more than one gene can
        # be affected by a hgvsg
        df = self.data[
            self.data.hgvsg.isin(hgvsgs) & self.data.Hugo_Symbol.isin(gene_symbols)
        ]
        return df

    def get_unique_mutations(
        self, gene_symbols: list = None, universe: str = "ensembl"
    ) -> list:
        """Get a unique list of mutations found in GENIE for a list of genes.

        For the *ensembl* universe, all variants affecting Ensembl transcripts
        that were included in the original GENIE mutation data file are
        returned. For the *mane* universe, only variants affecting MANE
        transcripts are returned.

        Args:
            gene_symbols: list of genes to query.
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.

        Returns:
            hgvsg for all unique mutations for the specified genes.
        """
        self._check_universe(universe)
        df = self.annot[universe]
        hgvsgs = set(df.loc[df.hgnc.isin(gene_symbols), "hgvsg"])
        return list(hgvsgs)

    def get_sample_mutation_profiles(
        self,
        gene_symbols: list,
        sample_ids: list,
        sample_info: SampleInfo,
        tested_mutations: TestedMutations,
        universe: str = "ensembl",
        hgvsgs: list = None,
    ) -> pd.DataFrame:
        """Get mutation profiles of genes for all tested samples.

        While the function `get_detected_mutations` returns
        mutation-sample-pairs only for those samples where a mutation was
        actually detected, this function adds all samples where a mutation was
        also tested but not found. The returned dataframe contains a column
        `mutated` that is either *True* or *False*.

        Args:
            gene_symbols: genes to include in the result
            sample_ids: sample identifiers of samples to include
            sample_info: annnotation of all samples
            tested_mutations: cache providing hgvsg vs panel matrix
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.
            hgvsgs: keep only these hvsgs, i.e., exclude all other hgvsgs for
                the specified gene_symbols. If None, keep all hgvsgs.
        """
        self._check_universe(universe)

        # Get HGVSG for all GENIE mutations for the specified genes
        all_hgvsgs = self.get_unique_mutations(gene_symbols, universe)
        if hgvsgs:
            hgvsgs = [x for x in all_hgvsgs if x in hgvsgs]
        else:
            hgvsgs = all_hgvsgs

        # Get mutations per sample for these HGVSGs
        muts = self.get_detected_mutations(gene_symbols, universe).rename(
            columns={"Tumor_Sample_Barcode": "SAMPLE_ID"}
        )
        muts = muts[muts.hgvsg.isin(hgvsgs)]
        muts = muts[muts.SAMPLE_ID.isin(sample_ids)]
        muts = muts[["hgvsg", "SAMPLE_ID"]].drop_duplicates()

        # Reduce list of unique mutations to those found in the specified
        # subset of samples
        hgvsgs = list(set(muts.hgvsg))
        muts = muts.assign(mutated=True).set_index(["hgvsg", "SAMPLE_ID"], drop=True)

        # Get info which of the hgvsgs was tested by which panel
        tested_muts = (
            tested_mutations.is_tested_matrix(hgvsgs)
            .melt(ignore_index=False, var_name="SEQ_ASSAY_ID")
            .query("value")  # Keep only the "True"s
            .reset_index()
            .drop(columns="value")
        )

        # Get a mapping from sample id to the panel id used for the sample
        sample_to_assay = sample_info.data.loc[sample_ids, "SEQ_ASSAY_ID"].reset_index()

        # Merge tested mutations with samples and assign "False" to everything
        tested_muts = (
            tested_muts.merge(sample_to_assay, on="SEQ_ASSAY_ID")
            .drop(columns="SEQ_ASSAY_ID")
            .assign(mutated=False)
            .set_index(["hgvsg", "SAMPLE_ID"], drop=True)
        )

        # Update "False" with "True" where a mutation was found
        tested_muts.update(muts)
        return tested_muts

    def get_mutation_annotations(
        self, hgvsgs: list, universe: str = "ensembl"
    ) -> pd.DataFrame:
        """Get annotations for a list of mutations.

        Get ICA re-annotations for a list of mutations specified by HGVSG.
        The returned dataframe can include more than one row per unique genomic
        variant. For example, a locus can be a *downstream\\_gene\\_variant* for
        one gene and a *3\\_prime\\_UTR\\_variant* for another gene.
        Furthermore, if the *mane* universe is specified, there can be a *MANE
        Select* and one or more *MANE Plus Clinical* transcript variants
        covering the genomic location, and there can be more than one gene
        covering that location. There can be up to 15 different MANE transcripts
        for a genomic locus. It depends on the use case which if these
        transcripts shall be included in downstream analyses. Therefore, the
        user of this function needs to add appropriate filtering to the returned
        annotations.

        Args:
            hgvsgs: list of HGVSG specifications of mutations.
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.

        Returns:
            Annotations for specified mutations.
        """
        self._check_universe(universe)
        df = self.annot[universe]
        df = df[df.hgvsg.isin(hgvsgs)]
        return df

    def hgvsp3_to_hgvsp1(self, hgvsp: str) -> str:
        """Translate 3-letter amino acid codes to 1-letter amino acid codes.

        Args:
            hgvsp: amino acid change with 3-letter amino acid codes

        Returns:
            hgvsp_short with 1-letter amino acid codes
        """
        if hgvsp is None:
            return None
        match = re.match("^.*p\\.\\(?([^)]*)\\)?", hgvsp)
        if match:
            hgvsp = match.group(1)
        aa_pattern = re.compile(r"([A-Z][a-z]{2})")
        hgvsp = re.sub(aa_pattern, lambda match: seq1(match.group(0)), hgvsp)
        return hgvsp


class Genie:
    """This is the main class for this module, holding all GENIE data.

    An object of this class holds and provides all data from a Genie release.
    This includes all panels (assays) and the genomic regions tested by each
    panel, patient information, sample information, copy number data, mutation
    data.
    """

    def __init__(self, genie_dir: str, config_file: str = None, verbose: bool = True):
        """Create a new `Genie` object.

        Args:
            genie_dir: path name of the directory with Genie data files.
            config_file: path name of an optional user config file, defining
                alternative names for GENIE data files if names should change
                for future versions of GENIE.

        Attributes:
            panel_set (PanelSet): all panels (assays) included in GENIE.
            patient_info (PatientInfo): clinical information about patients.
            sample_info (SampleInfo): sample annotations.
            cna (CNA): copy number data.
            mutations (Mutations): mutation data.
            version (str): Genie version number.
        """
        self.config = Configuration(genie_dir, config_file)
        if verbose:
            print("Loading panels", end=" ... ", flush=True)
        self.panel_set = PanelSet(self.config)
        if verbose:
            print("done.", flush=True)
            print("Loading patients", end=" ... ", flush=True)
        self.patient_info = PatientInfo(self.config)
        if verbose:
            print("done.", flush=True)
            print("Loading samples", end=" ... ", flush=True)
        self.sample_info = SampleInfo(self.config)
        if verbose:
            print("done.", flush=True)
            print("Loading copy number data", end=" ... ", flush=True)
        self.cna = CNA(self.config)
        if verbose:
            print("done.", flush=True)
            print("Loading mutations", end=" ... ", flush=True)
        self.mutations = Mutations(self.config)
        if verbose:
            print("done.", flush=True)
        self.version = self.config.get_genie_version()

    def summary(self) -> str:
        """Get a summary of data in this GENIE release.

        The summary includes the number of panels, number of patiens, number of
        samples, number of genes tested by at least one panel, the number of
        genes with copy number alterations identfied at least once.
        """
        version = self.config.get_genie_version()
        num_samples_with_mut_data = len(
            self.mutations.data.Tumor_Sample_Barcode.unique()
        )
        num_patients_with_mut_data = len(
            self.sample_info.data.loc[
                self.mutations.data.Tumor_Sample_Barcode.unique(), "PATIENT_ID"
            ].unique()
        )
        num_samples_with_cna_data = self.cna.data.shape[1]
        num_patients_with_cna_data = len(
            self.sample_info.data.loc[self.cna.data.columns, "PATIENT_ID"].unique()
        )
        num_unique_muts = len(self.mutations.data.hgvsg.unique())
        s = "\n".join(
            [
                f"{version}",
                f"{'=' * len(version)}",
                f"Number of panels: {len(self.panel_set.panels)}",
                f"Number of patients: {len(self.patient_info.data)}",
                f"Number of samples: {len(self.sample_info.data)}",
                f"Number of genes with mutations (tested by at least one panel): "
                f"{len(set(self.mutations.data.Hugo_Symbol))}",
                f"Number of genes with CNA (tested by at least one panel): "
                f"{len(self.cna.data)})",
                f"Number of samples with mutation data: {num_samples_with_mut_data}",
                f"Number of patients with mutation data: {num_patients_with_mut_data}",
                f"Number of samples with CNA data: {num_samples_with_cna_data}",
                f"Number of patients with CNA data: {num_patients_with_cna_data}",
                f"Number of unique mutations: {num_unique_muts}",
            ]
        )
        return s

    def __str__(self) -> str:
        return self.summary()

    def __repr__(self) -> str:
        s = f"genie.{type(self).__name__}" f'("{self.config.get_genie_dir()}")'
        return s

    def _check_cancer_types(self, cancer_types: list = None) -> None:
        """Check if cancer types are contained in Genie.

        This function checks whether the cancer types from a list are all
        contained in Genie. If some of them are not, a `ValueError` exception
        is raised. Otherwise, nothing happens. For the sake of convenience,
        `None` can also be provided as the list of cancer types to simplify
        coding in the calling functions.

        Args:
            cancer_types: list of cancer types to check.
        """
        if cancer_types:
            wrong_types = set(cancer_types) - set(self.sample_info.data.CANCER_TYPE)
            if wrong_types:
                raise ValueError(f"Unknown cancer types: {wrong_types}")

    def _check_cancer_subtypes(self, cancer_subtypes: list = None) -> None:
        """Check if cancer subtypes are contained in Genie.

        This function checks whether the cancer subtypes from a list are all
        contained in Genie. If some of them are not, a `ValueError` exception
        is raised. Otherwise, nothing happens. For the sake of convenience,
        `None` can also be provided as the list of cancer subtypes to simplify
        coding in the calling functions.

        Args:
            cancer_subtypes: list of cancer subtypes to check (ONCOTREE codes).
        """
        if cancer_subtypes:
            wrong_subtypes = set(cancer_subtypes) - set(
                self.sample_info.data.ONCOTREE_CODE
            )
            if wrong_subtypes:
                raise ValueError(f"Unknown cancer types: {wrong_subtypes}")

    def get_cancer_types(self) -> list:
        """Get a list of cancer type names as used in Genie.

        Returns:
            all cancer types used in Genie.
        """
        return sorted(self.sample_info.data.CANCER_TYPE.unique())

    def get_cancer_subtypes(self, cancer_types: list = None) -> pd.DataFrame:
        """Get cancer subtypes for all or selected cancer types.

        Args:
            cancer_types: optional list of cancer types for which to return
                cancer subtypes. If not specified, all cancer types are
                returned.

        Returns:
            cancer subtypes for all or specified cancer types. Cancer types as
                index of the returned dataframe, Oncotree codes and and names of
                cancer subtypes as columns.
        """
        df = self.sample_info.data
        if cancer_types:
            self._check_cancer_types(cancer_types)
            df = df[df.CANCER_TYPE.isin(cancer_types)]
        df = (
            df[["CANCER_TYPE", "ONCOTREE_CODE", "CANCER_TYPE_DETAILED"]]
            .drop_duplicates()
            .sort_values(["CANCER_TYPE", "ONCOTREE_CODE"])
            .set_index("CANCER_TYPE", drop=True)
        )
        return df

    def get_cancer_type_patient_counts(
        self, cancer_types_to_keep: list = None
    ) -> pd.Series:
        """Get cancer types and their patient numbers.

        Get the number of patients for each cancer type included in Genie. If
        the optional argument `cancer_types_to_keep`is specified, all other
        cancer types are summarized as *other*.

        The counts returned are for patients, not samples. Patients can have
        more than one sample and even more than one cancer type. Each patient
        is counted only once per cancer type, independent of the number of
        samples for that patient. However, if a patient has samples of more
        than one cancer type, the patient will be counted multiple times, once
        per cancer type. As a consequence, the sum of all counts returned by
        this function is larger than the total number of patients in Genie.

        Args:
            cancer_types_to_keep: cancer types not to summarized in the *other*
                category. If not specified, all cancer types will be returned.

        Returns:
            number of patients per cancer type
        """
        counts = (
            self.sample_info.data[["PATIENT_ID", "CANCER_TYPE"]]
            .drop_duplicates()
            .value_counts("CANCER_TYPE")
        )
        if cancer_types_to_keep:
            self._check_cancer_types(cancer_types_to_keep)
            other_types = list(set(counts.index) - set(cancer_types_to_keep))
            other_count = counts[other_types].sum()
            counts = counts[cancer_types_to_keep]
            counts["other"] = other_count
        return counts

    def get_cancer_subtype_sample_counts(
        self, cancer_types: list = None, cancer_subtypes_to_keep: list = None
    ):
        """Get cancer subtypes and their sample numbers.

        Get the number of samples for each cancer subtype included in Genie.
        This function returns sample numbers and not patient numbers because
        quite often patients have multiple samples of the same cancer type but
        slightly different cancer subtype annotations.
        """
        df = self.sample_info.data
        if cancer_types:
            self._check_cancer_types(cancer_types)
            df = df[df.CANCER_TYPE.isin(cancer_types)]
        counts = (
            df.reset_index()[["SAMPLE_ID", "CANCER_TYPE", "ONCOTREE_CODE"]]
            .drop_duplicates()
            .value_counts(["CANCER_TYPE", "ONCOTREE_CODE"])
        )
        if cancer_subtypes_to_keep:
            self._check_cancer_subtypes(cancer_subtypes_to_keep)
            other_subtypes = list(
                set(counts.index.get_level_values("ONCOTREE_CODE"))
                - set(cancer_subtypes_to_keep)
            )
            other_count = counts.loc[(slice(None), other_subtypes)].sum()
            counts = counts.loc[(slice(None), cancer_subtypes_to_keep)]
            counts["other"] = other_count
        return counts

    def get_sample_mutation_profiles(
        self,
        gene_symbols: list,
        hgvsgs: list = None,
        cancer_types: list = None,
        cancer_subtypes: list = None,
        sample_ids: list = None,
        universe: str = "ensembl",
    ) -> pd.DataFrame:
        """Get a mutation profile of selected genes across samples.

        This function returns the mutation profile of the specified genes across
        either all samples or across selected samples, based on sample ids or
        cancer types or cancer subtypes. If more than one sample selection
        criterion is specified, then the intersection of the specified criteria
        is used to determine the final set of samples. The number of genes
        should be small, because otherwise the volume of results is very high
        and it will take very long for the function to return.

        For each mutation detected in at least one sample, it is checked for
        each sample whether the mutation was tested by the panel used for this
        sample. If this is the case, `mutated` is returned as "False", unless
        the mutation was in fact detected for this sample, in which case "True"
        is returned. Samples that have not been tested for a mutation are not
        included in the result.

        Args:
            gene_symbols: genes to include in the result.
            hgvsgs: keep only these hvsgs, i.e., exclude all other hgvsgs for
                the specified gene_symbols. If None, keep all hgvsgs.
            cancer_types: cancer types to include in the result.
            cancer_subtypes: cancer subtypes to include in the result (Oncotree
                codes).
            sample_ids: sample identifiers to include in the result.
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.

        Returns:
            Mutation status (MUT or WT) for all tested mutations.
        """
        self.mutations._check_universe(universe)
        if sample_ids:
            sample_ids = set(sample_ids).intersection(self.sample_info.data.index)
        else:
            sample_ids = set(self.sample_info.data.index)
        if cancer_types:
            self._check_cancer_types(cancer_types)
            sample_ids = sample_ids.intersection(
                self.sample_info.data.loc[
                    self.sample_info.data.CANCER_TYPE.isin(cancer_types)
                ].index
            )
        if cancer_subtypes:
            self._check_cancer_subtypes(cancer_subtypes)
            sample_ids = sample_ids.intersection(
                self.sample_info.loc[
                    self.sample_info.data.ONCOTREE_CODE.isin(cancer_subtypes)
                ].index.to_list()
            )
        sample_ids = list(sample_ids)
        return self.mutations.get_sample_mutation_profiles(
            gene_symbols=gene_symbols,
            sample_ids=sample_ids,
            sample_info=self.sample_info,
            tested_mutations=self.panel_set.tested_mutations,
            universe=universe,
            hgvsgs=hgvsgs,
        )

    def get_nucleic_acid_level_frequencies(
        self,
        gene_symbols: list,
        hgvsgs: list = None,
        cancer_types: list = None,
        cancer_types_to_keep: list = None,
        cancer_subtypes: list = None,
        cancer_subtypes_to_keep: list = None,
        sample_ids: list = None,
        cancer_subtype_resolution: bool = False,
        extra_group_columns: dict = None,
        universe: str = "ensembl",
        precision: int = 1,
    ) -> pd.DataFrame:
        """Get a nucleic acid level mutation frequencies.

        This function returns the mutation frequencies of the specified genes
        across either all samples or across selected samples after filtering by
        sample ids or cancer types or cancer subtypes. Counts and frequencies
        are for patient numbers, not sample numbers. A patient is considered
        having a particular mutation if at least one sample of this patient has
        that mutation.

        The number of genes should be small, because otherwise the volume of
        results is very high and it will take very long for the function to
        return.

        For each mutation detected in at least one sample, it is checked for
        each sample whether the mutation was tested by the panel used for this
        sample. Samples not tested for a mutation are not included in the
        returned counts and frequencies.

        If the `cancer_types_to_keep` argument is specified, all cancer types
        not included in this list are summarized in the `other` category.

        If the `cancer_subtypes_to_keep` argument is specified, all cancer
        subtypes not included in this list are summarized in the `other`
        category.

        If `extra_group_columns` is specified, counts are provided for each
        value of these extra columns. One example for this is `PRIMARY_RACE`.
        The argument `extra_group_columns` needs to be provided as a dictionary,
        where the dict key is the name of the column, and the dict value is a
        list of all values of that column for which an extra row in the count
        table will be provided, while other values will be summarized as
        `other`. For example, if counts shall be returned separately for the
        "White", "Black" and "Asian" population and the remaining patients shall
        be summarized as "other", the `extra_group_columns` argument needs to
        be set to `{"PRIMARY_RACE": ["White", "Black", "Asian"]}`. If no
        aggregation to an `other` category is desired, set the dict value to
        `None`. For example, `{"PRIMARY_RACE": None}` would return separate
        counts for each race.

        Args:
            gene_symbols: genes to include in the result.
            hgvsgs: keep only these hvsgs, i.e., exclude all other hgvsgs for
                the specified gene_symbols. If None, keep all hgvsgs.
            cancer_types: cancer types to include in the result. All if None.
            cancer_types_to_keep: cancer types not to summarize in the *other*
                category. No summarization if None.
            cancer_subtypes: cancer subtypes to include in the result (Oncotree
                codes). All if None.
            cancer_subtypes_to_keep: cancer subtypes not to summarize in the
                *other* category (Oncotree codes). No summarization if None.
            sample_ids: Keep only these samples. No additional filtering if
                None.
            cancer_subtype_resolution: if True, provide frequencies by cancer
                subtype, otherwise return frequencies summarized by cancer type.
                If cancer_subtypes_to_keep is specified, this is automatically
                set to True, otherwise it defaults to False.
            extra_group_columns: by default, counts are returned for each
                combination of hgvsg, cancer type (CANCER_TYPE), and possibly
                cancer subtype (ONCOTREE_CODE). If counts should be split by
                additional factors, such as PRIMARY_RACE, it can be added here.
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.
            precision: number of fractional digits of formatted allele frequency
                percentages.
        """
        # Get the mutation profiles
        sample_mutation_profiles = self.get_sample_mutation_profiles(
            gene_symbols=gene_symbols,
            cancer_types=cancer_types,
            cancer_subtypes=cancer_subtypes,
            sample_ids=sample_ids,
            universe=universe,
            hgvsgs=hgvsgs,
        )
        freqs = self._get_mutation_frequencies(
            sample_mutation_profiles=sample_mutation_profiles,
            cancer_types_to_keep=cancer_types_to_keep,
            cancer_subtypes_to_keep=cancer_subtypes_to_keep,
            cancer_subtype_resolution=cancer_subtype_resolution,
            extra_group_columns=extra_group_columns,
            precision=precision,
        )
        return freqs

    def get_amino_acid_level_frequencies(
        self,
        gene_symbols: list,
        hgvsgs: list = None,
        hgvsps: list = None,
        cancer_types: list = None,
        cancer_types_to_keep: list = None,
        cancer_subtypes: list = None,
        cancer_subtypes_to_keep: list = None,
        sample_ids: list = None,
        cancer_subtype_resolution: bool = False,
        extra_group_columns: dict = None,
        universe: str = "ensembl",
        precision: int = 1,
        panel_coverage_threshold: float = 0.8,
        impute: bool = False,
    ) -> pd.DataFrame:
        """Get amino acid level mutation frequencies.

        This function returns the amino acid level mutation frequencies of the
        specified genes across either all samples or across selected samples
        after filtering based on sample ids or cancer types or cancer subtypes.
        Counts and frequencies are for patient numbers, not sample numbers. A
        patient is considered having a particular mutation if at least one
        sample of this patient has that mutation.

        There may be several different nucleic acid mutations leading to the
        same protein sequence change. At amino acid level, these different
        nucleic acid variants are integrated to a an amino acid level, that is,
        a sample is considered mutated for a particular amino acid change if it
        has any of the nucleic acid changes leading to that amino acid variant,
        and it is considered wild type if it has none of these mutations.

        The number of genes should be small, because otherwise the volume of
        results is very high and it will take very long for the function to
        return.

        For each mutation detected in at least one sample, it is checked for
        each sample whether the mutation was tested by the panel used for this
        sample. Samples not tested for the majority of mutation (less than
        `panel_coverage_threshold`) are not included in the returned counts and
        frequencies. For the remaining samples, missing values are imputed or
        set to wild type, depending on the `impute` argument. See
        `get_imputed_sample_mutation_profiles` for details.

        If the `cancer_types_to_keep` argument is specified, all cancer types
        not included in this list are summarized in the `other` category.

        If the `cancer_subtypes_to_keep` argument is specified, all cancer
        subtypes not included in this list are summarized in the `other`
        category.

        If `extra_group_columns` is specified, counts are provided for each
        value of these extra columns. One example for this is `PRIMARY_RACE`.
        The argument `extra_group_columns` needs to be provided as a dictionary,
        where the dict key is the name of the column, and the dict value is a
        list of all values of that column for which an extra row in the count
        table will be provided, while other values will be summarized as
        `other`. For example, if counts shall be returned separately for the
        "White", "Black" and "Asian" population and the remaining patients shall
        be summarized as "other", the `extra_group_columns` argument needs to be
        set to `{"PRIMARY_RACE": ["White", "Black", "Asian"]}`. If no
        aggregation to an `other` category is desired, set the dict value to
        `None`. For example, `{"PRIMARY_RACE": None}` would return separate
        counts for each race.

        Args:
            gene_symbols: genes to include in the result.
            hgvsgs: keep only these hvsgs, i.e., exclude all other hgvsgs for
                the specified gene_symbols. If None, keep all hgvsgs.
            hgvsps: keep only these hgvsps, i.e., exclude all other hgvsps for
                the specified gene symbols. If None, keep all hgvsps.
            cancer_types: cancer types to include in the result. All if None.
            cancer_types_to_keep: cancer types not to summarize in the *other*
                category. No summarization if None.
            cancer_subtypes: cancer subtypes to include in the result (Oncotree
                codes). All if None.
            cancer_subtypes_to_keep: cancer subtypes not to summarize in the
                *other* category (Oncotree codes). No summarization if None.
            sample_ids: Keep only these samples. No additional filtering if
                None.
            cancer_subtype_resolution: if True, provide frequencies by cancer
                subtype, otherwise return frequencies summarized by cancer type.
                If cancer_subtypes_to_keep is specified, this is automatically
                set to True, otherwise it defaults to False.
            extra_group_columns: by default, counts are returned for each
                combination of gene symbol, cancer type (CANCER_TYPE), and
                possibly cancer subtype (ONCOTREE_CODE). If counts should be
                split by additional factors, such as PRIMARY_RACE, it can be
                added here.
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.
            precision: number of fractional digits of formatted allele frequency
                percentages.
            panel_coverage_threshold: Exclude panels and samples profiled with
                these panels if the fraction of mutations tested by these panels
                is below that threshold. If this is set to 1, only panels that
                test all mutations are included. In that case, no imputation is
                needed. See `get_imputed_sample_muation_profiles` for details.
            impute: Whether to impute missing values (mutations not tested). If
                this is set to False (the default), then missing values are
                replaced with wild type. See
                `get_imputed_sample_muation_profiles` for details.

        Returns:
            Amino acid level counts and frequencies.
        """
        # Get the mutation profiles
        sample_mutation_profiles = self.get_sample_mutation_profiles(
            gene_symbols=gene_symbols,
            cancer_types=cancer_types,
            cancer_subtypes=cancer_subtypes,
            sample_ids=sample_ids,
            universe=universe,
            hgvsgs=hgvsgs,
        )
        # Exclude samples with panels that don't test enough mutations, and
        # impute missing values
        sample_mutation_profiles = self.get_imputed_sample_mutation_profiles(
            sample_mutation_profiles=sample_mutation_profiles,
            panel_coverage_threshold=panel_coverage_threshold,
            impute=impute,
        )
        # Aggregate to amino acid level
        sample_mutation_profiles = self.aggregate_to_amino_acid_level(
            sample_mutation_profiles=sample_mutation_profiles,
            gene_symbols=gene_symbols,
            universe=universe,
        )
        # Filter for hgvsps if specified
        if hgvsps:
            sample_mutation_profiles = sample_mutation_profiles.query(
                "hgvsp.isin(@hgvsps)"
            )
        # Add counts and frequencies
        freqs = self._get_mutation_frequencies(
            sample_mutation_profiles=sample_mutation_profiles,
            cancer_types_to_keep=cancer_types_to_keep,
            cancer_subtypes_to_keep=cancer_subtypes_to_keep,
            cancer_subtype_resolution=cancer_subtype_resolution,
            extra_group_columns=extra_group_columns,
            precision=precision,
        )
        return freqs

    def get_gene_level_frequencies(
        self,
        gene_symbols: list,
        hgvsgs: list = None,
        cancer_types: list = None,
        cancer_types_to_keep: list = None,
        cancer_subtypes: list = None,
        cancer_subtypes_to_keep: list = None,
        sample_ids: list = None,
        cancer_subtype_resolution: bool = False,
        extra_group_columns: dict = None,
        universe: str = "ensembl",
        min_mutations: int = 1,
        precision: int = 1,
        panel_coverage_threshold: float = 0.8,
        impute: bool = False,
    ) -> pd.DataFrame:
        """Get a mutation frequencies of selected genes across samples.

        This function returns the gene level mutation frequencies of the
        specified genes across either all samples or across selected samples
        after filtering by sample ids or cancer types or cancer subtypes. Counts
        and frequencies are for patient numbers, not sample numbers. A patient
        is considered having a particular mutation if at least one sample of
        this patient has that mutation.

        The number of genes should be small, because otherwise the volume of
        results is very high and it will take very long for the function to
        return.

        For each mutation detected in at least one sample, it is checked for
        each sample whether the mutation was tested by the panel used for this
        sample. Samples not tested for the majority of mutations (less than
        `panel_coverage_threshold`) are not included in the returned counts and
        frequencies. For the remaining samples, missing values are imputed or
        set to wild type, depending on the `impute` argument. See
        `get_imputed_sample_mutation_profiles` for details. Imputation is very
        time consuming, and for large enough `panel_coverage_threshold`s the
        differences in the results are marginal, which is why imputation is
        switched off by default.

        If the `cancer_types_to_keep` argument is specified, all cancer types
        not included in this list are summarized in the `other` category.

        If the `cancer_subtypes_to_keep` argument is specified, all cancer
        subtypes not included in this list are summarized in the `other`
        category.

        If `extra_group_columns` is specified, counts are provided for each
        value of these extra columns from sample or patient annotations. One
        example for this is `PRIMARY_RACE`, another one is `SEX`. The argument
        `extra_group_columns` needs to be provided as a dictionary, where the
        dict key is the name of the column, and the dict value is a list of all
        values of that column for which an extra row in the count table will be
        provided, while other values will be summarized as `other`. For example,
        if counts shall be returned separately for the "White", "Black" and
        "Asian" population and the remaining patients shall be summarized as
        "other", the `extra_group_columns` argument needs to be set to
        `{"PRIMARY_RACE": ["White", "Black", "Asian"]}`. If no aggregation to an
        `other` category is desired, set the dict value to `None`. For example,
        `{"PRIMARY_RACE": None}` would return separate counts for each race.

        Args:
            gene_symbols: genes to include in the result.
            hgvsgs: keep only these hvsgs, i.e., exclude all other hgvsgs for
                the specified gene_symbols. If None, keep all hgvsgs.
            cancer_types: cancer types to include in the result. All if None.
            cancer_types_to_keep: cancer types not to summarize in the *other*
                category. No summarization if None.
            cancer_subtypes: cancer subtypes to include in the result (Oncotree
                codes). All if None.
            cancer_subtypes_to_keep: cancer subtypes not to summarize in the
                *other* category (Oncotree codes). No summarization if None.
            sample_ids: Keep only these samples. No additional filtering if
                None.
            cancer_subtype_resolution: if True, provide frequencies by cancer
                subtype, otherwise return frequencies summarized by cancer type.
                If cancer_subtypes_to_keep is specified, this is automatically
                set to True, otherwise it defaults to False.
            extra_group_columns: by default, counts are returned for each
                combination of gene symbol, cancer type (CANCER_TYPE), and
                possibly cancer subtype (ONCOTREE_CODE). If counts should be
                split by additional factors, such as PRIMARY_RACE, it can be
                added here.
            universe: one of "ensembl" or "mane"; use "ensembl" (the default) to
                include the original Ensembl transcripts from GENIE, and use
                "mane" to use RefSeq MANE transcripts instead.
            min_mutations: Call a gene mutated if it has at least that many
                mutations. Should be 1 for oncogenes and may be set to 2 for
                tumor suppressor genes.
            precision: number of fractional digits of formatted allele frequency
                percentages.
            panel_coverage_threshold: Exclude panels and samples profiled with
                these panels if the fraction of mutations tested by these panels
                is below that threshold. If this is set to 1, only panels that
                test all mutations are included. In that case, no imputation is
                needed. See `get_imputed_sample_muation_profiles` for details.
            impute: Whether to impute missing values (mutations not tested). If
                this is set to False (the default), then missing values are
                replaced with wild type. See
                `get_imputed_sample_muation_profiles` for details.

        Returns:
            gene level counts and frequencies.
        """
        # Get the mutation profiles
        sample_mutation_profiles = self.get_sample_mutation_profiles(
            gene_symbols=gene_symbols,
            cancer_types=cancer_types,
            cancer_subtypes=cancer_subtypes,
            sample_ids=sample_ids,
            universe=universe,
            hgvsgs=hgvsgs,
        )
        # Filter for hgvsgs if specified
        if hgvsgs:
            sample_mutation_profiles = sample_mutation_profiles.query(
                "hgvsg.isin(@hgvsgs)"
            )
        # Exclude samples with panels that don't test enough mutations, and
        # impute missing values
        sample_mutation_profiles = self.get_imputed_sample_mutation_profiles(
            sample_mutation_profiles=sample_mutation_profiles,
            panel_coverage_threshold=panel_coverage_threshold,
            impute=impute,
        )
        # Aggregate to gene level
        sample_mutation_profiles = self.aggregate_to_gene_level(
            sample_mutation_profiles=sample_mutation_profiles,
            gene_symbols=gene_symbols,
            universe=universe,
            min_mutations=min_mutations,
        )
        # Add counts and frequencies
        freqs = self._get_mutation_frequencies(
            sample_mutation_profiles=sample_mutation_profiles,
            cancer_types_to_keep=cancer_types_to_keep,
            cancer_subtypes_to_keep=cancer_subtypes_to_keep,
            cancer_subtype_resolution=cancer_subtype_resolution,
            extra_group_columns=extra_group_columns,
            precision=precision,
        )
        return freqs

    def _get_mutation_frequencies(
        self,
        sample_mutation_profiles: pd.DataFrame,
        cancer_types_to_keep: list = None,
        cancer_subtypes_to_keep: list = None,
        cancer_subtype_resolution: bool = False,
        extra_group_columns: dict = None,
        precision: int = 1,
    ) -> pd.DataFrame:
        """Get a mutation counts and frequencies across samples.

        If the `cancer_types_to_keep` argument is specified, all cancer types
        not included in this list are summarized in the `other` category.

        If the `cancer_subtypes_to_keep` argument is specified, all cancer
        subtypes not included in this list are summarized in the `other`
        category.

        If `extra_group_columns` is specified, counts are provided for each
        value of these extra columns. One example for this is `PRIMARY_RACE`.
        The argument `extra_group_columns` needs to be provided as a dictionary,
        where the dict key is the name of the column, and the dict value is a
        list of all values of that column for which an extra row in the count
        table will be provided, while other values will be summarized as
        `other`. For example, if counts shall be returned separately for the
        "White", "Black" and "Asian" population and the remaining patients shall
        be summarized as "other", the `extra_group_columns` argument needs to
        be set to `{"PRIMARY_RACE": ["White", "Black", "Asian"]}`. If no
        aggregation to an `other` category is desired, set the dict value to
        `None`. For example, `{"PRIMARY_RACE": None}` would return separate
        counts for each race.

        Args:
            cancer_types_to_keep: cancer types not to summarize in the *other*
                category. No summarization if None.
            cancer_subtypes_to_keep: cancer subtypes not to summarize in the
                *other* category (Oncotree codes). No summarization if None.
            cancer_subtype_resolution: if True, provide frequencies by cancer
                subtype, otherwise return frequencies summarized by cancer type.
                If cancer_subtypes_to_keep is specified, this is automatically
                set to True, otherwise it defaults to False.
            extra_group_columns: by default, counts are returned for each
                combination of hgvsg, cancer type (CANCER_TYPE), and possibly
                cancer subtype (ONCOTREE_CODE). If counts should be split by
                additional factors, such as PRIMARY_RACE, it can be added here.
            precision: number of fractional digits of formatted allele frequency
                percentages.

        """
        muts = sample_mutation_profiles.copy()
        main_index = muts.index.names[0]  # can be hgvsg or hgvsp or gene_symbol
        # Add sample annotation
        sample_columns = ["PATIENT_ID", "CANCER_TYPE", "ONCOTREE_CODE"]
        if extra_group_columns:
            sample_columns += [
                x for x in extra_group_columns if x in self.sample_info.data.columns
            ]
            patient_columns = [
                x for x in extra_group_columns if x in self.patient_info.data.columns
            ]
            muts = muts.join(self.sample_info.data[sample_columns]).reset_index()
            if patient_columns:
                muts = muts.merge(
                    self.patient_info.data[patient_columns],
                    on="PATIENT_ID",
                )
            muts = muts.set_index([main_index, "SAMPLE_ID"], drop=True)
        else:
            extra_group_columns = {}
            muts = muts.join(self.sample_info.data[sample_columns])
        # Aggregate "other" cancer types if requested
        if cancer_types_to_keep:
            muts["CANCER_TYPE"] = [
                x if x in cancer_types_to_keep else "other" for x in muts.CANCER_TYPE
            ]
        # Aggregate "other" cancer subtypes if requested
        if cancer_subtypes_to_keep:
            cancer_subtype_resolution = True
            muts["ONCOTREE_CODE"] = [
                x if x in cancer_subtypes_to_keep else "other"
                for x in muts.ONCOTREE_CODE
            ]
            # Keep order of subtypes as specified by cancer_subtypes_to_keep
            subtype_ranks = dict(
                zip(
                    cancer_subtypes_to_keep + ["other"],
                    range(len(cancer_subtypes_to_keep) + 1, 0, -1),
                )
            )
        else:
            # Sort subtypes by their frequency
            subtype_ranks = muts.value_counts("ONCOTREE_CODE").to_dict()
        # Assign "True" to the mutation status of a patient and cancer (sub)type
        # if at least one sample for this (sub)type is "True"
        group_columns = [main_index, "CANCER_TYPE", "PATIENT_ID"] + list(
            extra_group_columns.keys()
        )
        if cancer_subtype_resolution:
            group_columns.append("ONCOTREE_CODE")
        muts = muts.groupby(group_columns, dropna=False).agg({"mutated": "any"})

        # If we have cancer subtype resolution, keep only the sample with the
        # highest subtype rank for patients with multiple samples per cancer
        if cancer_subtype_resolution:
            muts["subtype_rank"] = muts.index.get_level_values("ONCOTREE_CODE").map(
                subtype_ranks
            )
            muts = muts.sort_values(["mutated", "subtype_rank"], ascending=False)
            muts = (
                muts.reset_index(level="ONCOTREE_CODE")
                .groupby(
                    [main_index, "CANCER_TYPE", "PATIENT_ID"]
                    + list(extra_group_columns.keys()),
                    dropna=False,
                )
                .first()
                .drop(columns="subtype_rank")
            )

        # Create a table of counts
        group_columns = [main_index, "CANCER_TYPE"] + list(extra_group_columns.keys())
        if cancer_subtype_resolution:
            group_columns.append("ONCOTREE_CODE")
        counts = (
            muts.reset_index()
            .groupby(group_columns, dropna=False)["mutated"]
            .value_counts(dropna=False)
            .to_frame()
            .reset_index()
            .pivot(index=group_columns, columns="mutated", values="count")
            .fillna(0)
            .astype(int)
            .rename(columns={True: "MUT", False: "WT"})
        )
        for column, values_to_keep in extra_group_columns.items():
            counts = self._merge_counts_for_others(counts, column, values_to_keep)
        counts = self._add_counts_for_all(counts)
        counts = self._add_frequencies_to_counts(counts, precision=precision)
        return counts

    def _format_confidence_interval(self, row: pd.Series, precision: int = 1) -> str:
        """Summarize frequency and confidence interval in a single string.

        Creates a string of the format `freq% (ci_low% ... ci_high%)` from
        the three values (frequency and confidence interval limits).

        Args:
            row: dataframe row with frequency and confidence interval limits.
            precision: number of decimal digits.

        Returns:
            string aggregating frequency and confidence interval limits.
        """
        width = 3 + precision  # two leading digits, the dot, precision digits
        f = f"{{:>{width}.{precision}f}}"
        s = (
            f"{f.format(row['AF_PERC'])}% ({f.format(row['AF_PERC_CI_LOWER'])}%"
            + f" ... {f.format(row['AF_PERC_CI_UPPER'])}%)"
        )
        return s

    def _add_frequencies_to_counts(
        self, counts: pd.DataFrame, precision: int = 1
    ) -> pd.DataFrame:
        """Add total number of patients and allele frequencies.

        This function expects a dataframe with the two columns "MUT" and "WT"
        which contain the number of patients who have a mutation or are of wilde
        type. The function adds these columns:

        * "N" -  the total number of patients that were tested for this
          mutation.
        * "AF_PERC" - allele frequency of mutation (percentage)
        * "AF_PERC_CI_LOWER" - lower limit of Clopper-Pearson confidence
          interval
        * "AF_PERC_CI_UPPER" - upper limit of Clopper-Pearson confidence
          interval
        * "AF_FORMATTED" - frequency formatted as "XX% (LO% ... HI%)"

        If `counts` already contains some of these columns, they will be
        replaced.

        Args:
            counts: dataframe with columns "MUT" and "WT" containing the respective
                number of patients
            precision: the number of fractional digits to be used for the
                AF_FORMATTED columns
        """
        # Remove any possibly existing frequency columns
        drop_cols = [
            x
            for x in [
                "N",
                "AF_PERC",
                "AF_PERC_CI_LOWER",
                "AF_PERC_CI_UPPER",
                "AF_FORMATTED",
            ]
            if x in counts.columns
        ]
        if drop_cols:
            counts = counts.drop(columns=drop_cols)
        # Calculate the allele frequencies and confidence intervals
        counts["N"] = counts["MUT"] + counts["WT"]
        counts["AF_PERC"] = counts["MUT"] / counts["N"] * 100.0
        counts = pd.concat(
            [
                counts,
                pd.DataFrame(ci(counts["MUT"], counts["N"], method="beta"))
                .transpose()
                .rename(columns={0: "AF_PERC_CI_LOWER", 1: "AF_PERC_CI_UPPER"})
                * 100.0,
            ],
            axis=1,
        )
        # Add a formatted allele frequency of the type "XX% (LO% ... HI%)"
        counts["AF_FORMATTED"] = counts.apply(
            self._format_confidence_interval, axis=1, precision=precision
        )
        return counts

    def _merge_counts_for_others(
        self, counts: pd.DataFrame, column: str, values_to_keep: list = None
    ) -> pd.DataFrame:
        """Create an "other" category for column and sum up counts.

        This function replaces all values in `column` with "other" unless the
        value is in `values_to_keep`. Subsequently, all counts for `other` are
        summed up. For example, `column` might be "PRIMARY_RACE", and
        `values_to_keep` might be `["White", "Black", "Asian"]`, then all other
        races, including "Unknown", are summarized in the "other" category.

        If `values_to_keep" is not provided or None, the counts dataframe is
        returned without change.

        Args:
            counts: dataframe with columns "WT" and "MUT".
            column: name of the column with values to be merged as "other".
            values_to_keep: values in `column` not to summarize in "other".

        Returns:
            dataframe with modied `column` and aggregated rows for "other"
        """
        if not values_to_keep:
            return counts
        old_index = counts.index.names
        counts = counts[["WT", "MUT"]].reset_index()
        counts[column] = [x if x in values_to_keep else "other" for x in counts[column]]
        counts = counts.groupby(old_index, dropna=False).sum()
        return counts

    def _add_counts_for_all(self, counts: pd.DataFrame) -> pd.DataFrame:
        """Add counts for the "all" columns of the frequency tables.

        For each of the factors, such as "CANCER_TYPE" or "PRIMARY_RACE", the
        counts data frame contains the counts for each value of those factors,
        unless they are summarized in "other", and for "other". This function
        here adds total counts "all", which includes not just "all" for each
        single factor, but also for all possible combinations of factors, with
        the most extreme case of total counts across all patients and all factor
        values.

        Args:
            counts: a count table without "all" values.

        Returns:
            input counts dataframe with added columns for "all" total sums.

        """
        # Drop all columns except wildtype and mutation counts
        counts = counts[["WT", "MUT"]].copy()
        # Initialize some variables
        levels = counts.index.names
        main_level = levels[0]  # hgvsg or gene_symbol
        levels_to_aggregate = [x for x in levels if x != main_level]
        counts_all = [counts]
        # Get all possible combinations of factors
        all_combinations = []
        for r in range(1, len(levels_to_aggregate) + 1):
            combinations = itertools.combinations(levels_to_aggregate, r)
            all_combinations.extend(combinations)
        # Get the count sums for all of these combinations
        for levels_in_sum in all_combinations:
            other_levels = [x for x in levels if x not in levels_in_sum]
            df = counts.groupby(other_levels, dropna=False).sum()
            for level in levels_in_sum:
                df[level] = "all"
            df = df.reset_index().set_index(levels, drop=True)
            counts_all.append(df)
        # Join all "all" sums with the original count table and sort the index
        counts = pd.concat(counts_all, axis=0).sort_index(
            key=lambda x: x.str.replace("other", "zzother").str.replace("all", "zzzall")
        )
        return counts

    def get_annotated_unique_mutations(
        self, gene_symbols: list = None, universe: str = "ensembl"
    ) -> pd.DataFrame:
        """Get annotated mutations for all or selected genes.

        This function returns a dataframe with all unique mutations found in
        at least one sample in Genie for the specified genes. If no genes are
        specified, return annotations for all genes.

        If the universe _mane_ is specified, more than one transcript can be
        returned for a genomic variant, for example a _MANE Select_ and a
        _MANE Plus clinical_ variant.

        Args:
            gene_symbols: return mutations for these genes. If not specified,
                mutations for all genes will be returned.
            universe: _ensembl_ or _mane_. Return annotations for Ensembl
                transcripts or for MANE transcripts.

        Returns:
            dataframe with annotated unique mutations.
        """
        self.mutations._check_universe(universe)
        annot = self.mutations.annot[universe]
        if gene_symbols is None:
            return annot
        else:
            return annot[annot.hgnc.isin(gene_symbols)]

    def get_imputed_sample_mutation_profiles(
        self,
        sample_mutation_profiles: pd.DataFrame,
        panel_coverage_threshold: float = 0.8,
        impute: bool = False,
    ) -> pd.DataFrame:
        """Create a sample vs mutation matrix with imputation.

        This function accepts sample mutation profiles as obtained by the
        function `get_sample_mutation_profiles`. The dataframe returned by
        `get_sample_mutation_profiles` contains only values (True or False for
        MUT or WT) for those sample-mutation combinations that were actually
        tested by the panel used for a sample. For gene level aggregation, a
        full sample-vs-mutation matrix without missing values is needed.

        If `impute` is set to _Yes_, then this function calculates such a matrix
        using MICE (Multiple Imputation by Chained Equations) for imputing
        values. MUT is encoded as 1, WT as 0, and the imputation will result in
        a fractional number between 0 and 1 for each missing value. Random
        numbers (uniform distribution between 0 and 1) are then used to assign
        “MUT” or “WT” depending on the value of the random variable and on the
        imputed value from MICE. To be more precise, the originally missing
        value gets a “MUT” if the random number is larger than the imputed value
        from MICE, and “WT” otherwise. This way we get a full
        mutation-versus-sample matrix without missing values.

        Imputation is a time consuming process. For example, imputation for KRAS
        mutations in NSCLC and CRC with the default panel coverage threshold
        takes about 15 minutes. If `impute` is _False_, all missing values are
        replaced with "WT". For a large `panel_coverage_threshold`, imputation
        changes frequencies only marginally. Therefore, imputation is switched
        off by default.

        Please be careful when calling this function - the memory requirements
        for an all genes versus all samples matrix would most likely exceed what
        is available in the compute environment. Therefore, always work with
        subsets of genes and maybe also indications.

        Args:
            sample_mutation_profiles: the mutation profiles as obtained by
                `get_sample_mutation_profiles`.
            panel_coverage_threshold: Exclude panels and samples profiled with
                these panels if the fraction of mutations tested by these panels
                is below that threshold. If this is set to 1, only panels that
                test all mutations are included. In that case, no imputation is
                needed.
            impute: Whether to impute missing values (mutations not tested). If
                this is set to False (the default), then missing values are
                replaced with wild type. Imputation takes very long (maybe
                hours), and for large `panel_coverage_threshold`s, frequencies
                don't change much, which is why imputation is switched off by
                default.

        Returns:
            Sample-mutation profiles with no missing values.
        """
        if panel_coverage_threshold < 0.0 or panel_coverage_threshold > 1.0:
            raise ValueError("panel_coverage_threshold should be between 0 and 1")
        # We do not need to impute if panel_coverage_threshold is 1, because
        # then we won't have any missing values
        if panel_coverage_threshold == 1:
            impute = False
        # The first step is to get the fraction of mutations tested by each
        # panel
        panel_counts = self.panel_set.tested_mutations.is_tested_matrix(
            hgvsgs=sample_mutation_profiles.index.get_level_values(
                "hgvsg"
            ).drop_duplicates()
        )
        panel_coverage = panel_counts.sum() / len(panel_counts)
        # Panels are included only if they test a large enough fraction of
        # mutations
        panel_ids_to_include = panel_coverage[  # noqa: F841
            panel_coverage >= panel_coverage_threshold
        ].index.to_list()
        # Get the samples that are profiled with panels to include
        sample_ids = sample_mutation_profiles.index.get_level_values(
            "SAMPLE_ID"
        ).drop_duplicates()
        samples_ids_to_include = (
            self.sample_info.data.loc[sample_ids, ["SEQ_ASSAY_ID"]]
            .query("SEQ_ASSAY_ID.isin(@panel_ids_to_include)")
            .index.to_list()
        )
        # Filter for samples to include and transform to a sample-vs-mutation
        # matrix
        smp_matrix = sample_mutation_profiles.loc[
            sample_mutation_profiles.index.get_level_values("SAMPLE_ID").isin(
                samples_ids_to_include
            )
        ].unstack(level="hgvsg")
        smp_matrix.columns = smp_matrix.columns.droplevel(0)
        # If we don't impute, fill all NAs with False and return
        if impute:
            imputer = IterativeImputer(random_state=0, min_value=0, max_value=1)
            smp_matrix_imputed = imputer.fit_transform(smp_matrix)
            # Replace imputed values with True or False based on random numbers
            smp_matrix_imputed = pd.DataFrame(smp_matrix_imputed)
            smp_matrix_imputed.index = smp_matrix.index
            smp_matrix_imputed.columns = smp_matrix.columns
            rng = np.random.default_rng(seed=42)
            smp_matrix_imputed = smp_matrix_imputed.map(
                lambda x: x >= rng.uniform(0, 1)
            )
        else:
            smp_matrix = smp_matrix.fillna(False)
        smp = smp_matrix.stack().swaplevel().to_frame()
        smp.columns = ["mutated"]
        return smp

    def aggregate_to_amino_acid_level(
        self,
        sample_mutation_profiles: pd.DataFrame,
        gene_symbols: list,
        universe: str = "ensembl",
    ) -> pd.DataFrame:
        """Aggregate mutations to gene level.

        This function aggregates the sample mutation level profiles to amino
        acid level. Several different nucleic acid changes can lead to the same
        amino acid change. These are aggregated by this function.

        The `gene_symbols` argument needs to be specified because a mutation can
        affect multiple genes, so the mapping from hgvsg to gene symbol is not
        unambigious.

        Args:
            sample_mutation_profiles: a single column dataframe with hgvsg and
                sample_id as row index and the boolean mutation status as
                column.
            gene_symbols: genes of interest - the hgvsgs from the the mutation
                profiles may map to multiple genes, and only the gene symbols
                from this list will be kept.
            universe: one of _ensembl_ and _mane_, defines which annotations
                will be used to create the mapping from hgvsg to gene symbols.

        Returns:
            dataframe with two-dimensional row index (hgvsp, sample_id)
                and single column with the mutation status
        """
        hgvsg_to_hgvsp = (
            self.get_annotated_unique_mutations(
                gene_symbols=gene_symbols, universe=universe
            )
            .loc[:, ["hgvsg", "hgvsp"]]
            .set_index("hgvsg", drop=True)
        )
        smp = (
            sample_mutation_profiles.reset_index()
            .merge(hgvsg_to_hgvsp, on="hgvsg")
            .groupby(["hgvsp", "SAMPLE_ID"], dropna=False)["mutated"]
            .sum()
            >= 1
        )
        smp = smp.to_frame()
        smp.columns = ["mutated"]
        return smp

    def aggregate_to_gene_level(
        self,
        sample_mutation_profiles: pd.DataFrame,
        gene_symbols: list,
        universe: str = "ensembl",
        min_mutations: int = 1,
    ) -> pd.DataFrame:
        """Aggregate mutations to gene level.

        This function aggregates the sample mutation level profiles to gene
        level. A sample is called mutated for a gene if it has at least
        `min_mutations` of that gene. For oncogenes, `min_mutations` should  be
        kept at the default value of 1. For tumor suppressor genes, it may be
        desired to require at least two hits to call a sample mutated, assuming
        that the two hits would affect the two copies of that gene. Any
        filtering for the functional relevance of mutations needs to be done
        prior to calling this function.

        The `gene_symbols` argument needs to be specified because a mutation can
        affect multiple genes, so the mapping from hgvsg to gene symbol is not
        unambigious.

        Args:
            sample_mutation_profiles: a single column dataframe with hgvsg and
                sample_id as row index and the boolean mutation status as
                column.
            gene_symbols: genes of interest - the hgvsgs from the the mutation
                profiles may map to multiple genes, and only the gene symbols
                from this list will be kept.
            universe: one of _ensembl_ and _mane_, defines which annotations
                will be used to create the mapping from hgvsg to gene symbols.
            min_mutations: Call a gene mutated if it has at least that many
                mutations. Should be 1 for oncogenes and may be set to 2 for
                tumor suppressor genes.

        Returns:
            dataframe with two-dimensional row index (gene, sample_id) and
                single column with the mutation status
        """
        self.mutations._check_universe(universe)
        hgvsg_to_gene = (
            self.get_annotated_unique_mutations(universe=universe)
            .loc[:, ["hgvsg", "hgnc"]]
            .rename(columns={"hgnc": "gene_symbol"})
            .set_index("hgvsg", drop=True)
        )
        smp = (
            sample_mutation_profiles.reset_index()
            .merge(hgvsg_to_gene, on="hgvsg")
            .query("gene_symbol.isin(@gene_symbols)")
            .groupby(["gene_symbol", "SAMPLE_ID"], dropna=False)["mutated"]
            .sum()
            >= min_mutations
        )
        smp = smp.to_frame()
        smp.columns = ["mutated"]
        return smp

    def get_tmb(
        self,
        sample_ids: list = None,
        min_genomic_range: int = 10000,
        tmb_intermediate_threshold: int = 6,
        tmb_high_threshold: int = 20,
    ) -> pd.DataFrame:
        """Get the tumor mutational burdon (TMB) for all specified samples.

        The TMB is calculated by dividing the number of mutations reported for a
        sample by the total genomic range covered by the panel that is used to
        profile that sample. If the total genomic range of a panel is smaller
        than `min_genomic_range`, NA is returned.

        Args:
            sample_ids: samples for which the TMB is to be returned. If None
                (the default), TMB is returned for all GENIE samples.
            min_genomic_range: if the total genomic range covered by a panel
                is smaller that this value, TMB is returned as NA for all
                samples tested with such a panel.
            tmb_intermediate_threshold: samples are classified as "TMB
                intermediate" if the TMB is >= this threshold but smaller than
                `tmb_high_threshold`. If it is smaller than this threshold, it
                is classified as "TMB low".
            tmb_high_threshold: samples are classified as "TMB high" if the TMB
                is >= this threshold.

        Returns:
            Table with the mutation count, size of genomic range covered, TMB,
                and TMB_class for each sample, with SAMPLE_ID as index.
        """
        assay_ranges = pd.DataFrame(
            {
                "SEQ_ASSAY_ID": list(self.panel_set.panels.keys()),
                "assay_genomic_range": [
                    panel.get_total_range() for panel in self.panel_set.panels.values()
                ],
            }
        )
        assay_ranges["assay_genomic_range"] = assay_ranges["assay_genomic_range"].apply(
            lambda x: pd.NA if x < min_genomic_range else x
        )
        if sample_ids is None:
            mut_data = self.mutations.data
        else:
            mut_data = self.mutations.data.query(
                "Tumor_Sample_Barcode.isin(@sample_ids)"
            )
        mut_counts = (
            mut_data.groupby("Tumor_Sample_Barcode", dropna=False)
            .apply(len)
            .to_frame()
            .rename(columns={0: "mut_count"})
            .rename_axis(index="SAMPLE_ID")
        )
        tmb = (
            mut_counts.merge(self.sample_info.data[["SEQ_ASSAY_ID"]], on="SAMPLE_ID")
            .reset_index()
            .merge(assay_ranges, on="SEQ_ASSAY_ID")
            .set_index("SAMPLE_ID", drop=True)
            .drop(columns="SEQ_ASSAY_ID")
        )
        tmb["TMB"] = tmb["mut_count"] / tmb["assay_genomic_range"] * 1e6
        tmb["TMB_class"] = tmb["TMB"].apply(
            lambda x: x
            if pd.isna(x)
            else "low"
            if x < tmb_intermediate_threshold
            else ("intermediate" if x < tmb_high_threshold else "high")
        )
        return tmb

    def append_patient_info(self, more_info: pd.DataFrame) -> None:
        """Add more columns to the patient information data.

        Args:
            more_info: table with additional columns to be added to the patient
                information. This dataframe must have the `PATIENT_ID` as index.

        Returns:
            nothing
        """
        self.patient_info.append_info(more_info)

    def append_sample_info(self, more_info: pd.DataFrame) -> None:
        """Add more columns to the sample information data.

        Args:
            more_info: table with additional columns to be added to the sample
                information. This dataframe must have the `SAMPLE_ID` as index.

        Returns:
            nothing
        """
        self.sample_info.append_info(more_info)
