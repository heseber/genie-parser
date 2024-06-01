import os

import pandas as pd
from tqdm.auto import tqdm

import genie.genie as gd

# Where to find the GENIE data
genie_dir = os.path.join(
    os.environ.get("HOME"),
    "genie/Data/Original/consortium/Main_GENIE_cBioPortal_Releases/16.2-consortium/",
)
g = gd.Genie(genie_dir)

print("Loading re-annotations", end=" ... ", flush=True)
positions = (
    pd.read_table(
        g.config.get_aux_file_name("annot"),
        low_memory=False,
        usecols=["hgvsg", "chromosome", "begin", "end"],
    )
    .drop_duplicates()
    .set_index("hgvsg", drop=True)
)
print("done.", flush=True)

print("Checking mutations versus panels.", flush=True)
cols = []
panels = g.panel_set.panels
for panel_id, panel in tqdm(panels.items(), total=len(panels)):
    cols.append(
        pd.Series(
            positions.apply(
                lambda x: panel.range_is_tested(x["chromosome"], x["begin"], x["end"]),
                axis=1,
            ),
            name=panel_id,
        )
    )
tested = pd.concat(cols, axis=1)

print("Writing parquet file", end=" ... ", flush=True)
file_name = g.config.get_aux_file_name("panel_mut_tested")
tested.to_parquet(file_name)
print("done.", flush=True)
