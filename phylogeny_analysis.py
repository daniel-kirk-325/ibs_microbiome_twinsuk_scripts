from Bio import Phylo # 1.87
import pandas as pd
import os
import re

#### CODE EXPLANATION ####

## Obtains the closest known relative of unknown IBS-associated bacterial 
## species.

## The code was run as part of a bash script (phylogeny_bash.sh) in the CREATE HPC 
## environment. The results from each species were then loaded into R and the nearest 
## known relative was identified from each results dataframe

#### #### #### #### #### #

## Import array from bash
I = int(os.environ["SLURM_ARRAY_TASK_ID"])
print(I)

# Load tree
tree = Phylo.read("mpa_vJan21_CHOCOPhlAnSGB_202103.nwk", "newick")

# List unknown SGBs
unknown_sgbs = [
    "GGB3571 SGB4778", "GGB41458 SGB58520", "GGB9522 SGB14921",
    "GGB9559 SGB14969", "GGB9602 SGB15031", "GGB9614 SGB15049",
    "GGB9635 SGB15106", "GGB9699 SGB15216", "GGB9705 SGB15224",
    "GGB9705 SGB15225", "GGB9707 SGB15229", "GGB9712 SGB15244",
    "GGB9787 SGB15410", "Ruminococcaceae unclassified SGB15234",
    "Ruminococcaceae unclassified SGB15236", 'SGB15346'
]

sgb_numbers = [item.split("SGB")[-1] for item in unknown_sgbs if "SGB" in item]

species_name = sgb_numbers[I]

# Try exact match on full clade name first
target_clades = list(tree.find_clades({"name": species_name}))

if not target_clades:
    # Exact match on the SGB number as a whole token
    pattern = re.compile(rf"(?<!\d)SGB{re.escape(species_name)}(?!\d)")
    target_clades = [
        clade for clade in tree.find_clades()
        if clade.name and pattern.search(clade.name)
    ]

if not target_clades:
    print(f"Species '{species_name}' not found in the tree.")

target = target_clades[0]

distances = []
for clade in tree.get_terminals():
    if clade.name == target.name or not clade.name:
        continue
    dist = tree.distance(target, clade)
    distances.append({"SGB": clade.name, "distance": dist})

distances.sort(key=lambda x: x["distance"])
df = pd.DataFrame(distances[:150])


# Load mapping file and add names to SGB numbers in results
txt_file = "mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt"
mapping = {}
with open(txt_file, encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Handle common formats flexibly
        if '\t' in line:
            key, name = line.split('\t', 1)
        elif ' ' in line:
            parts = line.split(maxsplit=1)
            key = parts[0]
            name = parts[1] if len(parts) > 1 else key
        else:
            key = line
            name = line

        # Clean key (SGB12345 or just the number)
        key = key.strip()
        if key.startswith("SGB"):
            mapping[key] = name.strip()
        else:
            mapping[f"SGB{key}"] = name.strip()   # in case file has only numbers
            mapping[key] = name.strip()


# ====================== MAP NAMES IN RESULTS ======================

df['full_name'] = ('SGB' + df['SGB'].astype(str)).map(mapping)

# Account for blank entries (those with "_group" suffix)
df.loc[df['full_name'].isna(), 'full_name'] = ('SGB' + df.SGB[df['full_name'].isna()] + '_group').map(mapping)


# Reorder columns nicely
cols = ["SGB", 'full_name', "distance"]
df = df[cols]

df.to_csv(f'SGB{sgb_numbers[I]}.csv')