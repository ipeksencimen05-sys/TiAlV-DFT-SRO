
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

from ase.io import iread
from ase.neighborlist import neighbor_list


# Project folders
ROOT = Path(__file__).resolve().parent.parent

SOURCE = ROOT / "raw" / "ti64_combined_data_Apr08.extxyz"

OUTPUT = ROOT / "processed" / "TiAlV_HW1_cleaned.csv"


# Calculate atomic descriptors
def calculate_descriptors(atoms):

    symbols = np.array(atoms.get_chemical_symbols())

    i_atoms, j_atoms, distances = neighbor_list(
        "ijd", atoms, cutoff=5.0
    )

    pairs = []
    first_distances = []
    shell_gaps = []

    for i in range(len(atoms)):

        mask = i_atoms == i

        neighbor_ids = j_atoms[mask]
        neighbor_distances = distances[mask]

        order = np.argsort(neighbor_distances)

        neighbor_ids = neighbor_ids[order]
        neighbor_distances = neighbor_distances[order]

        if len(neighbor_ids) < 13:
            raise ValueError("Not enough neighbors")

        first12 = neighbor_ids[:12]

        first_distances.extend(neighbor_distances[:12])

        shell_gaps.append(
            neighbor_distances[12] - neighbor_distances[11]
        )

        for j in first12:
            pairs.append((symbols[i], symbols[j]))

    def probability(center, neighbor):

        selected_pairs = [
            p for p in pairs if p[0] == center
        ]

        return (
            sum(p[1] == neighbor for p in selected_pairs)
            / len(selected_pairs)
        )

    c_Al = np.mean(symbols == "Al")
    c_V = np.mean(symbols == "V")

    return {
        "SRO_Ti_Al": 1 - probability("Ti", "Al") / c_Al,
        "SRO_Ti_V": 1 - probability("Ti", "V") / c_V,
        "SRO_Al_V": 1 - probability("Al", "V") / c_V,

        "mean_neighbor_distance": np.mean(first_distances),
        "std_neighbor_distance": np.std(first_distances),
        "min_shell_gap": np.min(shell_gaps)
    }


# Step 1: Find HCP Ti13Al2V configurations
rows = []

for i, atoms in enumerate(iread(SOURCE, index=":")):

    counts = Counter(atoms.get_chemical_symbols())
    config_type = atoms.info.get("config_type", "")

    if (
        config_type.startswith("hcp_")
        and len(atoms) == 16
        and counts.get("Ti", 0) == 13
        and counts.get("Al", 0) == 2
        and counts.get("V", 0) == 1
    ):

        rows.append({
            "structure_id": i,
            "config_type": config_type,
            "n_atoms": len(atoms),
            "n_Ti": 13,
            "n_Al": 2,
            "n_V": 1,
            "volume_per_atom": atoms.get_volume() / len(atoms),
            "DFT_energy_per_atom":
                atoms.info["QM_energy"] / len(atoms)
        })


df = pd.DataFrame(rows)

# Step 2: Select comparable volumes
df["volume_group"] = df["volume_per_atom"].round(3)

group_counts = (
    df.groupby(["config_type", "volume_group"])
    .size()
    .sort_values(ascending=False)
)

best_type, best_volume = group_counts.index[0]

df_group = df[
    (df["config_type"] == best_type) &
    (df["volume_group"] == best_volume)
]

# Step 3: Select 50 structures
df_sample = (
    df_group.sample(n=50, random_state=42)
    .sort_values("structure_id")
    .reset_index(drop=True)
)

selected_ids = set(df_sample["structure_id"])

# Step 4: Calculate descriptors
features = []

for i, atoms in enumerate(iread(SOURCE, index=":")):

    if i in selected_ids:

        result = calculate_descriptors(atoms)
        result["structure_id"] = i

        features.append(result)

    if len(features) == len(selected_ids):
        break

df_features = pd.DataFrame(features)

df_final = df_sample.merge(
    df_features,
    on="structure_id",
    validate="one_to_one"
)

# Step 5: Calculate relative energies
E_min = df_final["DFT_energy_per_atom"].min()

df_final["delta_E_meV_atom"] = (
    df_final["DFT_energy_per_atom"] - E_min
) * 1000

# Step 6: Validate and save
assert len(df_final) == 50
assert df_final.isnull().sum().sum() == 0

OUTPUT.parent.mkdir(exist_ok=True)

df_final.to_csv(OUTPUT, index=False)

print("Data preparation completed!")
print("Selected volume group:", best_volume)
print("Final dataset shape:", df_final.shape)
print("Saved to:", OUTPUT)
