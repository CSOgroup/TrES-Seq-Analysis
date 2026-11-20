"""Paths and cell IDs for processed datasets."""

import re


MEX_FILES = ("matrix.mtx", "features.tsv", "barcodes.tsv")


STUDY_DATASET_NAMES = {
    "TrES": "Karpas422_TrES",
    "MTRseq": "Karpas422_scMTRseq_1h",
    "MTRseq_OverNight": "Karpas422_scMTRseq_12h",
    "uCoTargetX": "Karpas422_uCoTargetX_1h",
    "uCoTargetX_OverNight": "Karpas422_uCoTargetX_12h",
    "TrES_K422_K9_DT": "Karpas422_TrES_DualTagmentation",
}


def input_paths(root, dataset, suffixes):
    release = STUDY_DATASET_NAMES.get(dataset, dataset)
    candidates = (
        (root, release),
        (root / release, release),
        (root / dataset, dataset),
    )
    for directory, prefix in candidates:
        paths = tuple(directory / f"{prefix}_{suffix}" for suffix in suffixes)
        if all(path.is_file() for path in paths):
            return paths
    return tuple(root / f"{release}_{suffix}" for suffix in suffixes)


def rna_paths(root, dataset, *, for_input=False):
    suffixes = tuple(f"RNA_{name}.gz" for name in MEX_FILES)
    if for_input:
        return input_paths(root, dataset, suffixes)
    return tuple(root / dataset / f"{dataset}_{suffix}" for suffix in suffixes)


def fragment_path(root, dataset, mark, *, for_input=False):
    suffix = f"{mark}.fragments.tsv.gz"
    if for_input:
        return input_paths(root, dataset, (suffix,))[0]
    return root / dataset / f"{dataset}_{suffix}"


def library_pattern(libraries):
    names = "|".join(re.escape(name) for name in sorted(libraries, key=len, reverse=True))
    return re.compile(rf"({names})_(\S+)")


def split_cell_id(value, pattern):
    match = pattern.fullmatch(str(value))
    if not match:
        raise ValueError(f"Unexpected cell ID {value!r}. Expected a configured library prefix.")
    return match.groups()
