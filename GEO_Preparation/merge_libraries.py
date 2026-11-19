#!/usr/bin/env python3
"""Merge RNA matrices and histone fragments into the study's GEO input layout."""

import argparse
import gzip
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory


MEX_FILES = ("matrix.mtx", "features.tsv", "barcodes.tsv")
HISTONE_MARKS = ("H3K27ac", "H3K27me3", "H3K9me3")


def rna_paths(root, dataset):
    return tuple(root / dataset / f"{dataset}_RNA_{name}.gz" for name in MEX_FILES)


def fragment_path(root, dataset, mark):
    return root / dataset / f"{dataset}_{mark}.fragments.tsv.gz"


def library_pattern(libraries):
    names = "|".join(re.escape(name) for name in sorted(libraries, key=len, reverse=True))
    return re.compile(rf"({names})_(\S+)")


def split_cell_id(value, pattern):
    match = pattern.fullmatch(str(value))
    if not match:
        raise ValueError(f"Unexpected cell ID {value!r}. Expected a configured library prefix.")
    return match.groups()


def open_text(path):
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding="utf-8")


def fragment_rows(path):
    """Read an unshifted, 0-based, five-column fragment file."""
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            try:
                chrom, start, end, barcode, count = line.rstrip("\r\n").split("\t")
                start, end, count = int(start), int(end), int(count)
                if (not chrom or not barcode or any(c.isspace() for c in barcode)
                        or not 0 <= start < end or count < 1):
                    raise ValueError
            except ValueError:
                raise ValueError(f"{path}:{line_no}: expected chrom, start, end, barcode, count") from None
            yield chrom, start, end, barcode, count


def prefixed_cell_id(library, barcode):
    if not barcode or any(c.isspace() for c in barcode):
        raise ValueError(f"{library}: invalid cell barcode {barcode!r}")
    prefix = f"{library}_"
    return prefix + barcode.removeprefix(prefix)


def mex_member(directory, *names):
    candidates = [directory / f"{name}{suffix}" for name in names for suffix in (".gz", "")]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"{directory}: missing {' or '.join(names)}")


def read_mex(directory, library):
    import numpy as np
    from scipy import sparse
    from scipy.io import mmread

    features = {}
    with open_text(mex_member(directory, "features.tsv", "genes.tsv")) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            feature_id = fields[0]
            if not feature_id or feature_id in features:
                raise ValueError(f"{directory}: missing or duplicate feature ID {feature_id!r}")
            symbol = fields[1] if len(fields) > 1 and fields[1] else feature_id
            kind = fields[2] if len(fields) > 2 and fields[2] else "Gene Expression"
            features[feature_id] = (symbol, kind)
    with open_text(mex_member(directory, "barcodes.tsv")) as handle:
        cells = [prefixed_cell_id(library, line.rstrip("\r\n")) for line in handle]
    if not cells or len(cells) != len(set(cells)):
        raise ValueError(f"{directory}: missing or duplicate RNA cell barcodes in {library}")
    with open_text(mex_member(directory, "matrix.mtx")) as handle:
        counts = sparse.csr_matrix(mmread(handle))
    if (counts.shape != (len(features), len(cells)) or not np.isfinite(counts.data).all()
            or np.any(counts.data < 0) or np.any(counts.data != np.rint(counts.data))):
        raise ValueError(f"{directory}: expected a nonnegative integer feature-by-cell matrix")
    return counts.astype(np.int64, copy=False), features, cells


def merge_rna(libraries, input_root, output_paths):
    from scipy import sparse
    from scipy.io import mmwrite

    matrices, cells = [], []
    reference = None
    for name, spec in libraries.items():
        if "rna" not in spec:
            continue
        counts, features, ids = read_mex(input_root / spec["rna"], name)
        if reference is None:
            reference = features
        elif features != reference:
            raise ValueError(f"{name}: RNA feature IDs or annotation differ within this dataset")
        if list(features) != list(reference):
            positions = {feature_id: i for i, feature_id in enumerate(features)}
            counts = counts[[positions[feature_id] for feature_id in reference], :]
        matrices.append(counts)
        cells.extend(ids)
    if not matrices:
        return
    if len(cells) != len(set(cells)):
        raise ValueError("Library prefixes produce duplicate merged RNA cell IDs")
    matrix_file, features_file, barcodes_file = output_paths
    with gzip.open(matrix_file, "wb") as handle:
        mmwrite(handle, sparse.hstack(matrices, format="coo"), field="integer")
    with gzip.open(features_file, "wt", encoding="utf-8") as handle:
        handle.writelines(f"{key}\t{symbol}\t{kind}\n" for key, (symbol, kind) in reference.items())
    with gzip.open(barcodes_file, "wt", encoding="utf-8") as handle:
        handle.writelines(f"{cell}\n" for cell in cells)


def append_fragments(source, library, handle, work, barcode_tag):
    if source.suffix == ".bam":
        import snapatac2 as snap

        fragments = work / "from_bam.fragments.tsv.gz"
        snap.pp.make_fragment_file(str(source), str(fragments), barcode_tag=barcode_tag,
                                  shift_left=0, shift_right=0, min_mapq=20, tempdir=str(work))
    else:
        fragments = source
    seen = False
    for chrom, start, end, barcode, count in fragment_rows(fragments):
        cell = prefixed_cell_id(library, barcode)
        handle.write(f"{chrom}\t{start}\t{end}\t{cell}\t{count}\n")
        seen = True
    if not seen:
        raise ValueError(f"{source}: no fragments found")
    if fragments != source:
        fragments.unlink()


def merge_dataset(dataset, libraries, input_root, output_root):
    for name in (dataset, *libraries):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
            raise ValueError(f"Invalid dataset or library name: {name!r}")
    if not libraries:
        raise ValueError(f"{dataset}: at least one library is required")
    # A prefix must identify exactly one library, including for DNA-only cells.
    if any(a != b and b.startswith(a + "_") for a in libraries for b in libraries):
        raise ValueError(f"{dataset}: one library name is a prefix of another")
    marks = {mark for spec in libraries.values() for mark in spec.get("histone", {})}
    if marks - set(HISTONE_MARKS):
        raise ValueError(f"{dataset}: unsupported histone marks {sorted(marks - set(HISTONE_MARKS))}")
    if not marks and not any("rna" in spec for spec in libraries.values()):
        raise ValueError(f"{dataset}: no RNA or histone files were specified")
    destination = output_root / dataset
    if destination.exists():
        raise FileExistsError(f"Choose a new dataset output directory: {destination}")
    with TemporaryDirectory(prefix=f".{dataset}-", dir=output_root) as scratch:
        work = Path(scratch)
        (work / dataset).mkdir()
        merge_rna(libraries, input_root, rna_paths(work, dataset))
        for mark in sorted(marks):
            with gzip.open(fragment_path(work, dataset, mark), "wt", encoding="utf-8") as handle:
                for library, spec in libraries.items():
                    if mark in spec.get("histone", {}):
                        append_fragments(input_root / spec["histone"][mark], library, handle, work,
                                         spec.get("barcode_tag", "CB"))
        (work / dataset).rename(destination)
    print(f"{dataset}: wrote {destination}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="JSON mapping datasets to their source libraries")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--dataset", action="append", help="Merge only the named dataset, repeatable")
    args = parser.parse_args()
    manifest = args.manifest.expanduser().resolve()
    groups = json.loads(manifest.read_text(encoding="utf-8"))
    chosen = args.dataset or list(groups)
    if not chosen or set(chosen) - groups.keys():
        parser.error("Select dataset names listed in the manifest.")
    output = args.outdir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    for name in chosen:
        merge_dataset(name, groups[name], manifest.parent, output)


if __name__ == "__main__":
    main()
