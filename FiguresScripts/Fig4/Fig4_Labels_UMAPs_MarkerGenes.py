import matplotlib.pyplot as plt
import snapatac2 as snap
import anndata as ad
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
from matplotlib_venn import venn3
from upsetplot import UpSet, plot, from_indicators
import os
import csv

if not os.path.exists('figures'):
    os.makedirs('figures')
    
workdir = "/mnt/dataFast/ahrmad/triseq_2025052_2/processed/"

#################################################
#################### RNA QC IMPORT #####################
#################################################

print(f'IMPORT AND QC PLOTS RNA')

sample='Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil'

#gzip everything
mtxPath = "/mnt/dataFast/ahrmad/triseq_2025052_2/processed/RNA_ONLY/STARsoloFILT/"
#mtxPath = "/mnt/dataFast/ahrmad/triseq_2025052_2/processed/RNA_ONLY/STARsoloRAW/"
adata = sc.read_10x_mtx(mtxPath)

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
    adata,
    ["pct_counts_mt", 'pct_counts_in_top_50_genes'],
    jitter=0.4,
    multi_panel=True
    ,rotation=0.0000001,save=f'{sample}_Pct_PreFilter.pdf'
)

# Create the violin plot
sc.pl.violin(
    adata,
    ["log1p_n_genes_by_counts", "log1p_total_counts"],
    jitter=0.4,
    multi_panel=True,
    rotation=0.0000001,
    show=False  # Prevent auto-display so we can edit it
)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]  # reverse log1p
    ax.set_yticklabels(new_labels)

# Save or show the updated figure
plt.savefig(f"./figures/{sample}_PreFilter.pdf")
plt.close()

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",title=sample+'\nColor: pct_counts_mt',save=f'{sample}_PreFilter.pdf')
sc.pl.scatter(adata, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",title=sample+'\nColor: pct_counts_mt',save=f'{sample}_PreFilterlog.pdf')

adata.write(f'{sample}.h5ad')

#################################################
############# RNA FILTERING ###############
#################################################

exp ='Sc_VTHumanTonsil'

exp_rna='Sc_VTr10_S2_Sc_VTr11_S3_VTHumanTonsil'

#RNA
rna = sc.read_h5ad(f'{exp_rna}.h5ad')

#RNA
minUMI = 500
maxUMI = 8000
max_pct_mt = 2.5
max_pct_50g = 25


# --- RNA filtering criteria ---
filtered_rna = rna[rna.obs['total_counts'] >= minUMI]
filtered_rna = filtered_rna[filtered_rna.obs['total_counts'] <= maxUMI]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_mt'] < max_pct_mt]
filtered_rna = filtered_rna[filtered_rna.obs['pct_counts_in_top_50_genes'] < max_pct_50g]

sc.pl.violin(
    filtered_rna,
    ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt", 'pct_counts_in_top_50_genes'],
    jitter=0.4,
    multi_panel=True
    ,rotation=0.0000001,save=f'{exp}_rna_PostFilter.png'
)

# Create the violin plot
sc.pl.violin(
    filtered_rna,
    ["log1p_n_genes_by_counts", "log1p_total_counts"],
    jitter=0.4,
    multi_panel=True,
    rotation=0.0000001,
    show=False  # Prevent auto-display so we can edit it
)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]  # reverse log1p
    ax.set_yticklabels(new_labels)

# Save or show the updated figure
plt.suptitle(f'Number of cells: {filtered_rna.n_obs}\nMedian UMIs per cell: {np.median(filtered_rna.obs["total_counts"])}')
plt.savefig(f"./figures/{exp}_rna_PostFilterNoLog.pdf")
plt.close()


sc.pl.scatter(filtered_rna, "total_counts", "n_genes_by_counts", color="pct_counts_mt",title=exp+'\nColor: pct_counts_mt',save=f'{exp}_rna_PostFilter.png')
sc.pl.scatter(filtered_rna, "log1p_total_counts", "log1p_n_genes_by_counts", color="pct_counts_mt",title=exp+'\nColor: pct_counts_mt',save=f'{exp}_rna_PostFilterlog.png')


filtered_rna.write(f'{exp}_rna_CountFiltered.h5ad')

#################################################
############# RNA ANALYSIS ###############
#################################################

cell_lines = [
    "Sc_VTHumanTonsil",
]

# List to hold AnnData objects for each sample
adatas = []
sample_names = []

#For R code
import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
rcb.logger.setLevel(logging.ERROR)
#ro.pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython



##Doublet Analysis
%R library(Seurat)
%R library(scater)
%R library(scDblFinder)
%R library(BiocParallel)
%R library(scran)
%R library(scry)

import numpy as np
import scipy.sparse as sp

for sample in cell_lines:
    exp = sample
    try:
        adata_sample = sc.read_h5ad(f"{exp}_rna_CountFiltered.h5ad")
        adata_sample.obs['sample'] = exp  # Add sample name as a column for later coloring

        # genes x cells (scDblFinder expects features x cells)
        data_mat = adata_sample.X.T

        # rpy2 cannot handle SciPy sparse matrices – convert to dense
        if sp.issparse(data_mat):
            data_mat = data_mat.toarray()   # or data_mat.A

        # optionally ensure numeric type
        data_mat = data_mat.astype(float)

        # send to R
        %R -i data_mat
        %R set.seed(1)
        %R sce = scDblFinder(SingleCellExperiment(list(counts=data_mat)),artificialDoublets=4000,dbr=0.2, cluster=TRUE)
        %R doublet_scoreR <- sce$scDblFinder.score
        %R doublet_classR <- sce$scDblFinder.class

        adata_sample.obs["scDblFinder_score"] = %Rget doublet_scoreR
        adata_sample.obs["scDblFinder_class"] = %Rget doublet_classR

        adatas.append(adata_sample)
        sample_names.append(exp)
    except Exception as e:
        print(f"Could not load {exp}_rna_CountFiltered.h5ad: {e}")


sc.pl.scatter(adatas[0], "scDblFinder_score", "log1p_total_counts", color="scDblFinder_class",
              title='Color: scDblFinder_class', save=f'{exp}_rna_scDBL.png')


for ad in adatas:
    score = ad.obs['scDblFinder_score']
    ad.obs['scDblFinder_class_manual'] = np.where(score < 0.5, 'singlet', 'doublet')

adata = adatas[0][adatas[0].obs['scDblFinder_score'] < 0.6].copy()
#adata = adatas[0][adatas[0].obs['scDblFinder_class']=="singlet"].copy()
#adata = adatas[0].copy()

# Store raw counts in a layer
adata.layers["counts"] = adata.X.copy()

### Shifted logarithm Normalization
print('Shifted logarithm Normalization')
#scales_counts = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)
adata.layers["log1p_norm"] = sc.pp.log1p(adata.X, copy=True)


adata.write(f"{exp}_rna_raw.h5ad")

#SCT - R - Deviance

adata = sc.read_h5ad(f"{exp}_rna_raw.h5ad")

# Normalize and log1p
### SCTRANSFORM Normalization
print('SCTransform Normalization')
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
anndata2ri.activate()
#pandas2ri.activate()

sc.pp.filter_genes(adata,min_cells=5,inplace=True)
mat = adata.X.copy()

# Set names for the input matrix
cell_names = adata.obs_names
gene_names = adata.var_names
r.assign('mat', mat.T)
r.assign('cell_names', cell_names)
r.assign('gene_names', gene_names)
r('colnames(mat) <- cell_names')
r('rownames(mat) <- gene_names')

seurat = importr('Seurat')
r('seurat_obj <- CreateSeuratObject(mat)')

# Run
r(f'seurat_obj <- SCTransform(seurat_obj,vst.flavor="v2")')

# Extract the SCT data and add it as a new layer in the original anndata object
sct_data = np.asarray(r['as.matrix'](r('seurat_obj@assays$SCT@data')))
adata.layers['SCT_data'] = sct_data.T
sct_data = np.asarray(r['as.matrix'](r('seurat_obj@assays$SCT@counts')))
adata.layers['SCT_counts'] = sct_data.T

# choose one of the layers
#adata.X = adata.layers['SCT_data']
adata.X = adata.layers["log1p_norm"]


# Cell cycle scoring
cell_cycle_genes = [x.strip() for x in open('/mnt/dataFast/ahrmad/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

print('Deviance')
xTopDeviantGenes = 2000
import rpy2.robjects as ro
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri, pandas2ri
import anndata2ri

# No activate() calls
# remove common plotting-only offenders
to_drop = [k for k in adata.uns if k.endswith("_colors") or k.endswith("_color")]
for k in to_drop:
    adata.uns.pop(k, None)

# also drop colors saved under per-category keys, e.g. adata.uns['louvain_colors'], etc.
for k, v in list(adata.uns.items()):
    if isinstance(v, dict):
        for kk in list(v):
            if kk.endswith("_colors") or kk.endswith("_color"):
                del v[kk]

with localconverter(
    default_converter
    + numpy2ri.converter
    + pandas2ri.converter
    + anndata2ri.converter
):
    ro.globalenv["adata"] = adata


# Feature selection 
%R sce = devianceFeatureSelection(adata, assay="X")

binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T
idx = binomial_deviance.argsort()[-xTopDeviantGenes:]
mask = np.zeros(adata.var_names.shape, dtype=bool)
mask[idx] = True

adata.var["highly_deviant"] = mask
adata.var["binomial_deviance"] = binomial_deviance


# Highly variable genes (across all cells)
sc.pp.highly_variable_genes(adata,flavor='seurat_v3')

sc.pl.highly_variable_genes(adata, show=False, save=f'{exp}_HVG_ScanpySeuratV3_SCT.pdf')
adata.var["highly_variableScanpySeuratV3"] = adata.var["highly_variable"].copy()
adata.var["highly_variable"] = adata.var["highly_deviant"].copy()
sc.pl.highly_variable_genes(adata, show=False, save=f'{exp}_HVG_R_SCT.pdf')
adata.var["highly_variable"] = adata.var["highly_variableScanpySeuratV3"].copy()

# PCA
sc.tl.pca(adata,n_comps=20)
sc.pl.pca_variance_ratio(adata, n_pcs=20, log=False, show=False, save=f'{exp}_pca_variance_ratio.pdf')

# PCA plots
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase"],
    dimensions=[(0, 1), (0, 1), (0, 1), (0, 1)],
    ncols=3,
    show=False, save=f'{exp}_PCA1-2_QC.pdf'
)
sc.pl.pca(
    adata,
    color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase"],
    dimensions=[(2, 3), (2, 3), (2, 3), (2, 3)],
    ncols=3,
    show=False, save=f'{exp}_PCA3-4_QC.pdf'
)

# Neighbors/UMAP
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.umap(adata)

# UMAPs colored by QC and sample
sc.pl.umap(
    adata,
    color=["pct_counts_mt", "pct_counts_hb", "pct_counts_ribo", "pct_counts_in_top_50_genes", "log1p_total_counts"],
    ncols=3,
    show=False, save=f'{exp}_UMAP_QC.pdf'
)
sc.pl.umap(
    adata,
    color=["log1p_total_counts", "log1p_n_genes_by_counts","scDblFinder_score","scDblFinder_class"],
    ncols=3,
    show=False, save=f'{exp}_UMAP_QC_Log.pdf'
)
sc.pl.umap(
    adata,
    color=['S_score', 'G2M_score', 'phase', ],
    ncols=2,
    show=False, save=f'{exp}_UMAP_QC_CellCycle.pdf'
)

leiR = 0.15
leiRes = f'leiden_res{leiR}'
sc.tl.leiden(adata, key_added=leiRes, resolution=leiR)
sc.pl.umap(adata,color=leiRes,legend_loc="on data",save=f"{exp_rna}_{leiRes}",title=f"{exp_rna} {leiRes}",show=False)

adata.write('Tonsil_RNAOnly.h5ad')


#################################################
############# LABEL TRANSFER ###############
#################################################
#scvi ENV

import os
import scanpy as sc
import scvi
import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === Settings ===
#REF_PATH = "/mnt/svraw1/uporicchio-raw/divyanshu/server_backup_29112024/mnt_data1/exvivo/exvivo_home_backup/exvivo/visium/data/sc.h5ad"
REF_PATH = "sc.h5ad"

QUERY_PATH = "Tonsil_RNAOnly.h5ad"
LABEL_KEY = "new_celltype"  # column in adata_ref.obs
BATCH_KEY = "dataset"    # column to mark batch (you can rename it)
SAVE_DIR = "figures"
os.makedirs(SAVE_DIR, exist_ok=True)

# === Load Data ===
adata_ref = sc.read(REF_PATH)
adata_query = sc.read(QUERY_PATH)

#adata_ref.obs['total_counts'] = np.array(adata_ref.X.sum(axis=1)).flatten()
#sc.pl.violin(adata_ref, 'total_counts', show=False, save='_ref.png')
adata_query.X = adata_query.layers["counts"].copy()

sc.pl.umap(adata_ref, color=["new_celltype"],ncols=1, show=False,save=f"FullRef_CellType_OriUMAP.pdf")
#sc.pl.umap(adata_ref[adata_ref.obs["Method"]=="5GEX"], color=["new_celltype"],ncols=1, show=False,save=f"5GEXRef_CellType_OriUMAP.pdf")

#adata_ref = adata_ref[adata_ref.obs["Method"]=="5GEX"].copy()


import pandas as pd

# --- 0. Extract gene names ---
ref_genes = pd.Index(adata_ref.var_names)
query_genes = pd.Index(adata_query.var_names)

# === Ensure common genes ===
shared_genes = adata_ref.var_names.intersection(adata_query.var_names)
adata_ref = adata_ref[:, shared_genes].copy()
adata_query = adata_query[:, shared_genes].copy()

# === Add batch & tech labels ===
adata_ref.obs[BATCH_KEY] = "ref"
adata_query.obs[BATCH_KEY] = "query"

# === Add labels for scanvi ===
SCANVI_LABEL_KEY = "labels"
adata_ref.obs[SCANVI_LABEL_KEY] = adata_ref.obs[LABEL_KEY]
adata_query.obs[SCANVI_LABEL_KEY] = "Unknown"

# === Concatenate ===
adata = anndata.concat([adata_ref, adata_query], label="source", keys=["ref", "query"])

# === Store raw counts ===
adata.layers["counts"] = adata.X.copy()

# === HVG selection ===
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key=BATCH_KEY,
    subset=True,
)

# === Train scVI ===
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH_KEY)
scvi_model = scvi.model.SCVI(adata)
scvi_model.train()

# === Save UMAP (scVI latent) ===
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.3)

sc.pl.umap(adata, color=[BATCH_KEY], show=False)
plt.savefig(f"{SAVE_DIR}/umap_scvi_batch.png", bbox_inches="tight")
plt.close()

sc.pl.umap(adata, color=[SCANVI_LABEL_KEY], show=False)
plt.savefig(f"{SAVE_DIR}/umap_scvi_labels.png", bbox_inches="tight")
plt.close()

# Save the AnnData object
adata.write("my_adata.h5ad")

# Save the SCVI model
scvi_model.save("my_scvi_model/", overwrite=True)


import scvi
import scanpy as sc

# Load AnnData
adata = sc.read("my_adata.h5ad")

# Load the model
scvi_model = scvi.model.SCVI.load("my_scvi_model/", adata=adata)

# === Train scANVI for label transfer ===
scvi.model.SCANVI.setup_anndata(adata, labels_key=SCANVI_LABEL_KEY, unlabeled_category="Unknown")
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown", labels_key=SCANVI_LABEL_KEY)
scanvi_model.train(max_epochs=20, n_samples_per_label=100)

# === Predict labels and get latent space ===
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)

# Get label probabilities for each cell
probs = scanvi_model.predict(adata, soft=True)

# Add max confidence score and most likely label
adata.obs["predicted_labels"] = probs.idxmax(axis=1)   # already done if you used scanvi.predict()
adata.obs["label_confidence"] = probs.max(axis=1)
probs.to_csv("label_probabilities.csv")

# === UMAP visualization after scANVI ===
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)

sc.pl.umap(adata, color=["predicted_labels","label_confidence","dataset"],ncols=1, show=False,save=f"_scanvi_predicted_labels.pdf")

# === Save query with predicted labels ===
adata_query.obs["predicted_labels"] = adata[adata.obs["source"] == "query"].obs["predicted_labels"].values
adata_query.obs["label_confidence"] = adata[adata.obs["source"] == "query"].obs["label_confidence"].values
sc.pl.umap(adata_query, color=["predicted_labels","label_confidence","scDblFinder_score"],ncols=1, show=False,save=f"Annot_Filt_LT1.pdf")

# 1. Compute and sort cell counts
label_counts = adata_query.obs["predicted_labels"].value_counts()
label_counts = label_counts.sort_values(ascending=False)

# 2. Convert to DataFrame
label_counts_df = pd.DataFrame({
    "cell_type": label_counts.index,
    "cell_count": label_counts.values
})

# 3. Plot sorted counts
plt.figure(figsize=(10, 5))
plt.bar(label_counts_df["cell_type"], label_counts_df["cell_count"])
plt.xticks(rotation=90)
plt.ylabel("Number of cells")
plt.title("Cell counts per predicted label (sorted)")
plt.tight_layout()
plt.savefig("figures/query_celltype_distribution_sorted1.png", bbox_inches="tight")
plt.close()

adata_query.write("TonsilRNA_Annot_Filt_LT.h5ad")




#################################################
############# RNA Downstream ###############
#################################################



adata_query = sc.read_h5ad("./TonsilRNA_Annot_Filt_LT.h5ad")

sc.pl.umap(adata_query, color=["predicted_labels"],ncols=1, show=False,save=f"LabelTransfer_FULLRNA.pdf")

adata = adata_query.copy()

# Quick summary stats
print(adata.obs[['label_confidence', 'scDblFinder_score']].describe())

# Optional: visualize
sc.pl.violin(
    adata,
    ['label_confidence', 'scDblFinder_score'],
    groupby=None,
    jitter=0.4,
    multi_panel=True,
    show=False,save=f"LC_Dblt_FULLRNA.pdf"
)

# --- 2. Choose data-driven thresholds using quantiles ---
low_conf_q = 0.20   # bottom 20% of label_confidence
high_dbl_q = 0.80   # top 20% of scDblFinder_score

low_conf_thr = adata.obs['label_confidence'].quantile(low_conf_q)
high_dbl_thr = adata.obs['scDblFinder_score'].quantile(high_dbl_q)

print(f"Low label_confidence threshold (q={low_conf_q}): {low_conf_thr:.3f}")
print(f"High scDblFinder_score threshold (q={high_dbl_q}): {high_dbl_thr:.3f}")

# --- 3. Define cells to remove ---
to_remove_mask = (
    (adata.obs['label_confidence'] <= low_conf_thr) &
    (adata.obs['scDblFinder_score'] >= high_dbl_thr)
)

adata.obs['low_conf_high_doublet'] = to_remove_mask

sc.pl.umap(
    adata,
    color="low_conf_high_doublet",
    save=f"_low_conf_high_doublet_ALL.pdf",
)
print(f"Cells to remove: {to_remove_mask.sum()} / {adata.n_obs}")

# --- 4. Filter them out ---
adata_query = adata[~to_remove_mask].copy()

print(f"Remaining cells after filtering: {adata_query.n_obs}")


sc.pp.highly_variable_genes(
    adata_query,
    layer='counts',
    flavor='seurat_v3',
    n_top_genes=3000,
)
sc.pp.scale(adata_query,zero_center=False, max_value=10) # optional clipping
sc.pp.pca(adata_query,n_comps=15)
sc.pl.pca_variance_ratio(adata_query, n_pcs=15, log=False, show=False, save=f'Tonsils_RNA_PCAVar.pdf')


sc.pp.neighbors(adata_query, use_rep="X_pca", n_neighbors=30,metric='cosine')  # or adjust as needed
sc.tl.umap(adata_query)

sc.pl.umap(
    adata_query,
    color=["predicted_labels"],
    ncols=1,
    show=False,
    save="LabelTransferFiltered_FULLRNA.pdf"
)




LT_res=20
sc.tl.leiden(adata_query, resolution=LT_res, key_added=f'overcluster_{LT_res}')

# ----------------------- Params -----------------------
cluster_key = f'overcluster_{LT_res}'
score_key   = "scDblFinder_score"
score_key   = "label_confidence"
median_key  = f"{cluster_key}_medianLC"   # kept name for compatibility; now holds MEDIAN
mad_key     = f"{cluster_key}_madLC"    # new column
threshold   = 0.90

# ------------------ 1) Cluster stats ------------------
cluster_stats = (
    adata_query.obs
    .groupby(cluster_key, observed=True)
    .agg(
        median_label_confidence=(score_key, "median"),
        n_cells=(score_key, "size"),
    )
    .sort_values("median_label_confidence", ascending=False)
)

print("=== Per-cluster median label_confidence and sizes ===")
print(cluster_stats)

# ------------------ 2) Broadcast median ----------------
# Map the cluster median back to each cell and ensure float (not categorical)
median_map = cluster_stats["median_label_confidence"]
adata_query.obs[median_key] = adata_query.obs[cluster_key].map(median_map).astype(float)

# Safety: if something coerced it to categorical, force back to float
if pd.api.types.is_categorical_dtype(adata_query.obs[median_key]):
    adata_query.obs[median_key] = adata_query.obs[median_key].astype(float)

# Quick peek (unique cluster -> median pairs)
print("\n=== Unique cluster -> median (first few) ===")
print(
    adata_query.obs[[cluster_key, median_key]]
    .drop_duplicates()
    .sort_values(median_key)
    .head(20)
)

# Replace df and "column_name" with your DataFrame and column
plt.figure()
adata_query.obs[median_key].hist(bins=50)    # or df['column_name'].plot(kind='hist')
plt.title("Histogram of column_name")
plt.xlabel("column_name")
plt.ylabel("Frequency")

plt.savefig("figures/histogram_column_name.png", dpi=300, bbox_inches='tight')
plt.close()

# ------------------ 3) Identify removals ----------------
removed = cluster_stats[cluster_stats["median_label_confidence"] < threshold]
kept    = cluster_stats[cluster_stats["median_label_confidence"] >= threshold]

# Only print how many clusters were removed and how many cells that is
cells_removed = int(removed["n_cells"].sum())
cells_kept = int(kept["n_cells"].sum())
print(f"Removed clusters: {len(removed)}; cells removed: {cells_removed}")
print(f"Kept clusters: {len(kept)}; cells : {cells_kept}")

# ------------------ 4) Subset & plot -------------------
# All cells: UMAP colored by per-cluster median
sc.pl.umap(
    adata_query,
    color=median_key,
    cmap="viridis",
    vmin=0, vmax=1,         # adjust if your confidence isn't 0–1
    save=f"_{median_key}_ALL.pdf",
)

# Filter to kept clusters
clusters_to_keep = kept.index
adata_query = adata_query[adata_query.obs[cluster_key].isin(clusters_to_keep)].copy()

# UMAP of kept clusters only
sc.pl.umap(
    adata_query,
    color=median_key,
    cmap="viridis",
    vmin=0, vmax=1,
    save=f"_{median_key}_gt{threshold}.pdf",
)

import seaborn as sns
# Count cells per predicted_labels
# Prepare data
label_counts = (
    adata_query.obs["predicted_labels"]
    .value_counts()
    .sort_values(ascending=False)
    .reset_index()
)
label_counts.columns = ["predicted_labels", "cell_count"]
# Plot
plt.figure(figsize=(10, 6))
sns.barplot(
    data=label_counts,
    x="predicted_labels",
    y="cell_count",
    palette="tab20",
    order=label_counts["predicted_labels"],  # enforce sorted order
)

# Add count labels on top
for i, v in enumerate(label_counts["cell_count"]):
    plt.text(i, v + (v * 0.01), str(v), ha='center', fontsize=9)

plt.xticks(rotation=45, ha="right")
plt.ylabel("Number of Cells")
plt.xlabel("Predicted Labels")
plt.title("Cell Counts per Predicted Label")
plt.tight_layout()

# Save figure
fig_path = "figures/cell_counts_per_Predicted_label_sorted_FULLRNA.png"
plt.savefig(fig_path, dpi=300)
plt.close()

adata_q = adata_query.copy()

adata_query = adata_q.copy()


adata_query.obs["merged_type"] = adata_query.obs["predicted_labels"].copy()
#adata_query = adata_query[~np.isin(adata_query.obs["merged_type"],['DC_pDC','DC_CCR7+','T_Treg','NK','B_GC_prePB','Macrophages_M1'])].copy()
adata_query = adata_query[~np.isin(adata_query.obs["merged_type"],['ILC','NK','NKT','Monocytes','T_TIM3+'])].copy()

adata_query.obs["merged_type"] = adata_query.obs["merged_type"].cat.add_categories(['B_mem+activated','DC','Macrophages','T_CD8+','T_CD4+_FH'])
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['B_activated,4','B_activated,3','B_activated,0','B_activated,2','B_activated,1','B_mem','B_IFN']), "merged_type"] = 'B_mem+activated'
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['T_CD4+','T_TfR','T_Treg']), "merged_type"] = 'T_CD4+'
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['T_CD4+_TfH','T_CD4+_TfH_GC']), "merged_type"] = 'T_CD4+_FH'
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['B_GC_prePB']), "merged_type"] = 'B_GC_LZ'
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['T_CD8+_CD161+','T_CD8+_cytotoxic','T_CD8+_naive']), "merged_type"] = 'T_CD8+'
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['VSMC','Endo']), "merged_type"] = 'FDC'
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['Macrophages_M2','Macrophages_M1']), "merged_type"] = 'Macrophages'
adata_query.obs.loc[adata_query.obs["merged_type"].isin(['DC_cDC1','DC_cDC2','DC_pDC','DC_CCR7+']), "merged_type"] = 'DC'

adata_query.obs["merged_type"] = adata_query.obs["merged_type"].cat.remove_unused_categories()
adata_query.obs["merged_type"].cat.categories.tolist()


# plot
sc.pl.umap(
    adata_query,
    color=["merged_type","predicted_labels","scDblFinder_score","label_confidence","pct_counts_mt","log1p_total_counts","pct_counts_in_top_50_genes"],
    ncols=1,
    show=False,
    save="MergedAnnot1_FULLRNA.pdf",
)

# plot
sc.pl.umap(
    adata_query,
    color=["merged_type"],
    ncols=1,
    show=False,
    save="MergedAnnot1_FULLRNA_OnlyMerged.pdf",
)

# --- Leiden clustering ---
sc.tl.leiden(adata_query, key_added="leiden", resolution=5)

# --- Assign clusters to cell types (majority vote) ---
cluster_assign = (
    adata_query.obs.groupby("leiden")["merged_type"]
    .agg(lambda x: x.value_counts().idxmax())
)

adata_query.obs["leiden_merged_type"] = adata_query.obs["leiden"].map(cluster_assign)
adata_query.uns[f"leiden_merged_type_colors"] = adata_query.uns[f"{ct}_colors"]

sc.pl.umap(
    adata_query,
    color=["leiden_merged_type"],
    ncols=1,
    show=False,
    save=f"MergedAnnotClust_FULLRNA.pdf"
)

adata_query.write("Tonsil_FULL_RNA_ANNOTATED.h5ad")


sc.pl.violin(adata_query, 'log1p_total_counts',rotation=0.0000001, show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {adata_query.n_obs}\nMedian UMI: {np.median(adata_query.obs['total_counts'])}")


# Save or show the updated figure
plt.savefig(f'./figures/Tonsil1_rna_nUMIPostFilter.png')
plt.close()


import seaborn as sns
# Count cells per predicted_labels
# Prepare data
label_counts = (
    adata_query.obs["leiden_merged_type"]
    .value_counts()
    .sort_values(ascending=False)
    .reset_index()
)
label_counts.columns = ["leiden_merged_type", "cell_count"]
# Plot
plt.figure(figsize=(10, 6))
sns.barplot(
    data=label_counts,
    x="leiden_merged_type",
    y="cell_count",
    palette="tab20",
    order=label_counts["leiden_merged_type"],  # enforce sorted order
)

# Add count labels on top
for i, v in enumerate(label_counts["cell_count"]):
    plt.text(i, v + (v * 0.01), str(v), ha='center', fontsize=9)

plt.xticks(rotation=45, ha="right")
plt.ylabel("Number of Cells")
plt.xlabel("leiden_merged_type")
plt.title("Cell Counts per leiden_merged_type")
plt.tight_layout()

# Save figure
fig_path = "figures/cell_counts_per_leiden_merged_type_sorted_FULLRNA.png"
plt.savefig(fig_path, dpi=300)
plt.close()







#################################################
############# DNA MODALITY ###############
#################################################

exp ='Sc_VTHumanTonsil'

exp_ac = 'Sc_VTD9_S1_Sc_VTD11_S2_VTHumanTonsil_H3K27ac'
exp_me3 = 'Sc_VTD9_S1_Sc_VTD11_S2_VTHumanTonsil_H3K27me3'

#DNA
ac = sc.read_h5ad(f'../{exp_ac}.h5ad')
me3 = sc.read_h5ad(f'../{exp_me3}.h5ad')

#DNA
minFrags_me3 = 400
maxFrags_me3 = 10000
minFrags_ac = 400
maxFrags_ac = 10000

# --- DNA filtering criteria ---
filtered_ac = ac[snap.pp.filter_cells(ac, min_counts=minFrags_ac,max_counts=maxFrags_ac, min_tsse=0, inplace=False)]
filtered_me3 = me3[snap.pp.filter_cells(me3, min_counts=minFrags_me3,max_counts=maxFrags_me3, min_tsse=0, inplace=False)]


sc.pl.violin(filtered_me3, 'log10_n_fragment',rotation=0.0000001, show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {filtered_me3.n_obs}\nMedian n_fragment: {np.median(filtered_me3.obs['n_fragment'])}")


# Save or show the updated figure
plt.savefig(f'./figures/{exp}_me3_nfragPostFilter.png')
plt.close()


sc.pl.violin(filtered_ac, 'log10_n_fragment',rotation=0.0000001, show=False)

# Get the figure and axes
fig = plt.gcf()
axes = fig.axes

# Loop through the axes to update y-axis tick labels
for ax in axes:
    yticks = ax.get_yticks()
    new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
    ax.set_yticklabels(new_labels)
    ax.set_title(f"n_cells: {filtered_ac.n_obs}\nMedian n_fragment: {np.median(filtered_ac.obs['n_fragment'])}")


# Save or show the updated figure
plt.savefig(f'./figures/{exp}_ac_nfragPostFilter.png')
plt.close()

me3 = filtered_me3.copy()
ac = filtered_ac.copy()

rna = sc.read(f'./Tonsil_FULL_RNA_ANNOTATED.h5ad',backed=None)

ct_Cat = "leiden_merged_type"
#ct_Cat = "merged_type"


# ensure "NotInRNA" is part of categories
rna.obs[ct_Cat] = rna.obs[ct_Cat].cat.add_categories(["NotInRNA"])

def transfer_or_notinrna(target, ref, field):
    # cells shared by target and ref
    common = target.obs_names.intersection(ref.obs_names)

    # build category list = ref categories + "NotInRNA"
    ref_cats = pd.Categorical(ref.obs[field]).categories.tolist()
    if "NotInRNA" not in ref_cats:
        ref_cats.append("NotInRNA")

    # start with all "NotInRNA" but make sure that value is in categories
    vals = pd.Series(
        pd.Categorical(
            ["NotInRNA"] * target.n_obs,
            categories=ref_cats
        ),
        index=target.obs_names,
        dtype="category"
    )

    # assign PLAIN (object) values from ref to avoid categorical-vs-categorical mismatch
    src = ref.obs.loc[common, field].astype(object)
    vals.loc[common] = src.values

    # store as categorical in target.obs
    target.obs[field] = vals.astype("category")
    return target


# apply to both
me3 = transfer_or_notinrna(me3, rna, ct_Cat)
ac  = transfer_or_notinrna(ac,  rna, ct_Cat)


BL = '/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz'
cell_type_colors = {
    # B lineage — cool blues/teals/greens
    "B_Cycling":        "#00BFC4",  # bright cyan-teal (cycling pops)
    "B_GC_DZ":          "#98DF8A",  # light green (DZ)
    "B_GC_LZ":          "#2CA02C",  # medium green (LZ)
    "B_mem+activated":  "#2B8C9E",  # desaturated teal
    "B_naive":          "#1F77B4",  # deep blue
    "B_plasma":         "#76C1FF",  # sky blue

    # APC / Stromal
    "DC":               "#B39B00",  # mustard/gold
    "FDC":              "#7B1FA2",  # violet
    "Macrophages":      "#8C564B",  # brown

    # T lineage — warm reds/oranges/magenta
    "T_CD4+":           "#D62728",  # strong red (CD4 effector/memory)
    "T_CD4+_FH":        "#E377C2",  # magenta (TFH stands out near GCs)
    "T_CD4+_naive":     "#FF9896",  # light salmon (naïve = lighter shade)
    "T_CD8+":           "#FF7F0E",  # orange
    "NotInRNA":         (0.827, 0.827, 0.827, 0.5),
}

# Ensure correct order & mapping
cats = list(cell_type_colors.keys())
palette = [cell_type_colors[c] for c in cats]

for ad,ty in [(me3,'me3'),(ac,'ac')]:
    if ty == 'ac':
        n=150000
        bs=600
    elif ty == 'me3':
        n=120000
        bs=5000

    snap.pp.add_tile_matrix(ad, bin_size=bs ,counting_strategy='paired-insertion', chunk_size=500000, inplace=True)
    snap.pp.select_features(ad, n_features=n,inplace=True,blacklist=BL,filter_lower_quantile=0.001,filter_upper_quantile=0.001)
    #Perform dimension reduction using the spectrum of the normalized graph Laplacian defined by pairwise similarity between cells
    snap.tl.spectral(ad, n_comps=10, weighted_by_sd=False, chunk_size=80000, features='selected', distance_metric='cosine', inplace=True)
    

    n_dims = ad.obsm['X_spectral'].shape[1]
    snap.tl.umap(ad, use_rep='X_spectral', key_added='umap',min_dist=1, random_state=None,use_dims=list(range(1, n_dims)))
    #neighborhood graph of observations stored in data using the method specified by method.
    snap.pp.knn(ad, n_neighbors=30, use_rep='X_spectral', method='kdtree')

    #Plot
    snap.pl.spectral_eigenvalues(ad, width=600, height=400, show=False, interactive=False, out_file=f'./figures/{exp}_{ty}_Eigenvalues.png')
    snap.tl.leiden(ad, resolution=1)
    sc.pl.umap(ad, color=["leiden"],save=f'{exp}_{ty}_UMAP_Leiden.png',show=False)
    
    ad.obs[ct_Cat] = ad.obs[ct_Cat].astype("category")
    ad.obs[ct_Cat] = ad.obs[ct_Cat].cat.set_categories(cats)
    sc.pl.umap(
        ad,
        color=ct_Cat,
        palette=palette,
        title=f"Cell Types (from RNA) - Grey cells are not present in RNA",
        save=f'{exp}_{ty}_UMAP_TypesFromRNA.pdf',show=False
    )
    #sc.pl.umap(ad, color=["n_fragment","log10_n_fragment","FRiP","tsse"],save=f'{exp}_{ty}_UMAP_QC.png',show=False)
    sc.pl.umap(ad, color=["n_fragment","log10_n_fragment","tsse"],save=f'{exp}_{ty}_UMAP_QC.png',show=False)
    ad.write(f'{exp}_{ty}_processed.h5ad')
    del ad






#################################################
############# TRI MODALITY ###############
#################################################

from matplotlib_venn import venn3
# KEEP ONLY overlap CELLS in the 3 mods

exp = 'TonsilRNA'
Experiment = 'TonsilRNA_H2R_202505'
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

# Load adatas

adata1 = sc.read_h5ad("./Tonsil_FULL_RNA_ANNOTATED.h5ad")
adata2 = sc.read_h5ad('./Sc_VTHumanTonsil_ac_processed.h5ad')
adata3 = sc.read_h5ad('./Sc_VTHumanTonsil_me3_processed.h5ad')

# Define the set of barcodes for each dataset after filtering
rna_barcodes = set(adata1.obs.index)
ac_barcodes = set(adata2.obs.index)
me3_barcodes = set(adata3.obs.index)

# Get the intersection of barcodes across all datasets
shared_barcodes = rna_barcodes & ac_barcodes & me3_barcodes

print(f'Ac:{len(ac_barcodes)}\nMe3:{len(me3_barcodes)}\nRNA:{len(rna_barcodes)}\nShared:{len(shared_barcodes)}')

# Create a Venn diagram with the sets
plt.figure(figsize=(8, 8))
venn = venn3([rna_barcodes, ac_barcodes, me3_barcodes],
             set_labels=(f'RNA\n{minUMI}<UMIs<{maxUMI}\n<{max_pct_mt}% MT\n<{max_pct_50g}% 50TopGenes', 
                         f'H3K27ac\n{minFrags_ac}<fragments<{maxFrags_ac}', 
                         f'H3K27me3\n{minFrags_me3}<fragments<{maxFrags_me3}'))
plt.suptitle('Venn Diagram of Passing Cells\nfor RNA, H3K27ac, and H3K27me3', 
          fontsize=16, fontweight='bold')
plt.tight_layout()

# Show the plot
plt.savefig(f'./figures/{exp}_VennPassingCells_Proc.pdf')
plt.close()


# Subset adata2 and adata3 to only the cells that are present in adata1
adata1 = adata1[adata1.obs_names.isin(shared_barcodes)].copy()
adata2 = adata2[adata2.obs_names.isin(shared_barcodes)].copy()
adata3 = adata3[adata3.obs_names.isin(shared_barcodes)].copy()

adata1.write("./Tonsil_rna_Overlap.h5ad")
adata2.write("./Tonsil_ac_Overlap.h5ad")
adata3.write("./Tonsil_me3_Overlap.h5ad")


adata_query = adata1.copy()


adata_query = sc.read_h5ad("./Tonsil_rna_Overlap_Reclust.h5ad")

# Reclust RNA mod

#adata_query.X = adata_query.layers['SCT_counts'].copy()
#sc.pp.normalize_total(adata_query, target_sum=1e4, inplace=True)
#sc.pp.log1p(adata_query)
sc.pp.highly_variable_genes(
    adata_query,
    layer='counts',
    flavor='seurat_v3',
    n_top_genes=3000,
)

sc.pp.scale(adata_query,zero_center=False, max_value=4) # optional clipping
sc.pp.pca(adata_query,n_comps=20)
sc.pl.pca_variance_ratio(adata_query, n_pcs=30, log=False, show=False, save=f'Tonsils_RNA_PCAVarRNA_ALLMODS.pdf')

sc.pp.neighbors(adata_query, use_rep="X_pca", n_neighbors=30,metric='cosine')  # or adjust as needed
sc.tl.umap(adata_query)

sc.pl.umap(
    adata_query,
    color=["leiden_merged_type"],
    ncols=1,
    show=False,
    save="Annot_RNA_ALLMODS.pdf"
)


adata_query.write("./Tonsil_rna_Overlap_Reclust.h5ad")





#################################################
############# DNA TRI MODALITY ###############
#################################################
exp ='Tonsil'

me3 = snap.read(f'{exp}_me3_Overlap.h5ad',backed=None)
ac = snap.read(f'{exp}_ac_Overlap.h5ad',backed=None)


BL = '/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz'

for ad,ty in [(me3,'me3'),(ac,'ac')]:
    if ty == 'ac':
        n=100000
        bs=600
    elif ty == 'me3':
        n=80000
        bs=5000

    snap.pp.select_features(ad, n_features=n,inplace=True,blacklist=BL,filter_lower_quantile=0.001,filter_upper_quantile=0.001)

    #Perform dimension reduction using the spectrum of the normalized graph Laplacian defined by pairwise similarity between cells
    snap.tl.spectral(ad, n_comps=10, weighted_by_sd=False, chunk_size=80000, features='selected', distance_metric='cosine', inplace=True)
    n_dims = ad.obsm['X_spectral'].shape[1]
    snap.tl.umap(ad, use_rep='X_spectral', key_added='umap', random_state=None,use_dims=list(range(1, n_dims)))
    #neighborhood graph of observations stored in data using the method specified by method. The distance metric used is Euclidean.
    snap.pp.knn(ad, n_neighbors=20, use_rep='X_spectral', method='kdtree')

    #Plot
    snap.pl.spectral_eigenvalues(ad, width=600, height=400, show=False, interactive=False, out_file=f'./figures/{exp}_{ty}_Eigenvalues_overlap.png')
    snap.tl.leiden(ad, resolution=0.6)
    sc.pl.umap(ad, color=["leiden_merged_type"],save=f'{exp}_{ty}_UMAP_LeidenTYPE_overlap.pdf',show=False)
    sc.pl.umap(ad, color=["leiden","leiden_merged_type"],save=f'{exp}_{ty}_UMAP_Leiden_overlap.pdf',show=False)
    sc.pl.umap(ad, color=["n_fragment","log10_n_fragment","tsse"],save=f'{exp}_{ty}_UMAP_QC_overlap.png',show=False)
    ad.write(f'{exp}_{ty}_Overlap_Reclust.h5ad')
    del ad


ct_type = "leiden_merged_type"


#################################################
############# Linked UMAPs (scaled correctly) ###
#################################################
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import colors as mcolors
from anndata import AnnData
from typing import Union

def plot_linked_umaps_procrustes(
    adata_rna: Union[str, AnnData],
    adata_ac:  Union[str, AnnData],
    adata_me3: Union[str, AnnData],
    ct_type: str,
    output_prefix: str = "./figures/LinkedUMAPs_RNA_ac_me3_procrustes_lifted2",
    fig_size=(18, 6),
    fig_dpi=300,
    lift_rna_frac: float = 0.30,
    lift_me3_frac: float = 0.30,
    gap_frac: float = 0.30,
):
    """
    Create linked UMAP plots: RNA → H3K27ac → H3K27me3 with Procrustes alignment.
    Everything matches the original script's behavior and styling.

    Parameters
    ----------
    adata_rna, adata_ac, adata_me3 : str | anndata.AnnData
        File paths to .h5ad or preloaded AnnData objects for RNA, H3K27ac, H3K27me3.
    ct_type : str
        Column in RNA .obs that stores the cell type labels.
    output_prefix : str
        Output path prefix (PDF and PNG will be saved).
    fig_size : tuple
        Matplotlib figure size (width, height).
    fig_dpi : int
        Matplotlib DPI.
    lift_rna_frac, lift_me3_frac : float
        Vertical lifts (fractions of max panel height) applied to RNA and H3K27me3.
    gap_frac : float
        Horizontal gap as a fraction of max panel width.
    """

    # ---------- helpers ----------
    def _ensure_adata(x):
        if isinstance(x, AnnData):
            return x
        elif isinstance(x, str):
            return sc.read_h5ad(x)
        else:
            raise TypeError("adata inputs must be AnnData or file paths to .h5ad")

    def procrustes_align_to_base(target, base):
        """
        Align 'target' (n x 2) to 'base' (n x 2) in the *base's* scale and position.
        Returns the aligned target (n x 2).
        """
        X = np.asarray(base, dtype=float)
        Y = np.asarray(target, dtype=float)

        Xmean = X.mean(0, keepdims=True)
        Ymean = Y.mean(0, keepdims=True)
        Xc = X - Xmean
        Yc = Y - Ymean

        # Orthogonal Procrustes (rotation/reflection + *global* scale)
        M = Yc.T @ Xc
        U, S, Vt = np.linalg.svd(M, full_matrices=False)
        R = U @ Vt
        s = S.sum() / (Yc**2).sum()

        # Apply transform + translate into base coordinates
        Y_aligned = s * (Yc @ R) + Xmean
        return Y_aligned

    def x_extent(arr):  # width helper
        return float(arr[:, 0].max() - arr[:, 0].min())

    def y_extent(arr):  # height helper
        return float(arr[:, 1].max() - arr[:, 1].min())

    # ---------- I/O ----------
    a1 = _ensure_adata(adata_rna)  # base (RNA)
    a2 = _ensure_adata(adata_ac)   # H3K27ac
    a3 = _ensure_adata(adata_me3)  # H3K27me3

    # ---------- intersect cells ----------
    common = a1.obs_names.intersection(a2.obs_names).intersection(a3.obs_names)
    a1 = a1[common]
    a2 = a2[common]
    a3 = a3[common]

    umap1 = pd.DataFrame(a1.obsm['X_umap'], index=common, columns=['x1','y1'])
    umap2 = pd.DataFrame(a2.obsm['X_umap'], index=common, columns=['x2','y2'])
    umap3 = pd.DataFrame(a3.obsm['X_umap'], index=common, columns=['x3','y3'])

    # ---------- align ac & me3 to RNA (same size now) ----------
    umap2[['x2','y2']] = procrustes_align_to_base(
        umap2[['x2','y2']].values, umap1[['x1','y1']].values
    )
    umap3[['x3','y3']] = procrustes_align_to_base(
        umap3[['x3','y3']].values, umap1[['x1','y1']].values
    )

    # ---------- spread panels horizontally with consistent gap ----------
    w = max(
        x_extent(umap1[['x1','y1']].values),
        x_extent(umap2[['x2','y2']].values),
        x_extent(umap3[['x3','y3']].values),
    )
    gap = gap_frac * w  # consistent visual gap

    # shift panels
    x1min = umap1['x1'].min()
    x2min = umap2['x2'].min()
    x3min = umap3['x3'].min()

    umap2['x2'] += (umap1['x1'].max() - x2min) + gap
    umap3['x3'] += (umap2['x2'].max() - x3min) + gap

    # ---------- vertical lift to reduce line overlap ----------
    H = max(
        y_extent(umap1[['x1','y1']].values),
        y_extent(umap2[['x2','y2']].values),
        y_extent(umap3[['x3','y3']].values),
    )
    umap1['y1'] += lift_rna_frac * H   # move RNA up
    umap3['y3'] += lift_me3_frac * H   # move me3 up

    # ---------- combine + colors ----------
    df = pd.concat([umap1, umap2, umap3], axis=1)
    df['cell_type'] = a1.obs[ct_type].astype('category')

    cell_type_colors = {
        # B lineage — cool blues/teals/greens
        "B_Cycling":        "#00BFC4",
        "B_GC_DZ":          "#98DF8A",
        "B_GC_LZ":          "#2CA02C",
        "B_mem+activated":  "#2B8C9E",
        "B_naive":          "#1F77B4",
        "B_plasma":         "#76C1FF",

        # APC / Stromal
        "DC":               "#B39B00",
        "FDC":              "#7B1FA2",
        "Macrophages":      "#8C564B",

        # T lineage — warm reds/oranges/magenta
        "T_CD4+":           "#D62728",
        "T_CD4+_FH":        "#E377C2",
        "T_CD4+_naive":     "#FF9896",
        "T_CD8+":           "#FF7F0E",
    }

    to_rgba = {k: mcolors.to_rgba(v, alpha=0.22) for k,v in cell_type_colors.items()}
    fallback = mcolors.to_rgba("#BDBDBD", alpha=0.22)
    colors_arr = np.vstack([to_rgba.get(ct, fallback) for ct in df['cell_type'].astype(object)])

    # ---------- plot ----------
    fig, ax = plt.subplots(figsize=fig_size, dpi=fig_dpi)
    ax.set_title("RNA → H3K27ac → H3K27me3", fontsize=20, weight='bold')

    # lines behind points
    seg12 = np.stack([df[['x1','y1']].values, df[['x2','y2']].values], axis=1)
    seg23 = np.stack([df[['x2','y2']].values, df[['x3','y3']].values], axis=1)
    ax.add_collection(LineCollection(seg12, colors=colors_arr, linewidths=0.25, zorder=0))
    ax.add_collection(LineCollection(seg23, colors=colors_arr, linewidths=0.25, zorder=0))

    # points
    for (x,y,label_panel) in [('x1','y1','RNA'), ('x2','y2','H3K27ac'), ('x3','y3','H3K27me3')]:
        for ct, col in cell_type_colors.items():
            sub = df[df['cell_type'] == ct]
            if not len(sub):
                continue
            ax.scatter(sub[x], sub[y], s=30, c=col, edgecolors='none',
                       alpha=0.98, label=ct if label_panel=='RNA' else None,
                       zorder=1, rasterized=False)

    # keep geometry faithful and tidy
    ax.set_aspect('equal', adjustable='datalim')
    ax.margins(0.03)
    ax.axis('off')

    # one legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    leg = ax.legend(by_label.values(), by_label.keys(), title='Cell Type',
                    bbox_to_anchor=(1.01, 1), loc='upper left', markerscale=1.6, fontsize=9)
    leg.get_title().set_fontsize(10)

    plt.tight_layout()

    # ---------- save ----------
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pdf_path = f"{output_prefix}.pdf"
    png_path = f"{output_prefix}.png"
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(png_path, bbox_inches='tight')
    plt.close()

    return {"pdf": pdf_path, "png": png_path, "coords": df}

# ---- Example call (same files as your script) ----
results = plot_linked_umaps_procrustes(
     "Tonsil_rna_Overlap_Reclust.h5ad",
     "Tonsil_ac_Overlap_Reclust.h5ad",
     "Tonsil_me3_Overlap_Reclust.h5ad",
     ct_type="leiden_merged_type"  # or your column name
)
print(results["pdf"], results["png"])



def conditional_props_rowsum1(adata, subset_name, Experiment="TonsilRNA_H2R_202505"):
    """
    Heatmap of conditional proportions P(cluster | cell type).
    Rows (cell types) sum to 1.

    IMPORTANT: Sorting is matched to `jaccard_types`:
      - Columns sorted by column max of the Jaccard matrix (desc)
      - Rows: desired_order first, then remaining by row max of the Jaccard matrix (desc)
    """
    print(f'Plotting {Experiment} (row-conditional proportions, Jaccard-matched order)')
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import pandas as pd
    import numpy as np
    import os

    # Desired y-axis order (cell types)
    desired_order = [
        "B_Cycling","B_GC_DZ","B_GC_LZ","B_mem+activated","B_naive","B_plasma",
        "DC","FDC","Macrophages",
        "T_CD4+","T_CD4+_FH","T_CD4+_naive","T_CD8+",
    ]

    # Dataframes for types
    GEXdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': adata.obs_names, 'type': adata.obs['leiden']})
    ATACdf = pd.DataFrame({'Experiment': Experiment, 'obs_names': adata.obs_names, 'type': adata.obs['leiden_merged_type']})

    # Inner join on cells present in both
    merged_df = ATACdf.merge(GEXdf, on='obs_names', how='inner')

    # Drop "NotInRNA" rows before plotting
    merged_df = merged_df[merged_df['type_x'] != 'NotInRNA'].copy()
    if merged_df.empty:
        print("No cells remain after removing 'NotInRNA'; skipping plot.")
        return

    groups1 = np.asarray(merged_df['type_x'])  # cell types (y-axis)
    groups2 = np.asarray(merged_df['type_y'])  # clusters  (x-axis)

    cl1 = np.unique(groups1)  # unique cell types
    cl2 = np.unique(groups2)  # unique clusters

    # --- Precompute indices and sizes for speed ---
    idx_by_g1 = {g1: np.where(groups1 == g1)[0] for g1 in cl1}
    idx_by_g2 = {g2: np.where(groups2 == g2)[0] for g2 in cl2}
    size1_vec = np.array([idx_by_g1[g1].size for g1 in cl1], dtype=float)
    size2_vec = np.array([idx_by_g2[g2].size for g2 in cl2], dtype=float)

    # --- Build intersection COUNT matrix + Jaccard matrix ---
    count_mat = np.zeros((len(cl1), len(cl2)), dtype=float)
    jaccard_mat = np.zeros_like(count_mat)
    for i, g1 in enumerate(cl1):
        idx1 = idx_by_g1[g1]
        s1 = size1_vec[i]
        if s1 == 0:
            continue
        for j, g2 in enumerate(cl2):
            idx2 = idx_by_g2[g2]
            s2 = size2_vec[j]
            inter = np.intersect1d(idx1, idx2, assume_unique=False).size
            count_mat[i, j] = inter
            denom = s1 + s2 - inter
            jaccard_mat[i, j] = (inter / denom) if denom > 0 else 0.0

    # --- Convert to conditional proportions P(g2 | g1) (row-normalize counts) ---
    with np.errstate(invalid='ignore', divide='ignore'):
        prop_mat = (count_mat.T / np.where(size1_vec > 0, size1_vec, 1)).T
    prop_mat[np.isnan(prop_mat)] = 0.0

    # === ORDERING (match jaccard_types) ===
    # Columns: by column max of Jaccard (desc)
    col_max = np.max(jaccard_mat, axis=0)
    sorted_indices_x = np.argsort(-col_max)
    cl2_sorted = cl2[sorted_indices_x]

    # Rows: desired_order first (only those present), then remaining by row max of Jaccard (desc)
    row_max = np.max(jaccard_mat, axis=1)
    label_to_row = {label: idx for idx, label in enumerate(cl1)}
    present_desired = [lbl for lbl in desired_order if lbl in label_to_row]
    remaining = [lbl for lbl in cl1 if lbl not in present_desired]
    remaining_sorted = sorted(remaining, key=lambda lbl: -row_max[label_to_row[lbl]])
    cl1_sorted = np.array(present_desired + remaining_sorted, dtype=object)

    # Reindex matrices to the chosen orders
    row_idx = np.array([label_to_row[lbl] for lbl in cl1_sorted], dtype=int)
    prop_mat_ord = prop_mat[row_idx][:, sorted_indices_x]
    count_mat_ord = count_mat[row_idx][:, sorted_indices_x]  # optional counts for annotations

    # Counts for labels (match new orders)
    groups2_counts = [int(size2_vec[np.where(cl2 == g)[0][0]]) for g in cl2_sorted]
    groups1_counts = [int(size1_vec[np.where(cl1 == g)[0][0]]) for g in cl1_sorted]

    # Labels with counts
    x_labels = [f"{label} ({count})" for label, count in zip(cl2_sorted, groups2_counts)]
    y_labels = [f"{label} ({count})" for label, count in zip(cl1_sorted, groups1_counts)]

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(prop_mat_ord, cmap='coolwarm', vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(cl2_sorted)))
    ax.set_yticks(np.arange(len(cl1_sorted)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    ax.set_yticklabels(y_labels)

    # Annotate with proportions
    for i in range(prop_mat_ord.shape[0]):
        for j in range(prop_mat_ord.shape[1]):
            ax.text(j, i, f"{prop_mat_ord[i, j]:.2f}", ha="center", va="center",
                    color="black", fontsize=8)

    ax.set_xlabel(f'{subset_name} Leiden Clusters')
    ax.set_ylabel('Cell Type (# cells)')
    title = f'{subset_name[:-3]}' if len(subset_name) > 0 else f'{Experiment}'
    ax.set_title(f'{title}', pad=18)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    plt.colorbar(im, cax=cax)
    plt.tight_layout()

    outdir = './figures'
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(f'{outdir}/{Experiment}_{subset_name}_CondProp_rowsum1.pdf')
    plt.close()



adata1 = sc.read_h5ad("Tonsil_rna_Overlap_Reclust.h5ad")       # base (RNA)
adata2 = sc.read_h5ad("Tonsil_ac_Overlap_Reclust.h5ad")      # H3K27ac
adata3 = sc.read_h5ad("Tonsil_me3_Overlap_Reclust.h5ad")     # H3K27me3


from sklearn.metrics import normalized_mutual_info_score, silhouette_score

def _pick_rep(adata, use_rep):
    if use_rep is not None:
        if use_rep not in adata.obsm:
            raise ValueError(f"{use_rep} not found in adata.obsm")
        return use_rep
    for rep in ("X_pca", "X_lsi", "X_svd", "X_spectral"):
        if rep in adata.obsm:
            return rep
    raise ValueError("No suitable embedding found in .obsm. Compute PCA/LSI and retry.")

def _ensure_neighbors(adata, use_rep):
    metric = "cosine" if "spectral" in use_rep.lower() else "euclidean"
    sc.pp.neighbors(adata, use_rep=use_rep, metric=metric)

def _labels_key(r, res):
    return f"_tmp_leiden_r{r}_res{str(res).replace('.','_')}"

def _pairwise_nmi_values(label_matrix):
    n = len(label_matrix)
    vals = []
    for i in range(n):
        for j in range(i+1, n):
            vals.append(normalized_mutual_info_score(label_matrix[i], label_matrix[j]))
    return vals

def auto_res_candidates(
    adata,
    resolutions = np.round(np.linspace(0.1, 2, 60), 2),
    use_rep = None,
    n_repeats = 10,
    random_state = 0,
    min_cluster_size = 10,
    sample_for_silhouette = 10000
):
    """
    Returns:
      results_df, recommendations, stability_long
        - results_df: per-resolution metrics
        - recommendations: list of (resolution, n_clusters)
        - stability_long: long-form df with columns ['N. clusters','Stability'] for plotting
    """
    rep = _pick_rep(adata, use_rep)
    _ensure_neighbors(adata, rep)

    label_bank = {res: [] for res in resolutions}
    n_cells = adata.n_obs

    for r in range(n_repeats):
        seed = int(random_state + r)
        for res in resolutions:
            key = _labels_key(r, res)
            sc.tl.leiden(adata, resolution=float(res), key_added=key, random_state=seed)
            label_bank[res].append(adata.obs[key].to_numpy(dtype=str))

    # silhouette subsample
    if n_cells > sample_for_silhouette:
        rng = np.random.default_rng(random_state)
        idx_sil = np.sort(rng.choice(n_cells, size=sample_for_silhouette, replace=False))
    else:
        idx_sil = slice(None)
    X = adata.obsm[rep]
    X_for_sil = X[idx_sil]

    rows = []
    stab_rows = []  # for plotting like your autok_stability
    for res in resolutions:
        labels0 = label_bank[res][0]
        _, counts = np.unique(labels0, return_counts=True)
        n_clusters = int(len(counts))
        min_size = int(counts.min()) if n_clusters > 0 else 0

        if n_clusters >= 2:
            labs_for_sil = labels0[idx_sil] if isinstance(idx_sil, np.ndarray) else labels0
            try:
                sil = float(silhouette_score(X_for_sil, labs_for_sil, metric="euclidean"))
            except Exception:
                sil = np.nan
        else:
            sil = np.nan

        pair_vals = _pairwise_nmi_values(label_bank[res])
        mean_nmi = float(np.mean(pair_vals)) if len(pair_vals) else np.nan

        rows.append({
            "resolution": float(res),
            "n_clusters": n_clusters,
            "stability_nmi": mean_nmi,
            "silhouette": sil,
            "min_cluster_size": min_size
        })
        # build long df for stability plot (same spirit as your seaborn lineplot)
        stab_rows.extend([{"N. clusters": n_clusters, "Stability": v} for v in pair_vals])

    df = pd.DataFrame(rows).sort_values(["n_clusters", "resolution"]).reset_index(drop=True)
    stability_long = pd.DataFrame(stab_rows)

    # combined score & suggestions
    df["rank_stability"] = df["stability_nmi"].rank(ascending=False, method="min")
    df["rank_silhouette"] = df["silhouette"].rank(ascending=False, method="min")
    penalty = (df["min_cluster_size"] < min_cluster_size).astype(int) * 2.0
    df["score"] = (df["rank_stability"] + df["rank_silhouette"]) / 2.0 + penalty

    best = (
        df.sort_values("score")
          .groupby("n_clusters", as_index=False)
          .first()
          .sort_values("score")
          .head(5)
    )
    recommendations = list(best[["resolution", "n_clusters"]].itertuples(index=False, name=None))
    return df, recommendations, stability_long

def save_autok_plots(results_df, stability_long, outdir, prefix="autoK",
                     highlight_k=None, dot_size=160, dot_color="red"):
    """
    Save stability-vs-K and silhouette-vs-K plots.
    If highlight_k is provided, draw a larger red dot at that K.
    """
    os.makedirs(outdir, exist_ok=True)

    # ---- Stability (mean NMI per K) ----
    if not stability_long.empty:
        means = (stability_long
                 .groupby("N. clusters", as_index=False)["Stability"].mean()
                 .rename(columns={"N. clusters":"K"}))
        fig, ax = plt.subplots(figsize=(5.0, 3.2))
        ax.plot(means["K"], means["Stability"], marker="o")
        ax.set_xlabel("N. clusters")
        ax.set_ylabel("Stability (NMI)")
        ax.set_title("Clustering stability")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        if highlight_k is not None and highlight_k in set(means["K"]):
            y = float(means.loc[means["K"] == highlight_k, "Stability"].iloc[0])
            ax.scatter([highlight_k], [y], s=dot_size, color=dot_color, zorder=5)
        elif highlight_k is not None:
            print(f"[warn] highlight_k={highlight_k} not found in stability K values: {sorted(means['K'].unique())}")

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"{prefix}_stability_vs_K.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(outdir, f"{prefix}_stability_vs_K.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ---- Silhouette (mean per K across rows in results_df) ----
    sil = (results_df
           .dropna(subset=["silhouette"])
           .groupby("n_clusters", as_index=False)["silhouette"].mean()
           .rename(columns={"n_clusters":"K"}))
    if not sil.empty:
        fig, ax = plt.subplots(figsize=(5.0, 3.2))
        ax.plot(sil["K"], sil["silhouette"], marker="o")
        ax.set_xlabel("N. clusters")
        ax.set_ylabel("Silhouette")
        ax.set_title("Silhouette vs. number of clusters")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        if highlight_k is not None and highlight_k in set(sil["K"]):
            y = float(sil.loc[sil["K"] == highlight_k, "silhouette"].iloc[0])
            ax.scatter([highlight_k], [y], s=dot_size, color=dot_color, zorder=5)
        elif highlight_k is not None:
            print(f"[warn] highlight_k={highlight_k} not found in silhouette K values: {sorted(sil['K'].unique())}")

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"{prefix}_silhouette_vs_K.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(outdir, f"{prefix}_silhouette_vs_K.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)


# RNA
res1_df, res1_best, stab1 = auto_res_candidates(adata1, use_rep="X_pca")
# H3K27ac
res2_df, res2_best, stab2 = auto_res_candidates(adata2, use_rep="X_spectral")
# H3K27me3
res3_df, res3_best, stab3 = auto_res_candidates(adata3, use_rep="X_spectral")

# RNA
save_autok_plots(res1_df, stab1, outdir="Plots/RNA",     prefix="RNA",     highlight_k=7)

# H3K27ac
save_autok_plots(res2_df, stab2, outdir="Plots/H3K27ac", prefix="H3K27ac", highlight_k=5)

# H3K27me3
save_autok_plots(res3_df, stab3, outdir="Plots/H3K27me3",prefix="H3K27me3",highlight_k=5)

def pick_resolution_for_K(results_df, K):
    """Return the resolution that best yields K clusters."""
    dfK = results_df.loc[results_df["n_clusters"] == K].copy()
    if dfK.empty:
        raise ValueError(f"No rows with n_clusters == {K}. Try widening the resolution grid.")
    # Prefer lowest combined score if present, else maximize stability then silhouette, then larger min cluster
    if "score" in dfK:
        dfK = dfK.sort_values(["score", "resolution"], ascending=[True, True])
    else:
        dfK = dfK.sort_values(
            ["stability_nmi", "silhouette", "min_cluster_size", "resolution"],
            ascending=[False, False, False, True],
        )
    return float(dfK.iloc[0]["resolution"])

# Current K for each object (from your example)
K_current = {
    "RNA": 7,         # adata1 / res1_df / sc
    "H3K27ac": 5,     # adata2 / res2_df / snap
    "H3K27me3": 5     # adata3 / res3_df / snap
}

objs = {
    "RNA":     {"adata": adata1, "res_df": res1_df, "typo": "rna_Overlap3",  "engine": "scanpy"},
    "H3K27ac": {"adata": adata2, "res_df": res2_df, "typo": "ac_Overlap3",   "engine": "snap"},
    "H3K27me3":{"adata": adata3, "res_df": res3_df, "typo": "me3_Overlap3",  "engine": "snap"},
}

for name, meta in objs.items():
    ad = meta["adata"]
    res_df = meta["res_df"]
    typo = meta["typo"]
    engine = meta["engine"]
    K0 = K_current[name]

    # Iterate K = K0+1 ... 2*K0 (inclusive)
    for K in range(K0 + 1, 2 * K0 + 1):
        # Look up the resolution corresponding to this target K
        res = pick_resolution_for_K(res_df, K=K)
        if res is None:
            print(f"[{name}] Skipping K={K}: no matching resolution found.")
            continue

        # Run Leiden with the looked-up resolution (overwrite 'leiden' each time)
        if engine == "scanpy":
            sc.tl.leiden(ad, resolution=res, key_added="leiden")
        else:
            snap.tl.leiden(ad, resolution=res, key_added="leiden")
        try:
            res_str = f"{float(res):.3f}"
        except Exception:
            res_str = str(res)

        # Your overlap/summary helpers
        jaccard_types(ad, f"{typo}_K{K}_res{res_str}")
        conditional_props_rowsum1(ad, f"{typo}_K{K}_res{res_str}")

        # Save a UMAP for this K/res combo (file name encodes K and res)
        # Note: scanpy appends this to 'figures/umap{save}'
        sc.pl.umap(
            ad,
            color=["leiden"],
            save=f"_Tonsil_{typo}_K{K}_res{res_str}_UMAP_LeidenClusters.pdf",
            show=False
        )

        print(f"[{name}] Done: K={K}, resolution={res}")




adata = sc.read_h5ad("./Tonsil_FULL_RNA_ANNOTATED.h5ad")
adata1 = sc.read_h5ad("Tonsil_rna_Overlap_Reclust.h5ad")       # base (RNA)
adata2 = sc.read_h5ad("Tonsil_ac_Overlap_Reclust.h5ad")      # H3K27ac


adata3 = sc.read_h5ad("Tonsil_ac_Overlap_Reclust.h5ad")     # H3K27me3

snap.pp.select_features(adata3, n_features=100000,inplace=True,filter_lower_quantile=0.001,filter_upper_quantile=0.001)

snap.tl.spectral(adata3, n_comps=12, weighted_by_sd=False, chunk_size=80000, features='selected', distance_metric='cosine', inplace=True)

snap.pp.knn(adata3, n_neighbors=20, use_rep='X_spectral', method='kdtree')

snap.tl.leiden(adata3, resolution=1.1, key_added="leiden")
print(np.unique(adata3.obs['leiden']))

conditional_props_rowsum1(adata3, f"ac_overlap_K7")

sc.pl.umap(adata3,color=["leiden"],save=f"_Tonsil_me3_overlap_K7.pdf",show=False)


sc.pp.highly_variable_genes(
    adata,
    layer='counts',
    flavor='seurat_v3',
    n_top_genes=4000,
)
sc.pp.scale(adata,zero_center=False, max_value=10) # optional clipping
sc.pp.pca(adata,n_comps=20)

sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=20,metric='cosine')  # or adjust as needed
sc.tl.leiden(adata, resolution=0.54 , key_added="leiden")
print(adata.obs['merged_type'])

# --- pull resolutions for your three objects ---
# Assumes you already have: res1_df (RNA), res2_df (H3K27ac), res3_df (H3K27me3)
res_RNA     = pick_resolution_for_K(res1_df, K=7)   # adata1
res_H3K27ac = pick_resolution_for_K(res2_df, K=5)   # adata2
res_H3K27me3= pick_resolution_for_K(res3_df, K=5)   # adata3

print(f"Chosen resolutions → RNA: {res_RNA}, H3K27ac: {res_H3K27ac}, H3K27me3: {res_H3K27me3}")

# --- run Leiden with those exact resolutions
sc.tl.leiden(adata1, resolution=res_RNA, key_added="leiden")
snap.tl.leiden(adata2, resolution=res_H3K27ac, key_added="leiden")
snap.tl.leiden(adata3, resolution=res_H3K27me3, key_added="leiden")


for adaa,typo in [(adata1,"rna_Overlap3"),(adata2,"ac_Overlap3"),(adata3,"me3_Overlap3")]:
    conditional_props_rowsum1(adaa,typo)
    sc.pl.umap(adaa, color=["leiden"],save=f'Tonsil_{typo}_UMAP_LeidenClusters.pdf',show=False)




############################################
# --- Per-group & Full bigWig coverage ---
############################################
BL='../../../hg38-blacklist.bed.gz'
out_dir = "bw_tracks"
os.makedirs(out_dir, exist_ok=True)

for ty in ['me3','ac']:
    adata = sc.read_h5ad(f"Tonsil_{ty}_Overlap_Reclust.h5ad")
    # Build the group list, excluding a "NotInRNA" bucket if present
    group_col = "leiden_merged_type"
    groups = pd.Series(adata.obs[group_col].astype(str)).unique().tolist()
    selections = sorted([g for g in groups if g != "NotInRNA"])

    print(f"Exporting bigWig for {len(selections)} groups:", selections)

    snap.ex.export_coverage(
        adata,
        groupby=group_col,
        selections=selections,          # only these groups
        bin_size=20,                    # adjust if you want smoother/coarser tracks (e.g. 25/50/100)
        blacklist=BL,                 # e.g. Path("ENCFF356LFX_hg38_blacklist.bed")
        normalization="RPKM",           # RPKM/CPM/BPM or None
        include_for_norm=None,          # optionally restrict norm to promoters/peaks (BED or list like ["chr1:1-100"])
        exclude_for_norm=None,
        min_frag_length=None,
        max_frag_length=2000,
        counting_strategy="fragment",   # "fragment", "insertion", or "paired-insertion"
        smooth_base=None,               # e.g. 200 for light smoothing in the output track
        out_dir=out_dir,
        prefix=f"Overlap3Mods_{ty}_",                 # file name prefix, e.g. ATAC_B_naive.bw
        suffix=".bw",                   # ".bw" for bigWig, or ".bedgraph.gz" for bedGraph
        output_format="bigwig",         # or "bedgraph"
        compression=None,               # only relevant for bedGraph
        compression_level=None,
        tempdir=None,
        n_jobs=8,                       # parallelize across groups/chroms
    )

    snap.ex.export_coverage(
        adata,
        groupby='Exp',
        bin_size=20,                    # adjust if you want smoother/coarser tracks (e.g. 25/50/100)
        blacklist=BL,                 # e.g. Path("ENCFF356LFX_hg38_blacklist.bed")
        normalization="RPKM",           # RPKM/CPM/BPM or None
        include_for_norm=None,          # optionally restrict norm to promoters/peaks (BED or list like ["chr1:1-100"])
        exclude_for_norm=None,
        min_frag_length=None,
        max_frag_length=2000,
        counting_strategy="fragment",   # "fragment", "insertion", or "paired-insertion"
        smooth_base=None,               # e.g. 200 for light smoothing in the output track
        out_dir=out_dir,
        prefix=f"Overlap3Mods_Full",                 # file name prefix, e.g. ATAC_B_naive.bw
        suffix=".bw",                   # ".bw" for bigWig, or ".bedgraph.gz" for bedGraph
        output_format="bigwig",         # or "bedgraph"
        compression=None,               # only relevant for bedGraph
        compression_level=None,
        tempdir=None,
        n_jobs=8,                       # parallelize across groups/chroms
    )
    print("Done. Files are in:", os.path.abspath(out_dir))
















#!/usr/bin/env python3
# dge_markers_leiden_pseudobulk.py
# Pseudobulk cluster-marker analysis:
#   For each cell type in obs["leiden_merged_type"]:
#     cluster vs rest (logFC = cluster − rest)
# using edgeR QLF + DESeq2 with pseudobulk replicates.

import os
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from adjustText import adjust_text

# ------------------- Configuration -------------------
OUTDIR = "figures_cluster_markers_pseudobulk"
os.makedirs(OUTDIR, exist_ok=True)

# <<< CHANGE THIS TO YOUR H5AD FILE >>>
H5AD_PATH = "./Tonsil_FULL_RNA_ANNOTATED.h5ad"

# Cell-type / cluster annotation
CELLTYPE_KEY = "leiden_merged_type"

# Minimum cells in a cluster to consider it for marker calling
PB_MIN_CELLS_PER_GROUP = 50

# Minimum cells per condition (cluster / rest) for pseudobulk building
PB_MIN_CELLS_PER_CONDITION = 25

# How many pseudobulk replicates to create per condition
REPLICATES_PER_CONDITION = 3

# DE / plotting thresholds
FDR_THR = 0.01
LFC_THR = 1.0

# Annotation controls
TOP_ANNOTATE     = 80
ANNO_MIN_LOGCPM  = 1.0
ANNO_MIN_ABS_LFC = 1.0

# Toggle annotation strategy: "smart" or "simple"
ANNOTATION_MODE = "simple"

# Point-size scaling
POINTSIZE_MIN   = 6.0
POINTSIZE_MAX   = 60.0
POINTSIZE_POWER = 0.5  # sqrt-like scaling of detection rate

# Binary contrast labels: logFC = ALT - REF = cluster - rest
COND_REF = "rest"
COND_ALT = "cluster"


# ------------------- Helper Functions -------------------
def filter_genes_drop_malat1_mt(A):
    v  = pd.Index(A.var_names.astype(str))
    up = v.str.upper()
    keep = (
        (~up.str.contains("MALAT1"))
        & (~up.str.startswith(("MT-", "MT_", "MT.")))
        & (~up.str.contains("ENSG"))
    )
    if keep.sum() == 0:
        raise ValueError("After filtering MALAT1/MT/ENSG genes, no genes remain.")
    return A[:, keep].copy()


def pick_counts_layer(A):
    for cand in ("counts", "raw_counts"):
        if cand in A.layers:
            return cand
    X = A.X
    if sp.issparse(X):
        ok = np.all(np.equal(np.asarray(X.data), np.asarray(X.data).astype(int)))
    else:
        ok = np.all(np.equal(X, X.astype(int)))
    if not ok:
        raise ValueError(
            "No counts-like layer found (counts/SCT_counts/raw_counts) and .X is not integer."
        )
    return None


def _sum_over_cells(X):
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel()
    return X.sum(axis=0)


def _gene_detection(A, counts_layer):
    """
    Return DataFrame indexed by gene with:
      n_cells_detected: number of cells with counts > 0
      pct_cells_detected: proportion of cells with counts > 0
    """
    X = A.layers[counts_layer] if counts_layer is not None else A.X
    if sp.issparse(X):
        nnz = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        nnz = (X > 0).sum(axis=0)
    n_cells = A.n_obs
    out = pd.DataFrame(
        {
            "n_cells_detected": nnz.astype(float),
            "pct_cells_detected": nnz.astype(float) / float(n_cells),
        },
        index=A.var_names,
    )
    out.index.name = None
    return out, n_cells


def make_pseudobulk_binary(A, counts_layer, cond_key,
                           min_cells=30, reps_per_cond=4):
    """
    Build pseudobulk count matrices for a *binary* contrast:

    - cond_key must have exactly two levels: COND_REF and COND_ALT.
    - Cells from each condition are shuffled and split into `reps_per_cond`
      approximately equal-sized chunks.
    - Each chunk becomes one pseudobulk "sample".
    """
    if cond_key not in A.obs:
        raise KeyError(f"obs['{cond_key}'] missing")

    A = A.copy()
    A.obs[cond_key] = A.obs[cond_key].astype("category")

    if set(A.obs[cond_key].cat.categories) != {COND_REF, COND_ALT}:
        A = A[A.obs[cond_key].isin([COND_REF, COND_ALT])].copy()
        A.obs[cond_key] = A.obs[cond_key].astype("category")

    if A.obs[cond_key].nunique() != 2:
        raise ValueError(f"Need exactly two conditions: {COND_REF} and {COND_ALT}")

    Xmat = A.layers[counts_layer] if counts_layer is not None else A.X

    units = []
    metas = []

    for cond in (COND_REF, COND_ALT):
        mask = (A.obs[cond_key] == cond).to_numpy()
        idx = np.where(mask)[0]
        if idx.size < min_cells:
            print(
                f"  [WARN] condition '{cond}' has only {idx.size} cells "
                f"(< {min_cells}); skipping this contrast."
            )
            return None, None

        idx = np.random.permutation(idx)
        chunks = np.array_split(idx, reps_per_cond)

        for rep_i, chunk in enumerate(chunks, start=1):
            if chunk.size == 0:
                continue
            vec = _sum_over_cells(Xmat[chunk])
            units.append(vec)
            unit_name = f"{cond}__rep{rep_i}"
            metas.append(
                {
                    "unit": unit_name,
                    "condition": cond,
                    "sample": unit_name,  # each pseudobulk is its own "sample"
                }
            )

    counts = np.vstack(units)
    counts_df = pd.DataFrame(
        counts.T,
        index=A.var_names,
        columns=[m["unit"] for m in metas],
    )
    meta_df = pd.DataFrame(metas).set_index("unit")
    meta_df["condition"] = meta_df["condition"].astype("category")
    return counts_df, meta_df


def have_replicates(meta_df, min_units_per_condition=2):
    """
    Decide if we have enough pseudobulk samples to run edgeR/DESeq2 with p-values:

    Require >= min_units_per_condition pseudobulks per condition.
    """
    n_ref = (meta_df["condition"] == COND_REF).sum()
    n_alt = (meta_df["condition"] == COND_ALT).sum()
    return (n_ref >= min_units_per_condition) and (n_alt >= min_units_per_condition)


# ---------- Smart annotation ranking ----------
def _smart_annot_rank(df, top_n,
                      min_logcpm=ANNO_MIN_LOGCPM,
                      min_abs_lfc=ANNO_MIN_ABS_LFC):
    """
    df must contain: gene, logCPM, logFC (and optionally mlog10FDR).
    Returns top_n by equal-weight score of expression and |FC| past thresholds.
    """
    d = df.copy()
    d["abs_logFC"] = d["logFC"].abs()

    # Must exceed both minima
    cand = d[(d["logCPM"] > min_logcpm) & (d["abs_logFC"] > min_abs_lfc)].copy()
    if cand.empty:
        return cand

    # Excess over thresholds
    cand["expr_excess"] = cand["logCPM"] - min_logcpm
    cand["fc_excess"] = cand["abs_logFC"] - min_abs_lfc

    # Robust scaling by 90th percentile to balance terms
    def _scale(vals):
        vals = np.asarray(vals)
        if vals.size == 0:
            return vals
        p90 = np.percentile(vals, 90)
        denom = p90 if p90 > 0 else (vals.max() if vals.max() > 0 else 1.0)
        return vals / denom

    cand["score"] = 0.5 * _scale(cand["expr_excess"].values) + 0.5 * _scale(
        cand["fc_excess"].values
    )
    cand = cand.sort_values("score", ascending=False).head(top_n)
    return cand


def _sizes_from_detection(detected_prop):
    """
    Map detection proportion [0..1] -> point size.
    Use sqrt-like scaling to reduce dynamic range.
    """
    p = np.clip(np.asarray(detected_prop).astype(float), 0.0, 1.0)
    p = np.power(p, POINTSIZE_POWER)
    return POINTSIZE_MIN + (POINTSIZE_MAX - POINTSIZE_MIN) * p


# ---------- Plotting ----------
def volcano_plot(df, title, outbase,
                 fdr_thr=0.05, lfc_thr=1.0, annotate_top=200):
    """
    df must have: gene, logFC, FDR, logCPM, pct_cells_detected
    """
    d = (
        df[["gene", "logFC", "FDR", "logCPM", "pct_cells_detected"]]
        .dropna(subset=["gene", "logFC", "FDR"])
        .copy()
    )
    d["mlog10FDR"] = -np.log10(d["FDR"].clip(lower=np.finfo(float).tiny))

    sig_up = (d["FDR"] < fdr_thr) & (d["logFC"] >= lfc_thr)
    sig_dn = (d["FDR"] < fdr_thr) & (d["logFC"] <= -lfc_thr)
    ns = ~(sig_up | sig_dn)

    # For subtitle
    n_up = int(sig_up.sum())
    n_dn = int(sig_dn.sum())

    sizes = _sizes_from_detection(d["pct_cells_detected"].values)

    fig, ax = plt.subplots(figsize=(5.6, 5.6), dpi=220)
    ax.scatter(
        d.loc[ns, "logFC"],
        d.loc[ns, "mlog10FDR"],
        s=sizes[ns],
        c="#c7c7c7",
        alpha=0.45,
        linewidths=0,
        label="NS",
    )
    ax.scatter(
        d.loc[sig_dn, "logFC"],
        d.loc[sig_dn, "mlog10FDR"],
        s=sizes[sig_dn],
        c="#3b82f6",
        alpha=0.85,
        linewidths=0,
        label="Down",
    )
    ax.scatter(
        d.loc[sig_up, "logFC"],
        d.loc[sig_up, "mlog10FDR"],
        s=sizes[sig_up],
        c="#ef4444",
        alpha=0.85,
        linewidths=0,
        label="Up",
    )

    ax.axvline(lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axvline(-lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(-np.log10(fdr_thr), color="black", ls="--", lw=1, alpha=0.5)

    ax.set_xlabel(f"log2 fold-change ({COND_ALT} − {COND_REF})")
    ax.set_ylabel("-log10(FDR)")

    subtitle = f"Up: {n_up}   Down: {n_dn}   (FDR < {fdr_thr}, |log2FC| ≥ {lfc_thr})"
    ax.set_title(f"{title}\n{subtitle}", fontsize=9)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)

    # Gene labels
    texts = []
    if ANNOTATION_MODE == "smart":
        d_sig = d.loc[sig_up | sig_dn, ["gene", "logFC", "logCPM", "mlog10FDR"]].copy()
        top = _smart_annot_rank(d_sig, annotate_top)
    else:
        d_sig = d.loc[
            sig_up | sig_dn, ["gene", "logFC", "FDR", "logCPM", "mlog10FDR"]
        ].copy()
        if not d_sig.empty:
            d_sig["abs_logFC"] = d_sig["logFC"].abs()
            top = (
                d_sig.sort_values(["FDR", "abs_logFC"], ascending=[True, False])
                .head(annotate_top)
            )
        else:
            top = d_sig

    if not top.empty:
        y_lookup = d.set_index("gene")["mlog10FDR"].to_dict()
        for _, r in top.iterrows():
            x = r["logFC"]
            y = y_lookup.get(r["gene"], r.get("mlog10FDR", np.nan))
            texts.append(
                ax.text(x, y, r["gene"], fontsize=6, ha="center", va="bottom")
            )
        adjust_text(
            texts,
            ax=ax,
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
            arrowprops=dict(
                arrowstyle="-", lw=0.5, color="black", alpha=0.6
            ),
        )

    fig.tight_layout()
    fig.savefig(outbase + ".pdf")
    fig.savefig(outbase + ".png", dpi=220)
    plt.close(fig)


def ma_plot(df, title, outbase, lfc_thr=1.0, annotate_top=200):
    # df must have: gene, logCPM, logFC, pct_cells_detected
    x = df["logCPM"].values  # log2 CPM
    y = df["logFC"].values
    sizes = _sizes_from_detection(df["pct_cells_detected"].values)

    fig, ax = plt.subplots(figsize=(5.6, 5.6), dpi=220)
    ax.scatter(x, y, s=sizes, alpha=0.65, linewidths=0, c="#a8a8a8")
    ax.axhline(lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(-lfc_thr, color="black", ls="--", lw=1, alpha=0.5)
    ax.set_xlabel("mean expression (logCPM)")
    ax.set_ylabel(f"log2 fold-change ({COND_ALT} − {COND_REF})")
    ax.set_title(title, fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)

    # Gene labels
    texts = []
    if ANNOTATION_MODE == "smart":
        pick = df[["gene", "logCPM", "logFC"]].copy()
        top = _smart_annot_rank(pick, annotate_top)
    else:
        pick = df[["gene", "logCPM", "logFC"]].copy()
        if not pick.empty:
            pick["abs_logFC"] = pick["logFC"].abs()
            top = pick.sort_values("abs_logFC", ascending=False).head(annotate_top)
        else:
            top = pick

    if not top.empty:
        lut_x = df.set_index("gene")["logCPM"].to_dict()
        lut_y = df.set_index("gene")["logFC"].to_dict()
        for _, r in top.iterrows():
            gx = lut_x.get(r["gene"], r["logCPM"])
            gy = lut_y.get(r["gene"], r["logFC"])
            texts.append(
                ax.text(gx, gy, r["gene"], fontsize=6, ha="center", va="bottom")
            )
        adjust_text(
            texts,
            ax=ax,
            expand_points=(1.2, 1.4),
            expand_text=(1.1, 1.2),
            arrowprops=dict(
                arrowstyle="-", lw=0.5, color="black", alpha=0.6
            ),
        )

    fig.tight_layout()
    fig.savefig(outbase + ".pdf")
    fig.savefig(outbase + ".png", dpi=220)
    plt.close(fig)


def quick_pca_plot(counts_df, meta_df, title, out_png):
    try:
        from sklearn.decomposition import PCA
    except ImportError:
        return
    X = counts_df.T.values
    lib = X.sum(1, keepdims=True)
    Xnorm = np.log1p(1e6 * X / np.maximum(lib, 1))
    pc = PCA(n_components=2).fit_transform(Xnorm)
    fig, ax = plt.subplots(figsize=(4.5, 4))
    codes = pd.Categorical(meta_df["condition"]).codes
    ax.scatter(pc[:, 0], pc[:, 1], s=40, c=codes)
    for i, t in enumerate(meta_df.index):
        ax.text(pc[i, 0], pc[i, 1], t, fontsize=6, alpha=0.6)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


# ------------------- R-side via rpy2 -------------------
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

edgeR = importr("edgeR")
limma = importr("limma")
deseq2 = importr("DESeq2")

# Full edgeR pipeline for binary contrast: ALT - REF = cluster - rest
ro.r(f"""
run_edgeR_binary <- function(counts, group) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  if (nlevels(group) != 2) stop("Need exactly 2 groups: {COND_REF} and {COND_ALT}")
  library(edgeR); library(limma)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  y <- DGEList(counts=counts)
  keep <- filterByExpr(y, group=group)
  y <- y[keep,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  contrast <- makeContrasts(
    contrasts="{COND_ALT}-{COND_REF}",
    levels=design
  )
  qlf <- glmQLFTest(fit, contrast=contrast)
  tt  <- topTags(qlf, n=Inf)$table
  tt$gene <- rownames(tt)
  tt
}}
""")

# Effect-size-only path (no p-values), forced levels ALT - REF
ro.r(f"""
edgeR_effect_only_binary <- function(counts, group) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  if (nlevels(group) != 2) stop("Need exactly 2 groups: {COND_REF}, {COND_ALT}")
  y <- edgeR::DGEList(counts=counts)
  y <- edgeR::calcNormFactors(y)
  lcpms <- edgeR::cpm(y, log=TRUE, prior.count=0.5)
  mu_ref <- rowMeans(lcpms[, group=="{COND_REF}", drop=FALSE])
  mu_alt <- rowMeans(lcpms[, group=="{COND_ALT}", drop=FALSE])
  out  <- data.frame(
    gene   = rownames(counts),
    logFC  = mu_alt - mu_ref,           # cluster − rest
    logCPM = rowMeans(lcpms),
    PValue = NA_real_,
    FDR    = NA_real_
  )
  rownames(out) <- NULL
  out
}}
""")

# DESeq2 pipeline: same forced contrast ALT - REF
ro.r(f"""
run_DESeq2_binary <- function(counts, group) {{
  counts <- as.matrix(counts)
  group  <- factor(group, levels=c("{COND_REF}", "{COND_ALT}"))
  suppressPackageStartupMessages(library(DESeq2))
  coldata <- data.frame(group=group)
  rownames(coldata) <- colnames(counts)
  dds <- DESeqDataSetFromMatrix(countData=counts,
                                colData=coldata,
                                design=~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("group", "{COND_ALT}", "{COND_REF}"))
  logCPM <- log2(res$baseMean + 0.5)
  out <- data.frame(
    gene   = rownames(res),
    logFC  = res$log2FoldChange,
    logCPM = logCPM,
    PValue = res$pvalue,
    FDR    = res$padj
  )
  rownames(out) <- NULL
  out
}}
""")


def run_edgeR(counts_df, meta_df):
    # Ensure columns of counts_df are in the same order as rows of meta_df
    counts_df = counts_df.loc[:, meta_df.index]
    group = meta_df["condition"].astype(str).values

    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    rfun = ro.globalenv["run_edgeR_binary"]

    with localconverter(conv):
        dfR = rfun(counts_df, group)
        df = ro.conversion.rpy2py(dfR)

    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df


def effect_only_table(counts_df, meta_df):
    counts_df = counts_df.loc[:, meta_df.index]
    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    with localconverter(conv):
        dfR = ro.globalenv["edgeR_effect_only_binary"](
            counts_df, meta_df["condition"].astype(str).values
        )
        df = ro.conversion.rpy2py(dfR)
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df[["gene", "logFC", "logCPM", "PValue", "FDR"]]


def run_DESeq2(counts_df, meta_df):
    # Ensure columns of counts_df are in the same order as rows of meta_df
    counts_df = counts_df.loc[:, meta_df.index]
    group = meta_df["condition"].astype(str).values

    conv = default_converter + pandas2ri.converter + numpy2ri.converter
    rfun = ro.globalenv["run_DESeq2_binary"]

    with localconverter(conv):
        dfR = rfun(counts_df, group)
        df = ro.conversion.rpy2py(dfR)

    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    return df


# ------------------- Utility -------------------
def slugify(name):
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return s or "cluster"


# ------------------- Main Loop -------------------
def main():
    print(f"Reading AnnData from: {H5AD_PATH}")
    A0 = sc.read_h5ad(H5AD_PATH)
    if CELLTYPE_KEY not in A0.obs:
        raise KeyError(f"obs['{CELLTYPE_KEY}'] missing in AnnData.")

    print(f"Original shape: {A0.n_obs} cells × {A0.n_vars} genes")
    print("Cell-type counts (raw):")
    print(A0.obs[CELLTYPE_KEY].value_counts())

    A = filter_genes_drop_malat1_mt(A0)
    print(f"After gene filtering: {A.n_obs} cells × {A.n_vars} genes")

    counts_layer = pick_counts_layer(A)
    print(f"Using counts from layer = {counts_layer or 'X (integers)'}")

    # Detection stats over all cells (for point sizes)
    det_df, n_cells_total = _gene_detection(A, counts_layer)
    print(f"Computed detection stats over {n_cells_total} cells.")

    vc = A.obs[CELLTYPE_KEY].value_counts()
    print("\nWill process clusters with ≥ "
          f"{PB_MIN_CELLS_PER_GROUP} cells:\n{vc}")

    for ct, n_ct in vc.items():
        if n_ct < PB_MIN_CELLS_PER_GROUP:
            print(f"\n[SKIP] {ct}: only {n_ct} cells "
                  f"(< {PB_MIN_CELLS_PER_GROUP}).")
            continue

        print(f"\n=== {ct} ===  (n = {n_ct} cells)")

        cond_key = "__marker_cond__"
        cond = np.where(A.obs[CELLTYPE_KEY] == ct, COND_ALT, COND_REF)

        A_ct = A.copy()
        A_ct.obs[cond_key] = cond

        counts_df, meta_df = make_pseudobulk_binary(
            A_ct,
            counts_layer,
            cond_key=cond_key,
            min_cells=PB_MIN_CELLS_PER_CONDITION,
            reps_per_cond=REPLICATES_PER_CONDITION,
        )
        if counts_df is None:
            print(f"  -> Not enough cells for robust pseudobulk in {ct}; skipping.")
            continue

        ct_tag = slugify(ct)

        quick_pca_plot(
            counts_df,
            meta_df,
            title=f"{ct} vs rest pseudobulk PCA (log CPM)",
            out_png=os.path.join(OUTDIR, f"{ct_tag}_pseudobulk_PCA.png"),
        )

        base_edgeR = os.path.join(OUTDIR, f"{ct_tag}_pb_edgeR")
        base_DESeq2 = os.path.join(OUTDIR, f"{ct_tag}_pb_DESeq2")

        if have_replicates(meta_df):
            # ---------- edgeR QLF ----------
            title_edgeR = f"{ct} — {COND_ALT} vs {COND_REF} (edgeR QLF)"
            df_edgeR = run_edgeR(counts_df, meta_df)
            df_edgeR = df_edgeR.merge(
                det_df, left_on="gene", right_index=True, how="left"
            )
            keep_cols = [
                c
                for c in [
                    "gene",
                    "logFC",
                    "logCPM",
                    "F",
                    "PValue",
                    "FDR",
                    "n_cells_detected",
                    "pct_cells_detected",
                ]
                if c in df_edgeR.columns
            ]
            df_edgeR = df_edgeR[keep_cols]
            df_edgeR.to_csv(f"{base_edgeR}.csv", index=False)
            volcano_plot(
                df_edgeR,
                title=title_edgeR,
                outbase=f"{base_edgeR}__volcano",
                fdr_thr=FDR_THR,
                lfc_thr=LFC_THR,
                annotate_top=TOP_ANNOTATE,
            )

            # ---------- DESeq2 ----------
            title_DE = f"{ct} — {COND_ALT} vs {COND_REF} (DESeq2)"
            df_DE = run_DESeq2(counts_df, meta_df)
            df_DE = df_DE.merge(
                det_df, left_on="gene", right_index=True, how="left"
            )
            keep_cols_DE = [
                c
                for c in [
                    "gene",
                    "logFC",
                    "logCPM",
                    "PValue",
                    "FDR",
                    "n_cells_detected",
                    "pct_cells_detected",
                ]
                if c in df_DE.columns
            ]
            df_DE = df_DE[keep_cols_DE]
            df_DE.to_csv(f"{base_DESeq2}.csv", index=False)
            volcano_plot(
                df_DE,
                title=title_DE,
                outbase=f"{base_DESeq2}__volcano",
                fdr_thr=FDR_THR,
                lfc_thr=LFC_THR,
                annotate_top=TOP_ANNOTATE,
            )

        else:
            print(
                f"  -> {ct}: not enough pseudobulk replicates per group "
                "for edgeR/DESeq2; reporting effect sizes only."
            )
            df = effect_only_table(counts_df, meta_df)
            df = df.merge(det_df, left_on="gene", right_index=True, how="left")
            base = base_edgeR
            df.to_csv(base + "__effect_only.csv", index=False)
            ma_plot(
                df,
                title=f"{ct} — {COND_ALT} vs {COND_REF}",
                outbase=base + "__effect_only_MA",
                lfc_thr=LFC_THR,
                annotate_top=TOP_ANNOTATE,
            )

    print("\nDone. Outputs in:", OUTDIR)


if __name__ == "__main__":
    main()






# DOTPLOT


# ----------------- CONFIG -----------------
celltype_key       = "leiden_merged_type"
dge_dir            = "figures_cluster_markers_pseudobulk"  # where <slug>_pb_DESeq2.csv live
top_n              = 2          # top genes per cell type BEFORE dedup
max_fdr            = 0.00001       # FDR cutoff
min_abs_lfc        = 1.5        # |log2FC| cutoff
min_logCPM         = 0.5        # mean expression cutoff (DESeq2 logCPM)
min_pct_detected   = 0.02       # fraction of cells (0–1) with >0 counts


def slugify(name: str) -> str:
    """Must match how you named <cluster>_pb_DESeq2.csv."""
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return s or "cluster"


# ensure categorical order is defined
adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
clusters = list(adata.obs[celltype_key].cat.categories)

var_names = {}       # dict for scanpy.pl.dotplot: {cluster_label: [genes]}
used_genes = set()   # global set to de-duplicate genes across clusters

for ct in clusters:
    slug = slugify(ct)
    csv_path = os.path.join(dge_dir, f"{slug}_pb_DESeq2.csv")
    if not os.path.exists(csv_path):
        print(f"[skip] {ct}: {csv_path} not found")
        continue

    df = pd.read_csv(csv_path)
    if "gene" not in df.columns:
        print(f"[skip] {ct}: no 'gene' column in {csv_path}")
        continue

    # numeric columns + convenience abs_logFC
    for c in ["logFC", "logCPM", "PValue", "FDR", "pct_cells_detected"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["logFC"])
    df["abs_logFC"] = df["logFC"].abs()

    # ----------------- FILTERS (toggles) -----------------
    if max_fdr is not None and "FDR" in df.columns:
        df = df[df["FDR"] <= max_fdr]

    if min_abs_lfc is not None:
        df = df[df["abs_logFC"] >= min_abs_lfc]

    if min_logCPM is not None and "logCPM" in df.columns:
        df = df[df["logCPM"] >= min_logCPM]

    if min_pct_detected is not None and "pct_cells_detected" in df.columns:
        df = df[df["pct_cells_detected"] >= min_pct_detected]

    if df.empty:
        print(f"[skip] {ct}: all genes filtered out")
        continue

    # ----------------- RANKING (per cluster) -----------------
    # smaller FDR is better, larger abs_logFC is better
    # ranking within cluster: FDR then |logFC|
    if "FDR" not in df.columns:
        df["FDR"] = np.nan
    
    df = df.sort_values(
        ["FDR", "abs_logFC"],
        ascending=[True, False],
    )

 # pick TOP_N genes per cell type, scanning deeper to avoid global duplicates
    genes_for_ct = []
    dropped_dups = []

    for g in df["gene"].astype(str).tolist():
        if g not in adata.var_names:
            continue
        if g in used_genes:
            dropped_dups.append(g)
            continue
        genes_for_ct.append(g)
        used_genes.add(g)
        if len(genes_for_ct) == top_n:
            break

    if dropped_dups:
        print(f"[dedup] {ct}: skipped already-used genes (showing up to 10): {dropped_dups[:10]}")

    if len(genes_for_ct) < top_n:
        raise RuntimeError(
            f"{ct}: only found {len(genes_for_ct)}/{top_n} unique genes after scanning DESeq2 list. "
            f"Consider relaxing filters (max_fdr/min_abs_lfc/min_logCPM/min_pct_detected) or disabling global dedup."
        )

    var_names[ct] = genes_for_ct
    print(f"[OK] {ct}: using {len(genes_for_ct)} genes -> {genes_for_ct}")


print("\nFinal var_names dict going into scanpy.pl.dotplot:")
for ct, genes in var_names.items():
    print(f"  {ct}: {genes}")

sc.tl.dendrogram(
    adata,
    groupby=celltype_key,
    use_rep="X_pca",   # or None to use .X, or "X_umap" etc. if you prefer
)

# ---- the actual dotplot ----
sc.pl.dotplot(
    adata,
    var_names,               # dict: {cluster: [genes]}
    groupby=celltype_key,    # rows = your cell types
    standard_scale="var",    # scale per gene across groups (nice for markers)
    swap_axes=False,         # genes on x, clusters on y; flip if you prefer
    figsize=(8, 4),
    cmap="Reds",             # or e.g. "viridis"
    dendrogram=True,         # <---- this turns the tree on
    save="AllCellTypes_PseudobulkDESeq2.pdf",               # set to ".pdf" / ".png" to auto-save
)























#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.sparse as sp
import snapatac2 as snap
from matplotlib.gridspec import GridSpec

# ----------------- CONFIG -----------------

# Paths – adjust RNA_PATH if your RNA AnnData has a different name
RNA_PATH = "./Tonsil_FULL_RNA_ANNOTATED.h5ad"
AC_PATH  = "Sc_VTHumanTonsil_ac_processed.h5ad"   # H3K27ac fragments/peaks
ME3_PATH = "Sc_VTHumanTonsil_me3_processed.h5ad"  # H3K27me3 fragments/peaks
DGE_DIR  = "figures_cluster_markers_pseudobulk"   # <slug>_pb_DESeq2.csv directory

celltype_key       = "leiden_merged_type"
top_n              = 2          # top genes per cell type BEFORE dedup
max_fdr            = 0.00001       # FDR cutoff
min_abs_lfc        = 1.5        # |log2FC| cutoff
min_logCPM         = 0.5        # mean expression cutoff (DESeq2 logCPM)
min_pct_detected   = 0.02       # fraction of cells (0–1) with >0 counts

# Gene-activity window ±bp around gene TSS / body
PROMOTER_WINDOW = 5000

# Aggregation + normalization toggles
AGG_FUNC   = "mean"         # options: "mean", "median"
NORM_MODE  = "zscore_gene"         # options: "none", "minmax_gene", "zscore_gene"

OUTDIR = f"heatmaps_top{top_n}_markers_RNA_AC_ME3_groupmeans"
os.makedirs(OUTDIR, exist_ok=True)

# Genome annotation for gene activity (hg38)
try:
    GENE_ANNO = snap.genome.hg38
except Exception:
    GENE_ANNO = None  # will error later if not set

# ----------------- HELPERS -----------------

def slugify(name: str) -> str:
    """Must match how you named <cluster>_pb_DESeq2.csv."""
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return s or "cluster"

def _resolve_anno_for_snap(gene_anno):
    if gene_anno is None:
        raise ValueError(
            "GENE_ANNO is None. Set it to snap.genome.hg38 or a GTF/GFF path."
        )
    # snap.genome.hg38 has attribute .annotation, whereas a GTF path is used directly
    return getattr(gene_anno, "annotation", gene_anno)

def make_gene_activity(adata_mark, gene_anno, window_bp, layer_name="log1p"):
    """
    Build gene-activity matrix from a fragment/peak-level AnnData.

    - upstream/downstream = ±window_bp around genes
    - include_gene_body=True
    - id_type="gene" → gene symbols in var_names
    """
    ann = _resolve_anno_for_snap(gene_anno)
    ga = snap.pp.make_gene_matrix(
        adata_mark,
        ann,
        upstream=window_bp,
        downstream=window_bp,
        include_gene_body=True,
        id_type="gene",
    )
    # carry over useful obs annotations
    for col in ["sample", "leiden_merged_type", "leiden"]:
        if col in adata_mark.obs:
            ga.obs[col] = adata_mark.obs[col].copy()

    # add log1p layer
    X = ga.X
    if sp.issparse(X):
        Xlog = X.copy()
        Xlog.data = np.log1p(Xlog.data)
    else:
        Xlog = np.log1p(X)
    ga.layers[layer_name] = Xlog
    return ga

def compute_group_means(ad, genes, groupby, layer=None, group_order=None, agg_func="mean"):
    """
    Return DataFrame of aggregated expression per group:
      index   = groups (e.g. cell types)
      columns = genes (in the given order)

    agg_func: "mean" or "median"
    """
    df = sc.get.obs_df(ad, keys=genes, layer=layer)
    df[groupby] = ad.obs[groupby].astype(str).values

    if agg_func == "median":
        g = df.groupby(groupby).median(numeric_only=True)
    else:
        g = df.groupby(groupby).mean(numeric_only=True)

    # enforce ordering
    if group_order is not None:
        g = g.reindex(group_order)

    # ensure column order is exactly 'genes'
    g = g[genes]

    # replace any NaNs (e.g. empty groups) with 0 for plotting
    g = g.fillna(0.0)

    return g

def _normalize_matrix(M, mode="none", z_clip=1.5):
    """
    Normalize matrix M (genes × groups) according to mode.
    """
    if mode == "none":
        return M, "mean expression"

    if mode == "minmax_gene":
        M_scaled = np.zeros_like(M, dtype=float)
        for i in range(M.shape[0]):
            row = M[i, :]
            rmin = float(np.min(row))
            rmax = float(np.max(row))
            if rmax > rmin:
                M_scaled[i, :] = (row - rmin) / (rmax - rmin)
            else:
                M_scaled[i, :] = 0.0
        return M_scaled, "per-gene min–max scaled mean"

    if mode == "zscore_gene":
        M_scaled = np.zeros_like(M, dtype=float)
        for i in range(M.shape[0]):
            row = M[i, :]
            mu  = float(np.mean(row))
            sd  = float(np.std(row))
            if sd > 0:
                z = (row - mu) / sd
                if z_clip is not None:
                    z = np.clip(z, -z_clip, z_clip)
                M_scaled[i, :] = z
            else:
                M_scaled[i, :] = 0.0
        return M_scaled, f"per-gene z-score (clip ±{z_clip})"

    # fallback
    return M, "mean expression"

def plot_group_mean_heatmap(
    mean_df,
    title,
    fname,
    cmap="viridis",
    norm_mode="none",
    var_group_positions=None,
    var_group_labels=None,
    x_label=None,
    agg_label="mean",
):
    """
    mean_df: DataFrame (index = groups, columns = genes)

    Plots genes on rows, groups on columns.
    Supports several normalization modes.
    Uses separate colorbar axis so it does not overlap right-side labels.
    """
    groups = list(mean_df.index)      # e.g. clusters / cell types
    genes  = list(mean_df.columns)    # genes_ordered

    # matrix for plotting: genes × groups
    M = mean_df.to_numpy().T         # shape: (n_genes, n_groups)

    # normalization
    M_scaled, norm_label = _normalize_matrix(M, mode=norm_mode)

    n_genes, n_groups = M_scaled.shape

    fig_w = max(4.5, 0.5 * n_groups + 2.5)
    fig_h = max(4.5, 0.25 * n_genes + 2.5)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)
    gs = GridSpec(nrows=1, ncols=2, width_ratios=[20, 1], wspace=0.15, figure=fig)

    ax  = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    im = ax.imshow(M_scaled, aspect="auto", interpolation="nearest", cmap=cmap)

    # axes / labels
    ax.set_xticks(np.arange(n_groups))
    ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(genes, fontsize=7)

    if x_label is not None:
        ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel("marker genes", fontsize=9)
    ax.set_title(title, fontsize=10)

    # horizontal lines between marker blocks (per original cluster)
    if var_group_positions is not None and len(var_group_positions) > 1:
        for (_, end) in var_group_positions[:-1]:
            ax.axhline(end + 0.5, color="white", linewidth=0.5)

    # optional left-side group labels (cluster-of-origin for gene blocks) on a twin y-axis
    if var_group_positions is not None and var_group_labels is not None:
        centers = [0.5 * (s + e) for (s, e) in var_group_positions]
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(centers)
        ax2.set_yticklabels(var_group_labels, fontsize=8)
        ax2.tick_params(axis="y", length=0)
        ax2.set_ylabel("cluster of origin", fontsize=9)

    # colorbar in a dedicated axis (no overlap with right labels)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(f"{norm_label} (agg: {agg_label})", fontsize=8)

    pdf = os.path.join(OUTDIR, f"{fname}.pdf")
    png = os.path.join(OUTDIR, f"{fname}.png")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[saved] {pdf}")

# ----------------- LOAD DATA -----------------

print(f"[INFO] Loading RNA AnnData from: {RNA_PATH}")
adata_rna = sc.read_h5ad(RNA_PATH)

print(f"[INFO] Loading H3K27ac from:  {AC_PATH}")
print(f"[INFO] Loading H3K27me3 from: {ME3_PATH}")
adata_ac_raw  = sc.read_h5ad(AC_PATH)
adata_me3_raw = sc.read_h5ad(ME3_PATH)

# ----------------- ALIGN CELLS ACROSS MODALITIES -----------------

common_cells = (
    adata_rna.obs_names
    .intersection(adata_ac_raw.obs_names)
    .intersection(adata_me3_raw.obs_names)
)
print(f"[INFO] Common cells across RNA / H3K27ac / H3K27me3: {len(common_cells)}")

adata     = adata_rna[common_cells].copy()
adata_ac  = adata_ac_raw[common_cells].copy()
adata_me3 = adata_me3_raw[common_cells].copy()

# ensure cell type column exists in histone objects
if celltype_key not in adata.obs:
    raise KeyError(f"obs['{celltype_key}'] is missing from RNA AnnData.")

for mod in (adata_ac, adata_me3):
    if celltype_key not in mod.obs:
        mod.obs[celltype_key] = adata.obs[celltype_key].reindex(mod.obs_names)

# make celltype categorical and fix order
adata.obs[celltype_key] = adata.obs[celltype_key].astype("category")
clusters = list(adata.obs[celltype_key].cat.categories)
print("[INFO] cell types (order):", clusters)

# ----------------- TOP3 MARKERS FROM DESeq2 CSVs -----------------

var_names = {}       # {cluster_label: [genes]}
used_genes = set()   # global set to avoid duplicate genes across clusters

for ct in clusters:
    slug = slugify(ct)
    csv_path = os.path.join(DGE_DIR, f"{slug}_pb_DESeq2.csv")
    if not os.path.exists(csv_path):
        print(f"[skip] {ct}: {csv_path} not found")
        continue

    df = pd.read_csv(csv_path)
    if "gene" not in df.columns:
        print(f"[skip] {ct}: no 'gene' column in {csv_path}")
        continue

    # numeric columns + abs_logFC
    for c in ["logFC", "logCPM", "PValue", "FDR", "pct_cells_detected"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["logFC"])
    df["abs_logFC"] = df["logFC"].abs()

    # filters / toggles
    if max_fdr is not None and "FDR" in df.columns:
        df = df[df["FDR"] <= max_fdr]

    if min_abs_lfc is not None:
        df = df[df["abs_logFC"] >= min_abs_lfc]

    if min_logCPM is not None and "logCPM" in df.columns:
        df = df[df["logCPM"] >= min_logCPM]

    if min_pct_detected is not None and "pct_cells_detected" in df.columns:
        df = df[df["pct_cells_detected"] >= min_pct_detected]

    if df.empty:
        print(f"[skip] {ct}: all genes filtered out after thresholds")
        continue

    # ranking within cluster: FDR then |logFC|
    if "FDR" not in df.columns:
        df["FDR"] = np.nan
    
    df = df.sort_values(
        ["FDR", "abs_logFC"],
        ascending=[True, False],
    )

 # pick TOP_N genes per cell type, scanning deeper to avoid global duplicates
    genes_for_ct = []
    dropped_dups = []

    for g in df["gene"].astype(str).tolist():
        if g not in adata.var_names:
            continue
        if g in used_genes:
            dropped_dups.append(g)
            continue
        genes_for_ct.append(g)
        used_genes.add(g)
        if len(genes_for_ct) == top_n:
            break

    if dropped_dups:
        print(f"[dedup] {ct}: skipped already-used genes (showing up to 10): {dropped_dups[:10]}")

    if len(genes_for_ct) < top_n:
        raise RuntimeError(
            f"{ct}: only found {len(genes_for_ct)}/{top_n} unique genes after scanning DESeq2 list. "
            f"Consider relaxing filters (max_fdr/min_abs_lfc/min_logCPM/min_pct_detected) or disabling global dedup."
        )

    var_names[ct] = genes_for_ct
    print(f"[OK] {ct}: using {len(genes_for_ct)} genes -> {genes_for_ct}")


print("\nFinal var_names dict (top3 markers per cell type):")
for ct, genes in var_names.items():
    print(f"  {ct}: {genes}")

# ----------------- BUILD GENE-ACTIVITY MATRICES -----------------

print("\n[INFO] Building gene-activity matrices for H3K27ac / H3K27me3 ...")

ga_ac  = make_gene_activity(adata_ac,  GENE_ANNO, PROMOTER_WINDOW)
ga_me3 = make_gene_activity(adata_me3, GENE_ANNO, PROMOTER_WINDOW)

# Ensure celltype_key is present and categorical in GA matrices
for ga in (ga_ac, ga_me3):
    if celltype_key not in ga.obs:
        ga.obs[celltype_key] = adata.obs[celltype_key].reindex(ga.obs_names)
    ga.obs[celltype_key] = ga.obs[celltype_key].astype("category")
    ga.obs[celltype_key] = ga.obs[celltype_key].cat.set_categories(clusters)

# ----------------- GENE LIST & GROUPING FOR HEATMAPS -----------------

# flatten markers in cluster order
genes_flat = []
for ct in clusters:
    if ct not in var_names:
        continue
    genes_flat.extend(var_names[ct])

# uniqueness preserving order
seen = set()
genes_flat_unique = []
for g in genes_flat:
    if g not in seen:
        seen.add(g)
        genes_flat_unique.append(g)

print("\nTop-marker genes (before GA intersection):")
print(genes_flat_unique)

# intersect with genes present in RNA + gene-activity matrices
var_rna  = set(adata.var_names)
var_ac   = set(ga_ac.var_names)
var_me3  = set(ga_me3.var_names)

genes_shared = []
missing_ac = []
missing_me3 = []

for g in genes_flat_unique:
    in_rna = g in var_rna
    in_ac  = g in var_ac
    in_me3 = g in var_me3
    if in_rna and in_ac and in_me3:
        genes_shared.append(g)
    else:
        if not in_ac:
            missing_ac.append(g)
        if not in_me3:
            missing_me3.append(g)

if missing_ac:
    print(f"[WARN] genes missing in H3K27ac GA: {missing_ac}")
if missing_me3:
    print(f"[WARN] genes missing in H3K27me3 GA: {missing_me3}")

if not genes_shared:
    raise RuntimeError("No top marker genes are present in ALL of RNA/GA_ac/GA_me3.")

print("\nGenes present in ALL modalities (base list for heatmaps):")
print(genes_shared)

# Filter var_names per cluster to only keep genes that survived the shared filter
filtered_var_names = {}
for ct in clusters:
    if ct not in var_names:
        continue
    kept = [g for g in var_names[ct] if g in genes_shared]
    if kept:
        filtered_var_names[ct] = kept

# Rebuild final ordered gene list + var_group_positions/labels
genes_ordered = []
var_group_positions = []
var_group_labels    = []

start = 0
for ct in clusters:
    genes_ct = filtered_var_names.get(ct, [])
    if not genes_ct:
        continue
    genes_ordered.extend(genes_ct)
    end = start + len(genes_ct) - 1
    var_group_positions.append((start, end))
    var_group_labels.append(ct)
    start = end + 1

if not genes_ordered:
    raise RuntimeError("No genes left after intersecting with GA var_names.")

print("\nFinal gene order for heatmaps (per cluster of origin):")
for ct, (g_start, g_end) in zip(var_group_labels, var_group_positions):
    print(f"  {ct}: {genes_ordered[g_start:g_end+1]}")

# ----------------- LAYERS -----------------

# RNA: assume adata.X is already normalized / log1p
RNA_LAYER = None
# GA matrices use 'log1p' layer
AC_LAYER  = "log1p"
ME3_LAYER = "log1p"

# ----------------- GROUP-MEAN MATRICES -----------------

rna_mean  = compute_group_means(
    ad=adata,
    genes=genes_ordered,
    groupby=celltype_key,
    layer=RNA_LAYER,
    group_order=clusters,
    agg_func=AGG_FUNC,
)
ac_mean   = compute_group_means(
    ad=ga_ac,
    genes=genes_ordered,
    groupby=celltype_key,
    layer=AC_LAYER,
    group_order=clusters,
    agg_func=AGG_FUNC,
)
me3_mean  = compute_group_means(
    ad=ga_me3,
    genes=genes_ordered,
    groupby=celltype_key,
    layer=ME3_LAYER,
    group_order=clusters,
    agg_func=AGG_FUNC,
)

print("[INFO] rna_mean shape:",  rna_mean.shape)
print("[INFO] ac_mean shape:",   ac_mean.shape)
print("[INFO] me3_mean shape:",  me3_mean.shape)

# ----------------- GROUP-MEAN HEATMAPS -----------------

groupby_label = celltype_key  # e.g. "leiden_merged_type"

# RNA
plot_group_mean_heatmap(
    mean_df=rna_mean,
    title=f"RNA: {AGG_FUNC} expression per cell type (top3 DESeq2 markers)",
    fname="groupmean_RNA_top3_markers_by_celltype",
    cmap="viridis",
    norm_mode=NORM_MODE,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    x_label=groupby_label,
    agg_label=AGG_FUNC,
)

# H3K27ac GA
plot_group_mean_heatmap(
    mean_df=ac_mean,
    title=f"H3K27ac gene activity (±{PROMOTER_WINDOW//1000}kb): {AGG_FUNC} per cell type",
    fname=f"groupmean_H3K27ac_top{top_n}_markers_by_celltype",
    cmap="magma",
    norm_mode=NORM_MODE,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    x_label=groupby_label,
    agg_label=AGG_FUNC,
)

# H3K27me3 GA
plot_group_mean_heatmap(
    mean_df=me3_mean,
    title=f"H3K27me3 gene activity (±{PROMOTER_WINDOW//1000}kb): {AGG_FUNC} per cell type",
    fname=f"groupmean_H3K27me3_top{top_n}_markers_by_celltype",
    cmap="magma",
    norm_mode=NORM_MODE,
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    x_label=groupby_label,
    agg_label=AGG_FUNC,
)

print(f"\nDone. Group-mean heatmaps written to: {OUTDIR}")




#CheckGenesonUMAP

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import snapatac2 as snap
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap

# ------------------- Config -------------------
RNA_PATH = "Tonsil_rna_Overlap_Reclust.h5ad"
AC_PATH  = "Tonsil_ac_Overlap_Reclust.h5ad"
ME3_PATH = "Tonsil_me3_Overlap_Reclust.h5ad"

OUTDIR = "Tonsil_SpeGenes"
os.makedirs(OUTDIR, exist_ok=True)

# Genes of interest (uppercased to avoid common alias case issues)
GENES = ["BCL6", "IRF4", "XBP1", "SDC1","MYC","CCNB1","FOXM1"]

# Column to color/clump by if present (fallback handled below)
PREFERRED_GROUPBY = "leiden_merged_type"

# Window size for histone gene activity (in bp). Set to 5000 for ±5 kb.
PROMOTER_WINDOW = 5000

# --- UMAP color scaling ---
# smart vmax: 3× non-zero median, capped at 99th percentile; fallback 2× non-zero mean
MEDIAN_FACTOR   = 3.0
MEAN_FACTOR     = 2.0
CAP_PERCENTILE  = 99.0
# tiny epsilon so 0's render "under" (light gray)
ZERO_EPS = 1e-12


# ------------------- NEW: Analysis parameters -------------------
# Choose the RNA layer used for heatmaps/UMAPs:
#   - "auto": picks the first found from RNA_LAYER_PRIORITY
#   - or a specific layer name present in adata_rna.layers (e.g. "SCT_data", "SCT_counts", "counts", "log1p_norm")
RNA_LAYER_FOR_HEATMAP = "auto"
RNA_LAYER_PRIORITY    = ["counts", "SCT_data", "SCT_counts", "log1p_norm", None]  # None → use X
HTM_LAYER  = "counts" #"log1p"
rna_layer_hard = "counts"
# Binary call rules for ordering cells by overlap pattern in the heatmaps
# Allowed: "auto", "gt0", "ge1", "ge2", "z>1"
RNA_BIN_RULE  = "auto"        # auto = gt0 if layer ~log; otherwise ge1
AC_BIN_RULE   = "ge1"         # GA counts ≥1 → positive
ME3_BIN_RULE  = "ge1"

# How to display the overlap patterns (priority in sorting inside each group)
# (bits: RNA<<2 | AC<<1 | ME3<<0)
PATTERN_ORDER = [
    0b111,  # RNA+AC+ME3
    0b110,  # RNA+AC
    0b100,  # RNA only
    0b010,  # AC only
    0b001,  # ME3 only
    0b011,  # AC+ME3
    0b101,  # RNA+ME3
    0b000,  # none
]


# Genome annotation for gene activity (GTF/GFF).
try:
    GENE_ANNO = snap.genome.hg38  # requires that hg38 is available to SnapATAC2
except Exception:
    GENE_ANNO = None  # set manually if needed

# ------------------- I/O -------------------
adata_rna = sc.read_h5ad(RNA_PATH)
adata_ac  = sc.read_h5ad(AC_PATH)
adata_me3 = sc.read_h5ad(ME3_PATH)

# ------------------- Helpers -------------------
def resolve_rna_layer(adata, choice="auto", priority=None):
    """Return an existing layer name or None to use X."""
    if choice != "auto":
        if choice in adata.layers:
            return choice
        print(f"[WARN] RNA layer '{choice}' not found; falling back to auto.")
    priority = (priority or []) + [None]
    for lay in priority:
        if lay is None:
            return None
        if lay in adata.layers:
            return lay
    return None

def _binarize_vector(v: np.ndarray, mode: str, layer_name: str | None = None) -> np.ndarray:
    """Return boolean positives according to mode."""
    m = (mode or "auto").lower()
    if m == "auto":
        # treat 'SCT_*', 'log1p*', 'logcounts' as log→use >0; else counts→use ≥1
        if layer_name and (layer_name.lower().startswith("sct_")
                           or layer_name.lower().startswith("log1p")
                           or layer_name.lower().startswith("logcounts")):
            return v > 0
        return v >= 1
    if m == "gt0":  return v > 0
    if m == "ge1":  return v >= 1
    if m == "ge2":  return v >= 2
    if m == "z>1":
        mu = v.mean(); sd = v.std()
        if sd == 0: return np.zeros_like(v, dtype=bool)
        return (v - mu) / sd > 1.0
    raise ValueError(f"Unknown binarization mode: {mode}")

def _ga_gene_vectors(ga_adata, gene, order_index=None):
    """
    For a GA AnnData (counts in X, log1p in layers['log1p']):
    return counts_vec, log1p_vec (both aligned to RNA group order if order_index provided).
    """
    try:
        j = ga_adata.var_names.get_loc(gene)
    except KeyError:
        n = ga_adata.n_obs
        zeros = np.zeros(n, dtype=float)
        return (zeros if order_index is None else zeros[order_index],
                zeros if order_index is None else zeros[order_index])
    Xcnt = ga_adata.X
    v_cnt = Xcnt.getcol(j).toarray().ravel() if sp.issparse(Xcnt) else np.asarray(Xcnt[:, j]).ravel()
    Xlog = ga_adata.layers["log1p"]
    v_log = Xlog.getcol(j).toarray().ravel() if sp.issparse(Xlog) else np.asarray(Xlog[:, j]).ravel()
    if order_index is None:
        return v_cnt, v_log
    return v_cnt[order_index], v_log[order_index]

def _order_cells_by_pattern(
    rna_vals, ac_cnt, me3_cnt, group_index: np.ndarray,
    rna_layer_name: str | None,
    pattern_order=PATTERN_ORDER,
    rna_rule=RNA_BIN_RULE, ac_rule=AC_BIN_RULE, me3_rule=ME3_BIN_RULE
) -> np.ndarray:
    """
    Within a given group (indices into the global order), return the indices
    sorted by (pattern priority, then mean of the three values desc).
    """
    rna_g  = rna_vals[group_index]
    ac_g   = ac_cnt[group_index]
    me3_g  = me3_cnt[group_index]

    b_rna  = _binarize_vector(rna_g,  rna_rule, layer_name=rna_layer_name)
    b_ac   = _binarize_vector(ac_g,   ac_rule,  layer_name=None)      # counts by default
    b_me3  = _binarize_vector(me3_g,  me3_rule, layer_name=None)

    pattern = (b_rna.astype(int) << 2) | (b_ac.astype(int) << 1) | (b_me3.astype(int) << 0)
    pri_map = {p:i for i,p in enumerate(pattern_order)}
    pri = np.array([pri_map.get(int(x), len(pattern_order)) for x in pattern], dtype=int)

    # tie-breaker: stronger average signal first
    mean_sig = (rna_g + ac_g + me3_g) / 3.0

    # stable sort by priority then -mean_sig
    local_order = np.lexsort(( -mean_sig, pri ))
    return group_index[local_order]

def ensure_umap(adata, n_neighbors=15, n_pcs=50, key="X_umap"):
    """Compute UMAP if not present."""
    if key in adata.obsm:
        return
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=n_pcs)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca" if "X_pca" in adata.obsm else None)
    sc.tl.umap(adata)  # writes to .obsm['X_umap']

def copy_umap_from_reference(ref, target, key="X_umap"):
    """
    Copy UMAP coordinates from ref -> target (for shared barcodes).
    Adds a shallow copy of the embedding into target.obsm[key].
    """
    common = ref.obs_names.intersection(target.obs_names)
    if len(common) == 0:
        raise ValueError("No shared barcodes between reference and target.")
    target._inplace_subset_obs(common)
    ref_umap = pd.DataFrame(ref.obsm[key], index=ref.obs_names, columns=["UMAP1", "UMAP2"])
    target.obsm[key] = ref_umap.loc[common].to_numpy()

def add_obs_from_rna(reference_rna, target, columns=("leiden_merged_type", "leiden")):
    """Bring selected annotation columns from RNA to another AnnData by reindexing."""
    for col in columns:
        if col in reference_rna.obs:
            target.obs[col] = reference_rna.obs[col].reindex(target.obs_names)

def pick_groupby(adata):
    if PREFERRED_GROUPBY in adata.obs:
        return PREFERRED_GROUPBY
    if "leiden" in adata.obs:
        return "leiden"
    # otherwise make a dummy single group
    adata.obs["_all"] = "all"
    return "_all"

def ensure_counts_layer_for_rna(adata):
    """Prefer raw counts for RNA violins."""
    if "counts" in adata.layers:
        return "counts", None  # (layer, use_raw)
    if adata.raw is not None:
        return None, True
    return None, False

def upper_intersection(list_like, var_names):
    # case-insensitive intersection of genes of interest with var_names
    wanted = {g.upper() for g in list_like}
    mapping = {v.upper(): v for v in var_names}
    return [mapping[g] for g in wanted if g in mapping]

def _resolve_anno_for_snap(gene_anno):
    """Allow passing a Genome or a path-like/string."""
    if gene_anno is None:
        raise ValueError("Please set GENE_ANNO to your hg38 GTF/GFF path (or ensure snap.genome.hg38 is available).")
    return getattr(gene_anno, "annotation", gene_anno)

def make_gene_activity(adata_mark, gene_anno, window_bp):
    """
    Build gene activity matrix (counts in gene body ± window).
    Returns a NEW AnnData with genes as .var_names, cells unchanged.
    """
    ann = _resolve_anno_for_snap(gene_anno)
    ga = snap.pp.make_gene_matrix(
        adata_mark,
        ann,
        upstream=window_bp,
        downstream=window_bp,
        include_gene_body=True,
        id_type="gene",
    )
    # carry over group labels for plotting
    for col in [PREFERRED_GROUPBY, "leiden"]:
        if col in adata_mark.obs:
            ga.obs[col] = adata_mark.obs[col].copy()
    # add log1p layer for coloring UMAPs
    X = ga.X
    if sp.issparse(X):
        X_log = X.copy()
        X_log.data = np.log1p(X_log.data)
    else:
        X_log = np.log1p(X)
    ga.layers["log1p"] = X_log
    ga.layers["counts"] = X.copy()
    return ga

# ---- Zero=gray then viridis colormap ----
def make_zero_gray_viridis():
    base = mpl.cm.get_cmap('viridis', 256)
    colors = base(np.linspace(0, 1, 256))
    cmap = ListedColormap(colors)
    cmap.set_under('lightgray')  # zeros (values < vmin) appear light gray
    return cmap

ZERO_GRAY_VIRIDIS = make_zero_gray_viridis()

def _get_gene_values(adata, gene, layer=None):
    """Return 1D np.array of per-cell values for a gene from X or a given layer; empty array if missing."""
    try:
        idx = adata.var_names.get_loc(gene)
    except KeyError:
        return np.array([])
    Xsrc = adata.layers[layer] if layer is not None else adata.X
    if sp.issparse(Xsrc):
        return Xsrc.getcol(idx).toarray().ravel()
    else:
        return np.asarray(Xsrc[:, idx]).ravel()

def _smart_vmax(values, med_factor=MEDIAN_FACTOR, mean_factor=MEAN_FACTOR, cap_pct=CAP_PERCENTILE):
    """Compute vmax = min(cap_percentile, med_factor*median_nonzero) with sensible fallbacks."""
    if values.size == 0:
        return None
    nz = values[values > 0]
    if nz.size == 0:
        return None
    med = np.median(nz)
    if med > 0:
        vmax_candidate = med_factor * med
    else:
        m = np.mean(nz)
        vmax_candidate = mean_factor * m if m > 0 else None
    if vmax_candidate is None:
        return None
    cap = np.percentile(nz, cap_pct)
    vmax = float(max(1e-9, min(vmax_candidate, cap))) if np.isfinite(cap) else float(max(1e-9, vmax_candidate))
    return vmax

def save_umap(adata, color, fname, legend_loc="on data", layer=None,
              vmin=0, vmax=None, color_map=None):
    """Generic UMAP saver (supports 'layer' & custom colormap)."""
    sc.pl.umap(
        adata,
        color=color,
        layer=layer,
        legend_loc=legend_loc if (color in adata.obs and adata.obs[color].dtype.name == "category") else None,
        frameon=False,
        title=fname.replace(".pdf", "").replace(".png", ""),
        vmin=vmin,
        vmax=vmax,
        color_map=color_map,
        sort_order=True,        # <<< make high values (colored) plot on top of zeros (grey)
        show=False,
        save=None,
    )
    plt.savefig(os.path.join(OUTDIR, fname), bbox_inches="tight", dpi=200)
    plt.close()

def save_violin(adata, keys, groupby, fname, layer=None, use_raw=None, ylabel="counts"):
    # limit to keys that actually exist
    present = upper_intersection(keys, adata.var_names)
    if len(present) == 0:
        print(f"[WARN] None of {keys} found in var_names.")
        return
    sc.pl.violin(
        adata, keys=present, groupby=groupby,
        layer=layer, use_raw=use_raw, multi_panel=True,
        rotation=30, show=False
    )
    plt.suptitle(f"{fname.replace('.png','')} — {ylabel}")
    plt.savefig(os.path.join(OUTDIR, fname), bbox_inches="tight", dpi=200)
    plt.close()

def save_gene_umaps(adata, genes, layer, prefix, suffix_pdf=True):
    """
    Loop genes and save one UMAP per gene, manually plotting:
      - grey zero cells in the background
      - colored non-zero cells on top (zero→gray then viridis + smart vmax)
    Uses Scanpy's default figure size so it matches sc.pl.umap.
    """
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP coordinates not found in adata.obsm['X_umap'].")

    present = upper_intersection(genes, adata.var_names)
    missing = [g for g in genes if g.upper() not in {p.upper() for p in present}]
    if missing:
        print(f"[WARN] Missing genes in var_names: {missing}")

    umap = adata.obsm["X_umap"]
    x = umap[:, 0]
    y = umap[:, 1]

    # use the same default figure size as Scanpy's plotting
    sc_figsize = getattr(sc.settings, "figsize", None)

    for g in present:
        vals = _get_gene_values(adata, g, layer=layer)
        if vals.size == 0:
            continue

        vmax = _smart_vmax(vals)
        if vmax is None:
            vmax_use = ZERO_EPS       # all zeros → just show grey background
        else:
            vmax_use = vmax

        nonzero = vals > 0
        zero = ~nonzero

        if sc_figsize is not None:
            fig, ax = plt.subplots(figsize=sc_figsize)
        else:
            fig, ax = plt.subplots()

        # 1) grey zeros in the background
        if zero.any():
            ax.scatter(
                x[zero], y[zero],
                s=20,
                c="lightgray",
                linewidths=0,
                rasterized=False,
            )

        # 2) colored non-zero cells on top
        if nonzero.any():
            norm = mpl.colors.Normalize(vmin=ZERO_EPS, vmax=vmax_use)
            cvals = np.clip(vals[nonzero], ZERO_EPS, vmax_use)
            sc_hdl = ax.scatter(
                x[nonzero], y[nonzero],
                s=24,
                c=cvals,
                cmap=ZERO_GRAY_VIRIDIS,
                norm=norm,
                linewidths=0,
                rasterized=False,
            )
            cbar = plt.colorbar(sc_hdl, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label(f"{g} ({layer})")

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_title(g)

        plt.tight_layout()
        fn = (
            f"{prefix}_{g}_{('SCT' if layer=='SCT_data' else layer)}.pdf"
            if suffix_pdf else f"{prefix}_{g}_{layer}.png"
        )
        plt.savefig(os.path.join(OUTDIR, fn), bbox_inches="tight")
        plt.close(fig)


# ------------------- Full-cohort UMAPs (RNA base) -------------------
ensure_umap(adata_rna)  # compute if missing
groupby_rna = pick_groupby(adata_rna)
save_umap(adata_rna,   color=PREFERRED_GROUPBY if PREFERRED_GROUPBY in adata_rna.obs   else pick_groupby(adata_rna),   fname="UMAP_RNA_all_byType.pdf")

# Project RNA UMAP to histone modalities and carry labels
for ad_mod in (adata_ac, adata_me3):
    add_obs_from_rna(adata_rna, ad_mod, columns=(PREFERRED_GROUPBY, "leiden"))
    copy_umap_from_reference(adata_rna, ad_mod)  # use .obsm['X_umap'] from RNA

# ------------------- Histone gene activity (gene body ± window) + log1p layer -------------------
ga_ac  = make_gene_activity(adata_ac,  GENE_ANNO, PROMOTER_WINDOW)   # ±5 kb
ga_me3 = make_gene_activity(adata_me3, GENE_ANNO, PROMOTER_WINDOW)

# Copy RNA UMAP into gene activity matrices (for shared cells)
copy_umap_from_reference(adata_rna, ga_ac)
copy_umap_from_reference(adata_rna, ga_me3)

# ------------------- UMAPs colored by gene values -------------------
# RNA: color by gene using SCT_data layer (fallback to counts/raw/X if missing)
rna_color_layer = rna_layer_hard
save_gene_umaps(adata_rna, GENES, layer=rna_color_layer, prefix="UMAP_RNA_all", suffix_pdf=True)

# Histone marks: color by gene activity using log1p layer
save_gene_umaps(ga_ac,  GENES, layer=HTM_LAYER, prefix="UMAP_H3K27ac_all(RNA-embed)", suffix_pdf=True)
save_gene_umaps(ga_me3, GENES, layer=HTM_LAYER, prefix="UMAP_H3K27me3_all(RNA-embed)", suffix_pdf=True)

# ------------------- Violin plots: RNA expression (counts) -------------------
rna_layer_counts, rna_use_raw = ensure_counts_layer_for_rna(adata_rna)
save_violin(
    adata=adata_rna, keys=GENES, groupby=groupby_rna,
    layer=rna_layer_counts, use_raw=rna_use_raw,
    fname="Violin_RNA_counts_all.png", ylabel="RNA counts"
)

# ------------------- Violin plots: Histone gene activity (counts, NOT log1p) -------------------
save_violin(
    adata=ga_ac, keys=GENES, groupby=pick_groupby(ga_ac),
    fname=f"Violin_H3K27ac_geneActivity_pm{PROMOTER_WINDOW//1000}kb_all.png",
    ylabel=f"H3K27ac counts in gene body ±{PROMOTER_WINDOW//1000}kb"
)
save_violin(
    adata=ga_me3, keys=GENES, groupby=pick_groupby(ga_me3),
    fname=f"Violin_H3K27me3_geneActivity_pm{PROMOTER_WINDOW//1000}kb_all.png",
    ylabel=f"H3K27me3 counts in gene body ±{PROMOTER_WINDOW//1000}kb"
)

# ------------------- B cells only: recluster RNA, project embedding, repeat -------------------
B_TYPES = ["B_Cycling","B_GC_DZ","B_GC_LZ","B_mem+activated","B_naive","B_plasma"]
B_TYPES = ["B_Cycling","B_GC_DZ","B_GC_LZ","B_mem+activated"]

if groupby_rna not in adata_rna.obs:
    raise ValueError(f"Grouping column '{groupby_rna}' not found in RNA AnnData.")

# Select B cells using leiden_merged_type if available; otherwise fallback to startswith('B_')
if PREFERRED_GROUPBY in adata_rna.obs:
    mask_B = adata_rna.obs[PREFERRED_GROUPBY].astype(str).isin(B_TYPES)
else:
    mask_B = adata_rna.obs[groupby_rna].astype(str).str.startswith("B_")

adata_rna_B = adata_rna[mask_B].copy()

# Reclustering on RNA subset
sc.pp.pca(adata_rna_B, n_comps=10)
sc.pp.neighbors(adata_rna_B, n_neighbors=15)
sc.tl.umap(adata_rna_B)
sc.tl.leiden(adata_rna_B, key_added="leiden_B")
save_umap(
    adata_rna_B,
    color=PREFERRED_GROUPBY if PREFERRED_GROUPBY in adata_rna_B.obs else "leiden_B",
    fname="UMAP_RNA_Bcells_byType.pdf"
)

# Build B-only histone AnnDatas (by overlapping cells)
adata_ac_B  = adata_ac[adata_ac.obs_names.intersection(adata_rna_B.obs_names)].copy()
adata_me3_B = adata_me3[adata_me3.obs_names.intersection(adata_rna_B.obs_names)].copy()

# Bring annotations from RNA B into histone B objects
add_obs_from_rna(adata_rna_B, adata_ac_B,  columns=(PREFERRED_GROUPBY, "leiden", "leiden_B"))
add_obs_from_rna(adata_rna_B, adata_me3_B, columns=(PREFERRED_GROUPBY, "leiden", "leiden_B"))

# Project the B-only RNA UMAP to histone B objects
copy_umap_from_reference(adata_rna_B, adata_ac_B)
copy_umap_from_reference(adata_rna_B, adata_me3_B)

# Gene activity objects for B-only (subset from full GA; reuse log1p layer)
ga_ac_B  = ga_ac[ga_ac.obs_names.intersection(adata_rna_B.obs_names)].copy()
ga_me3_B = ga_me3[ga_me3.obs_names.intersection(adata_rna_B.obs_names)].copy()

# Copy the B-only RNA UMAP to GA B objects
copy_umap_from_reference(adata_rna_B, ga_ac_B)
copy_umap_from_reference(adata_rna_B, ga_me3_B)

# ------------------- B-only UMAPs colored by gene values (zero=gray then viridis) -------------------
rna_B_color_layer = rna_layer_hard

save_gene_umaps(adata_rna_B, GENES, layer=rna_B_color_layer, prefix="UMAP_RNA_Bcells", suffix_pdf=True)

save_gene_umaps(ga_ac_B,  GENES, layer=HTM_LAYER, prefix="UMAP_H3K27ac_Bcells(RNA-embed)", suffix_pdf=True)
save_gene_umaps(ga_me3_B, GENES, layer=HTM_LAYER, prefix="UMAP_H3K27me3_Bcells(RNA-embed)", suffix_pdf=True)

# ------------------- HEATMAPS (publication-style) -------------------
# Uses the same logic and helpers defined above:
# - GA from snap.pp.make_gene_matrix (gene body ± window) with log1p layer
# - zero→gray viridis via ZERO_GRAY_VIRIDIS
# - smart-vmax via _smart_vmax
HEATMAP_DIR = os.path.join(OUTDIR, "Heatmaps")
os.makedirs(HEATMAP_DIR, exist_ok=True)

def _group_order_and_index(adata, groupby):
    """Return ordered group names and column index (cells) sorted by group."""
    if groupby not in adata.obs:
        groupby = pick_groupby(adata)
    ser = adata.obs[groupby].astype("category") if groupby in adata.obs else pd.Series(index=adata.obs_names, data="all")
    if str(ser.dtype) == "category":
        order = list(ser.cat.categories)
    else:
        order = sorted(ser.unique().tolist())
    idx = []
    for g in order:
        idx.extend(list(np.where(adata.obs[groupby].astype(str).values == str(g))[0]))
    return groupby, order, np.array(idx, dtype=int)

def _get_gene_vector_aligned(adata, gene, layer, order_index=None):
    """1D values for a gene, optionally re-ordered by order_index; zeros if gene is absent."""
    try:
        j = adata.var_names.get_loc(gene)
    except KeyError:
        return np.zeros(adata.n_obs, dtype=float)
    Xsrc = adata.layers[layer] if layer is not None else adata.X
    if sp.issparse(Xsrc):
        v = Xsrc.getcol(j).toarray().ravel()
    else:
        v = np.asarray(Xsrc[:, j]).ravel()
    return v if order_index is None else v[order_index]

def _scale_rows_zero_to_one(rows):
    """
    rows: list of 1D arrays (already log1p where appropriate)
    Uses smart-vmax per row and divides by vmax (clipped 0..1).
    Returns 2D array shape (n_rows, n_cells) and list of per-row vmax used.
    """
    scaled = []
    vmax_used = []
    for v in rows:
        vmax = _smart_vmax(v)
        if vmax is None or not np.isfinite(vmax) or vmax <= 0:
            scaled.append(np.zeros_like(v))
            vmax_used.append(0.0)
        else:
            s = np.clip(v / vmax, 0.0, 1.0)
            scaled.append(s)
            vmax_used.append(float(vmax))
    return np.vstack(scaled), vmax_used

def plot_gene_heatmap_three_rows(
    gene, adata_rna, ga_ac, ga_me3, rna_layer_for_heatmap=None,
    groupby=None, fname_prefix="Heatmap_all"
):
    # --- grouping & base order from RNA ---
    if groupby is None:
        groupby = pick_groupby(adata_rna)
    groupby, group_order, base_index = _group_order_and_index(adata_rna, groupby)

    # --- choose RNA layer (param or auto) ---
    rna_layer = resolve_rna_layer(adata_rna, RNA_LAYER_FOR_HEATMAP, RNA_LAYER_PRIORITY) \
                if rna_layer_for_heatmap is None else rna_layer_for_heatmap

    # --- get vectors in RNA order, then refine order per group by pattern ---
    # 1) raw vectors in RNA order (not grouped yet)
    # RNA vector from chosen layer (or X if None)
    try:
        j_rna = adata_rna.var_names.get_loc(gene)
        Xr = adata_rna.layers[rna_layer] if rna_layer is not None else adata_rna.X
        rna_vec = Xr.getcol(j_rna).toarray().ravel() if sp.issparse(Xr) else np.asarray(Xr[:, j_rna]).ravel()
    except KeyError:
        rna_vec = np.zeros(adata_rna.n_obs, dtype=float)

    # GA vectors: counts AND log1p in RNA order
    # First map GA obs → RNA obs order
    pos = pd.Series(index=adata_rna.obs_names, data=np.arange(adata_rna.n_obs))
    ac_cnt_rna, ac_log_rna   = _ga_gene_vectors(ga_ac,  gene, order_index=pos.loc[ga_ac.obs_names].values)
    me3_cnt_rna, me3_log_rna = _ga_gene_vectors(ga_me3, gene, order_index=pos.loc[ga_me3.obs_names].values)

    # 2) take the RNA "grouped" order, but refine inside groups by pattern
    ser_grp = adata_rna.obs[groupby].astype(str).values
    group_names = [g for g in group_order if np.any(ser_grp == str(g))]

    refined_indices = []
    for g in group_names:
        idx_g = np.where(ser_grp == str(g))[0]
        # sort these idx_g by overlap pattern using COUNT vectors for binarization
        ord_g = _order_cells_by_pattern(
            rna_vals=rna_vec, ac_cnt=ac_cnt_rna, me3_cnt=me3_cnt_rna,
            group_index=idx_g, rna_layer_name=rna_layer,
            pattern_order=PATTERN_ORDER, rna_rule=RNA_BIN_RULE, ac_rule=AC_BIN_RULE, me3_rule=ME3_BIN_RULE
        )
        refined_indices.append(ord_g)
    order_index = np.concatenate(refined_indices, axis=0) if len(refined_indices) else base_index

    # 3) final per-row values (display): RNA uses chosen layer; GA uses log1p layer
    rna_v = rna_vec[order_index]
    ac_v  = ac_log_rna[order_index]
    me3_v = me3_log_rna[order_index]

    # --- scale 0..1 per row (smart vmax) ---
    M, vmaxs = _scale_rows_zero_to_one([rna_v, ac_v, me3_v])

    # --- plot ---
    n = M.shape[1]
    fig_w = max(6.0, min(18.0, n / 80.0 + 6.0))
    fig_h = 2.6

    plt.figure(figsize=(fig_w, fig_h), dpi=300)
    ax = plt.gca()
    im = ax.imshow(M, aspect="auto", interpolation="nearest",
                   cmap=ZERO_GRAY_VIRIDIS, vmin=ZERO_EPS, vmax=1.0)
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["RNA", "H3K27ac", "H3K27me3"])
    ax.set_xticks([]); ax.set_title(gene, pad=6)

    # --- draw group boundaries + labels (TOP axis) in the *refined* order ---
    counts, labels = [], []
    for g in group_names:
        counts.append(int((adata_rna.obs[groupby].astype(str).values[order_index] == str(g)).sum()))
        labels.append(str(g))
    edges = np.cumsum([0] + counts)
    for boundary in edges[1:-1]:
        ax.vlines(boundary - 0.5, -0.5, 2.5, color="white", lw=0.6, alpha=0.8)
    centers, left = [], 0
    for c in counts:
        centers.append(left + (c - 1) / 2.0); left += c
    top_ax = ax.secondary_xaxis('top')
    top_ax.set_xticks(centers); top_ax.set_xticklabels(labels)
    plt.setp(top_ax.get_xticklabels(), rotation=40, ha='left', va='bottom', fontsize=8)
    top_ax.tick_params(axis='x', pad=2, length=0)

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("scaled intensity (0..1)")

    plt.tight_layout()
    png = os.path.join(HEATMAP_DIR, f"{fname_prefix}_{gene}_RNA_AC_ME3.png")
    pdf = os.path.join(HEATMAP_DIR, f"{fname_prefix}_{gene}_RNA_AC_ME3.pdf")
    plt.savefig(pdf, bbox_inches="tight"); plt.close()
    print(f"[heatmap] {gene} — vmax: RNA={vmaxs[0]:.3g}, ac={vmaxs[1]:.3g}, me3={vmaxs[2]:.3g}  |  RNA layer='{rna_layer}'")

# ---- run heatmaps for all cells ----
rna_layer_for_heatmap = resolve_rna_layer(adata_rna, RNA_LAYER_FOR_HEATMAP, RNA_LAYER_PRIORITY)
if rna_layer_for_heatmap is None:
    print("[INFO] RNA heatmaps/UMAPs will use X (no layer).")
save_gene_umaps(adata_rna, GENES, layer=rna_layer_for_heatmap, prefix="UMAP_RNA_all", suffix_pdf=True)
groupby_rna = pick_groupby(adata_rna)

for g in upper_intersection(GENES, adata_rna.var_names):
    plot_gene_heatmap_three_rows(
        g, adata_rna=adata_rna, ga_ac=ga_ac, ga_me3=ga_me3,
        rna_layer_for_heatmap=rna_layer_for_heatmap, groupby=groupby_rna,
        fname_prefix="Heatmap_AllCells"
    )

# ---- run heatmaps for B-only subset ----
if 'adata_rna_B' in globals():
    rna_layer_B_for_heatmap = resolve_rna_layer(adata_rna_B, RNA_LAYER_FOR_HEATMAP, RNA_LAYER_PRIORITY)
    groupby_B = PREFERRED_GROUPBY if PREFERRED_GROUPBY in adata_rna_B.obs else "leiden_B"
    save_gene_umaps(adata_rna_B, GENES, layer=rna_layer_B_for_heatmap, prefix="UMAP_RNA_Bcells", suffix_pdf=True)

    for g in upper_intersection(GENES, adata_rna_B.var_names):
        plot_gene_heatmap_three_rows(
            g, adata_rna=adata_rna_B, ga_ac=ga_ac_B, ga_me3=ga_me3_B,
            rna_layer_for_heatmap=rna_layer_B_for_heatmap, groupby=groupby_B,
            fname_prefix="Heatmap_Bcells"
        )

print(f"Done. Figures in: {OUTDIR} (heatmaps in {os.path.join(OUTDIR, 'Heatmaps')})")



#QC Violin plots


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc

sns.set(style="whitegrid")

def _first_present(A: ad.AnnData, candidates):
    for c in candidates:
        if c in A.obs.columns:
            return c
    return None

def fmt_log_ticks(ax):
    """
    Set y-ticks for axes where y = log10(raw + 1), with
    nice raw values like 0, 100, 250, 500, 1000, 2500, ...
    """
    ymin, ymax = ax.get_ylim()
    # convert current limits back to raw scale
    raw_min = max(0.0, 10**ymin - 1.0)
    raw_max = max(0.0, 10**ymax - 1.0)

    # candidate raw tick values
    candidates = np.array([
        0, 10, 25, 50,
        100, 250, 500,
        1000, 2500, 5000,
        10000, 25000, 50000, 100000
    ], dtype=float)

    ticks_raw = candidates[(candidates >= raw_min) & (candidates <= raw_max)]

    # if range is small / weird, fall back to something reasonable
    if ticks_raw.size < 2:
        # simple fallback: min, mid, max in raw
        rmax = max(raw_max, 1.0)
        ticks_raw = np.array([0, rmax / 2, rmax])

    ticks_raw = np.unique(ticks_raw)
    ticks_log = np.log10(ticks_raw + 1.0)

    ax.set_yticks(ticks_log)
    ax.set_yticklabels([f"{int(tr):,}" for tr in ticks_raw])

def _annotate_global_count_and_median(ax, values, y_is_log=False):
    """Add 'N cells / median' label above a single violin."""
    ymin, ymax = ax.get_ylim()
    headroom = 0.12 * (ymax - ymin)
    ax.set_ylim(ymin, ymax + headroom)
    ypos = ymax + 0.04 * (ymax - ymin)

    v = np.asarray(values)
    if v.size == 0:
        return

    if y_is_log:
        med_raw = int(np.median(10**v - 1))
        txt = f"{v.size} cells\n{med_raw:,} median"
    else:
        med = float(np.median(v))
        txt = f"{v.size} cells\n{med:.2f} median"

    ax.text(0, ypos, txt, ha="center", va="bottom", fontsize=10)

def violin_global(A: ad.AnnData, value_col: str, ylabel: str, outpath: str,
                  title: str, log10_plus1: bool):
    """Single global violin (all cells) with box overlay and nice annotation."""
    if value_col not in A.obs.columns:
        return

    vals = pd.to_numeric(A.obs[value_col], errors="coerce").fillna(0)
    y = np.log10(vals.to_numpy() + 1) if log10_plus1 else vals.to_numpy()

    df = pd.DataFrame({"group": ["all"] * len(y), "value": y})

    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    # Violin
    sns.violinplot(
        data=df, x="group", y="value",
        inner=None, ax=ax
    )
    # Box on top
    sns.boxplot(
        data=df, x="group", y="value",
        width=0.25, showcaps=True, dodge=False,
        boxprops={'facecolor': 'none', 'linewidth': 1.5},
        whiskerprops={'linewidth': 1.5},
        medianprops={'color': 'black'},
        flierprops={'marker': '', 'markersize': 0},
        ax=ax
    )

    if log10_plus1:
        # set ticks/labels to nice raw values (0, 100, 250, 500, ...)
        fmt_log_ticks(ax)

    _annotate_global_count_and_median(ax, y, y_is_log=log10_plus1)

    ax.set_xlabel("")
    ax.set_xticklabels(["All cells"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath)
    plt.close(fig)

# ------------------- Public entrypoint -------------------
def make_qc_plots_simple(rna_ad: ad.AnnData, ac_ad: ad.AnnData, me3_ad: ad.AnnData,
                         outdir: str, prefix: str = "Showcase"):
    """
    Create global violins (no conditions) and scatters for:
      - RNA: UMI per cell (log10) and nGenes per cell (log10)
      - AC: Fragments per cell (log10), TSSE violin, TSSE vs log10(Fragments) scatter
      - ME3: Fragments per cell (log10), TSSE violin (if present), scatter (if present)
    """
    os.makedirs(outdir, exist_ok=True)

    # ----- RNA -----
    if rna_ad is not None:
        umi_col = _first_present(rna_ad, ["total_counts", "n_counts", "nCount_RNA", "umi", "counts", "SCT_counts"])
        genes_col = _first_present(rna_ad, ["n_genes_by_counts", "nFeature_RNA", "genes"])

        if umi_col is not None:
            violin_global(
                rna_ad, umi_col, "UMIs per cell",
                os.path.join(outdir, f"{prefix}_RNA_UMI_global.pdf"),
                title=f"RNA – UMIs per Cell\nCells: {rna_ad.n_obs:,}",
                log10_plus1=True
            )

        if genes_col is not None:
            violin_global(
                rna_ad, genes_col, "Genes per cell",
                os.path.join(outdir, f"{prefix}_RNA_nGenes_global.pdf"),
                title=f"RNA – Genes per Cell\nCells: {rna_ad.n_obs:,}",
                log10_plus1=True
            )

    # ----- H3K27ac -----
    if ac_ad is not None:
        frag_col = _first_present(ac_ad, ["n_fragment", "passed_filters", "nCount_ATAC"])
        if frag_col is not None:
            violin_global(
                ac_ad, frag_col, "Fragments per cell",
                os.path.join(outdir, f"{prefix}_H3K27ac_Fragments_global.pdf"),
                title=f"H3K27ac – Fragments per Cell\nCells: {ac_ad.n_obs:,}",
                log10_plus1=True
            )

        tsse_col = _first_present(ac_ad, ["tsse", "tss_enrichment", "TSSE"])
        if tsse_col is not None:
            violin_global(
                ac_ad, tsse_col, "TSSE",
                os.path.join(outdir, f"{prefix}_H3K27ac_TSSE_global.pdf"),
                title=f"H3K27ac – TSSE\nCells: {ac_ad.n_obs:,}",
                log10_plus1=False
            )

    # ----- H3K27me3 -----
    if me3_ad is not None:
        frag_col = _first_present(me3_ad, ["n_fragment", "passed_filters", "nCount_ATAC"])
        if frag_col is not None:
            violin_global(
                me3_ad, frag_col, "Fragments per cell",
                os.path.join(outdir, f"{prefix}_H3K27me3_Fragments_global.pdf"),
                title=f"H3K27me3 – Fragments per Cell\nCells: {me3_ad.n_obs:,}",
                log10_plus1=True
            )

        tsse_col = _first_present(me3_ad, ["tsse", "tss_enrichment", "TSSE"])
        if tsse_col is not None:
            violin_global(
                me3_ad, tsse_col, "TSSE",
                os.path.join(outdir, f"{prefix}_H3K27me3_TSSE_global.pdf"),
                title=f"H3K27me3 – TSSE\nCells: {me3_ad.n_obs:,}",
                log10_plus1=False
            )




adataRNA = sc.read_h5ad("Tonsil_rna_Overlap_Reclust.h5ad")
adataMe3 = sc.read_h5ad("Tonsil_me3_Overlap_Reclust.h5ad")
adataAc  = sc.read_h5ad("Tonsil_ac_Overlap_Reclust.h5ad")

make_qc_plots_simple(
    rna_ad=adataRNA,
    ac_ad=adataAc,
    me3_ad=adataMe3,
    outdir="QC_Overlap3",
    prefix="Overlap3"
)


adataRNA = sc.read_h5ad("Tonsil_FULL_RNA_ANNOTATED.h5ad")
adataMe3 = sc.read_h5ad("Sc_VTHumanTonsil_me3_processed.h5ad")
adataAc  = sc.read_h5ad("Sc_VTHumanTonsil_ac_processed.h5ad")


make_qc_plots_simple(
    rna_ad=adataRNA,
    ac_ad=adataAc,
    me3_ad=adataMe3,
    outdir="QC_Full",
    prefix="Full"
)




