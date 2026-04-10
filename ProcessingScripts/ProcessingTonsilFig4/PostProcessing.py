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
import scipy.sparse as sp
import anndata
import scvi

if not os.path.exists('figures'):
    os.makedirs('figures')

#################################################
############# RNA ANALYSIS ###############
#################################################

# Use sample stems only (without the .h5ad suffix)
cell_lines = [
    "Sc_TP3r_TP3_A",
    "Sc_TP3r_TP3_B",
    "Sc_TP4r_TP4_A",
    "Sc_TP4r_TP4_B",
    "Sc_HumanTonsil2", 
]

# For R code
import anndata2ri
import logging
import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.robjects import default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri, pandas2ri

rcb.logger.setLevel(logging.ERROR)
anndata2ri.activate()
%load_ext rpy2.ipython

## Doublet Analysis
%R library(Seurat)
%R library(scater)
%R library(scDblFinder)
%R library(BiocParallel)
%R library(scran)
%R library(scry)

for exp in cell_lines:
    print(f"\n{'='*80}")
    print(f"Processing {exp}")
    print(f"{'='*80}")

    try:
        #################################################
        ############# LOAD SAMPLE ###############
        #################################################
        adata_sample = sc.read_h5ad(f"{exp}_rna_CountFiltered.h5ad")
        adata_sample.obs['sample'] = exp

        #################################################
        ############# DOUBLET ANALYSIS ###############
        #################################################
        data_mat = adata_sample.X.T

        # rpy2 cannot handle SciPy sparse matrices – convert to dense
        if sp.issparse(data_mat):
            data_mat = data_mat.toarray()

        data_mat = data_mat.astype(float)

        %R -i data_mat
        %R set.seed(1)
        %R sce = scDblFinder(SingleCellExperiment(list(counts=data_mat)), artificialDoublets=4000, dbr=0.1, cluster=TRUE)
        %R doublet_scoreR <- sce$scDblFinder.score
        %R doublet_classR <- sce$scDblFinder.class

        adata_sample.obs["scDblFinder_score"] = %Rget doublet_scoreR
        adata_sample.obs["scDblFinder_class"] = %Rget doublet_classR

        sc.pl.scatter(
            adata_sample,
            "scDblFinder_score",
            "log1p_total_counts",
            color="scDblFinder_class",
            title='Color: scDblFinder_class',
            show=False,
            save=f'_{exp}_rna_scDBL.png'
        )

        score = adata_sample.obs['scDblFinder_score']
        adata_sample.obs['scDblFinder_class_manual'] = np.where(score < 0.5, 'singlet', 'doublet')

        adata = adata_sample[adata_sample.obs['scDblFinder_score'] < 0.6].copy()
        # adata = adata_sample[adata_sample.obs['scDblFinder_class']=="singlet"].copy()
        # adata = adata_sample.copy()

        # Store raw counts in a layer
        adata.layers["counts"] = adata.X.copy()

        #################################################
        ############# NORMALIZATION ###############
        #################################################
        print('Shifted logarithm Normalization')
        adata.layers["log1p_norm"] = sc.pp.log1p(adata.X, copy=True)

        adata.write(f"{exp}_rna_raw.h5ad")
        adata = sc.read_h5ad(f"{exp}_rna_raw.h5ad")

        print('SCTransform Normalization')
        sc.pp.filter_genes(adata, min_cells=5, inplace=True)
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
        r('seurat_obj <- SCTransform(seurat_obj, vst.flavor="v2")')

        # Extract the SCT data and add it as a new layer
        sct_data = np.asarray(r['as.matrix'](r('seurat_obj@assays$SCT@data')))
        adata.layers['SCT_data'] = sct_data.T
        sct_counts = np.asarray(r['as.matrix'](r('seurat_obj@assays$SCT@counts')))
        adata.layers['SCT_counts'] = sct_counts.T

        # choose one of the layers
        # adata.X = adata.layers['SCT_data']
        adata.X = adata.layers["log1p_norm"]

        #################################################
        ############# CELL CYCLE / FEATURE SELECTION ###############
        #################################################
        cell_cycle_genes = [x.strip() for x in open('/mnt/dataFast/ahrmad/regev_lab_cell_cycle_genes.txt')]
        s_genes = cell_cycle_genes[:43]
        g2m_genes = cell_cycle_genes[43:]
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

        print('Deviance')
        xTopDeviantGenes = 2000

        # remove plotting-only offenders before conversion
        to_drop = [k for k in adata.uns if k.endswith("_colors") or k.endswith("_color")]
        for k in to_drop:
            adata.uns.pop(k, None)

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

        %R sce = devianceFeatureSelection(adata, assay="X")

        binomial_deviance = ro.r("rowData(sce)$binomial_deviance").T
        idx = binomial_deviance.argsort()[-xTopDeviantGenes:]
        mask = np.zeros(adata.var_names.shape, dtype=bool)
        mask[idx] = True

        adata.var["highly_deviant"] = mask
        adata.var["binomial_deviance"] = binomial_deviance

        # Highly variable genes
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3')

        sc.pl.highly_variable_genes(adata, show=False, save=f'_{exp}_HVG_ScanpySeuratV3_SCT.pdf')
        adata.var["highly_variableScanpySeuratV3"] = adata.var["highly_variable"].copy()
        adata.var["highly_variable"] = adata.var["highly_deviant"].copy()
        sc.pl.highly_variable_genes(adata, show=False, save=f'_{exp}_HVG_R_SCT.pdf')
        adata.var["highly_variable"] = adata.var["highly_variableScanpySeuratV3"].copy()

        #################################################
        ############# PCA / UMAP / LEIDEN ###############
        #################################################
        sc.tl.pca(adata, n_comps=20)
        sc.pl.pca_variance_ratio(adata, n_pcs=20, log=False, show=False, save=f'_{exp}_pca_variance_ratio.pdf')

        sc.pl.pca(
            adata,
            color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase"],
            dimensions=[(0, 1), (0, 1), (0, 1), (0, 1)],
            ncols=3,
            show=False, save=f'_{exp}_PCA1-2_QC.pdf'
        )
        sc.pl.pca(
            adata,
            color=["pct_counts_in_top_50_genes", "pct_counts_mt", "log1p_total_counts", "phase"],
            dimensions=[(2, 3), (2, 3), (2, 3), (2, 3)],
            ncols=3,
            show=False, save=f'_{exp}_PCA3-4_QC.pdf'
        )

        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.umap(adata)

        sc.pl.umap(
            adata,
            color=["pct_counts_mt", "pct_counts_hb", "pct_counts_ribo", "pct_counts_in_top_50_genes", "log1p_total_counts"],
            ncols=3,
            show=False, save=f'_{exp}_UMAP_QC.pdf'
        )
        sc.pl.umap(
            adata,
            color=["log1p_total_counts", "log1p_n_genes_by_counts", "scDblFinder_score", "scDblFinder_class"],
            ncols=3,
            show=False, save=f'_{exp}_UMAP_QC_Log.pdf'
        )
        sc.pl.umap(
            adata,
            color=['S_score', 'G2M_score', 'phase'],
            ncols=2,
            show=False, save=f'_{exp}_UMAP_QC_CellCycle.pdf'
        )

        leiR = 0.15
        leiRes = f'leiden_res{leiR}'
        sc.tl.leiden(adata, key_added=leiRes, resolution=leiR)
        sc.pl.umap(
            adata,
            color=leiRes,
            legend_loc="on data",
            save=f"_{exp}_{leiRes}.pdf",
            title=f"{exp} {leiRes}",
            show=False
        )

        #################################################
        ############# BCELL SIGS ###############
        #################################################
        b_cells_dict = {
            "B_cells_general": ["MS4A1", "CD79A", "CD79B", "VPREB3", "CD19"],
            "B_cells_naive": ["IGHM", "IGHD", "CCR7", "SELL"],
            "B_cells_activated": ["CD69", "CD83", "JUN", "IRF4", "MYC"],
            "B_cells_memory": ["CD27", "XBP1", "CLEC2D", "SSPN"],
            "B_cells_unswitched memory": ["IGHM", "IGHD", "CD27", "FCRL5"],
            "B_cells_proliferative": ["TCL1A", "CD79B", "CD79A", "MS4A1", "MKI67"],
            "Plasma_B_cells": [
                "JCHAIN", "IGHGP", "IGKC", "IGHM", "IGHG3",
                "IGKV1D-39", "IGHG2", "IGLC3", "IGHG1",
                "IGHG4", "IGHV4-4", "XBP1", "CD79A", "CD38",
                "MZB1", "SSR4", "DERL3", "FKBP11", "PRDX4"
            ]
        }

        for cl, genelist in b_cells_dict.items():
            sc.tl.score_genes(adata, genelist, score_name=f'BCellSubSig_{cl}')

        sc.pl.umap(
            adata,
            color=[f'BCellSubSig_{clu}' for clu in list(b_cells_dict.keys())],
            ncols=4,
            color_map='viridis',
            vmax=0.3, vmin=0,
            show=False,
            save=f'_{exp}_UMAP_BCellSubSig.pdf'
        )

        marker_genes_dict = {
            "B cells": ["CD79A"],
            "T cells": ["CD3D"],
            "NK cells": ["NKG7"],
            "Myeloid": ["CD68"],
            "Fibroblasts": ["COL1A1"],
            "Epithelial": ["EPCAM", "CDH1"],
            "Endothelial": ["PTPRB"]
        }

        sc.tl.dendrogram(adata, groupby=leiRes)

        sc.set_figure_params(figsize=(10, 10))
        pl_fig = sc.pl.dotplot(
            adata,
            marker_genes_dict,
            leiRes,
            dendrogram=True,
            return_fig=True
        )
        pl_fig.add_totals(color='lightpink').savefig(f'./figures/{exp}_RNA_Imm.pdf')
        plt.tight_layout()
        plt.close()

        sc.tl.rank_genes_groups(adata, groupby=leiRes, method="wilcoxon")
        pl_fig = sc.pl.rank_genes_groups_dotplot(
            adata, groupby=leiRes, n_genes=10, show=False, return_fig=True
        )
        pl_fig.add_totals(color='lightpink').savefig(f'./figures/{exp}_RNA_ImmDotPlot10.pdf')
        plt.close()

        #################################################
        ############# SAVE SAMPLE-LEVEL OUTPUT ###############
        #################################################
        query_path = f"{exp}_Tonsil_RNAOnly.h5ad"
        adata.write(query_path)

        #################################################
        ############# LABEL TRANSFER ###############
        #################################################
        REF_PATH = "sc.h5ad"
        QUERY_PATH = query_path
        LABEL_KEY = "new_celltype"
        BATCH_KEY = "dataset"
        SAVE_DIR = "figures"
        os.makedirs(SAVE_DIR, exist_ok=True)

        adata_ref = sc.read(REF_PATH)
        adata_query = sc.read(QUERY_PATH)

        adata_query.X = adata_query.layers["counts"].copy()
        sc.pl.umap(adata_ref, color=["new_celltype"], ncols=1, show=False, save=f"_{exp}_FullRef_CellType_OriUMAP.pdf")

        shared_genes = adata_ref.var_names.intersection(adata_query.var_names)
        print(f"{exp}: Number of shared genes: {len(shared_genes)}")

        # --- 0. Extract gene names ---
        ref_genes = pd.Index(adata_ref.var_names)
        query_genes = pd.Index(adata_query.var_names)

        # --- 1. Standardize gene names ---
        def clean_gene_name(g):
            g = g.upper()
            g = g.split(".")[0]
            return g

        ref_clean = ref_genes.map(clean_gene_name)
        query_clean = query_genes.map(clean_gene_name)

        # --- 2. Compute smart shared genes ---
        shared_genes_clean = pd.Index(ref_clean).intersection(query_clean)

        # Map back to original names
        shared_ref = ref_genes[ref_clean.isin(shared_genes_clean)]
        shared_query = query_genes[query_clean.isin(shared_genes_clean)]

        print(f"{exp}: Number of shared genes WITH UPPERCASE: {len(shared_genes_clean)}")

        # Keeping your original behavior here
        adata_ref = adata_ref[:, shared_genes].copy()
        adata_query = adata_query[:, shared_genes].copy()

        adata_ref.obs[BATCH_KEY] = "ref"
        adata_query.obs[BATCH_KEY] = "query"

        SCANVI_LABEL_KEY = "labels"
        adata_ref.obs[SCANVI_LABEL_KEY] = adata_ref.obs[LABEL_KEY]
        adata_query.obs[SCANVI_LABEL_KEY] = "Unknown"

        adata_combined = anndata.concat([adata_ref, adata_query], label="source", keys=["ref", "query"])
        adata_combined.layers["counts"] = adata_combined.X.copy()

        sc.pp.normalize_total(adata_combined, target_sum=1e4)
        sc.pp.log1p(adata_combined)
        adata_combined.raw = adata_combined
        sc.pp.highly_variable_genes(
            adata_combined,
            flavor="seurat_v3",
            n_top_genes=2000,
            layer="counts",
            batch_key=BATCH_KEY,
            subset=True,
        )

        # === Train scVI ===
        scvi.model.SCVI.setup_anndata(adata_combined, layer="counts", batch_key=BATCH_KEY)
        scvi_model = scvi.model.SCVI(adata_combined)
        scvi_model.train()

        adata_combined.obsm["X_scVI"] = scvi_model.get_latent_representation()
        sc.pp.neighbors(adata_combined, use_rep="X_scVI")
        sc.tl.umap(adata_combined, min_dist=0.3)

        sc.pl.umap(adata_combined, color=[BATCH_KEY], show=False)
        plt.savefig(f"{SAVE_DIR}/{exp}_umap_scvi_batch.png", bbox_inches="tight")
        plt.close()

        sc.pl.umap(adata_combined, color=[SCANVI_LABEL_KEY], show=False)
        plt.savefig(f"{SAVE_DIR}/{exp}_umap_scvi_labels.png", bbox_inches="tight")
        plt.close()

        # Save sample-specific AnnData and model
        adata_combined.write(f"{exp}_my_adata.h5ad")
        scvi_model.save(f"{exp}_my_scvi_model/", overwrite=True)

        # Reload
        adata_combined = sc.read(f"{exp}_my_adata.h5ad")
        scvi_model = scvi.model.SCVI.load(f"{exp}_my_scvi_model/", adata=adata_combined)

        # === Train scANVI for label transfer ===
        scvi.model.SCANVI.setup_anndata(
            adata_combined,
            labels_key=SCANVI_LABEL_KEY,
            unlabeled_category="Unknown"
        )
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            unlabeled_category="Unknown",
            labels_key=SCANVI_LABEL_KEY
        )
        scanvi_model.train(max_epochs=20, n_samples_per_label=100)

        # === Predict labels and get latent space ===
        adata_combined.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata_combined)

        probs = scanvi_model.predict(adata_combined, soft=True)
        adata_combined.obs["predicted_labels"] = probs.idxmax(axis=1)
        adata_combined.obs["label_confidence"] = probs.max(axis=1)
        probs.to_csv(f"{exp}_label_probabilities.csv")

        sc.pp.neighbors(adata_combined, use_rep="X_scANVI")
        sc.tl.umap(adata_combined)

        sc.pl.umap(
            adata_combined,
            color=["predicted_labels", "label_confidence", "dataset"],
            ncols=1,
            show=False,
            save=f"_{exp}_scanvi_predicted_labels.pdf"
        )

        # === Save query with predicted labels ===
        query_mask = adata_combined.obs["source"] == "query"
        adata_query.obs["predicted_labels"] = adata_combined[query_mask].obs["predicted_labels"].values
        adata_query.obs["label_confidence"] = adata_combined[query_mask].obs["label_confidence"].values

        sc.pl.umap(
            adata_query,
            color=["predicted_labels", "label_confidence", "scDblFinder_score"],
            ncols=1,
            show=False,
            save=f"_{exp}_Annot_Filt_LT1.pdf"
        )

        # Cell counts by predicted label
        label_counts = adata_query.obs["predicted_labels"].value_counts().sort_values(ascending=False)
        label_counts_df = pd.DataFrame({
            "cell_type": label_counts.index,
            "cell_count": label_counts.values
        })

        plt.figure(figsize=(10, 5))
        plt.bar(label_counts_df["cell_type"], label_counts_df["cell_count"])
        plt.xticks(rotation=90)
        plt.ylabel("Number of cells")
        plt.title(f"{exp}: Cell counts per predicted label (sorted)")
        plt.tight_layout()
        plt.savefig(f"figures/{exp}_query_celltype_distribution_sorted1.png", bbox_inches="tight")
        plt.close()

        adata_query.write(f"{exp}_TonsilRNA_Annot_Filt_LT.h5ad")

        print(f"Finished {exp}")

    except Exception as e:
        print(f"Could not process {exp}: {e}")
















import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp

os.makedirs("figures", exist_ok=True)
sc.settings.figdir = "figures"

#################################################
############# SETTINGS ###############
#################################################

# Auto-discover all sample-specific outputs from your prior label-transfer loop
input_files = sorted(glob.glob("*_TonsilRNA_Annot_Filt_LT.h5ad"))

# If you want to hardcode instead, use this:
# input_files = [
#     "Sc_TP3r_TP3_A_TonsilRNA_Annot_Filt_LT.h5ad",
#     "Sc_TP3r_TP3_B_TonsilRNA_Annot_Filt_LT.h5ad",
#     "Sc_TP4r_TP4_A_TonsilRNA_Annot_Filt_LT.h5ad",
#     "Sc_TP4r_TP4_B_TonsilRNA_Annot_Filt_LT.h5ad",
# ]

if len(input_files) == 0:
    raise FileNotFoundError(
        "No '*_TonsilRNA_Annot_Filt_LT.h5ad' files found. "
        "Run the previous label-transfer loop first, or point input_files to the correct sample files."
    )

print("Detected input files:")
for f in input_files:
    print("  ", f)

# Cell-level filter
LOW_CONF_Q = 0.20
HIGH_DBL_Q = 0.80

# Cluster-level filter
OVERCLUSTER_RES = 20
CLUSTER_CONF_THRESHOLD = 0.90

# Neighbor agreement filter
NEIGHBOR_AGREEMENT_THRESHOLD = 0.10
NEIGHBOR_WEIGHTED = True

# Final clustering
FINAL_LEIDEN_RES = 5

# Use normalized data for PCA/UMAP if available
USE_LOG1P_LAYER_FOR_EMBEDDING = True

#################################################
############# LABEL MERGING RULES ###############
#################################################

REMOVE_LABELS = ['ILC', 'NK', 'NKT', 'Monocytes', 'T_TIM3+','Mast']

MERGE_RULES = {
    'B_mem': ['B_mem', 'B_IFN'],
    'B_naive': ['B_activated,4', 'B_activated,3', 'B_activated,0', 'B_activated,2', 'B_activated,1', 'B_naive'],
    'T_CD4+': ['T_CD4+', 'T_TfR', 'T_Treg'],
    'B_GC_LZ': ['B_GC_prePB'],
    'T_CD8+': ['T_CD8+_CD161+', 'T_CD8+_cytotoxic', 'T_CD8+_naive'],
    'FDC': ['VSMC', 'Endo'],
    'Macrophages': ['Macrophages_M2', 'Macrophages_M1'],
    'DC': ['DC_cDC1', 'DC_cDC2', 'DC_pDC', 'DC_CCR7+'],
}

CELL_TYPE_COLORS = {
    "B_Cycling":        "#00BFC4",
    "B_GC_DZ":          "#98DF8A",
    "B_GC_LZ":          "#2CA02C",
    "B_mem":  "#2B8C9E",
    "B_naive":          "#1F77B4",
    "B_plasma":         "#76C1FF",
    "DC":               "#B39B00",
    "FDC":              "#7B1FA2",
    "Macrophages":      "#8C564B",
    "T_CD4+":           "#D62728",
    "T_CD4+_TfH": "#E377C2",
    "T_CD4+_TfH_GC": "#C05AA0",
    "T_CD4+_naive":     "#FF9896",
    "T_CD8+":           "#FF7F0E",
}

DESIRED_ORDER = [
    "B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem+activated", "B_naive", "B_plasma",
    "DC", "FDC", "Macrophages",
    "T_CD4+", "T_CD4+_FH", "T_CD4+_naive", "T_CD8+",
]

#################################################
############# HELPERS ###############
#################################################

def get_sample_name(path):
    base = os.path.basename(path)
    suffix = "_TonsilRNA_Annot_Filt_LT.h5ad"
    if base.endswith(suffix):
        return base[:-len(suffix)]
    return os.path.splitext(base)[0]

def prepare_X_for_embedding(adata):
    """
    Keep your intended workflow, but make PCA/UMAP robust:
    - prefer log1p_norm layer if present
    - otherwise build log-normalized X from counts
    """
    if USE_LOG1P_LAYER_FOR_EMBEDDING and "log1p_norm" in adata.layers:
        adata.X = adata.layers["log1p_norm"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        if USE_LOG1P_LAYER_FOR_EMBEDDING:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
    return adata

def recompute_embedding(adata, sample, tag, n_top_genes=3000, n_pcs=15, n_neighbors=30, metric='cosine'):
    adata = adata.copy()
    adata = prepare_X_for_embedding(adata)

    hvg_kwargs = dict(
        flavor='seurat_v3',
        n_top_genes=n_top_genes,
    )
    if "counts" in adata.layers:
        hvg_kwargs["layer"] = "counts"

    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
    sc.pp.scale(adata, zero_center=False, max_value=10)
    sc.pp.pca(adata, n_comps=n_pcs)
    sc.pl.pca_variance_ratio(
        adata,
        n_pcs=n_pcs,
        log=False,
        show=False,
        save=f"_{sample}_{tag}_PCAVar.pdf"
    )

    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=n_neighbors, metric=metric)
    sc.tl.umap(adata)
    return adata

def plot_bar_counts(series, xlab, title, outpath):
    counts = series.value_counts().sort_values(ascending=False).reset_index()
    counts.columns = [xlab, "cell_count"]

    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=counts,
        x=xlab,
        y="cell_count",
        palette="tab20",
        order=counts[xlab],
    )

    for i, v in enumerate(counts["cell_count"]):
        plt.text(i, v + (v * 0.01), str(v), ha='center', fontsize=9)

    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Number of Cells")
    plt.xlabel(xlab)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

def filter_low_confidence_clusters(adata, sample, cluster_key, score_key="label_confidence", threshold=0.90):
    median_key = f"{cluster_key}_medianLC"

    cluster_stats = (
        adata.obs
        .groupby(cluster_key, observed=True)
        .agg(
            median_label_confidence=(score_key, "median"),
            n_cells=(score_key, "size"),
        )
        .sort_values("median_label_confidence", ascending=False)
    )

    print(f"\n=== {sample}: per-cluster median label_confidence and sizes ===")
    print(cluster_stats)

    adata.obs[median_key] = adata.obs[cluster_key].map(
        cluster_stats["median_label_confidence"]
    ).astype(float)

    plt.figure()
    adata.obs[median_key].hist(bins=50)
    plt.title(f"{sample}: histogram of {median_key}")
    plt.xlabel(median_key)
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(f"figures/{sample}_{median_key}_hist.png", dpi=300, bbox_inches='tight')
    plt.close()

    removed = cluster_stats[cluster_stats["median_label_confidence"] < threshold]
    kept = cluster_stats[cluster_stats["median_label_confidence"] >= threshold]

    cells_removed = int(removed["n_cells"].sum())
    cells_kept = int(kept["n_cells"].sum())

    print(f"{sample}: removed clusters = {len(removed)}; cells removed = {cells_removed}")
    print(f"{sample}: kept clusters = {len(kept)}; cells kept = {cells_kept}")

    sc.pl.umap(
        adata,
        color=median_key,
        cmap="viridis",
        vmin=0, vmax=1,
        show=False,
        save=f"_{sample}_{median_key}_ALL.pdf",
    )

    if len(kept) == 0:
        raise ValueError(f"{sample}: all clusters failed the median label_confidence threshold of {threshold}")

    adata_kept = adata[adata.obs[cluster_key].isin(kept.index)].copy()

    sc.pl.umap(
        adata_kept,
        color=median_key,
        cmap="viridis",
        vmin=0, vmax=1,
        show=False,
        save=f"_{sample}_{median_key}_gt{threshold}.pdf",
    )

    return adata_kept, cluster_stats, median_key, cells_removed

def apply_merged_labels(adata):
    adata = adata.copy()

    # safer in a loop than category juggling from the start
    adata.obs["merged_type"] = adata.obs["predicted_labels"].astype(str)

    # remove unwanted types
    adata = adata[~adata.obs["merged_type"].isin(REMOVE_LABELS)].copy()

    # merge fine labels into broader biology
    for merged_label, old_labels in MERGE_RULES.items():
        adata.obs.loc[adata.obs["merged_type"].isin(old_labels), "merged_type"] = merged_label

    return adata

def order_and_color_categories(adata, ct="merged_type"):
    vals = adata.obs[ct].astype(str)
    present = list(pd.Index(vals.unique()))
    extras = sorted([c for c in present if c not in DESIRED_ORDER])
    new_order = [c for c in DESIRED_ORDER if c in present] + extras

    adata.obs[ct] = pd.Categorical(vals, categories=new_order, ordered=True)
    fallback = "#BDBDBD"
    adata.uns[f"{ct}_colors"] = [CELL_TYPE_COLORS.get(c, fallback) for c in new_order]
    return adata

def filter_by_neighbor_agreement(adata, ct="merged_type", thr=0.10, weighted=True):
    adata = adata.copy()

    if "connectivities" not in adata.obsp:
        sc.pp.neighbors(adata, n_neighbors=15)

    if not pd.api.types.is_categorical_dtype(adata.obs[ct]):
        adata.obs[ct] = adata.obs[ct].astype("category")

    if adata.obs[ct].isna().any():
        adata.obs[ct] = adata.obs[ct].cat.add_categories(["Unknown"]).fillna("Unknown")

    W = adata.obsp["connectivities"]
    if not weighted:
        W = (W > 0).astype(np.float32)

    deg = np.asarray(W.sum(axis=1)).ravel()
    n = adata.n_obs

    codes = adata.obs[ct].cat.codes.to_numpy()
    n_cats = len(adata.obs[ct].cat.categories)
    rows = np.arange(n)

    L = sp.csr_matrix(
        (np.ones(n, dtype=np.float32), (rows, codes)),
        shape=(n, n_cats)
    )

    S = W @ L
    same_weight = np.asarray((S.multiply(L)).sum(axis=1)).ravel()

    agreement = np.divide(
        same_weight,
        deg,
        out=np.zeros_like(deg, dtype=float),
        where=deg > 0
    )

    adata.obs["neighbor_agreement"] = agreement

    to_drop = agreement < thr
    print(f"Cells to drop (<{thr:.0%} agreement): {to_drop.sum()} / {n}")

    adata_filtered = adata[~to_drop].copy()
    print(f"Remaining cells: {adata_filtered.n_obs}")

    return adata_filtered, int(to_drop.sum())

#################################################
############# MAIN LOOP ###############
#################################################

summary_rows = []

for in_file in input_files:
    sample = get_sample_name(in_file)
    print("\n" + "=" * 90)
    print(f"Processing sample: {sample}")
    print("=" * 90)

    try:
        #################################################
        ############# LOAD ###############
        #################################################
        adata_query = sc.read_h5ad(in_file)

        if "sample" not in adata_query.obs.columns:
            adata_query.obs["sample"] = sample

        n_start = adata_query.n_obs

        sc.pl.umap(
            adata_query,
            color=["predicted_labels"],
            ncols=1,
            show=False,
            save=f"_{sample}_LabelTransfer_FULLRNA.pdf"
        )

        #################################################
        ############# STEP 1: CELL-LEVEL FILTER ###############
        #################################################
        print(adata_query.obs[['label_confidence', 'scDblFinder_score']].describe())

        sc.pl.violin(
            adata_query,
            ['label_confidence', 'scDblFinder_score'],
            groupby=None,
            jitter=0.4,
            multi_panel=True,
            show=False,
            save=f"_{sample}_LC_Dblt_FULLRNA.pdf"
        )

        low_conf_thr = adata_query.obs['label_confidence'].quantile(LOW_CONF_Q)
        high_dbl_thr = adata_query.obs['scDblFinder_score'].quantile(HIGH_DBL_Q)

        print(f"{sample}: low label_confidence threshold (q={LOW_CONF_Q}) = {low_conf_thr:.3f}")
        print(f"{sample}: high scDblFinder_score threshold (q={HIGH_DBL_Q}) = {high_dbl_thr:.3f}")

        to_remove_mask = (
            (adata_query.obs['label_confidence'] <= low_conf_thr) &
            (adata_query.obs['scDblFinder_score'] >= high_dbl_thr)
        )

        adata_query.obs['low_conf_high_doublet'] = to_remove_mask

        sc.pl.umap(
            adata_query,
            color="low_conf_high_doublet",
            show=False,
            save=f"_{sample}_low_conf_high_doublet_ALL.pdf",
        )

        n_removed_step1 = int(to_remove_mask.sum())
        print(f"{sample}: cells to remove in step 1 = {n_removed_step1} / {adata_query.n_obs}")

        adata_query = adata_query[~to_remove_mask].copy()
        n_after_step1 = adata_query.n_obs
        print(f"{sample}: remaining cells after step 1 = {n_after_step1}")

        #################################################
        ############# STEP 2: RE-EMBED AFTER STEP 1 ###############
        #################################################
        adata_query = recompute_embedding(
            adata_query,
            sample=sample,
            tag="postCellFilter",
            n_top_genes=3000,
            n_pcs=15,
            n_neighbors=30,
            metric='cosine'
        )

        sc.pl.umap(
            adata_query,
            color=["predicted_labels"],
            ncols=1,
            show=False,
            save=f"_{sample}_LabelTransferFiltered_step1_FULLRNA.pdf"
        )

        #################################################
        ############# STEP 3: OVERCLUSTER + MEDIAN LC FILTER ###############
        #################################################
        cluster_key = f'overcluster_{OVERCLUSTER_RES}'
        sc.tl.leiden(adata_query, resolution=OVERCLUSTER_RES, key_added=cluster_key)

        adata_query, cluster_stats, median_key, n_removed_step2 = filter_low_confidence_clusters(
            adata_query,
            sample=sample,
            cluster_key=cluster_key,
            score_key="label_confidence",
            threshold=CLUSTER_CONF_THRESHOLD
        )

        n_after_step2 = adata_query.n_obs

        plot_bar_counts(
            adata_query.obs["predicted_labels"],
            xlab="predicted_labels",
            title=f"{sample}: Cell Counts per Predicted Label",
            outpath=f"figures/{sample}_cell_counts_per_Predicted_label_sorted_FULLRNA.png"
        )

        #################################################
        ############# STEP 4: MERGE ANNOTATIONS ###############
        #################################################
        adata_query = apply_merged_labels(adata_query)
        adata_query = order_and_color_categories(adata_query, ct="merged_type")

        sc.pl.umap(
            adata_query,
            color=[
                "merged_type",
                "predicted_labels",
                "scDblFinder_score",
                "label_confidence",
                "pct_counts_mt",
                "log1p_total_counts",
                "pct_counts_in_top_50_genes"
            ],
            ncols=1,
            show=False,
            save=f"_{sample}_MergedAnnot1_FULLRNA.pdf",
        )

        sc.pl.umap(
            adata_query,
            color=["merged_type"],
            ncols=1,
            show=False,
            save=f"_{sample}_MergedAnnot1_FULLRNA_OnlyMerged.pdf",
        )

        #################################################
        ############# STEP 5: NEIGHBOR AGREEMENT FILTER ###############
        #################################################
        adata_query, n_removed_step3 = filter_by_neighbor_agreement(
            adata_query,
            ct="merged_type",
            thr=NEIGHBOR_AGREEMENT_THRESHOLD,
            weighted=NEIGHBOR_WEIGHTED
        )

        #################################################
        ############# STEP 6: RE-EMBED AFTER STEP 5 ###############
        #################################################
        adata_query = recompute_embedding(
            adata_query,
            sample=sample,
            tag="postNeighborFilter",
            n_top_genes=4000,
            n_pcs=20,
            n_neighbors=20,
            metric='cosine'
        )

        adata_query = order_and_color_categories(adata_query, ct="merged_type")

        sc.pl.umap(
            adata_query,
            color=["predicted_labels"],
            ncols=1,
            show=False,
            save=f"_{sample}_LabelTransferFiltered_step2_FULLRNA.pdf"
        )

        sc.pl.umap(
            adata_query,
            color=["merged_type"],
            ncols=1,
            show=False,
            save=f"_{sample}_MergedAnnot3_FULLRNA.pdf",
        )

        sc.pl.umap(
            adata_query,
            color=["neighbor_agreement", "label_confidence"],
            ncols=1,
            show=False,
            save=f"_{sample}_MergedAnnot3_neighbor_agreement_FULLRNA.pdf",
        )

        sc.pl.umap(
            adata_query,
            color=["merged_type"],
            ncols=1,
            show=False,
            save=f"_{sample}_MergedAnnot_FULLRNA.pdf",
        )

        #################################################
        ############# STEP 7: FINAL CLUSTER-LEVEL LABEL ###############
        #################################################
        sc.tl.leiden(adata_query, key_added="leiden", resolution=FINAL_LEIDEN_RES)

        cluster_assign = (
            adata_query.obs
            .groupby("leiden")["merged_type"]
            .agg(lambda x: x.astype(str).value_counts().idxmax())
        )

        adata_query.obs["leiden_merged_type"] = adata_query.obs["leiden"].map(cluster_assign)

        # match category order/colors to merged_type
        merged_categories = list(adata_query.obs["merged_type"].cat.categories)
        adata_query.obs["leiden_merged_type"] = pd.Categorical(
            adata_query.obs["leiden_merged_type"].astype(str),
            categories=merged_categories,
            ordered=True
        )
        adata_query.uns["leiden_merged_type_colors"] = adata_query.uns["merged_type_colors"].copy()

        sc.pl.umap(
            adata_query,
            color=["leiden_merged_type"],
            ncols=1,
            show=False,
            save=f"_{sample}_MergedAnnotClust_FULLRNA.pdf"
        )

        #################################################
        ############# STEP 8: SAVE FINAL H5AD ###############
        #################################################
        out_h5ad = f"{sample}_Tonsil_FULL_RNA_ANNOTATED.h5ad"
        adata_query.write(out_h5ad)

        #################################################
        ############# STEP 9: QC PLOTS ###############
        #################################################
        sc.pl.violin(adata_query, 'log1p_total_counts', rotation=0.0000001, show=False)

        fig = plt.gcf()
        axes = fig.axes

        for ax in axes:
            yticks = ax.get_yticks()
            new_labels = [f"{np.expm1(tick):.0f}" for tick in yticks]
            ax.set_yticklabels(new_labels)
            ax.set_title(
                f"n_cells: {adata_query.n_obs}\n"
                f"Median UMI: {np.median(adata_query.obs['total_counts']):.0f}"
            )

        plt.savefig(f'./figures/{sample}_rna_nUMIPostFilter.png', dpi=300, bbox_inches='tight')
        plt.close()

        plot_bar_counts(
            adata_query.obs["leiden_merged_type"].astype(str),
            xlab="leiden_merged_type",
            title=f"{sample}: Cell Counts per leiden_merged_type",
            outpath=f"figures/{sample}_cell_counts_per_leiden_merged_type_sorted_FULLRNA.png"
        )

        #################################################
        ############# SUMMARY ###############
        #################################################
        summary_rows.append({
            "sample": sample,
            "input_file": in_file,
            "n_start": n_start,
            "low_conf_thr": low_conf_thr,
            "high_dbl_thr": high_dbl_thr,
            "removed_lowconf_highdbl": n_removed_step1,
            "n_after_step1": n_after_step1,
            "removed_low_medianLC_clusters": n_removed_step2,
            "n_after_step2": n_after_step2,
            "removed_low_neighbor_agreement": n_removed_step3,
            "n_final": adata_query.n_obs,
            "status": "OK",
        })

        print(f"{sample}: finished successfully -> {out_h5ad}")

    except Exception as e:
        print(f"{sample}: FAILED -> {e}")
        summary_rows.append({
            "sample": sample,
            "input_file": in_file,
            "n_start": np.nan,
            "low_conf_thr": np.nan,
            "high_dbl_thr": np.nan,
            "removed_lowconf_highdbl": np.nan,
            "n_after_step1": np.nan,
            "removed_low_medianLC_clusters": np.nan,
            "n_after_step2": np.nan,
            "removed_low_neighbor_agreement": np.nan,
            "n_final": np.nan,
            "status": f"FAILED: {e}",
        })

#################################################
############# SAVE GLOBAL SUMMARY ###############
#################################################

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("figures/RNA_postLT_cleanup_summary.csv", index=False)

print("\nFinal summary:")
print(summary_df)





import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import scanpy.external as sce  # for bbknn

os.makedirs("figures", exist_ok=True)
sc.settings.figdir = "figures"

#################################################
############# INPUTS / SAMPLE RENAMING ##########
#################################################

file_map = {
    "Tonsil_FULL_RNA_ANNOTATED.h5ad": "Tonsil1",
    "Sc_HumanTonsil2_Tonsil_FULL_RNA_ANNOTATED.h5ad": "Tonsil2",
    "Sc_TP3r_TP3_A_Tonsil_FULL_RNA_ANNOTATED.h5ad": "Tonsil3_K27ac",
    "Sc_TP3r_TP3_B_Tonsil_FULL_RNA_ANNOTATED.h5ad": "Tonsil3_K9me3",
    "Sc_TP4r_TP4_A_Tonsil_FULL_RNA_ANNOTATED.h5ad": "Tonsil4_K27ac",
    "Sc_TP4r_TP4_B_Tonsil_FULL_RNA_ANNOTATED.h5ad": "Tonsil4_K9me3",
}

sample_order = [
    "Tonsil1",
    "Tonsil2",
    "Tonsil3_K27ac",
    "Tonsil3_K9me3",
    "Tonsil4_K27ac",
    "Tonsil4_K9me3",
]

cell_type_colors = {
    "B_Cycling":        "#00BFC4",
    "B_GC_DZ":          "#98DF8A",
    "B_GC_LZ":          "#2CA02C",
    "B_mem+activated":  "#2B8C9E",
    "B_naive":          "#1F77B4",
    "B_plasma":         "#76C1FF",
    "DC":               "#B39B00",
    "FDC":              "#7B1FA2",
    "Macrophages":      "#8C564B",
    "T_CD4+":           "#D62728",
    "T_CD4+_FH":        "#E377C2",
    "T_CD4+_naive":     "#FF9896",
    "T_CD8+":           "#FF7F0E",
}

desired_ct_order = [
    "B_Cycling","B_GC_DZ","B_GC_LZ","B_mem+activated","B_naive","B_plasma",
    "DC","FDC","Macrophages",
    "T_CD4+","T_CD4+_FH","T_CD4+_naive","T_CD8+",
]

#################################################
############# LOAD + STANDARDIZE ################
#################################################

adatas = []
for f, sample_name in file_map.items():
    print(f"Loading {f} -> {sample_name}")
    a = sc.read_h5ad(f)

    # make obs names unique across samples
    a.obs_names = [f"{sample_name}_{x}" for x in a.obs_names]
    a.obs_names_make_unique()

    # sample metadata
    a.obs["sample"] = sample_name
    a.obs["tonsil_id"] = sample_name.split("_")[0]  # Tonsil1/2/3/4

    if "K27ac" in sample_name:
        a.obs["mark"] = "K27ac"
    elif "K9me3" in sample_name:
        a.obs["mark"] = "K9me3"
    else:
        a.obs["mark"] = "none"

    # choose a consistent cell-type annotation column
    if "merged_type" not in a.obs.columns:
        if "leiden_merged_type" in a.obs.columns:
            a.obs["merged_type"] = a.obs["leiden_merged_type"].copy()
        elif "predicted_labels" in a.obs.columns:
            a.obs["merged_type"] = a.obs["predicted_labels"].copy()
        else:
            raise ValueError(f"{f} has no merged_type / leiden_merged_type / predicted_labels column")

    # make sure counts layer exists
    if "counts" not in a.layers:
        a.layers["counts"] = a.X.copy()

    adatas.append(a)

#################################################
############# MERGE ON SHARED GENES #############
#################################################

shared_genes = sorted(set(adatas[0].var_names).intersection(*[set(x.var_names) for x in adatas[1:]]))
print(f"Shared genes across all samples: {len(shared_genes)}")

adatas = [x[:, shared_genes].copy() for x in adatas]

adata_merged = ad.concat(
    adatas,
    join="inner",
    merge="same",
    label="orig_sample_file",
    keys=list(file_map.keys()),
    index_unique=None,
)

print(adata_merged)
print(adata_merged.obs["sample"].value_counts())

adata_merged.obs["sample"] = pd.Categorical(
    adata_merged.obs["sample"],
    categories=sample_order,
    ordered=True
)

#################################################
############# ORDER / COLOR CELL TYPES ##########
#################################################

adata_merged.obs["merged_type"] = adata_merged.obs["merged_type"].astype(str)
present_ct = list(pd.Index(adata_merged.obs["merged_type"].unique()))
extra_ct = sorted([x for x in present_ct if x not in desired_ct_order])
new_ct_order = [x for x in desired_ct_order if x in present_ct] + extra_ct

adata_merged.obs["merged_type"] = pd.Categorical(
    adata_merged.obs["merged_type"],
    categories=new_ct_order,
    ordered=True
)

fallback = "#BDBDBD"
adata_merged.uns["merged_type_colors"] = [cell_type_colors.get(c, fallback) for c in new_ct_order]

#################################################
############# PREPROCESS MERGED #################
#################################################

# start from raw counts for joint preprocessing
adata_merged.X = adata_merged.layers["counts"].copy()

sc.pp.normalize_total(adata_merged, target_sum=1e4)
sc.pp.log1p(adata_merged)
adata_merged.layers["log1p_norm"] = adata_merged.X.copy()

# keep raw log-normalized values
adata_merged.raw = adata_merged

# batch-aware HVGs
sc.pp.highly_variable_genes(
    adata_merged,
    layer="counts",
    flavor="seurat_v3",
    n_top_genes=3000,
    batch_key="sample",
    subset=False,
)

print(f"HVGs selected: {adata_merged.var['highly_variable'].sum()}")

sc.pp.scale(adata_merged, zero_center=False, max_value=10)
sc.tl.pca(adata_merged, n_comps=30, use_highly_variable=True)

sc.pl.pca_variance_ratio(
    adata_merged,
    n_pcs=30,
    log=False,
    show=False,
    save="_MergedTonsils_uncorrected_PCAVar.pdf"
)

#################################################
############# UNCORRECTED UMAP / CLUSTERING #####
#################################################

adata_uncorrected = adata_merged.copy()

sc.pp.neighbors(
    adata_uncorrected,
    use_rep="X_pca",
    n_neighbors=20,
    metric="cosine",
)

sc.tl.umap(adata_uncorrected)
sc.tl.leiden(adata_uncorrected, resolution=2, key_added="leiden_uncorrected")

sc.pl.umap(
    adata_uncorrected,
    color=["merged_type"],
    ncols=1,
    show=False,
    save="_MergedTonsils_uncorrected_celltype.pdf"
)

sc.pl.umap(
    adata_uncorrected,
    color=["sample"],
    ncols=1,
    show=False,
    save="_MergedTonsils_uncorrected_sample.pdf"
)

sc.pl.umap(
    adata_uncorrected,
    color=["leiden_uncorrected"],
    ncols=1,
    legend_loc="on data",
    show=False,
    save="_MergedTonsils_uncorrected_leiden.pdf"
)

sc.pl.umap(
    adata_uncorrected,
    color=["merged_type", "sample"],
    ncols=1,
    show=False,
    save="_MergedTonsils_uncorrected_celltype_sample_panel.pdf"
)

adata_uncorrected.write("MergedTonsils_uncorrected.h5ad")






import os
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sp

os.makedirs("figures", exist_ok=True)
sc.settings.figdir = "figures"

adata_filtered = sc.read("MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad")
adata_filtered.X = adata_filtered.layers["SCT_data"]

sc.pp.highly_variable_genes(
    adata_filtered,
    layer='counts',
    flavor='seurat_v3',
    n_top_genes=3000,
)
sc.pp.scale(adata_filtered,zero_center=True, max_value=4) # optional clipping
sc.pp.pca(adata_filtered,n_comps=20)
sc.pl.pca_variance_ratio(adata_filtered, n_pcs=20, log=True, show=False, save=f'Tonsils_RNA_PCAVar.pdf')


sc.pp.neighbors(adata_filtered, use_rep="X_pca", n_neighbors=20,metric='cosine')  # or adjust as needed
sc.tl.umap(adata_filtered)


sc.pl.umap(
    adata_filtered,
    color=["predicted_labels","merged_type","tonsil_id"],
    ncols=1,
    show=False,
    save="LabelTransferFiltered_FULLRNA.pdf"
)

sc.tl.leiden(adata_filtered, key_added="leiden", resolution=4)

# --- Assign clusters to cell types (majority vote) ---
cluster_assign = (
    adata_filtered.obs.groupby("leiden")["merged_type"]
    .agg(lambda x: x.value_counts().idxmax())
)

adata_filtered.obs["leiden_merged_type"] = adata_filtered.obs["leiden"].map(cluster_assign)
adata_filtered.uns[f"leiden_merged_type_colors"] = adata_filtered.uns[f"merged_type_colors"]

sc.pl.umap(
    adata_filtered,
    color=["merged_type","leiden_merged_type"],
    ncols=1,
    show=False,
    save=f"MergedAnnotClust_FULLRNA.pdf"
)

adata_filtered_path = f"MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"
adata_filtered.write(adata_filtered_path)




########################
#### FIGURES FOR PAPER 
########################
import pandas as pd
import matplotlib.pyplot as plt

def sorted_barplot_from_obs(
    adata,
    obs_col,
    color_map,
    title,
    outpath,
    pretty_names=None,
):
    counts = adata.obs[obs_col].astype(str).value_counts().sort_values(ascending=False)
    bar_colors = [color_map.get(x, "#BDBDBD") for x in counts.index]
    plot_labels = [pretty_names.get(x, x) for x in counts.index] if pretty_names else counts.index

    plt.figure(figsize=(8, 6))
    bars = plt.bar(
        plot_labels,
        counts.values,
        color=bar_colors,
        edgecolor="black",
        linewidth=0.5
    )

    offset = max(counts.max() * 0.01, 5)
    for bar, val in zip(bars, counts.values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + offset,
            str(val),
            ha="center",
            va="bottom",
            fontsize=11
        )

    plt.ylabel("Number of cells", fontsize=16)
    plt.xlabel("")
    plt.xticks(rotation=50, ha="right", fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(title, fontsize=18)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


ct_col = "leiden_merged_type"
color_key = f"{ct_col}_colors"

if color_key not in adata_filtered.uns:
    raise KeyError(f"{color_key} not found in adata_filtered.uns")

# map category -> color using the categorical order stored in adata
cat_order = list(adata_filtered.obs[ct_col].cat.categories)
ct_colors = dict(zip(cat_order, adata_filtered.uns[color_key]))

pretty_names = {
    "B_mem": "Memory",
    "B_naive": "Naive",
    "B_GC_LZ": "GC Light Zone",
    "B_GC_DZ": "GC Dark Zone",
    "B_Cycling": "Cycling",
    "T_CD4+": "CD4+",
    "T_CD4+_naive": "CD4+/Naive",
    "T_CD4+_TfH": "CD4+/Follicular Helper",
    "T_CD4+_TfH_GC": "CD4+/Follicular Helper GC",
    "T_CD8+": "CD8+",
    "B_plasma": "Plasma",
    "DC": "Dendritic Cells",
    "FDC": "Follicular Dendritic Cells",
    "Macrophages": "Macrophages",
}

sorted_barplot_from_obs(
    adata=adata_filtered,
    obs_col="leiden_merged_type",
    color_map=ct_colors,
    title="Cells per cell_type",
    outpath="figures/Supp4c.pdf",
    pretty_names=pretty_names,
)

tonsil_colors = {
    "Tonsil1": "#1f77b4",
    "Tonsil2": "#ff7f0e",
    "Tonsil3": "#2ca02c",
    "Tonsil4": "#d62728",
}

sorted_barplot_from_obs(
    adata=adata_filtered,
    obs_col="tonsil_id",
    color_map=tonsil_colors,
    title="Cells per tonsil_id",
    outpath="figures/Supp4_tonsilID_BARPLOT.pdf",
)

sc.pl.umap(
    adata_filtered,
    color=["leiden_merged_type"],
    ncols=1,
    show=False,
    save=f"Supp4d.pdf"
)

sc.pl.umap(
    adata_filtered,
    color=["tonsil_id"],
    ncols=1,
    show=False,
    save=f"Supp4_tonsilID_UMAP.pdf"
)


##############################
############### DNA
##############################




import os
import glob
import csv
import numpy as np
import scanpy as sc
import snapatac2 as snap
import matplotlib.pyplot as plt

os.makedirs("figures", exist_ok=True)
os.makedirs("DNA_processed", exist_ok=True)

srcdir = "DNA"
workdir = "DNA_processed"
maxFrags = 80000
n_jobs = 64

print("IMPORT AND QC PLOTS DNA")

bam_files = sorted(glob.glob(os.path.join(srcdir, "*_NoDup.bam")))

if len(bam_files) == 0:
    raise FileNotFoundError(f"No *_NoDup.bam files found in {srcdir}")

for bam_path in bam_files:
    bam_name = os.path.basename(bam_path)
    exp = bam_name.replace("_NoDup.bam", "")
    frag_file = os.path.join(workdir, f"{exp}.tsv.gz")
    h5ad_file = os.path.join(workdir, f"{exp}.h5ad")
    csv_file = os.path.join(workdir, f"{exp}.csv")

    print(f"\nProcessing {exp}")

    # Make fragment file
    frag_info = snap.pp.make_fragment_file(
        bam_path,
        frag_file,
        is_paired=True,
        barcode_tag="CB",
        shift_left=4,
        shift_right=-5,
        min_mapq=20,
        chunk_size=500000000,
        chrM=["chrM", "M"],
        tempdir=workdir,
        compression_level=1,compression="gzip",
    )

    with open(csv_file, mode="w", newline="") as file:
        writer = csv.writer(file)
        for key, value in frag_info.items():
            writer.writerow([key, value])

    # Import fragments
    adata = snap.pp.import_fragments(
        frag_file,
        snap.genome.GRCh38,
        min_num_fragments=0,
        sorted_by_barcode=True,
        chrM=["chrM", "M"],
        shift_left=0,
        shift_right=0,
        chunk_size=50000,
        tempdir=workdir,
        backend="hdf5",
        n_jobs=n_jobs,
    )

    # Metrics + filtering
    snap.metrics.tsse(adata, gene_anno=snap.genome.GRCh38, inplace=True)
    snap.pp.filter_cells(adata, min_counts=0, max_counts=maxFrags, min_tsse=0, inplace=True)
    snap.metrics.frag_size_distr(adata, add_key="frag_size_distr", inplace=True)

    adata.obs["log10_n_fragment"] = np.log10(adata.obs["n_fragment"])

    # Violin
    sc.pl.violin(adata, "log10_n_fragment", rotation=0.0000001, show=False)
    fig = plt.gcf()
    axes = fig.axes
    for ax in axes:
        yticks = ax.get_yticks()
        new_labels = [f"{10**tick - 1:.0f}" for tick in yticks]
        ax.set_yticklabels(new_labels)
        ax.set_title(f"n_cells: {adata.n_obs}\nMedian n_fragment: {np.median(adata.obs['n_fragment']):.0f}")
    plt.savefig(f"./figures/{exp}_FragViolin.pdf")
    plt.close()

    # Scatter
    sc.pl.scatter(
        adata,
        "log10_n_fragment",
        "tsse",
        title=exp,
        show=False,
        save=f"_{exp}_tsse.pdf",
    )

    # SnapATAC2 QC plots
    snap.pl.tsse(
        adata,
        min_fragment=0,
        width=750,
        height=600,
        interactive=False,
        show=False,
        out_file=f"./figures/{exp}_tsseDensity_PreFilter.pdf",
    )

    snap.pl.frag_size_distr(
        adata,
        width=750,
        height=600,
        interactive=False,
        show=False,
        out_file=f"./figures/{exp}_frag_size_distr_PreFilter.pdf",
    )

    # Add experiment name
    adata.obs["Exp"] = exp

    # Export filtered fragments + save AnnData
    adata.write(h5ad_file)

    print(f"Finished {exp}")



import os
import re
import glob
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import snapatac2 as snap
import matplotlib.pyplot as plt

os.makedirs("figures", exist_ok=True)
os.makedirs("DNA_merged", exist_ok=True)

INPUT_DIR = "DNA_processed"
OUTPUT_DIR = "DNA_merged"

BL = "/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz"
GG = snap.genome.GRCh38

# Map DNA sample ROOTS to the RNA sample namespace.
# For TP3/TP4, the RNA labels K27ac/K9me3 are just names for A/B.
SAMPLE_ROOT_TO_RNA_NAME = {
    "Sc_VTD9_S1_Sc_VTD11_S2_VTHumanTonsil": "Tonsil1",
    "Sc_VTD_HumanTonsil2": "Tonsil2",
    "Sc_TP3D_TP3_A": "Tonsil3_K27ac",
    "Sc_TP3D_TP3_B": "Tonsil3_K9me3",
    "Sc_TP4D_TP4_A": "Tonsil4_K27ac",
    "Sc_TP4D_TP4_B": "Tonsil4_K9me3",
}

FILTER_PARAMS = {
    "H3K27ac":  {"min_counts": 500, "max_counts": 10000, "min_tsse": 0},
    "H3K27me3": {"min_counts": 500, "max_counts": 10000, "min_tsse": 0},
    "H3K9me3":  {"min_counts": 500, "max_counts": 10000, "min_tsse": 0},
}

ANALYSIS_PARAMS = {
    "H3K27ac":  {"n_features": 60000, "bin_size": 1000},
    "H3K27me3": {"n_features": 60000, "bin_size": 5000},
    "H3K9me3":  {"n_features": 60000, "bin_size": 20000},
}


def parse_name(stem):
    m = re.match(r"(.+)_(H3K27ac|H3K27me3|H3K9me3)$", stem)
    if m is None:
        return None, None
    return m.group(1), m.group(2)


def tonsil_id_from_sample(sample_name):
    if sample_name.startswith("Tonsil1"):
        return "Tonsil1"
    if sample_name.startswith("Tonsil2"):
        return "Tonsil2"
    if sample_name.startswith("Tonsil3"):
        return "Tonsil3"
    if sample_name.startswith("Tonsil4"):
        return "Tonsil4"
    return "Unknown"


all_h5ads = sorted(glob.glob(os.path.join(INPUT_DIR, "*.h5ad")))
mark_to_files = {k: [] for k in ANALYSIS_PARAMS.keys()}

for path in all_h5ads:
    stem = os.path.splitext(os.path.basename(path))[0]
    sample_root, mark = parse_name(stem)
    if mark is None:
        print(f"MARK NOT FOUND: {stem}")
        exit()
    if sample_root not in SAMPLE_ROOT_TO_RNA_NAME:
        print(f"Skipping (no RNA-matched sample root): {stem}")
        exit()
    mark_to_files[mark].append(path)

print(mark_to_files)

for mark, files in mark_to_files.items():
    if len(files) == 0:
        continue

    print("\n" + "=" * 80)
    print(f"Processing modality: {mark}")
    print("=" * 80)

    filt = FILTER_PARAMS[mark]
    pars = ANALYSIS_PARAMS[mark]
    experiment = f"MergedTonsils_{mark}"

    filtered_file_list = []
    concat_list = []
    sample_names = []

    # -----------------------------
    # per-sample filtering
    # -----------------------------
    for path in files:
        stem = os.path.splitext(os.path.basename(path))[0]
        sample_root, detected_mark = parse_name(stem)
        sample_name = SAMPLE_ROOT_TO_RNA_NAME[sample_root]

        print(f"Filtering {stem} -> {sample_name}")
        a = snap.read(path, backed=None)

        keep = snap.pp.filter_cells(
            a,
            min_counts=filt["min_counts"],
            max_counts=filt["max_counts"],
            min_tsse=filt["min_tsse"],
            inplace=False,
        )
        a = a[keep].copy()

        a.obs["sample"] = sample_name
        a.obs["sample_root"] = sample_root
        a.obs["tonsil_id"] = tonsil_id_from_sample(sample_name)
        a.obs["mark"] = detected_mark
        a.obs["log10_n_fragment"] = np.log10(a.obs["n_fragment"])

        # Prefix exactly like RNA:
        # Tonsil1_<barcode>, Tonsil3_K27ac_<barcode>, etc.
        a.obs_names = [f"{sample_name}_{x}" for x in a.obs_names.astype(str)]
        a.obs_names_make_unique()

        out_filt = os.path.join(OUTPUT_DIR, f"{stem}_CountFiltered.h5ad")
        a.write(out_filt)

        filtered_file_list.append((sample_name, out_filt))
        sample_names.append(sample_name)

        sc.pl.violin(a, "log10_n_fragment", rotation=0.0000001, show=False)
        fig = plt.gcf()
        for ax in fig.axes:
            yticks = ax.get_yticks()
            ax.set_yticklabels([f"{10**tick - 1:.0f}" for tick in yticks])
            ax.set_title(
                f"{sample_name}\n"
                f"n_cells: {a.n_obs}\n"
                f"Median n_fragment: {np.median(a.obs['n_fragment']):.0f}"
            )
        plt.savefig(f"./figures/{stem}_FragViolin_CountFiltered.pdf")
        plt.close()

        del a

    # -----------------------------
    # export merged fragments
    # -----------------------------
    ds_file = os.path.join(OUTPUT_DIR, f"{experiment}.h5ad")
    ds = snap.AnnDataSet(filtered_file_list, filename=ds_file, add_key="sample")
    ds.obs["Exp"] = [experiment] * ds.shape[0]
    snap.ex.export_fragments(
        ds,
        groupby="Exp",
        suffix=".bed.zst",
        compression_level=1,
    )
    ds.close()

    merged_fragment_file = os.path.join(f"{experiment}.bed.zst")

    # -----------------------------
    # concatenate filtered objects
    # -----------------------------
    for sample_name, filt_file in filtered_file_list:
        a = snap.read(filt_file, backed=None)
        a.obsm.clear()
        concat_list.append(a)

    adata = ad.concat(
        concat_list,
        join="inner",
        label="sample_concat",
        keys=sample_names,
        index_unique=None,
    )

    print(adata)
    print(adata.obs["sample"].value_counts())

    # -----------------------------
    # import merged fragments and sync cells
    # -----------------------------
    b = snap.pp.import_fragments(
        fragment_file=merged_fragment_file,
        chrom_sizes=GG,
        sorted_by_barcode=True,
        min_num_fragments=0,
        tempdir=OUTPUT_DIR,
    )

    common = adata.obs_names.intersection(b.obs_names)
    print(f"Cells kept after obs-name intersection: {len(common)} / {adata.n_obs}")

    adata = adata[common].copy()
    b = b[common].copy()

    adata.obsm = b.obsm.copy()
    adata.uns["reference_sequences"] = b.uns["reference_sequences"].copy()
    del b

    adata.write(os.path.join(OUTPUT_DIR, f"{experiment}_merged_UNprocessed.h5ad"))

    # -----------------------------
    # modality analysis
    # -----------------------------
    snap.pp.add_tile_matrix(
        adata,
        bin_size=pars["bin_size"],
        counting_strategy="paired-insertion",
        chunk_size=500000,
        inplace=True,
    )

    snap.pp.select_features(
        adata,
        n_features=pars["n_features"],
        inplace=True,
        blacklist=BL,
    )

    snap.tl.spectral(
        adata,
        n_comps=20,
        weighted_by_sd=False,
        chunk_size=80000,
        features="selected",
        distance_metric="cosine",
        inplace=True,
    )
    adata.obs["tagmentation"] = adata.obs["tonsil_id"].map({
    "Tonsil1": "single",
    "Tonsil2": "single",
    "Tonsil3": "dual",
    "Tonsil4": "dual",
    })
    import scanpy.external as sce
    #snap.pp.mnc_correct(adata, batch="tagmentation", n_clusters=12,use_rep='X_spectral', key_added='X_spectral')
    #sce.pp.bbknn(
    #adata,
    #batch_key="tagmentation",
    #use_rep="X_spectral",
    #metric="cosine",
    #n_pcs=10,
    #neighbors_within_batch=5,
    #)
    sce.pp.harmony_integrate(
    adata,
    key="tagmentation",
    basis="X_spectral",
    adjusted_basis="X_spectral",
    )
    snap.tl.umap(adata, use_rep="X_spectral", key_added="umap", random_state=None)
    snap.pp.knn(adata, n_neighbors=30, use_rep="X_spectral")
    snap.tl.leiden(adata, resolution=1)

    # -----------------------------
    # plots
    # -----------------------------
    snap.pl.spectral_eigenvalues(
        adata,
        width=600,
        height=400,
        show=False,
        interactive=False,
        out_file=f"./figures/{experiment}_Eigenvalues.png",
    )

    sc.pl.umap(
        adata,
        color=["n_fragment", "log10_n_fragment", "tsse", "tonsil_id"],
        ncols=2,
        show=False,
        save=f"_{experiment}_UMAP_QC.pdf",
    )

    sc.pl.umap(
        adata,
        color=["leiden", "tonsil_id"],
        legend_loc="on data",
        ncols=2,
        show=False,
        save=f"_{experiment}_UMAP_Leiden.pdf",
    )

    adata.write(os.path.join(OUTPUT_DIR, f"{experiment}_merged_processed.h5ad"))
    print(f"Finished {experiment}")










import os
from collections import defaultdict

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

os.makedirs("figures", exist_ok=True)
os.makedirs("overlap_objects", exist_ok=True)

sc.settings.figdir = "figures"

RNA_FILE = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"

DNA_FILES = {
    "H3K27ac":  "DNA_merged/MergedTonsils_H3K27ac_merged_processed.h5ad",
    "H3K27me3": "DNA_merged/MergedTonsils_H3K27me3_merged_processed.h5ad",
    "H3K9me3":  "DNA_merged/MergedTonsils_H3K9me3_merged_processed.h5ad",
}

ct_Cat = "leiden_merged_type"

CELL_TYPE_COLORS = {
    "B_Cycling": "#00BFC4",
    "B_GC_DZ": "#98DF8A",
    "B_GC_LZ": "#2CA02C",
    "B_mem": "#2B8C9E",
    "B_naive": "#1F77B4",
    "B_plasma": "#76C1FF",
    "DC": "#B39B00",
    "FDC": "#7B1FA2",
    "Macrophages": "#8C564B",
    "T_CD4+": "#D62728",
    "T_CD4+_TfH": "#E377C2",
    "T_CD4+_TfH_GC": "#C05AA0",
    "T_CD4+_naive": "#FF9896",
    "T_CD8+": "#FF7F0E",
    "NotInRNA": "#D3D3D3",
}

DESIRED_ORDER = [
    "B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem", "B_naive", "B_plasma",
    "DC", "FDC", "Macrophages",
    "T_CD4+", "T_CD4+_TfH", "T_CD4+_TfH_GC", "T_CD4+_naive", "T_CD8+",
    "NotInRNA"
]

# =============================================================================
# Inline SB rewrite: DT(3nt) -> RNA(4nt)
# Applied only to Tonsil3 / Tonsil4 DNA cells
# =============================================================================

IDX_DT = {
    "GAT": "1",
    "AGT": "2",
    "TCA": "3",
    "ACG": "4",
    "GTC": "5",
    "CTG": "6",
    "GAC": "7",
    "AGA": "8",
}

IDX_RNA = {
    "GCAT": "1",
    "CGAT": "2",
    "GACT": "3",
    "AGCT": "4",
    "CAGT": "5",
    "ACGT": "6",
    "GCTA": "7",
    "CGTA": "8",
}

NUM_TO_RNA = {v: k for k, v in IDX_RNA.items()}
DT_TO_RNA = {dt: NUM_TO_RNA[num] for dt, num in IDX_DT.items()}

INLINE_OFFSET = 1
LEN_DT = 3

for modality, f in DNA_FILES.items():
    if not os.path.exists(f):
        raise FileNotFoundError(f"Missing DNA file: {f}")

rna = sc.read_h5ad(RNA_FILE)

if ct_Cat not in rna.obs.columns:
    raise KeyError(f"{ct_Cat} not found in RNA obs")

rna.obs[ct_Cat] = rna.obs[ct_Cat].astype("category")
if "NotInRNA" not in rna.obs[ct_Cat].cat.categories:
    rna.obs[ct_Cat] = rna.obs[ct_Cat].cat.add_categories(["NotInRNA"])


def _uniquify_names(names):
    s = pd.Series([str(x) for x in names], dtype=object)
    if s.duplicated().any():
        dup = s.duplicated(keep=False)
        s.loc[dup] = s.loc[dup] + "_" + s.loc[dup].groupby(s.loc[dup]).cumcount().astype(str)
    return pd.Index(s.astype(str).to_numpy(dtype=object))


def split_full_name(x):
    sample, barcode = x.rsplit("_", 1)
    return sample, barcode


def overlap_by_sample(rna_names, dna_names):
    rna_by_sample = defaultdict(set)
    dna_by_sample = defaultdict(set)

    for x in rna_names:
        s, b = split_full_name(x)
        rna_by_sample[s].add(b)

    for x in dna_names:
        s, b = split_full_name(x)
        dna_by_sample[s].add(b)

    all_samples = sorted(set(rna_by_sample) | set(dna_by_sample))
    rows = []
    for s in all_samples:
        r = rna_by_sample.get(s, set())
        d = dna_by_sample.get(s, set())
        rows.append({
            "sample": s,
            "n_rna": len(r),
            "n_dna": len(d),
            "n_shared": len(r & d),
        })
    return pd.DataFrame(rows)


def rewrite_dna_inline_sb_to_rna_sb_tonsil34(
    adata,
    inline_offset=INLINE_OFFSET,
    min_recognized_frac=0.95,
    strict=True,
):
    """
    Rewrite obs_names only for DNA cells belonging to Tonsil3* or Tonsil4*.

    Full obs_name format:
      <sample_prefix>_<barcode>

    For Tonsil3/4 DNA barcodes:
      input barcode : <BASE><DT3><rest...>
      output barcode: <BASE><RNA4><rest...>

    Tonsil1/2 are left unchanged.
    """
    adata = adata.copy()

    old_full = pd.Series([str(x) for x in adata.obs_names], index=adata.obs_names, dtype=object)

    sample_prefix = old_full.str.rsplit("_", n=1).str[0]
    barcode = old_full.str.rsplit("_", n=1).str[-1]

    target_mask = sample_prefix.str.startswith(("Tonsil3", "Tonsil4"))

    dt_tag = pd.Series(pd.NA, index=old_full.index, dtype="object")
    rna_tag = pd.Series(pd.NA, index=old_full.index, dtype="object")
    barcode_rewritten = barcode.copy()

    if target_mask.any():
        target_barcode = barcode.loc[target_mask]
        dt_tag.loc[target_mask] = target_barcode.str.slice(inline_offset, inline_offset + LEN_DT)
        rna_tag.loc[target_mask] = dt_tag.loc[target_mask].map(DT_TO_RNA)

        recognized = rna_tag.loc[target_mask].notna()
        frac = float(recognized.mean()) if target_mask.sum() else 1.0

        if frac < min_recognized_frac:
            bad_idx = recognized.index[~recognized.to_numpy()]
            ex = None
            if len(bad_idx) > 0:
                i = bad_idx[0]
                ex = f"{old_full.loc[i]} (DT_tag='{dt_tag.loc[i]}')"
            msg = (
                f"[rewrite_dna_inline_sb_to_rna_sb_tonsil34] "
                f"recognized_frac={frac:.3f} "
                f"({int(recognized.sum())}/{int(target_mask.sum())}). "
                f"Example bad: {ex}"
            )
            if strict:
                raise ValueError(msg)
            else:
                print("WARNING:", msg)

        base_part = target_barcode.str.slice(0, inline_offset)
        rest_part = target_barcode.str.slice(inline_offset + LEN_DT, None)
        barcode_rewritten.loc[target_mask] = (base_part + rna_tag.loc[target_mask] + rest_part).astype(str)

    new_full = sample_prefix + "_" + barcode_rewritten.astype(str)

    adata.obs["sample_prefix"] = sample_prefix.values
    adata.obs["barcode_raw"] = barcode.values
    adata.obs["barcode_rewritten"] = barcode_rewritten.values
    adata.obs["inline_sb_dt"] = dt_tag.values
    adata.obs["inline_sb_rna"] = rna_tag.values
    adata.obs["barcode_was_rewritten"] = target_mask.values

    adata.obs_names = _uniquify_names(new_full.values)
    return adata


def transfer_or_notinrna(target, ref, field):
    common = target.obs_names.intersection(ref.obs_names)

    ref_cats = ref.obs[field].cat.categories.tolist()
    if "NotInRNA" not in ref_cats:
        ref_cats.append("NotInRNA")

    vals = pd.Series(
        pd.Categorical(["NotInRNA"] * target.n_obs, categories=ref_cats),
        index=target.obs_names,
        dtype="category"
    )

    src = ref.obs.loc[common, field].astype(object)
    vals.loc[common] = src.values
    target.obs[field] = vals.astype("category")

    return target, common


def apply_order_and_colors(ad, field, desired_order, color_dict):
    present = set(ad.obs[field].astype(str).unique())
    ordered = [x for x in desired_order if x in present]
    extras = sorted([x for x in present if x not in ordered])
    cats = ordered + extras

    ad.obs[field] = pd.Categorical(
        ad.obs[field].astype(str),
        categories=cats,
        ordered=True
    )
    ad.uns[f"{field}_colors"] = [str(color_dict.get(c, "#BDBDBD")) for c in cats]
    return ad


summary_rows = []

for modality, dna_file in DNA_FILES.items():
    print("\n" + "=" * 90)
    print(f"Processing {modality}")
    print("=" * 90)

    dna_orig = sc.read_h5ad(dna_file)

    rna_barcodes = set(rna.obs_names)
    dna_barcodes_before = set(dna_orig.obs_names)
    shared_before = rna_barcodes & dna_barcodes_before

    print("\nBefore barcode rewrite")
    print(f"RNA:    {len(rna_barcodes)}")
    print(f"DNA:    {len(dna_barcodes_before)}")
    print(f"Shared: {len(shared_before)}")
    print(overlap_by_sample(rna_barcodes, dna_barcodes_before).to_string(index=False))

    dna = rewrite_dna_inline_sb_to_rna_sb_tonsil34(
        dna_orig,
        inline_offset=INLINE_OFFSET,
        min_recognized_frac=0.95,
        strict=True,
    )

    dna_barcodes_after = set(dna.obs_names)
    shared_after = rna_barcodes & dna_barcodes_after

    print("\nAfter barcode rewrite")
    print(f"RNA:    {len(rna_barcodes)}")
    print(f"DNA:    {len(dna_barcodes_after)}")
    print(f"Shared: {len(shared_after)}")
    print(overlap_by_sample(rna_barcodes, dna_barcodes_after).to_string(index=False))

    if ct_Cat in dna.obs.columns:
        del dna.obs[ct_Cat]
    if f"{ct_Cat}_colors" in dna.uns:
        del dna.uns[f"{ct_Cat}_colors"]

    dna, common = transfer_or_notinrna(dna, rna, ct_Cat)
    dna = apply_order_and_colors(dna, ct_Cat, DESIRED_ORDER, CELL_TYPE_COLORS)

    print("\nLabel transfer after barcode rewrite")
    print(f"Matched to RNA: {len(common)}")
    print(f"Not in RNA: {(dna.obs[ct_Cat].astype(str) == 'NotInRNA').sum()}")

    if "X_umap" in dna.obsm:
        sc.pl.umap(
            dna,
            color=ct_Cat,
            title=f"{modality} - Cell types from RNA (barcode-rewritten)",
            show=False,
            save=f"_MergedTonsils_{modality}_UMAP_TypesFromRNA_barcodeRewritten.pdf"
        )

    dna_labeled_out = (
        f"DNA_merged/MergedTonsils_{modality}_merged_processed_barcodeRewritten_withRNAlabels.h5ad"
    )
    dna.write(dna_labeled_out)
    print(f"Saved labeled DNA: {dna_labeled_out}")

    # Pairwise Venn after rewrite
    plt.figure(figsize=(7, 7))
    venn2([rna_barcodes, dna_barcodes_after], set_labels=("RNA", modality))
    plt.title(f"Passing cells overlap: RNA vs {modality} (barcode-rewritten)")
    plt.tight_layout()
    plt.savefig(f"figures/MergedTonsils_RNA_vs_{modality}_Venn_barcodeRewritten.pdf")
    plt.close()

    # Pairwise overlap objects after rewrite
    shared_sorted = sorted(shared_after)

    rna_overlap = rna[shared_sorted].copy()
    dna_overlap = dna[shared_sorted].copy()

    rna_out = f"overlap_objects/MergedTonsils_RNA_overlap_{modality}_barcodeRewritten.h5ad"
    dna_out = f"overlap_objects/MergedTonsils_{modality}_overlap_RNA_barcodeRewritten.h5ad"

    rna_overlap.write(rna_out)
    dna_overlap.write(dna_out)

    summary_rows.append({
        "modality": modality,
        "n_rna": len(rna_barcodes),
        "n_dna_before": len(dna_barcodes_before),
        "n_shared_before": len(shared_before),
        "n_dna_after": len(dna_barcodes_after),
        "n_shared_after": len(shared_after),
        "dna_labeled_file": dna_labeled_out,
        "rna_overlap_file": rna_out,
        "dna_overlap_file": dna_out,
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("figures/MergedTonsils_pairwise_RNA_DNA_overlap_summary_barcodeRewritten.csv", index=False)

print("\nDone.")
print(summary_df)













import os
import scanpy as sc
import snapatac2 as snap

os.makedirs("figures", exist_ok=True)

BL = "/mnt/dataFast/ahrmad/hg38-blacklist.bed.gz"

PAIR_FILES = {
    "H3K27ac": {
        "rna": "overlap_objects/MergedTonsils_RNA_overlap_H3K27ac_barcodeRewritten.h5ad",
        "dna": "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_barcodeRewritten.h5ad",
        "dna_n_features": 100000,
        "dna_bin_size": 600,
    },
    "H3K27me3": {
        "rna": "overlap_objects/MergedTonsils_RNA_overlap_H3K27me3_barcodeRewritten.h5ad",
        "dna": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_barcodeRewritten.h5ad",
        "dna_n_features": 80000,
        "dna_bin_size": 40000,
    },
    "H3K9me3": {
        "rna": "overlap_objects/MergedTonsils_RNA_overlap_H3K9me3_barcodeRewritten.h5ad",
        "dna": "overlap_objects/MergedTonsils_H3K9me3_overlap_RNA_barcodeRewritten.h5ad",
        "dna_n_features": 80000,
        "dna_bin_size": 40000,
    },
}

for modality, cfg in PAIR_FILES.items():
    print("\n" + "=" * 80)
    print(f"Processing {modality}")
    print("=" * 80)

    # =================================================
    # RNA overlap reclustering
    # =================================================
    rna = sc.read_h5ad(cfg["rna"])

    sc.pp.highly_variable_genes(
        rna,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=2000,
    )
    sc.pp.scale(rna, zero_center=False, max_value=4)
    sc.pp.pca(rna, n_comps=10)

    sc.pl.pca_variance_ratio(
        rna,
        n_pcs=30,
        log=False,
        show=False,
        save=f"_MergedTonsils_RNA_{modality}_PCAVar_overlap.pdf",
    )

    sc.pp.neighbors(rna, use_rep="X_pca", n_neighbors=30, metric="cosine")
    sc.tl.umap(rna)

    sc.pl.umap(
        rna,
        color=["leiden_merged_type"],
        ncols=1,
        show=False,
        save=f"_MergedTonsils_RNA_{modality}_Annot_overlap.pdf",
    )

    rna_out = cfg["rna"].replace(".h5ad", "_Reclust.h5ad")
    rna.write(rna_out)
    print(f"Saved RNA: {rna_out}")

    # =================================================
    # DNA overlap reclustering
    # =================================================
    dna = snap.read(cfg["dna"], backed=None)

    # add tile matrix only if needed
    if "tile_matrix" not in dna.obsm and getattr(dna, "X", None) is None:
        snap.pp.add_tile_matrix(
            dna,
            bin_size=cfg["dna_bin_size"],
            counting_strategy="paired-insertion",
            chunk_size=500000,
            inplace=True,
        )

    snap.pp.select_features(
        dna,
        n_features=cfg["dna_n_features"],
        inplace=True,
        blacklist=BL,
        filter_lower_quantile=0.001,
        filter_upper_quantile=0.001,
    )

    snap.tl.spectral(
        dna,
        n_comps=10,
        weighted_by_sd=False,
        chunk_size=80000,
        features="selected",
        distance_metric="cosine",
        inplace=True,
    )

    n_dims = dna.obsm["X_spectral"].shape[1]
    snap.pp.mnc_correct(dna, batch="tagmentation", n_clusters=12,use_rep='X_spectral', key_added='X_spectral')
    #snap.pp.harmony(
    #dna,
    #batch="tagmentation",
    #use_rep="X_spectral",
    #key_added="X_spectral",
    #use_dims=list(range(1, n_dims-1)),
    #)
    #snap.pp.scanorama_integrate(
    #dna,
    #batch="tagmentation",
    #use_rep="X_spectral",
    #key_added="X_spectral",
    #use_dims=list(range(1, n_dims-1)),
    #batch_size=50000,
    #)
    snap.tl.umap(
        dna,
        use_rep="X_spectral",
        key_added="umap",
        random_state=None,)

#        use_dims=list(range(1, n_dims-1)),

    snap.pp.knn(dna, n_neighbors=20, use_rep="X_spectral", method="kdtree")
    snap.tl.leiden(dna, resolution=0.6)

    snap.pl.spectral_eigenvalues(
        dna,
        width=600,
        height=400,
        show=False,
        interactive=False,
        out_file=f"./figures/MergedTonsils_{modality}_DNA_Eigenvalues_overlap.png",
    )

    sc.pl.umap(
        dna,
        color=["leiden_merged_type"],
        show=False,
        save=f"_MergedTonsils_{modality}_DNA_UMAP_LeidenTYPE_overlap.pdf",
    )

    sc.pl.umap(
        dna,
        color=["leiden", "leiden_merged_type"],
        show=False,
        save=f"_MergedTonsils_{modality}_DNA_UMAP_Leiden_overlap.pdf",
    )

    sc.pl.umap(
        dna,
        color=["n_fragment", "log10_n_fragment", "tsse"],
        show=False,
        save=f"_MergedTonsils_{modality}_DNA_UMAP_QC_overlap.png",
    )

    dna_out = cfg["dna"].replace(".h5ad", "_Reclust.h5ad")
    dna.write(dna_out)
    print(f"Saved DNA: {dna_out}")










import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from anndata import AnnData
from typing import Union
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib import colors as mcolors

os.makedirs("figures", exist_ok=True)

# -----------------------------
# Inputs
# -----------------------------
RNA_FILE = "MergedTonsils_Tonsil_FULL_RNA_ANNOTATED.h5ad"

DNA_FILES = {
    "H3K27ac":  "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_barcodeRewritten_Reclust.h5ad",
    "H3K27me3": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
    "H3K9me3":  "overlap_objects/MergedTonsils_H3K9me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
}

CELL_TYPE_COLORS = {
    "B_Cycling": "#00BFC4",
    "B_GC_DZ": "#98DF8A",
    "B_GC_LZ": "#2CA02C",
    "B_mem": "#2B8C9E",
    "B_naive": "#1F77B4",
    "B_plasma": "#76C1FF",
    "DC": "#B39B00",
    "FDC": "#7B1FA2",
    "Macrophages": "#8C564B",
    "T_CD4+": "#D62728",
    "T_CD4+_TfH": "#E377C2",
    "T_CD4+_TfH_GC": "#C05AA0",
    "T_CD4+_naive": "#FF9896",
    "T_CD8+": "#FF7F0E",
}

DESIRED_ORDER = [
    "B_Cycling", "B_GC_DZ", "B_GC_LZ", "B_mem", "B_naive", "B_plasma",
    "DC", "FDC", "Macrophages",
    "T_CD4+", "T_CD4+_TfH", "T_CD4+_TfH_GC", "T_CD4+_naive", "T_CD8+",
]

TONSIL_COLORS = {
    "Tonsil1": "#1f77b4",
    "Tonsil2": "#ff7f0e",
    "Tonsil3": "#2ca02c",
    "Tonsil4": "#d62728",
}

TONSIL_ORDER = ["Tonsil1", "Tonsil2", "Tonsil3", "Tonsil4"]

# rasterized artists inside PDF
PDF_RASTER_DPI = 600


# -----------------------------
# Helpers
# -----------------------------
def _ensure_adata(x: Union[str, AnnData]) -> AnnData:
    if isinstance(x, AnnData):
        return x
    if isinstance(x, str):
        return sc.read_h5ad(x)
    raise TypeError("Input must be an AnnData or a path to .h5ad")


def procrustes_align_to_base(target_xy, base_xy):
    """
    Align target to base using rotation/reflection + global scale, then
    translate into base coordinates.
    """
    X = np.asarray(base_xy, dtype=float)
    Y = np.asarray(target_xy, dtype=float)

    Xmean = X.mean(0, keepdims=True)
    Ymean = Y.mean(0, keepdims=True)
    Xc = X - Xmean
    Yc = Y - Ymean

    M = Yc.T @ Xc
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    R = U @ Vt
    s = S.sum() / (Yc ** 2).sum()

    return s * (Yc @ R) + Xmean


def extent(arr):
    arr = np.asarray(arr)
    return float(arr[:, 0].max() - arr[:, 0].min()), float(arr[:, 1].max() - arr[:, 1].min())


def ordered_categories(series, desired_order):
    vals = pd.Index(series.astype(str).unique())
    ordered = [x for x in desired_order if x in vals]
    extras = sorted([x for x in vals if x not in ordered])
    return ordered + extras


def make_legend(ax, palette_dict, order, title):
    handles = [
        Line2D([0], [0], marker='o', linestyle='None',
               markerfacecolor=palette_dict[k], markeredgecolor='none',
               markersize=7, label=k)
        for k in order if k in palette_dict
    ]
    leg = ax.legend(
        handles=handles,
        title=title,
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        frameon=False,
        fontsize=9
    )
    leg.get_title().set_fontsize(10)


def build_linked_layout(
    adata_rna,
    dna_dict,
    center_field="leiden_merged_type",
):
    """
    Returns:
      rna_center : RNA subset used in center panel
      aligned    : dict of modality -> DataFrame with columns:
                   x_rna, y_rna, x_mod, y_mod, cell
    """
    rna = _ensure_adata(adata_rna).copy()
    dna_loaded = {k: _ensure_adata(v).copy() for k, v in dna_dict.items()}

    if "X_umap" not in rna.obsm:
        raise ValueError("RNA object has no X_umap")

    for k, ad in dna_loaded.items():
        if "X_umap" not in ad.obsm:
            raise ValueError(f"{k} object has no X_umap")

    # union of all RNA cells shared with at least one modality
    union_shared = set()
    common_per_mod = {}

    for mod, ad in dna_loaded.items():
        common = rna.obs_names.intersection(ad.obs_names)
        common_per_mod[mod] = common
        union_shared.update(common)

    union_shared = pd.Index(sorted(union_shared))
    if len(union_shared) == 0:
        raise ValueError("No RNA cells shared with any DNA modality")

    rna_center = rna[union_shared].copy()

    # central RNA coordinates
    rna_xy = pd.DataFrame(
        rna_center.obsm["X_umap"],
        index=rna_center.obs_names,
        columns=["x_rna", "y_rna"]
    )

    # align each modality to RNA shared cells
    aligned = {}
    for mod, ad in dna_loaded.items():
        common = common_per_mod[mod]
        if len(common) == 0:
            continue

        rna_sub = rna[common].copy()
        dna_sub = ad[common].copy()

        base_xy = rna_sub.obsm["X_umap"]
        mod_xy = dna_sub.obsm["X_umap"]

        mod_aligned = procrustes_align_to_base(mod_xy, base_xy)

        df = pd.DataFrame(
            {
                "x_rna": base_xy[:, 0],
                "y_rna": base_xy[:, 1],
                "x_mod": mod_aligned[:, 0],
                "y_mod": mod_aligned[:, 1],
            },
            index=common
        )
        df["cell"] = df.index.astype(str)
        aligned[mod] = df

    return rna_center, dna_loaded, aligned


def plot_linked_umaps_hub(
    adata_rna,
    dna_dict,
    color_field,
    palette_dict,
    desired_order,
    output_prefix,
    panel_positions=None,
    panel_scales=None,
    fig_size=(14, 12),
    fig_dpi=300,
    point_size=18,
    line_width=0.25,
    line_alpha=0.18,
):
    """
    Hub-and-spoke linked UMAP:
      RNA in center, DNA modalities placed around it.
    """
    rna_center, dna_loaded, aligned = build_linked_layout(adata_rna, dna_dict)

    if color_field not in rna_center.obs.columns:
        raise KeyError(f"{color_field} not found in RNA obs")

    for mod, ad in dna_loaded.items():
        if color_field not in ad.obs.columns:
            raise KeyError(f"{color_field} not found in {mod} obs")

    # RNA center, K27 panels a bit wider, K9 panel closer
    if panel_positions is None:
        panel_positions = {
            "RNA":      (0.0,  0.0),
            "H3K27ac":  (-1.30,  0.45),
            "H3K27me3": ( 1.30,  0.45),
            "H3K9me3":  ( 0.0, -1.05),
        }

    # make K27 panels slightly larger than K9
    if panel_scales is None:
        panel_scales = {
            "H3K27ac": 1.12,
            "H3K27me3": 1.12,
            "H3K9me3": 0.92,
        }

    # Determine reference width/height from RNA
    rna_xy = pd.DataFrame(
        rna_center.obsm["X_umap"],
        index=rna_center.obs_names,
        columns=["x", "y"]
    )
    w, h = extent(rna_xy[["x", "y"]].values)
    scale_x = max(w, 1e-6)
    scale_y = max(h, 1e-6)

    cats_present = ordered_categories(rna_center.obs[color_field], desired_order)
    fallback = "#BDBDBD"

    def get_colors(series):
        return [palette_dict.get(str(v), fallback) for v in series.astype(str)]

    # Central RNA shifted coordinates
    x0, y0 = panel_positions["RNA"]
    rna_plot = rna_xy.copy()
    rna_plot["x"] = rna_plot["x"] + x0 * scale_x
    rna_plot["y"] = rna_plot["y"] + y0 * scale_y

    fig, ax = plt.subplots(figsize=fig_size, dpi=fig_dpi)

    panel_titles = {"RNA": "RNA", "H3K27ac": "H3K27ac", "H3K27me3": "H3K27me3", "H3K9me3": "H3K9me3"}

    for mod, df in aligned.items():
        if mod not in panel_positions:
            continue

        dx, dy = panel_positions[mod]
        scale_mod = panel_scales.get(mod, 1.0)

        dfp = df.copy()

        # RNA anchor points
        dfp["x_rna_p"] = dfp["x_rna"] + panel_positions["RNA"][0] * scale_x
        dfp["y_rna_p"] = dfp["y_rna"] + panel_positions["RNA"][1] * scale_y

        # Scale modality panel around its own centroid
        x_center = dfp["x_mod"].mean()
        y_center = dfp["y_mod"].mean()
        dfp["x_mod_scaled"] = (dfp["x_mod"] - x_center) * scale_mod + x_center
        dfp["y_mod_scaled"] = (dfp["y_mod"] - y_center) * scale_mod + y_center

        # Shift scaled modality panel to target location
        dfp["x_mod_p"] = dfp["x_mod_scaled"] + dx * scale_x
        dfp["y_mod_p"] = dfp["y_mod_scaled"] + dy * scale_y

        line_colors = [
            mcolors.to_rgba(palette_dict.get(str(v), fallback), alpha=line_alpha)
            for v in rna_center.obs.loc[dfp.index, color_field].astype(str)
        ]

        seg = np.stack(
            [dfp[["x_rna_p", "y_rna_p"]].values, dfp[["x_mod_p", "y_mod_p"]].values],
            axis=1
        )
        lc = LineCollection(seg, colors=line_colors, linewidths=line_width, zorder=0)
        lc.set_rasterized(True)
        ax.add_collection(lc)

        mod_colors = get_colors(dna_loaded[mod].obs.loc[dfp.index, color_field])
        ax.scatter(
            dfp["x_mod_p"],
            dfp["y_mod_p"],
            s=point_size,
            c=mod_colors,
            edgecolors="none",
            alpha=0.95,
            zorder=2,
            rasterized=True
        )

        ax.text(
            dfp["x_mod_p"].median(),
            dfp["y_mod_p"].max() + 0.08 * scale_y,
            panel_titles.get(mod, mod),
            ha="center",
            va="bottom",
            fontsize=15,
            weight="bold"
        )

    # RNA center points
    rna_colors = get_colors(rna_center.obs[color_field])
    ax.scatter(
        rna_plot["x"],
        rna_plot["y"],
        s=point_size,
        c=rna_colors,
        edgecolors="none",
        alpha=0.98,
        zorder=3,
        rasterized=True
    )

    ax.text(
        rna_plot["x"].median(),
        rna_plot["y"].max() + 0.08 * scale_y,
        "RNA",
        ha="center",
        va="bottom",
        fontsize=16,
        weight="bold"
    )

    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.03)
    ax.axis("off")

    make_legend(ax, palette_dict, cats_present, title=color_field)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}.pdf", bbox_inches="tight", dpi=PDF_RASTER_DPI)
    plt.savefig(f"{output_prefix}.png", bbox_inches="tight", dpi=300)
    plt.close()

    return {
        "pdf": f"{output_prefix}.pdf",
        "png": f"{output_prefix}.png",
    }


# -----------------------------
# Make the figures
# -----------------------------
res_ct1 = plot_linked_umaps_hub(
    adata_rna=RNA_FILE,
    dna_dict=DNA_FILES,
    color_field="leiden_merged_type",
    palette_dict=CELL_TYPE_COLORS,
    desired_order=DESIRED_ORDER,
    output_prefix="figures/LinkedUMAPs_RNA_center_DNA_byCellType1",
)

res_ct2 = plot_linked_umaps_hub(
    adata_rna=RNA_FILE,
    dna_dict=DNA_FILES,
    color_field="leiden_merged_type",
    palette_dict=CELL_TYPE_COLORS,
    desired_order=DESIRED_ORDER,
    output_prefix="figures/LinkedUMAPs_RNA_center_DNA_byCellType2",
    panel_positions={
        "RNA":      (0.0,  0.0),
        "H3K27ac":  (-1.35,  0.45),
        "H3K27me3": ( 1.35,  0.45),
        "H3K9me3":  ( 0.0, -0.95),
    },
    panel_scales={
        "H3K27ac": 1.15,
        "H3K27me3": 1.15,
        "H3K9me3": 0.90,
    },
)

res_tonsil = plot_linked_umaps_hub(
    adata_rna=RNA_FILE,
    dna_dict=DNA_FILES,
    color_field="tonsil_id",
    palette_dict=TONSIL_COLORS,
    desired_order=TONSIL_ORDER,
    output_prefix="figures/LinkedUMAPs_RNA_center_DNA_byTonsil",
)

print(res_ct1["pdf"], res_ct1["png"])
print(res_ct2["pdf"], res_ct2["png"])
print(res_tonsil["pdf"], res_tonsil["png"])


import scanpy as sc

DNA_OVERLAP_FILES = {
    "H3K27ac":  "overlap_objects/MergedTonsils_H3K27ac_overlap_RNA_barcodeRewritten_Reclust.h5ad",
    "H3K27me3": "overlap_objects/MergedTonsils_H3K27me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
    "H3K9me3":  "overlap_objects/MergedTonsils_H3K9me3_overlap_RNA_barcodeRewritten_Reclust.h5ad",
}

for modality, path in DNA_OVERLAP_FILES.items():
    ad = sc.read_h5ad(path)

    sc.pl.umap(
        ad,
        color=["leiden_merged_type"],
        show=False,
        save=f"_MergedTonsils_{modality}_DNA_overlap_leidenMergedType.pdf",
    )

    sc.pl.umap(
        ad,
        color=["tonsil_id"],
        show=False,
        save=f"_MergedTonsils_{modality}_DNA_overlap_tonsilID.pdf",
    )

    sc.pl.umap(
        ad,
        color=["leiden_merged_type", "tonsil_id"],
        ncols=2,
        show=False,
        save=f"_MergedTonsils_{modality}_DNA_overlap_leidenMergedType_tonsilID.pdf",
    )




