#### Visium HD Colon Cancer Spatial Transcriptomics Analysis
This repository contains a complete Visium HD spatial transcriptomics workflow in Python, applied to colon cancer and matched normal tissues. It starts from Space Ranger outputs, converts them into a unified SpatialData object, performs quality control, normalization, dimensionality reduction, clustering, batch correction, and marker gene discovery, and generates interpretable spatial visualizations.

### 1. Project Overview
This notebook demonstrates how to:
- Ingest 10x Genomics Visium HD outputs (Space Ranger).
- Convert them into Zarr and construct a SpatialData object.
- Integrate multiple samples (tumor vs normal) into a single analysis.
- Apply QC filtering on segmented cells.
- Run normalization, PCA, neighborhood graph, UMAP, and clustering.
- Perform batch correction across samples.
- Identify cluster-specific marker genes.
- Visualize H&E images, tissue segmentation, and spatial gene/cluster maps.

## The workflow is built around the colon cancer dataset:
Tumor samples: Cancer_P1, Cancer_P2
Normal samples: Normal_P3, Normal_P5

### 2. Input Data & Space Ranger Outputs
The analysis assumes standard Space Ranger outputs are available for each sample, including:
- web_summary.html
- spatial/aligned_fiducials.jpg
- spatial/detected_tissue_image.jpg
- spatial/tissue_hires_image.png
- spatial/scalefactors_json.json
- filtered_feature_bc_matrix/
- Segmentation-derived count matrices / Zarr outputs
These are subsequently converted to Zarr and read into spatialdata.

### 3. Environment & Dependencies
python3.10 -m venv spatialenv
source spatialData/bin/activate
jupyter notebook
Core Python packages used in the notebook:
- numpy
- pandas
- matplotlib
- scanpy
- scanpy.external (for batch correction / integration)
- spatialdata
- spatialdata-io
- spatialdata_plot
- shapely
- geopandas
- geosketch
- Pillow
- pydeseq2 (optional downstream DE)

You can install a compatible environment (example):
pip install numpy pandas matplotlib scanpy spatialdata spatialdata-io shapely geopandas geosketch pydeseq2 pillow
(Adjust versions to match what was used when this notebook was created.)

### 4. Workflow Summary

## 4.1. Zarr Conversion & SpatialData Construction
For each sample:
Convert Space Ranger output to Zarr.
Load with:
import spatialdata as spd

sdata = spd.read_zarr("<sample_zarr_path>")
Ensure unique feature names and annotate sample identity in each AnnData table.
Concatenate all samples into a single SpatialData object:
concatenated_sdata = spd.concatenate(sdatas, concatenate_tables=True)
This produces a unified object with:
High-res H&E / tissue images
Segmentation-based count matrices
Sample metadata

## 4.2. Quality Control (QC)
The notebook applies QC at the segmented cell/bin level (on segmentation_counts), including:
Minimum counts / features per cell
Filtering extreme total counts / mitochondrial content
Removing low-quality cells to improve downstream clustering
Customize thresholds directly in the QC cells to match your dataset and platform metrics.

## 4.3. Normalization & Dimensionality Reduction
On the filtered object:
Library size normalization
Log1p transform
Highly variable gene selection
PCA for dimensionality reduction
Neighborhood graph construction
These steps are implemented with standard Scanpy functions on concatenated_sdata["segmentation_counts"].

## 4.4. Clustering & UMAP
Using the neighborhood graph:
Compute Leiden clustering to identify spatial cell states / niches.
Embed cells into UMAP for visualization.
Visualize clusters colored by:
Cluster ID
Sample (tumor vs normal)
Marker gene expression
This helps distinguish immune, stromal, tumor, and other niches in the tissue.

## 4.5. Batch Correction
If needed, batch / sample effects are corrected across samples with external integration tools (via scanpy.external). The corrected embeddings are then reused for:
Clustering
UMAP visualization
Downstream marker analyses
The README assumes moderate batch effects; you can toggle correction steps depending on your design.

## 4.6. Marker Gene Identification
Cluster-specific marker genes are derived using:
sc.tl.rank_genes_groups(
    adata=concatenated_sdata["segmentation_counts"],
    groupby="clusters",
    method="wilcoxon"
)
Top markers are exported to marker_genes_pval.csv and visualized via:
Dot plots
Heatmaps
Cluster-wise marker expression summaries

### 5. Expected Outputs
Running the notebook will generate:
- A concatenated SpatialData object for all samples
- QC-filtered AnnData tables
- UMAP plots colored by sample, cluster, and key marker genes
- Spatial feature maps overlaying expression/clusters on the tissue
- marker_genes_pval.csv: cluster-level marker statistics

### 6. Figures & Suggested Image Layout
To make the repository visually informative, add a figures/ folder and include key panels. Suggested filenames and markdown snippets (you can export these directly from the notebook’s plotting cells):

## 6.1. H&E Tissue Overview
A representative high-resolution H&E image from Visium HD:
![H&E-stained Visium HD tissue section](figures/01_he_tissue_overview.png)

*Example H&E image from a Visium HD colon cancer section used for spatial transcriptomic profiling.*

## 6.2. Tissue Detection / Segmentation
Overlay of detected tissue / segmented bins or cells:
![Tissue detection and segmentation](figures/02_tissue_detection_segmentation.png)

*Detected tissue mask and segmented units used for downstream quantification.*
## 6.3. Spatial Expression Map
Spatial feature plot of a canonical gene (e.g. epithelial, immune, or tumor marker) projected on the tissue:
![Spatial gene expression map](figures/03_spatial_gene_expression.png)

*Spatial distribution of a selected marker gene overlaid on the H&E image.*
## 6.4. Spatial Clusters / Niches
Clustered spatial map showing distinct niches:
![Spatial clustering of Visium HD data](figures/04_spatial_clusters_umap.png)

*Leiden clusters projected in spatial coordinates, highlighting distinct tumor–immune–stromal niches.*
Once you generate these panels in the notebook, save them into figures/ with the above filenames so the README renders them automatically.

### 7. How to Run
Clone/download this repository.
Place your Space Ranger-derived Zarr folders (e.g. Colon_Cancer_P1, etc.) in the project directory.
Create and activate a Python environment with required packages.
Open the notebook:
jupyter lab VisiumHD_python_analysis.ipynb
Run cells sequentially:
Inspect Space Ranger outputs & tissue detection
Load Zarr → build SpatialData
Perform QC, normalization, clustering
Run marker discovery & export results
Save figures into figures/ for README visualization
