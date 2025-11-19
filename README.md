# DGE-Analysis-of-GSE153089-Using-DESeq2-with-PCA-Heatmap-and-Volcano-Plot-Visualization.
DGE Analysis of GSE153089 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization.
# 🔬 DESeq2 DGE Pipeline: Hepatoblastoma vs. Normal Liver (GSE153089)

This R script automates a complete bioinformatics workflow for **Differential Gene Expression (DGE)** analysis of the **GSE153089** dataset, which investigates microRNA clusters in hepatoblastoma. The pipeline specifically compares **Hepatoblastoma Tumor** samples against **Normal Liver** tissue.

The pipeline utilizes the **`DESeq2`** package for count-based DGE and provides essential steps for automated data fetching, sample grouping, QC, and comprehensive visualization.

## 🚀 Key Features

* **Automated Data Fetching:** Downloads expression data and metadata directly from **GEO (GSE153089)** using `GEOquery`.
* **Robust Sample Grouping:** Correctly identifies and groups samples into **"Cancer"** and **"Normal"** categories based on the `characteristics_ch1` metadata column (`tissue: normal liver`).
* **DESeq2 Analysis:** Implements the standard `DESeq2` workflow, including count rounding and low-count gene filtering.
* **QC Transformation:** Applies **Variance Stabilizing Transformation (VST)** for accurate distance-based visualization.
* **Integrated Visualization:** Generates three essential plots for interpreting the results: **PCA**, **Heatmap**, and **Volcano Plot**.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE153089 | Study of the DLK1-DIO3 imprinted region microRNA cluster in hepatoblastoma. |
| **DGE Tool** | `DESeq2` | Statistical method optimized for count data (RNA-Seq). |
| **Comparison** | Cancer vs. Normal | Identifies the gene expression signature associated with hepatoblastoma development. |
| **Significance** | $\text{pCutoff} = 0.05$, $\text{FCcutoff} = 1.5$ | Used for highlighting significant findings in the Volcano Plot. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically checks for and installs the necessary Bioconductor and CRAN packages:
* `GEOquery` (For data download)
* `DESeq2` (For DGE analysis)
* `pheatmap` (For Heatmap visualization)
* `EnhancedVolcano` (For Volcano Plot visualization)
* `ggplot2` (For PCA visualization)

### ⚙️ Execution

1.  **Download** the `DGE Analysis of GSE153089 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization..R` file.
2.  **Optional:** Before running, you would typically define an output directory for saving plots (though this step is missing in the uploaded script, it's a best practice for a full pipeline).
3.  **Execute** the script in your R environment:
    ```R
    source("DGE Analysis of GSE153089 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization..R")
    ```

---

## 📊 Visualization and Output

The script generates the following plots, which are displayed in the R graphics device upon execution:

| Visualization | Analysis Stage | Description |
| :--- | :--- | :--- |
| **PCA Plot** | QC / Results | **Principal Component Analysis** plot demonstrating global clustering and separation of Cancer vs. Normal samples. |
| **Top 50 Variable Genes Heatmap** | QC | **Heatmap** of the 50 genes with the highest variance (VST-transformed data) to visualize sample grouping quality. |
| **Volcano Plot** | Results | **Volcano Plot** showing the $\log_2 \text{Fold Change}$ vs. $P_{\text{value}}$, highlighting significant and highly changed genes between the groups. |
