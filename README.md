# MGP API

This repository contains a Plumber API for performing Morph-Genetic-Profiling (MGP). MGP is a multivariate method used to analyze the relationship between high-dimensional genetic data (such as gene sets from biological pathways) and high-dimensional shape data (such as 3D landmark coordinates).

This API implements the `ddsPLS` (discriminant and de-sparsified PLS) method to find associations between genotype and phenotype. It provides endpoints for running analyses based on Gene Ontology (GO) terms or custom gene lists, for both Diversity Outbred (DO) mouse skull data and human facial data.

## 1. Setup & Dependencies

### R Packages

This API requires numerous R packages. You can install them from CRAN and Bioconductor.

**CRAN Packages:**
```R
install.packages(c(
  "plumber", "Morpho", "ddsPLS", "Jovid", "dplyr", "dbplyr", "future", 
  "promises", "abind", "devtools", "rgl", "ggplot2", "cmocean", "pals",
  "RColorBrewer", "plot3D", "Rvcg", "geomorph", "shapes", "foreach",
  "parallel", "doParallel", "biomaRt", "biomartr", "plinkFile", "LDlinkR",
  "paran", "gmodels", "DBI", "RSQLite"
))
```

**Bioconductor Packages:**
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GenomicFeatures", "org.Mm.eg.db", "org.Hs.eg.db", "GO.db", "GOSemSim"
))
```

### External Dependencies
*   For human data analysis, the **PLINK** command-line tool is required. The API script has a hardcoded path to it, which you may need to adjust.

## 2. Data Requirements

The API relies on several pre-processed data objects being loaded into the R environment. The script contains hardcoded paths to these files, which you will need to update to match your local file structure.

### For Mouse DO Skulls

1.  **Phenotype Data (`Y`, `A_lm_raw`, etc.)**
    *   **Format:** An R matrix loaded from an `.Rdata` file.
    *   **Dimensions:** `n_individuals` x `(n_landmarks * 3)`.
    *   **Content:** Rows represent individual mice, and columns contain the flattened 3D coordinates of landmarks (e.g., x1, y1, z1, x2, y2, z2, ...).

2.  **Genotype Data (`DO_probs_DB`, etc.)**
    *   **Format:** An SQLite database file (`.sqlite`).
    *   **Content:** The database should contain one table per genetic marker (e.g., SNP). Each table is named with the marker ID.
    *   **Table Schema:** Each marker table should have `n_individuals` rows and 8 columns, with each column containing the haplotype probability for one of the 8 founder strains (A/J, C57BL/6J, etc.).

3.  **Marker Information (`combined.markers`, etc.)**
    *   **Format:** An R data frame loaded from an `.Rdata` file.
    *   **Content:** Must contain information for each marker, including:
        *   A column with the marker ID (matching table names in the genotype SQLite database).
        *   `chr`: The chromosome.
        *   `Mbp_mm10`: The marker's position in megabase pairs.

4.  **Gene Annotations**
    *   `mmusculusEnsembl`: An `AnnotationDbi` SQLite database (`ensemble.sqlite`) containing mouse gene transcript information (chromosome, start, end).
    *   `DO.go`: An R data frame mapping GO terms to GO IDs.
    *   `bad_genes_mm.csv`, `chrom_two_genes_mm.csv`: CSV files listing genes to be excluded from analysis.

### For Human Faces

1.  **Phenotype Data (`YHuman`, `YHumanTANZ`, etc.)**
    *   **Format:** An R matrix loaded from an `.Rdata` file.
    *   **Dimensions:** `n_individuals` x `n_features`.
    *   **Content:**
        *   For **sparse landmarks** (`YHuman`, `YHumanTANZ`): `n_individuals` x `(n_landmarks * 3)`.
        *   For **dense landmarks** (`YHumanDense`, `YHumanTANZDense`): A matrix of principal components (`n_individuals` x `n_PCs`) representing the shape data. The reversal of PCA also requires `pca.eigvec` (eigenvectors) and `pheno.avg` (mean shape) to be loaded.
    *   An accompanying `pheno.id.human` vector of individual IDs is also required.

2.  **Genotype Data**
    *   **Format:** PLINK binary files (`.bed`, `.bim`, `.fam`). The API expects these files to be in a specific directory structure.

3.  **Gene Annotations**
    *   The `org.Hs.eg.db` Bioconductor package is used for human gene annotations.
    *   `bad_genes_human.csv`: A CSV file listing human genes to be excluded.

## 3. API Usage Guide

The API is served using Plumber. Once running, you can send requests to the following endpoints.

---

### `GET /mgp`

Runs the primary MGP analysis using a Gene Ontology (GO) term.

*   **Parameters:**
    *   `GO.term` (string): The GO term to analyze (e.g., `"chondrocyte differentiation"`).
    *   `pheno` (string): The phenotype dataset to use (e.g., `"Y"`, `"A_lm_gen"`, `"YHuman"`).
    *   `pheno_index` (string): A string representing the range of landmark columns to use (e.g., `"1:54"`).
    *   `lambda` (numeric): Regularization strength for ddsPLS. Defaults to `0.06`.
    *   `pls_axis` (integer): The number of PLS components to compute. Defaults to `1`.
    *   `remove_pc1` (boolean): If `TRUE`, removes the first principal component from the phenotype data before analysis.
    *   `use_standardized_PCA` (boolean): If `TRUE`, standardizes the phenotype data using PCA.
    *   `doPermutation` (boolean): If `TRUE`, runs a permutation test to calculate a p-value against shuffled subjects.
    *   `permutation` (integer): Number of permutations to run (max 200).
    *   `doRangenes` (boolean): If `TRUE`, runs a coherence test against random gene sets.
    *   `rangenes` (integer): Number of random gene sets to test against (max 200).

*   **Example:**
    ```bash
    curl "http://127.0.0.1:6234/mgp?GO.term=chondrocyte%20differentiation&pheno=Y&pls_axis=2"
    ```

---

### `GET /custom_mgp`

Runs an MGP analysis using a custom list of genes.

*   **Parameters:**
    *   `genelist` (string): A comma-separated list of gene symbols (e.g., `"Bmp7,Bmp2,Bmp4"`).
    *   Accepts the same `pheno`, `pheno_index`, `lambda`, `pls_axis`, and optional analysis parameters as `/mgp`.

*   **Example:**
    ```bash
    curl "http://127.0.0.1:6234/custom_mgp?genelist=Bmp7,Bmp2,Bmp4&pheno=Y"
    ```

---

### `POST /autoLambda`

Automatically explores a range of lambda values to find the optimal one based on the Bayesian Information Criterion (BIC).

*   **Body (form-data):**
    *   Takes the same parameters as `/mgp` or `/custom_mgp` (e.g., `GO.term` or `genelist`, `pheno`, etc.).
*   **Returns:** The best lambda value found.

*   **Example:**
    ```bash
    curl -X POST -F "GO.term=chondrocyte differentiation" -F "pheno=Y" http://127.0.0.1:6234/autoLambda
    ```

---

### `GET /mgi`

Queries gene information (location, chromosome) for a given biological process from the database.

*   **Parameters:**
    *   `process` (string): The name of the GO term process.

*   **Example:**
    ```bash
    curl "http://127.0.0.1:6234/mgi?process=chondrocyte%20differentiation"
    ```

---

### `GET /mutant_list`

Retrieves a list of available mouse mutants for comparison analyses.

*   **Parameters:** None.

*   **Example:**
    ```bash
    curl http://127.0.0.1:6234/mutant_list
    ```

---

### `GET /mutant_comparison`

Compares the shape effect of an MGP analysis with the shape change observed in a known mouse mutant.

*   **Parameters:**
    *   `MGP_pheno_loadings` (string): A comma-separated string of numeric phenotype loadings from a previous `/mgp` or `/custom_mgp` result.
    *   `mutant` (string): The name of the mutant to compare against (from `/mutant_list`).

*   **Example:**
    ```bash
    curl "http://127.0.0.1:6234/mutant_comparison?MGP_pheno_loadings=0.1,0.2,...&mutant=Bmp2"
    ```
---

### Other Endpoints

*   `GET /GO_network`: Returns the gene-GO term network for a given process or gene list.
*   `POST /mgp_cor`: Calculates a correlation matrix between the effects of several randomly chosen biological processes from a cached set of results.
