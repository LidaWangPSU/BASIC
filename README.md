# BASIC
Bulk And Single cell eQTL Integration across Cell states (BASIC)


## Table of contents
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Quick tutorial](#Quick_tutorial)
* [Input files](#Input_files)
* [Contact](#Contact)

## Introduction
To enhance the efficacy of sc-eQTL studies, we propose a novel method that integrates both bulk and single-cell RNA sequencing data to improve the detection of cell type-specific eQTLs. 

It is developed and maintained by Lida Wang at [Dajiang Liu's Group](https://dajiangliu.blog).

## Installation
BASIC is hosted on GitHub, making it easy to install and update. Before installing. The following R packages are required: MASS, dplyr, data.table, tidyr, R.utils, pracma, quadprog and ATAC.

Steps to Install:
```r
install.packages("devtools")
library(devtools)
```

Install ACAT and BASIC from the GitHub repository:

```
devtools::install_github("yaowuliu/ACAT")
library(ACAT)
devtools::install_github("LidaWangPSU/BASIC/BASIC")
library(BASIC)
```


## Quick tutorial
### 1. Prepare bulk and single cell eQTLs summary statistics (effect size and s.e.). 

#### Effect Size Matrix
This matrix contains eQTL effect sizes for both bulk and single-cell data.

- **Columns**:
  - **Column 1**: Gene-SNP pair identifiers.
  - **Column 2**: Bulk effect size.
  - **Columns 3+**: Cell type-specific eQTL effect sizes.
 
**Example**:

|    Gene-snp pair    |      Bulk     | Cell type 1  |  Cell type 2 | ...... |  Cell type k |
| ------------------- |      ----     | -----------  |  ----------- | ------ |  ----------- |
| ENSG00000XXXXX-snp1 |      0.65     |      0.55    |      0.30    | ...... |      0.75    |  
| ENSG00000XXXXX-snp2 |     -0.50     |     -0.62    |     -0.35    | ...... |      0.02    |


#### Standard Error Matrix (S.E.)
This matrix has the same dimensions as the effect size matrix and represents the standard errors associated with the eQTL effect sizes.

- **Columns**:
  - **Column 1**: Gene-SNP pair identifiers (should match the effect size matrix).
  - **Column 2**: Bulk effect size standard error.
  - **Columns 3+**: Standard errors for the cell type-specific eQTL effect sizes.

**Example**:

|    Gene-snp pair    |      Bulk     | Cell type 1  |  Cell type 2 | ...... |  Cell type k |
| ------------------- |      ----     | -----------  |  ----------- | ------ |  ----------- |
| ENSG00000XXXXX-snp1 |    0.0071     |    0.0316    |    0.0316    | ...... |    0.0316    |  
| ENSG00000XXXXX-snp2 |    0.0071     |    0.0316    |    0.0316    | ...... |    0.0316    |

- It is acceptable to have missing data in the input. We will analyze all gene-SNP pairs unless sc-eQTL data is entirely missing for all cell types. In other cases, we will simply ignore the missing values and keep "NA" in the output.
   
We have a small example input data [here](https://github.com/LidaWangPSU/BASIC/tree/main/example_data/). Data were subsetted from brain bulk and single cell eQTLs as an example to run the script.

```r
library(data.table)
beta <- as.data.frame(fread("~/example_beta_chr22.txt.gz"))
se <- as.data.frame(fread("~/example_beta_se_chr22.txt.gz"))
```
  
### 2. Run BASIC
#### BASIC Step 1: cell type weights
Here, we used Non-negative least squares to estimate cell type weights.
You can also specify the weights by yourself, e.g. estimated from scRNAseq data or other methods.
```r
weight <- nnls.weights(beta,se)
```


#### BASIC Step 2: project sc-eQTLs onto PCs 
- **nPC**: Number of PCs to use. Maximum is the number of cell types - 1 (K-1)
```r
K <- 8 #8 cell types in example 
meta <- meta_regression_fast(beta,se,nPC=(K-1))
MDS_ALL <- meta[[1]]
gamma_all <- meta[[2]]
```
The step 2 output includes:
* PC loadings: the loadings of PC 0 to PC (K-1)
* Projected eQTL by meta-regression: effect sizes and s.e. for K PCs

#### BASIC Step 3: Integrate Bulk eQTL and joint analyze
```r
basic_res <- basic_axisQTL(beta,se,gamma_all,MDS_ALL,weight,nPC=(K-1))
```

### 3. Output results
The output includes:
* K PCs axis-QTL p value and Cauchy combined p value
* K PCs axis-QTL effect size and s.e.

Example of EXPRESSO output can be found [here](https://github.com/LidaWangPSU/BASIC/tree/main/basic_output).

## Contact
Lida Wang [lida.wang.96@gmail.com](lida.wang.96@gmail.com)
