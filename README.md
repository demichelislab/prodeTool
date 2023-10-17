# PRODE: Neighborhood Informed Essential and Context Essential scores

PRODE is a tool that integrates Gene Effects data and Protein-Protein Interactions to
compute, for each gene, Neighborhood-Informed Essential (NIE) or Neighborhood-Informd Context Essential (NICE) scores. 

## Installation

PRODE is currently available as the R package prodeTool. 
It can be installed from github with: 

```R
devtools::install_github('cantorethomas/prodeTool')
```

## Usage 
Depending on the user-defined analysis setting, PRODE computes NIE or NICE scores. In both cases, prodeTool runs PRODE in two main steps. First, a `prodeInput` object is computed leveraging the `getProdeInput()` function. Second, PRODE workflow is ran through the `runProde()` function. Here follows a more detailed example when running PRODE to compute NIE or NICE scores. 

### Input data example 
Here follows a description of required input information from the `getProdeInput()` function documentation. 

* `design` this is the design formula required for NICE score computation. In the case of
   NIE score computation, it is ignored (default is set to `NULL`). `design` has to
   be a `formula` object. Last variable in formula will be considered as the group
   variable for the linear model fit, being 1 for case group and 0 for control group.

* `score_matrix` matrix of input scores with genes on rows
   and samples on columns. No `NA` values are allowed. Gene names are required
   as `rownames`. Sample names are required as `colnames` and need to be
   identical to `col_data` column names.

* `col_data` data.frame of column data, as described by `SummarizedExperiment`. Row names
   need to be identical to score-matrix column names and correspond to samples.
   Every sample-level included in design formula (if provided) should be included in `col_data`.
   The group variable, encoding the sample groups that need to be compared for NICE score
   computation, has to be a binary integer vector, encoding 1 for the case-group and 0 for the
   control-group.

* `edge_table` matrix of interactions between gene pairs. This matrix is composed by interactions on rows
   and two columns (one for each gene in the pair). Genes not present in `score_matrix` row names
   will be discarded when constructing `prodeInput` object.

```R
library(prodeTool)

set.seed(100)

N_GENES = 100
N_SAMPLES = 10

## This is an example of a score_matrix 
ds <-
    matrix(
        rnorm(N_GENES*N_SAMPLES),
        nrow=N_GENES
    )

rownames(ds) <- paste0("YY", 1:N_GENES)
colnames(ds) <- paste0("S", 1:ncol(ds))

## This is an example of a column_data table
dm <- data.frame(
    a = rep(c(1, 0), each=N_SAMPLES/2),
    b = sample(c("a", "b"), N_SAMPLES, replace=T)
)
rownames(dm) <- paste0("S", 1:ncol(ds))

## This is an example of an edge_table
gr <- data.frame(
    n1 =  sample(paste0("YY", 1:(N_GENES-1)), N_GENES*100, replace=T),
    n2 =  sample(paste0("YY", 1:N_GENES), N_GENES*100, replace=T)
)

```

### NIE scores  

`getProdeInput()` function allows to construct a `prodeInput` object. NIE or NICE analyses will be performed depending on design formula, as described by the function documentation. `runProde()` function 
runs PRODE workflow. During NIE scores computation, no groups are compared, hence, design formula is not necessary. 

#### Minimal run 
```R
prodeInputNIE <- getProdeInput(
    score_matrix   = ds,
    col_data       = dm,
    edge_table     = gr,
)

outputNIE <- runProde(
    prodeInput        = prodeInputNIE,
    scaledEst         = F # When computing NIE scores, scaling is suggested to be set as F
)

```
#### Output description 
Here follows a description of major output columns. Missing ones are better descibed in `runProde()` documentation. 

##### Gene-level results
* `Estimate` in case of NIE scores, it corresponds to the intercept of model fit,
        i.e., for each gene, the average values across samples. In case of NICE scores, it's
        the coefficient of to the variable encoding the condition of
        interest, as represented in the `prodeInput` object.,


##### Neighborhood-level results

* `rra_score` is the $\rho$ value computed by RRA algorithm for each gene.
* `rra_p` is the p-value corresponding to each $\rho$ value (computed according to neighborhood size).

##### Final Score results

* `u_gene` is the percentile of gene-level signal (`Estimate` or `t.value` columns, depeindin if `scaledEst=T`)
* `u_neigh` is the percentile of neighborhood-level signal (`rra_p`)
* `NIE_score`computed as $log(u_{gene} \times u_{neigh})$

### NICE scores  
NICE scores are obtained from the comparison of two sample groups, hence, `design` formula, need to include the grouping variable present in `col_data` object, on the left side as the last variable. In the design, other variables can be added, as confounding factors in the analysis. In the following example, in the $~b+a$, $a$ is a binary vector encoding the group variable, while $b$ is a covariate. lower NICE scores correspond to genes which display neighborhood-informed context essential  signal in samples part of the case group, compared to the control group. 

```R
prodeInputNICE <- getProdeInput(
    score_matrix   = ds,
    col_data       = dm,
    edge_table     = gr,
    design         = as.formula("~b+a")
)

outputNICE <- runProde(
    prodeInput        = prodeInputNICE,
    filterCtrl        = T,
    extendedNICEStats = T,
    scaledEst         = T
)

outputNICE
```
#### Output description  
When running PRODE with `extendedNICEStats=T`, additional columns are computed: 

*`ctrl_mean` is the average value, for each gene, of control samples.
*`case_mean` is the average value, for each gene, of case samples.
*`ctrl_sd` is the standard deviation of each gene values in control samples.
*`case_sd` is the standard deviation of each gene values in case samples.
*`ctrl_n` is the number of samples in the control group.
*`case_n` is the number of samples in the case group.
