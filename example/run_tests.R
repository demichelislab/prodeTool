# This is a quick example of a possible run when performing essentiality or
# context-essentiality analysis
set.seed(100)

library(prodeTool)

# TOY EXAMPLE ==================================================================

N_GENES = 100
N_SAMPLES = 10

## This is an example of an edge table
gr <- data.frame(
    n1 =  sample(paste0("YY", 1:(N_GENES-1)), N_GENES*100, replace=T),
    n2 =  sample(paste0("YY", 1:N_GENES), N_GENES*100, replace=T)
)

## This is an example of a gene effects matrix
ds <-
    matrix(
        rnorm(N_GENES*N_SAMPLES),
        nrow=N_GENES
    )

rownames(ds) <- paste0("YY", 1:N_GENES)
colnames(ds) <- paste0("S", 1:ncol(ds))

## This is an example of a column-data table
dm <- data.frame(
    a = rep(c(1, 0), each=N_SAMPLES/2),
    b = sample(c("a", "b"), N_SAMPLES, replace=T)
)
rownames(dm) <- paste0("S", 1:ncol(ds))

# FOR NIE SCORES (ESSENTIALITY ANALYSIS) .......................................

prodeInputNIE <- getProdeInput(
    score_matrix   = ds,
    col_data       = dm,
    edge_table     = gr,
)

outputNIE <- runProde(
    prodeInput        = prodeInputNIE,
    scaledEst         = F # When computing NIE scores, scaling is suggested to be set as F
)

outputNIE

# FOR NICE SCORES (CONTEXT ESSENTIALITY ANALYSIS) ..............................
# with pre-computed back-ground distribution (suggested option)

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

# # FOR NICE SCORES (CONTEXT ESSENTIALITY ANALYSIS)
# # computing back-ground distribution (takes more time)
#
# outputNICE.v2 <- runProde(
#     prodeInput        = prodeInputNICE,
#     filterCtrl        = T,
#     computeBack       = T,
#     extendedNICEStats = T,
#     n_iter            = 10^5,
#     cores             = 2
# )
#
# cor(outputNICE@results$NICE_score, outputNICE.v2@results$NICE_score)
