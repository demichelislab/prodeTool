# Tests for functions in utils file --------------------------------------------

# Prepare sample input for tests -----------------------------------------------

{
    set.seed(100)

    N_GENES = 10
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
        n1 =  sample(paste0("YY", 1:(N_GENES-1)), N_GENES*2, replace=T),
        n2 =  sample(paste0("YY", 1:N_GENES), N_GENES*2, replace=T)
    )
}

# run tests --------------------------------------------------------------------

# Tests for .getAdjMatr ........................................................

test_that('.getAdjMatr creates a matrix with all genes in ds', {

    gr <- data.frame(
        n1 =  sample(paste0("YY", 1:(N_GENES-2)), N_GENES*2, replace=T),
        n2 =  sample(paste0("YY", 1:(N_GENES-2)), N_GENES*2, replace=T)
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))

    expect_equal(nrow(adj_m), nrow(ds))

})

test_that('.getAdjMatr creates a matrix correct degree', {

    gr <- data.frame(
        n1 =  c(paste0("YY", 1:(N_GENES)), paste0("YY", (N_GENES):1)),
        n2 =  c(paste0("YY", (N_GENES):1), paste0("YY", 1:(N_GENES)))
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))

    expect_true(all(colSums(as.matrix(adj_m)) == 2))

})

test_that('.getAdjMatr creates a matrix correct degree', {

    gr <- data.frame(
        n1 =  c(paste0("YY", 1:(N_GENES))),
        n2 =  c(paste0("YY", 1:(N_GENES)))
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))

    expect_true(all(colSums(as.matrix(adj_m)) == 1))

})

test_that('.getAdjMatr creates a matrix correct degree', {

    gr <- data.frame(
        n1 =  c(paste0("YY", seq(1, N_GENES, 2))),
        n2 =  c(paste0("YY", seq(2, N_GENES, 2)))
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))

    expect_equal({
        tmp <- as.data.frame(t(apply(as.matrix(adj_m), 1, function(x){
          colnames(adj_m)[which(x == 1)]
        }))[seq(1, 10, 2),])
        colnames(tmp) <- c('n1', 'n2')
        rownames(tmp) <- NULL
        tmp
    }, gr)

})

test_that('.getAdjMatr creates a squared matrix', {

    gr <- data.frame(
        n1 =  c(paste0("YY", seq(1, N_GENES, 2))),
        n2 =  c(paste0("YY", seq(2, N_GENES, 2)))
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))

    expect_equal(rownames(adj_m), rownames(adj_m))

})

test_that('.getAdjMatr creates a matrix that is symmetric', {

    gr <- data.frame(
        n1 =  sample(paste0("YY", 1:(N_GENES-1)), N_GENES*2, replace=T),
        n2 =  sample(paste0("YY", 1:N_GENES), N_GENES*2, replace=T)
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))

    expect_true(Matrix::isSymmetric(adj_m))

})

# Tests for .checkAdjMatr ......................................................

test_that('.checkAdjMatr raises error if genes are not concordant', {

    adj_m <- .getAdjMatr(gr, rownames(ds))

    expect_error(.checkAdjMatr(colnames(adj_m), letters))

})

# Tests for .filterAdjMatrix ...................................................

test_that('.filterAdjMatrix output contains 1 degree genes', {

    gr <- data.frame(
        n1 =  c(paste0('YY', 1:N_GENES), 'YY1'),
        n2 =  c(paste0('YY', 1:N_GENES), 'YY2')
    )

    prInput <- suppressMessages(
        getProdeInput(score_matrix = ds, col_data = dm, edge_table = gr)
    )

    fit_tab <- fitLms(designMatrix(prInput), assay(prInput))

    filtMat <- suppressWarnings(.filterAdjMatrix(F, prInput, fit_tab = fit_tab))

    expect_equal(nrow(filtMat[['filtered']]), 8)

})

test_that('.filterAdjMatrix output contains high ctrl mean gns', {

    ds <- rbind(
        matrix(
            rnorm((N_GENES*N_SAMPLES)/2, 1, 0.1),
            nrow=N_GENES/2
        ),
        matrix(
            rnorm((N_GENES*N_SAMPLES)/2, -1, 0.1),
            nrow=N_GENES/2
        )
    )

    rownames(ds) <- paste0("YY", 1:N_GENES)
    colnames(ds) <- paste0("S", 1:ncol(ds))

    gr <- data.frame(
        n1 =  c(paste0("YY", seq(1, N_GENES, 2))),
        n2 =  c(paste0("YY", seq(2, N_GENES, 2)))
    )

    prInput <- suppressMessages(
        getProdeInput(score_matrix = ds, col_data = dm, edge_table = gr, as.formula('~a'))
    )

    fit_tab <- fitLms(designMatrix(prInput), assay(prInput), extendedNICEStats = T)

    filtMat <- suppressWarnings(
        .filterAdjMatrix(filterCtrl = T, prodeInput = prInput, fit_tab = fit_tab)
    )

    oo <- rowMeans(ds[,which(dm$a == 0)])
    rem <- names(oo)[which(oo > 0)]

    expect_true(all(rem %in% rownames(filtMat[['filtered']])))

})

test_that('.filterAdjMatrix removes 1 degree nodes', {

    gr <- data.frame(
        n1 =  c(paste0('YY', 1:N_GENES), 'YY1'),
        n2 =  c(paste0('YY', 1:N_GENES), 'YY2')
    )

    prInput <- suppressMessages(
        getProdeInput(score_matrix = ds, col_data = dm, edge_table = gr)
    )

    filtMat <- suppressWarnings(.filterAdjMatrix(F, prInput, fit_tab = ds))

    expect_equal(nrow(filtMat[['adjMatrix']]), 2)

})

test_that('.filterAdjMatrix raises warning if removing more than 50% of genes', {

    gr <- data.frame(
        n1 =  c(paste0('YY', 1:N_GENES), 'YY1'),
        n2 =  c(paste0('YY', 1:N_GENES), 'YY2')
    )

    prInput <- suppressMessages(
        getProdeInput(score_matrix = ds, col_data = dm, edge_table = gr)
    )

    expect_warning(.filterAdjMatrix(F, prInput, fit_tab = ds))

})

test_that('.filterAdjMatrix raises warning if removing more than 50% of genes', {

    gr <- data.frame(
        n1 =  c(paste0('YY', 1:N_GENES)),
        n2 =  c(paste0('YY', 1:N_GENES))
    )

    prInput <- suppressMessages(
        getProdeInput(score_matrix = ds, col_data = dm, edge_table = gr)
    )

    expect_error(.filterAdjMatrix(F, prInput, fit_tab = ds))

})

# Tests for .filter0DegAdj .....................................................

test_that('.filter0DegAdj  removes 0 degrees',{

    gr <- data.frame(
        n1 =  sample(paste0("YY", 1:(N_GENES-1)), N_GENES*2, replace=T),
        n2 =  sample(paste0("YY", 1:N_GENES), N_GENES*2, replace=T)
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))

    dds <- colSums(as.matrix(adj_m))

    expect_equal(nrow(.filter0DegAdj(adj_m, dds)), 0)

})


# Tests for .inputCheck ........................................................

test_that('.inputCheck raises error in case of wrong input (1)',{

    expect_error(
        .inputCheck(
            score_matrix = ds,
            col_data = dm[-1,],
            design = as.formula('~1'),
            edge_table = gr
        )
    )

})

test_that('.inputCheck raises error in case of wrong input (2)',{

    rownames(dm) <- NULL

    expect_error(
        .inputCheck(
            score_matrix = ds,
            col_data = dm,
            design = as.formula('~1'),
            edge_table = gr
        )
    )

})

test_that('.inputCheck raises error in case of wrong input (3)',{

    expect_error(
        .inputCheck(
            score_matrix = ds,
            col_data = dm,
            design = as.formula('~aa'),
            edge_table = gr
        )
    )

})

test_that('.inputCheck raises error in case of wrong input (4)',{

    expect_error(
        .inputCheck(
            score_matrix = as.data.frame(ds),
            col_data = dm,
            design = as.formula('~1'),
            edge_table = gr
        )
    )

})

test_that('.inputCheck raises error in case of wrong input (5)',{

    ds[1] <- Inf

    expect_error(
        .inputCheck(
            score_matrix = ds,
            col_data = dm,
            design = as.formula('~1'),
            edge_table = gr
        )
    )

})


test_that('.inputCheck raises error in case of wrong input (5)',{

    mode(ds) <- "character"

    expect_error(
        .inputCheck(
            score_matrix = ds,
            col_data = dm,
            design = as.formula('~1'),
            edge_table = gr
        )
    )

})




