# Tests for allClasses.R -------------------------------------------------------

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

    adj_m <- .getAdjMatr(gr, rownames(ds))
    wws <- 1:nrow(adj_m)

    design <- stats::model.matrix.default(
        as.formula('~a'),
        as.data.frame(dm)
    )
    
    obj <- newProdeInput(ds, dm, design, adj_m, 'NIE', wws, gr)

    fit_tab <- fitLms(
        x = stats::model.matrix.default(as.formula('~a'), dm),
        y = ds
    )

    rra_tab <-  getRealRhosFitDistr(
        bet_tab  = fit_tab,
        adj_m    = adj_m,
        back_par = ww,
        scaledEst = F
    )

    obj2 <- newProdeResults(fit_tab, rra_tab, adj_m, 'NIE',
                            filtered = DataFrame())

}

# Run tests --------------------------------------------------------------------

# Test newProdeInput ...........................................................

test_that("newProdeInput deals with wrong input", {
    expect_error(newProdeInput(ds, dm, design, adj_m, NA, wws))
})

test_that("newProdeInput assay returns correct data", {
    expect_equal(assay(obj), ds)
})

test_that("newProdeInput colData returns correct data", {
    expect_equal(colData(obj), DataFrame(dm))
})

test_that("newProdeInput adjMatrix returns correct data", {
    expect_equal(adjMatrix(obj), adj_m)
})

test_that("newProdeInput modality returns correct data", {
    expect_equal(modality(obj), 'NIE')
})

test_that("newProdeInput subsetAdjMatrix creates correct output (1)", {
    expect_equal(nrow(subsetAdjMatrix(obj, rownames(assay(obj))[1:2])), 2)
})

test_that("newProdeInput subsetAdjMatrix creates correct output (2)", {
    expect_equal(rownames(subsetAdjMatrix(obj, rownames(assay(obj))[1:2])),
                 rownames(assay(obj))[1:2])
})

test_that("newProdeInput class is correct", {
    expect_equal(class(obj)[1], 'prodeInput')
})

# Test newProdeResults .........................................................

test_that("newProdeResults deals with wrong input", {
    expect_equal(class(obj2)[1], 'prodeResults')
})

test_that("newProdeInput assay returns correct data (1)", {
    expect_equal(filteredGenes(obj2), DataFrame())
})

test_that("newProdeInput assay returns correct data (2)", {
    expect_equal(results(obj2), {
        oo <- DataFrame('gene' = rownames(fit_tab), fit_tab, rra_tab)
        colnames(oo)[ncol(oo)] <- modality(obj)
        oo
    })
})

test_that("newProdeInput assay returns correct data (3)", {
    expect_equal(adjMatrixResults(obj2), adj_m)
})

test_that("newProdeInput assay returns correct data (4)", {

    obj2 <- newProdeResults(fit_tab, rra_tab, adj_m, 'NIE',
                            filtered = DataFrame())

    expect_length(rownames(adjMatrixResults(obj2)), nrow(fit_tab))

})


test_that("newProdeInput assay returns correct data (5)", {

    obj2 <- newProdeResults(fit_tab, rra_tab, adj_m, 'NIE',
                            filtered = DataFrame())

    expect_equal(rownames(adjMatrixResults(obj2)), rownames(fit_tab))

})


test_that("newProdeInput assay returns correct data (5)", {

    obj2 <- newProdeResults(fit_tab[1:3,], rra_tab[1:3,], adj_m, 'NIE',
                            filtered = DataFrame())

    expect_gt(nrow(adjMatrixResults(obj2)), nrow(results(obj2)))

})


test_that("newProdeInput assay returns correct data (6)", {

    expect_error(results(obj))

})


















