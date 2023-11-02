# Tests for fitLms functions

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

    # Fit tab output
    fit_tab <- fitLms(
        x = stats::model.matrix.default(as.formula('~a'), dm),
        y = ds
    )

    adj_m <- .getAdjMatr(gr, rownames(ds))
}

# Run tests --------------------------------------------------------------------

# Test getRealRhosFitDistr .....................................................

test_that("getRealRhosFitDistr ranking is correct with scaledEst = F", {

    rra_tab <-  getRealRhosFitDistr(
        bet_tab  = fit_tab,
        adj_m    = adj_m,
        back_par = ww,
        scaledEst = F
    )

    expect_equal(rra_tab[,'u_gene'], rank(fit_tab[,'Estimate'])/nrow(fit_tab))

})

test_that("getRealRhosFitDistr ranking is correct with scaledEst = T", {

    rra_tab <-  getRealRhosFitDistr(
        bet_tab  = fit_tab,
        adj_m    = adj_m,
        back_par = ww,
        scaledEst = T
    )

    expect_equal(rra_tab[,'u_gene'], rank(fit_tab[,'t value'])/nrow(fit_tab))

})

rra_tab0 <-  getRealRhosFitDistr(
    bet_tab  = fit_tab,
    adj_m    = adj_m,
    back_par = ww,
    scaledEst = F
)

test_that("getRealRhosFitDistr ranking works for u_neigh", {

    expect_equal(rra_tab0[,'u_neig'], rank(rra_tab0[,'rra_p'])/nrow(rra_tab0))

})

test_that("getRealRhosFitDistr rra_p is comprised between 1 and 0", {

    expect_true(all(rra_tab0[,'rra_p'] <= 1))

})

test_that("getRealRhosFitDistr rra_p is comprised between 1 and 0", {

    expect_true(all(rra_tab0[,'rra_p'] >= 0))

})

test_that("getRealRhosFitDistr score computed as formula", {

    expect_equal(rra_tab0[,'score'], log(rra_tab0[,"u_gene"]*rra_tab0[,"u_neig"]))

})

test_that("getRealRhosFitDistr gets error if adjMatrix misses any gene", {

    expect_error(
        getRealRhosFitDistr(
            bet_tab  = fit_tab,
            adj_m    = adj_m[-1,],
            back_par = ww,
            scaledEst = F
        )
    )

})

## Test getRandomRhos ..........................................................

test_that("getRandomRhos output columns corresponds to degree distribution", {

    back_dis <-  getRandomRhos(
        adj_m,
        n_iter = 10
    )

    expect_equal(as.integer(colnames(back_dis)),
                 sort(unique(colSums(as.matrix(adj_m)))))

})

test_that("getRandomRhos output corresponds to n.iterations", {

    back_dis <-  getRandomRhos(
        adj_m,
        n_iter = 10
    )

    expect_equal(nrow(back_dis), 10)

})

## Test getRandomPar ...........................................................

test_that("getRandomRhosPar raises error if cores is set to 0", {

    expect_error(getRandomRhosPar(
        adj_m,
        n_iter = 10,
        cores = 0
    ))

})

## Test getRealRhos ............................................................

back_dis <- getRandomRhos(adj_m, 10)

test_that("getRealRhos ranking is correct with scaledEst = F", {

    rra_tab <-  getRealRhos(
        bet_tab  = fit_tab,
        adj_m    = adj_m,
        back_dis = back_dis,
        scaledEst = F
    )

    expect_equal(rra_tab[,'u_gene'], rank(fit_tab[,'Estimate'])/nrow(fit_tab))

})

test_that("getRealRhos ranking is correct with scaledEst = T", {

    rra_tab <-  getRealRhos(
        bet_tab  = fit_tab,
        adj_m    = adj_m,
        back_dis = back_dis,
        scaledEst = T
    )

    expect_equal(rra_tab[,'u_gene'], rank(fit_tab[,'t value'])/nrow(fit_tab))

})








