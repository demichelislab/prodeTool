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

    design <- stats::model.matrix.default(
        as.formula('~a'),
        as.data.frame(dm)
    )

    obj <- newProdeInput(ds, dm, design, adj_m, 'NIE', weights = NULL)

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

# Tests for runProde -----------------------------------------------------------

test_that('runProde works (1).', {

    oo <- suppressMessages(runProde(obj, scaledEst = F, extendedNICEStats = F))
    expect_true(!'ctrl_mean' %in% colnames(results(oo)))

})


test_that('runProde extededNICE stats works correctly.', {

    oo <- suppressMessages(runProde(obj, scaledEst = F, extendedNICEStats = T))
    expect_true('ctrl_mean' %in% colnames(results(oo)))

})

test_that('runProde output dimension is as expected.', {

    oo <- suppressMessages(runProde(obj, scaledEst = F, extendedNICEStats = T))
    expect_equal(nrow(results(oo)), 10)

})


test_that('runProde raises error When recieving NA as input.', {

    expect_error(suppressMessages(runProde(obj, scaledEst = NA, extendedNICEStats = T)))

})

test_that('runProde scaledEst works as expected compared to ranks.', {

    oo <- suppressMessages(runProde(obj, scaledEst = T, extendedNICEStats = T))
    expect_equal(results(oo)[,'u_gene'], rank(results(oo)[['t.value']])/nrow(results(oo)))

})

test_that('runProde scaledEst works as expected compared to ranks.', {

    oo <- suppressMessages(runProde(obj, scaledEst = F, extendedNICEStats = T))
    expect_equal(results(oo)[,'u_gene'], rank(results(oo)[['Estimate']])/nrow(results(oo)))

})

test_that('runProde filters out compared to ctrl_mean.', {

    oo <- suppressWarnings(suppressMessages(
        runProde(obj, scaledEst = F, extendedNICEStats = T, filterCtrl = T)))
    expect_true(all(results(oo)[['ctrl_mean']] < 0))

})


test_that('runProde filters out 1 degree', {

    ## This is an example of an edge_table
    gr <- data.frame(
        n1 =  sample(paste0("YY", 1:(N_GENES-1)), N_GENES*2, replace=T),
        n2 =  sample(paste0("YY", 1:N_GENES), N_GENES*2, replace=T)
    )

    excl <- gr[,'n1'][1]
    gr <- gr[-which(gr[,'n1'] == gr[,'n1'][1] | gr[,'n2'] == gr[,'n1'][1]),]

    obj <- suppressMessages(getProdeInput(ds, dm, gr))

    oo <- suppressWarnings(suppressMessages(
        runProde(obj, scaledEst = F, extendedNICEStats = F, filterCtrl = F)))

    expect_true(excl %in% rownames(filteredGenes(oo)))

})















