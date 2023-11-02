# Tests for fitLms functions

# Prepare sample input for tests -----------------------------------------------

y <- matrix(c(1:50, 50:1), 10, 10)

col_data <- data.frame(
    'gr' = rep.int(c(1,0), c(5,5)),
    'z'  = 1:10,
    'k'  = letters[rep(1:2, each=5)]
)

x1 <- stats::model.matrix.default(
    stats::as.formula("~gr"),
    as.data.frame(col_data)
)

x2 <- stats::model.matrix.default(
    stats::as.formula("~1"),
    as.data.frame(col_data)
)

# Run tests --------------------------------------------------------------------

test_that("fitLms produces the estimate for each row (check first fit)", {

    expect_equal(
        fitLms(x1, y)[1,],
        summary(
            lm(y[1,]~col_data[['gr']])
        )[[4]][2,]
    )

})

test_that("fitLms produces the estimate for each row (check last fit)", {

    expect_equal(
        fitLms(x1, y)[nrow(y),],
        summary(
            lm(y[nrow(y),]~col_data[['gr']])
        )[[4]][2,]
    )

})


test_that("fitLms computes average when formula ~1.", {

    expect_equal(
        fitLms(x2, y)[,1],
        base::rowMeans(y)
    )

})


test_that("Directionality of estimate is coherent with groups average", {

        expect_true({
            oo <-  fitLms(x1, y, extendedNICEStats=T)
            all(
                sign(oo[,'Estimate']) ==
                    sign(oo[,'case_mean'] - oo[,'ctrl_mean'])
            )
        })

})

test_that("Directionality of estimate is coherent with groups average", {

    expect_warning(
        fitLms(x2, y, extendedNICEStats=T)
    )

})


test_that("Input matrix contains NAs raises error", {

    y[seq(1, 100, 25)] <- NA

    expect_error(
        fitLms(x1, y, extendedNICEStats=T)
    )

})

test_that("Input group variable contains NAs raises error", {

    x1[,2][1] <- NA

    expect_warning(
        fitLms(x1, y, extendedNICEStats=T)
    )

})







