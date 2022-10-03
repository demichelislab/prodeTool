test_that("splitting works works", {
  expect_equal(
      strsplit1("a,2,3", ","),  c("a", "2", "3")
  )
})

