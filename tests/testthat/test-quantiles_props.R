test_that("quantiles_props() returns a vector of the correct length", {
  expect_equal(length(quantiles_props(x = rnorm(1e2), quants = c(0.2, 0.4, 0.6, 0.8) ) ), length(c(0.2, 0.4, 0.6, 0.8) ) + 1)
})
