test_that("results of ptweibull match base", {
  expect_equal(ptweibull(q = 1.5, shape = 3.2, scale = 1.75), pweibull(q = 1.5, shape = 3.2, scale = 1.75))
  expect_equal(ptweibull(q = 1.5, shape = 3.2, scale = 1.75, log.p = TRUE), pweibull(q = 1.5, shape = 3.2, scale = 1.75, log.p = TRUE))
  expect_equal(ptweibull(q = 1.5, shape = 3.2, scale = 1.75, lower.tail = FALSE), pweibull(q = 1.5, shape = 3.2, scale = 1.75, lower.tail = FALSE))
  expect_equal(ptweibull(q = 2.3, shape = 1.5, scale = 5.2, lower.tail = FALSE, log.p = TRUE), pweibull(q = 2.3, shape = 1.5, scale = 5.2, lower.tail = FALSE, log.p = TRUE))
})

test_that("results of ptweibull match base with truncation", {
  pbase = function(q, shape, scale = 1, a = 0, b = Inf, lower.tail = TRUE){
    pw = local({shape = shape; scale = scale; function(q) pweibull(q, shape = shape, scale = scale)})
    if(isTRUE(lower.tail)){
      (pw(q) - pw(a)) / (pw(b) - pw(a))
    } else {
      (pw(b) - pw(q)) / (pw(b) - pw(a))
    }
  }

  expect_equal(
    ptweibull(q = 7.8, shape = 2.5, scale = 5.6, a = 2),
    pbase(q = 7.8, shape = 2.5, scale = 5.6, a = 2)
  )
  expect_equal(
    ptweibull(q = 5.8, shape = 2, scale = 2.5, b = 10),
    pbase(q = 5.8, shape = 2, scale = 2.5, b = 10)
  )
  expect_equal(
    ptweibull(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10),
    pbase(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10)
  )

  expect_equal(
    ptweibull(q = 7.8, shape = 2.5, scale = 5.6, a = 2, lower.tail = FALSE),
    pbase(q = 7.8, shape = 2.5, scale = 5.6, a = 2, lower.tail = FALSE)
  )
  expect_equal(
    ptweibull(q = 5.8, shape = 2, scale = 2.5, b = 10, lower.tail = FALSE),
    pbase(q = 5.8, shape = 2, scale = 2.5, b = 10, lower.tail = FALSE)
  )
  expect_equal(
    ptweibull(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10, lower.tail = FALSE),
    pbase(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10, lower.tail = FALSE)
  )

  expect_equal(
    ptweibull(q = 7.8, shape = 2.5, scale = 5.6, a = 2, log.p = TRUE),
    log(pbase(q = 7.8, shape = 2.5, scale = 5.6, a = 2))
  )
  expect_equal(
    ptweibull(q = 5.8, shape = 2, scale = 2.5, b = 10, log.p = TRUE),
    log(pbase(q = 5.8, shape = 2, scale = 2.5, b = 10))
  )
  expect_equal(
    ptweibull(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10, log.p=TRUE),
    log(pbase(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10))
  )

  expect_equal(
    ptweibull(q = 7.8, shape = 2.5, scale = 5.6, a = 2, lower.tail = FALSE, log.p = TRUE),
    log(pbase(q = 7.8, shape = 2.5, scale = 5.6, a = 2, lower.tail = FALSE))
  )
  expect_equal(
    ptweibull(q = 5.8, shape = 2, scale = 2.5, b = 10, log.p = TRUE, lower.tail = FALSE),
    log(pbase(q = 5.8, shape = 2, scale = 2.5, b = 10, lower.tail = FALSE))
  )
  expect_equal(
    ptweibull(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10, log.p=TRUE, lower.tail = FALSE),
    log(pbase(q = 2.2, shape = 5.2, scale = 1.8, a = 2, b = 10, lower.tail = FALSE))
  )
})
