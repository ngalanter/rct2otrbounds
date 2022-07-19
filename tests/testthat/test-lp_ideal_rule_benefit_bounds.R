test_that("transform range works for positive ranges with higher", {

  test <- transform_range(y_range = 1:10,mean_worse = 4,mean_better = 6,scale = "higher")

  expect_setequal(0:9, test$new_y_range)

  expect_equal(3, test$new_mean_worse-min(test$new_y_range))

  expect_equal(4, max(test$new_y_range)-test$new_mean_better)

  expect_equal(2, test$new_mean_better - test$new_mean_worse)

})

test_that("transform range works for non-integer values of means with higher", {

  test <- transform_range(y_range = 1:10,mean_worse = 4.4,mean_better = 6.1,scale = "higher")

  expect_setequal(0:9, test$new_y_range)

  expect_equal(3.4, test$new_mean_worse-min(test$new_y_range))

  expect_equal(3.9, max(test$new_y_range)-test$new_mean_better)

  expect_equal(1.7, test$new_mean_better - test$new_mean_worse)

})

test_that("transform range works for negative values with higher", {

  test <- transform_range(y_range = -1*(1:10),mean_worse = -7,mean_better = -5,scale = "higher")

  expect_setequal(0:9, test$new_y_range)

  expect_equal(3, test$new_mean_worse-min(test$new_y_range))

  expect_equal(4, max(test$new_y_range)-test$new_mean_better)

  expect_equal(2, test$new_mean_better - test$new_mean_worse)
})

test_that("transform range works for positive ranges with lower", {

  test <- transform_range(y_range = 1:10,mean_better = 5,mean_worse = 7,scale = "lower")

  expect_setequal(0:9, test$new_y_range)

  expect_equal(3, test$new_mean_worse-min(test$new_y_range))

  expect_equal(4, max(test$new_y_range)-test$new_mean_better)

  expect_equal(2, test$new_mean_better - test$new_mean_worse)

})

test_that("transform range works for negative values with lower", {

  test <- transform_range(y_range = -1*(1:10),mean_better = -6,mean_worse = -4,scale = "lower")

  expect_setequal(0:9, test$new_y_range)

  expect_equal(3, test$new_mean_worse-min(test$new_y_range))

  expect_equal(4, max(test$new_y_range)-test$new_mean_better)

  expect_equal(2, test$new_mean_better - test$new_mean_worse)
})
