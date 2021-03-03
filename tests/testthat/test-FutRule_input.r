# Test the input of the FutRule function

# Normal case
normalCase = list(
  # Endpoint type
  endpoint_type = "Normal",

  # Number of enrolled patients (control, two treatments)
  sample_size = c(100, 100, 100),

  # Patient dropout rate
  dropout_rate = 0.05,

  # Mean and SD in the control arm 
  control_mean = 0,
  control_sd = 1,

  # Means and SDs in the treatment arms 
  treatment_mean = c(0.25, 0.30),
  treatment_sd = c(1, 1),

  # Information fraction
  info_frac = 0.5,

  # One-sided alpha level
  alpha = 0.025,

  # Number of simulations
  nsims = 1000
)

test_that("Success run FutRule", {

  # Success run
  results = FutRule(normalCase)
  expect_is(results, "FutRuleResults")
  expect_named(results, c("parameters", "sim_summary"))

  sim_summary = results$sim_summary
  expect_named(sim_summary, c('sensitivity', 'specificity', 'accuracy', 'cp_threshold'))
  expect_length(sim_summary$sensitivity, 101)
  expect_length(sim_summary$specificity, 101)
  expect_length(sim_summary$accuracy, 101)
  expect_length(sim_summary$cp_threshold, 101)

  accuracy = sim_summary$accuracy
  cp_threshold = sim_summary$cp_threshold

  index = which.max(accuracy)
  optimal_point = round(100 * cp_threshold[index], 1)
  level = 0.95 * max(accuracy)
  zone = cp_threshold[accuracy >= level]
  optimal_lower = round(100 * min(zone), 1)
  optimal_upper = round(100 * max(zone), 1)

  expect_true(
    abs(optimal_lower - 6) < 4, 
    info = paste("optimal_lower(",optimal_lower,") is out of range 6±4"))
  expect_true(
    abs(optimal_upper - 65) < 15,
    info = paste("optimal_upper(",optimal_upper,") is out of range 65±15"))

  # Check for report generation
  FutRuleReport(results, tempfile(fileext = ".docx"))
})

test_that("Input parameters errors check FutRule", {
  # Errors check
  expect_error(
    FutRule(
      c("Not a list")
    ),
    info = "Checking for wrong parameters collection type"
  )

  expect_error(
    FutRule(
      list(
        # missing
        #endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing endpoint type"
  )

  expect_error(
    FutRule(
      list(
        # wrong
        endpoint_type = "SomeWrongType", # normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong endpoint type"
  )

  expect_error(
    FutRule(
      list(
        endpoint_type = normalCase$endpoint_type,
        # missing
        #sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing number of enrolled patients (sample_size)"
  )

  expect_error(
    FutRule(
      list(
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # missing
        #info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing information fraction"
  )

  expect_error(
    FutRule(
      list(
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        # missing
        #control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effect in the control arm"
  )

  expect_error(
    FutRule(
      list(
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        # missing
        #treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effects in the treatment arms"
  )

  expect_error(
    FutRule(
      list(
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        # missing
        #control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviation in the control arm"
  )

  expect_error(
    FutRule(
      list(
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        # missing
        #treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviations in the treatment arms"
  )

})
