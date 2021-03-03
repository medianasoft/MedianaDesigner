# Test the input of the ADTreatSel function

# Normal case
normalCase = list(

  # Number of trial arms
  narms = 3,

  # Sample size
  sample_size = c(120, 120, 120),

  # Primary endpoint's type
  endpoint_type = "Normal",

  # Mean and SD in the control arm 
  control_mean = 0,
  control_sd = 1,

  # Mean and SD in the treatment arm 
  treatment_mean = c(0.3, 0.3),
  treatment_sd = c(1, 1),

  # Information fractions at IA1, IA2, FA
  info_frac = c(0.4, 0.6, 1.0),

  # Futility threshold for conditional power at IA1
  futility_threshold = 0.1,

  # Dropout rate at the end of the treatment period
  dropout_rate = 0.05,

  # One-sided Type I error rate
  alpha = 0.025,

  # Number of simulations
  nsims = 1000
)

# Run simulations to compute key characteristics

test_that("Success run ADTreatSel", {

  # Success run
  results = ADTreatSel(normalCase)
  expect_is(results, "ADTreatSelResults")
  expect_type(  results$sim_results, "double")
  expect_length(results$sim_results, 13000)
  
  expect_type(    results$sim_summary, "list")
  expect_true(abs(results$sim_summary$ad_power - 0.5) < 0.2)
  expect_true(abs(results$sim_summary$overall_futility - 0.1) < 0.1)
  
  expect_is(      results$sim_summary$select, "numeric")
  expect_length(  results$sim_summary$select, 3)
  expect_true(abs(results$sim_summary$select[1] - 0.4) < 0.1)
  expect_true(abs(results$sim_summary$select[2] - 0.4) < 0.1)
  expect_true(abs(results$sim_summary$select[3] - 0.1) < 0.1)

  expect_is(      results$sim_summary$futility, "numeric")
  expect_length(  results$sim_summary$futility, 2)
  expect_true(abs(results$sim_summary$futility[1] - 0.2) < 0.1)
  expect_true(abs(results$sim_summary$futility[2] - 0.2) < 0.1)

  expect_is(      results$sim_summary$trad_power, "numeric")
  expect_length(  results$sim_summary$trad_power, 2)
  expect_true(abs(results$sim_summary$trad_power[1] - 0.5) < 0.2)
  expect_true(abs(results$sim_summary$trad_power[2] - 0.5) < 0.2)

  # Check for report generation
  ADTreatSelReport(results, tempfile(fileext = ".docx"))
})

test_that("Input parameters errors check ADTreatSel", {
  # Errors check
  expect_error(
    ADTreatSel(
      c("Not a list")
    ),
    info = "Checking for wrong parameters collection type"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        # missing
        #endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing endpoint type"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        # wrong
        endpoint_type = "SomeType", #normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing endpoint type"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        # missing
        #sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing sample size"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # missing
        #info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing information fractions at IA1, IA2, FA"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong (IA1 >= IA2)
        info_frac = c(0.8, 0.6, 1.0),  # normalCase$info_frac == c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong information fraction (IA1 must be < IA2)"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong (IA2 >= FA)
        info_frac = c(0.4, 1.6, 1.0),  # normalCase$info_frac == c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong information fraction (IA2 must be < FA)"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong (FA != 1)
        info_frac = c(0.4, 0.6, 1.1),  # normalCase$info_frac == c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong information fraction (FA must be 1)"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        # missing
        #futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing futility threshold at IA1"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        # missing
        #control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effect in the control arm"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        # missing
        #treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effects in the treatment arms"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        # missing
        #control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviation in the control arm"
  )

  expect_error(
    ADTreatSel(
      list(
        narms = normalCase$narms,
        sample_size = normalCase$sample_size,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        # missing
        #treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviations in the treatment arms"
  )

})
