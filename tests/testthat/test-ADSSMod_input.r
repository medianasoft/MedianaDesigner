# Test the input of the ADSSMod function

# Normal case
normalCase = list(

  # Endpoint type
  endpoint_type = "Normal",

  # Number of enrolled patients (control, treatment) 
  sample_size = c(120, 120),

  # Annual patient dropout rate based on an exponential dropout distribution
  dropout_rate = 0.05,

  # Mean and SD in the control arm 
  control_mean = 0,
  control_sd = 1,

  # Mean and SD in the treatment arm 
  treatment_mean = 0.3,
  treatment_sd = 1,

  # Information fractions at IA1, IA2, FA (before sample size adjustment) and FA (after sample size adjustment)
  info_frac = c(0.4, 0.6, 1, 1.3),

  # Futility threshold for conditional power at IA1
  futility_threshold = 0.1,

  # Promising interval for conditional power at IA2
  promising_interval = c(0.5, 0.9),

  # Target conditional power for increasing the number of events at IA2
  target_power = 0.9,

  # One-sided alpha level
  alpha = 0.025,

  # Number of simulations
  nsims = 100
)

# Run simulations to compute key characteristics

test_that("Check endpoint_type", {

  # Success run
  results = ADSSMod(normalCase)
  expect_type(results$sim_results, "double")
  expect_type(results$sim_summary, "list")
  expect_length(results$sim_results, 1900)
  expect_true(abs(results$sim_summary$futility - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$increase - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$ad_power - 0.5) < 0.2)
  expect_true(abs(results$sim_summary$ad_under - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$ad_prom - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$ad_over - 0.9) < 0.1)
  expect_true(abs(results$sim_summary$trad_under - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$trad_prom - 0.7) < 0.3)
  expect_true(abs(results$sim_summary$trad_over - 0.9) < 0.1)

  # Check for report generation
  ADSSModReport(results, tempfile(fileext = ".docx"))

  # Errors check
  expect_error(
    ADSSMod(
      c("Not a list")
    ),
    info = "Checking for wrong parameters collection type"
  )

  expect_error(
    ADSSMod(
      list(
        # miss
        # endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing endpoint type"
  )

  expect_error(
    ADSSMod(
      list(
        # wrong
        endpoint_type      = "SomeStrangeType",
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for incorrect value of endpoint type"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        # missing
        #sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing sample size"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        # missing
        #info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing information fractions"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        # invalid value
        info_frac          = c(0.4, 0.6, 1),      # Right value is c(0.4, 0.6, 1, 1.3)
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for invalid information fractions"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        # invalid value
        info_frac          = c(0.6, 0.4, 1, 1.3),     # Right value is c(0.4, 0.6, 1, 1.3)
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for inconsistent information fractions (first value greater than second)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        # invalid value
        info_frac          = c(0.4, 1.1, 1, 1.3),     # Right value is c(0.4, 0.6, 1, 1.3)
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for inconsistent information fractions (second value greater than third)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        # invalid value
        info_frac          = c(0.4, 0.6, 0.9, 1.3),     # Right value is c(0.4, 0.6, 1, 1.3)
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for inconsistent information fractions (third value is not 1)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        # invalid value
        info_frac          = c(0.4, 0.6, 1, 0.8),     # Right value is c(0.4, 0.6, 1, 1.3)
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for inconsistent information fractions (third value greater than forth)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        # missing
        #futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing futility threshold"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        # missing
        #promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing promising interval"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        # wrong
        promising_interval = c(0.5), #normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for wrong promising interval (short vector)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        # wrong
        promising_interval = c(0.9, 0.5), #normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for wrong promising interval (first value greater than second)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        # missing
        #target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing target conditional power"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        # missing
        #control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effect in the control arm"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        # missing
        #treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effect in the treatment arm"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        # missing
        #control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviation in the control arm"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type      = normalCase$endpoint_type,
        sample_size        = normalCase$sample_size,
        dropout_rate       = normalCase$dropout_rate,
        control_mean       = normalCase$control_mean,
        control_sd         = normalCase$control_sd,
        treatment_mean     = normalCase$treatment_mean,
        # missing
        #treatment_sd       = normalCase$treatment_sd,
        info_frac          = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power       = normalCase$target_power,
        alpha              = normalCase$alpha,
        nsims              = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviation in the treatment arm"
  )

})
