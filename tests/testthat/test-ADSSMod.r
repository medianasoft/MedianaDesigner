# Test the input of the ADSSMod function

# Normal case parameters
normalCase = list(

  # Endpoint type
  endpoint_type = "Normal",

  # Direction of favorable outcome (Higher or Lower)
  # Default: Higher
  #direction = "Higher",

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

# Binary case parameters
binaryCase = list(

  # Endpoint type
  endpoint_type = "Binary",

  # Direction of favorable outcome (Higher or Lower)
  # Default: Higher
  direction = "Higher",

  # Number of enrolled patients (control, treatment) 
  sample_size = c(100, 200),

  # Patient dropout rate
  dropout_rate = 0.1,

  # Response rate in the control arm 
  control_rate = 0.1,

  # Response rate in the treatment arm 
  treatment_rate = 0.25,

  # Information fractions at IA1, IA2, FA (before sample size adjustment) 
  # and FA (after sample size adjustment)
  info_frac = c(0.4, 0.6, 1, 1.4),

  # Futility threshold for conditional power at IA1
  futility_threshold = 0.2,

  # Promising interval for conditional power at IA2
  promising_interval = c(0.5, 0.9),

  # Target conditional power for increasing the sample size at IA2
  target_power = 0.9,

  # One-sided alpha level
  alpha = 0.025,

  # Number of simulations
  nsims = 100
)

# Time-to-event case parameters
timeToEventCase = list(

  # Endpoint type
  endpoint_type = "Time-to-event",

  # Direction of favorable outcome (Higher or Lower)
  # Default: Higher
  direction = "Higher",

  # Number of enrolled patients (control, treatment) 
  sample_size = c(220, 220),

  # Annual patient dropout rate
  dropout_rate = 0.05,

  # Median time in the control arm 
  control_time = 7.5,

  # Median time in the treatment arm 
  treatment_time = 10.5,

  # Target event count at FA (before event count adjustment)
  event_count = 300,

  # Information fractions at IA1, IA2, FA (before event count adjustment) 
  # and FA (after event count adjustment)
  info_frac = c(0.4, 0.6, 1, 1.3),

  # Futility threshold for conditional power at IA1
  futility_threshold = 0.1,

  # Promising interval for conditional power at IA2
  promising_interval = c(0.5, 0.9),

  # Target conditional power for increasing the number of events at IA2
  target_power = 0.9,

  # Enrollment period
  enrollment_period = 12,

  # Median enrollment time
  enrollment_parameter = 8,

  # One-sided alpha level
  alpha = 0.025,

  # Number of simulations
  nsims = 100
)

context("ADSSMod - Success runs")

test_that("Success run ADSSMod with Normal case", {

  # Set the seed of R‘s random number generator.
  # It also takes effect to Rcpp randome generation functions.
  # https://stackoverflow.com/questions/60119621/get-the-same-sample-of-integers-from-rcpp-as-base-r
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  set.seed(5)

  # Success run with check default values
  results = ADSSMod(
	  list(
	    endpoint_type      = normalCase$endpoint_type,
	    sample_size        = normalCase$sample_size,
	    # check defaults:
	    #dropout_rate       = normalCase$dropout_rate,
	    control_mean       = normalCase$control_mean,
	    control_sd         = normalCase$control_sd,
	    treatment_mean     = normalCase$treatment_mean,
	    treatment_sd       = normalCase$treatment_sd,
	    info_frac          = normalCase$info_frac,
	    futility_threshold = normalCase$futility_threshold,
	    promising_interval = normalCase$promising_interval,
	    target_power       = normalCase$target_power #,
	    # check defaults:
	    #alpha              = normalCase$alpha,
	    #nsims              = normalCase$nsims
	  )
	)
  expect_type(results$sim_results, "double")
  expect_type(results$sim_summary, "list")
  expect_length(results$sim_results, 19000)
  expect_true(abs(results$sim_summary$futility - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$increase - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$ad_power - 0.5) < 0.3)
  expect_true(abs(results$sim_summary$ad_under - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$ad_prom - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$ad_over - 0.9) < 0.1)
  expect_true(abs(results$sim_summary$trad_under - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$trad_prom - 0.7) < 0.3)
  expect_true(abs(results$sim_summary$trad_over - 0.9) < 0.1)

  # Check for report generation
  ADSSModReportDoc(results)
  GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Success run ADSSMod with Binary case", {

  # Success run
  results = ADSSMod(binaryCase)
  expect_type(results$sim_results, "double")
  expect_type(results$sim_summary, "list")
  expect_length(results$sim_results, 1900)

  expect_true(abs(results$sim_summary$futility - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$increase - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$ad_power - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$ad_prom - 0.9) < 0.2)
  expect_true(abs(results$sim_summary$ad_over - 0.9) < 0.2)
  expect_true(abs(results$sim_summary$trad_prom - 0.8) < 0.3)
  expect_true(abs(results$sim_summary$trad_over - 0.9) < 0.2)

  # Check for report generation
  ADSSModReportDoc(results)
})

test_that("Success run ADSSMod with Time-to-event case", {

  # Success run
  results = ADSSMod(timeToEventCase)
  expect_type(results$sim_results, "double")
  expect_type(results$sim_summary, "list")
  expect_length(results$sim_results, 1900)

  expect_true(abs(results$sim_summary$futility - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$increase - 0.2) < 0.2)
  expect_true(abs(results$sim_summary$ad_power - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$ad_prom - 0.9) < 0.2)
  expect_true(abs(results$sim_summary$ad_over - 0.9) < 0.2)
  expect_true(abs(results$sim_summary$trad_prom - 0.8) < 0.3)
  expect_true(abs(results$sim_summary$trad_over - 0.9) < 0.2)

  # Check for report generation
  ADSSModReportDoc(results)
})

test_that("Success run ADSSMod with Time-to-event case with Lower direction", {

  # Success run
  results = ADSSMod(
    list(
      endpoint_type = timeToEventCase$endpoint_type,
      direction = "Lower",
      sample_size = timeToEventCase$sample_size,
      dropout_rate = timeToEventCase$dropout_rate,
      control_time = timeToEventCase$control_time,
      treatment_time = timeToEventCase$treatment_time,
      event_count = timeToEventCase$event_count,
      info_frac = timeToEventCase$info_frac,
      futility_threshold = timeToEventCase$futility_threshold,
      promising_interval = timeToEventCase$promising_interval,
      target_power = timeToEventCase$target_power,
      enrollment_period = timeToEventCase$enrollment_period,
      enrollment_parameter = timeToEventCase$enrollment_parameter,
      alpha = timeToEventCase$alpha,
      nsims = timeToEventCase$nsims
    )
  )
  expect_type(results$sim_results, "double")
  expect_type(results$sim_summary, "list")
  expect_length(results$sim_results, 1900)

  # Check for report generation
  ADSSModReportDoc(results)
})


context("ADSSMod - Error checks")

test_that("Input parameters errors check ADSSMod", {


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
        # wrong
        direction          = "Wrong", # normalCase$direction,
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
    info = "Checking for incorrect value of direction"
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
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong
        info_frac = c(0.4, 0.6, 1, 1.3, 1.5),  # c(0.4, 0.6, 1, 1.3) = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power = normalCase$target_power,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong Information fractions at IA1, IA2, FA (incorrect value size)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong
        info_frac = c(0, 0.6, 1, 1.3),  # c(0.4, 0.6, 1, 1.3) = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power = normalCase$target_power,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong Information fractions at IA1, IA2, FA (value <= 0)"
  )

  expect_error(
    ADSSMod(
      list(
        endpoint_type = normalCase$endpoint_type,
        sample_size = normalCase$sample_size,
        dropout_rate = normalCase$dropout_rate,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong
        info_frac = c(0.4, 0.6, 1, 3.001),  # c(0.4, 0.6, 1, 1.3) = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        promising_interval = normalCase$promising_interval,
        target_power = normalCase$target_power,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong Information fractions at IA1, IA2, FA (value > 3)"
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

  testParameterErrors = function(params, paramName, paramDesc, 
    checkMissing = TRUE, checkSize = TRUE, checkMin = NA, checkMax = NA) {

    func = ADSSMod

    paramDesc = paste0(paramDesc, " (", paramName, ")")
    if (!is.null(params$endpoint_type)) 
      paramDesc = paste0(params$endpoint_type, ": ", paramDesc)
    # Missing
    if (checkMissing) {
      testParams = params
      testParams[paramName] <- NULL
      expect_error(func(testParams), 
        info = paste0("Checking for missing ", paramDesc))
    }
    # Check size
    if (checkSize) {
      testParams = params
      testParams[[paramName]] <- append(testParams[[paramName]], testParams[[paramName]][1])
      expect_error(func(testParams), 
        info = paste0("Checking for wrong ", paramDesc, " (incorrect value size)"))
    }
    # Check below min value
    if (!is.null(checkMin) && !is.na(checkMin)) {
      testParams = params
      testParams[[paramName]][1] <- checkMin
      expect_error(func(testParams), 
        info = paste0("Checking for wrong ", paramDesc, " (incorrect value < min)"))
    }
    # Check under max value
    if (!is.null(checkMax) && !is.na(checkMax)) {
      testParams = params
      testParams[[paramName]][length(testParams[[paramName]])] <- checkMax
      expect_error(func(testParams), 
        info = paste0("Checking for wrong ", paramDesc, " (incorrect value > max)"))
    }
  }

  testParameterErrors(normalCase, 
    'sample_size', 
    'Number of enrolled patients',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = 1001)

  testParameterErrors(normalCase, 
    'futility_threshold', 
    'Futility threshold at IA1',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = -0.001,
    checkMax = 0.999)

  testParameterErrors(normalCase, 
    'promising_interval', 
    'Promising interval at IA2',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0.001,
    checkMax = 0.999)

  testParameterErrors(normalCase, 
    'target_power', 
    'Target conditional power',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0.001,
    checkMax = 0.999)

  testParameterErrors(normalCase, 
    'dropout_rate', 
    'Patient dropout rate',
    checkMissing = FALSE,
    checkSize = TRUE,
    checkMin = -0.001,
    checkMax = 1)

  testParameterErrors(normalCase, 
    'nsims', 
    'Number of simulations',
    checkMissing = FALSE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = 10001)

  testParameterErrors(normalCase, 
    'alpha', 
    'One-sided Type I error rate',
    checkMissing = FALSE,
    checkSize = TRUE,
    checkMin = 0.001,
    checkMax = 0.5)

  testParameterErrors(normalCase, 
    'control_mean', 
    'Mean effect in the control arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = NA,
    checkMax = NA)

  testParameterErrors(normalCase, 
    'treatment_mean', 
    'Mean effect in the treatment arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = NA,
    checkMax = NA)

  testParameterErrors(normalCase, 
    'control_sd', 
    'Standard deviation in the control arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(normalCase, 
    'treatment_sd', 
    'Standard deviation in the treatment arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(binaryCase, 
    'control_rate', 
    'Response rate in the control arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = 1)

  testParameterErrors(binaryCase, 
    'treatment_rate', 
    'Response rate in the treatment arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = 1)

  testParameterErrors(timeToEventCase, 
    'control_time', 
    'Median time in the control arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(timeToEventCase, 
    'treatment_time', 
    'Median time in the treatment arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(timeToEventCase, 
    'enrollment_period', 
    'Patient enrollment period',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(timeToEventCase, 
    'enrollment_parameter', 
    'Median enrollment time',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = timeToEventCase$enrollment_period)

  testParameterErrors(timeToEventCase, 
    'event_count', 
    'Target number of events at FA',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = sum(timeToEventCase$sample_size))

})

test_that("Input parameters errors check ADSSModReportDoc", {

  expect_error(
    ADSSModReportDoc(""),
    info = "Checking for wrong parameter type for report generator"
  )

})
