# Test the input of the ADTreatSel function

# Normal case parameters
normalCase = list(

  # Number of trial arms
  narms = 3,

  # Sample size
  sample_size = c(120, 120, 120),

  # Primary endpoint's type
  endpoint_type = "Normal",

  # Direction of favorable outcome (Higher or Lower)
  # Default: Higher
  #direction = "Higher",

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

  # Number of selected treatments
  treatment_count = 1,

  # Multiple testing procedure (Bonferroni, Holm or Hochberg)
  mult_test = "Bonferroni",

  # Dropout rate at the end of the treatment period
  dropout_rate = 0.05,

  # One-sided Type I error rate
  alpha = 0.025,

  # Number of simulations
  nsims = 1000
)

# Binary case parameters
binaryCase = list(

  # Endpoint type
  endpoint_type = "Binary",

  # Direction of favorable outcome (Higher or Lower)
  # Default: Higher
  direction = "Higher",

  # Number of enrolled patients (control, three treatments)
  sample_size = c(150, 150, 150, 150),

  # Patient dropout rate
  dropout_rate = 0.1,

  # Response rate in the control arm 
  control_rate = 0.1,

  # Response rates in the treatment arms 
  treatment_rate = c(0.25, 0.25, 0.25),

  # Information fractions at IA1, IA2, FA
  info_frac = c(0.4, 0.6, 1),

  # Futility threshold for conditional power at IA1
  futility_threshold = 0.15,

  # Number of selected treatments
  treatment_count = 1,

  # Multiple testing procedure (Bonferroni, Holm or Hochberg)
  mult_test = "Bonferroni",

  # One-sided alpha level
  alpha = 0.025,

  # Number of simulations
  nsims = 1000
)

# Time-to-event case parameters
timeToEventCase = list(

  # Endpoint type
  endpoint_type = "Time-to-event",

  # Direction of favorable outcome (Higher or Lower)
  # Default: Higher
  direction = "Higher",

  # Number of enrolled patients (control, two treatments)
  sample_size = c(180, 180, 180),

  # Annual patient dropout rate
  dropout_rate = 0.05,

  # Median time in the control arm 
  control_time = 7.5,

  # Median times in the treatment arms
  treatment_time = c(10.5, 10.5),

  # Target event count at FA
  event_count = 450,

  # Information fractions at IA1, IA2, FA
  info_frac = c(0.4, 0.6, 1),

  # Futility threshold for conditional power at IA1
  futility_threshold = 0.1,

  # Number of selected treatments
  treatment_count = 1,

  # Multiple testing procedure (Bonferroni, Holm or Hochberg)
  mult_test = "Bonferroni",

  # Enrollment period
  enrollment_period = 24,

  # Median enrollment time
  enrollment_parameter = 18,

  # One-sided alpha level
  alpha = 0.025,

  # Number of simulations
  nsims = 1000
)

context("ADTreatSel - Success runs")

test_that("Success run ADTreatSel with Normal case", {

  # Success run
  results = ADTreatSel(
    list(
      narms = normalCase$narms,
      sample_size = normalCase$sample_size,
      endpoint_type = normalCase$endpoint_type,
      control_mean = normalCase$control_mean,
      control_sd = normalCase$control_sd,
      treatment_mean = normalCase$treatment_mean,
      treatment_sd = normalCase$treatment_sd,
      info_frac = normalCase$info_frac,
      futility_threshold = normalCase$futility_threshold,
      treatment_count = normalCase$treatment_count,
      mult_test = normalCase$mult_test,
      # use default
      #dropout_rate = normalCase$dropout_rate,
      alpha = normalCase$alpha #,
      # use default
      #nsims = normalCase$nsims
    )
  )
  expect_is(results, "ADTreatSelResults")
  expect_type(  results$sim_results, "double")
  expect_length(results$sim_results, 7000)
  
  expect_type(    results$sim_summary, "list")
  expect_true(abs(results$sim_summary$ad_power - 0.5) < 0.3)
  expect_true(abs(results$sim_summary$overall_futility - 0.1) < 0.1)
  
  expect_is(      results$sim_summary$select, "numeric")
  expect_length(  results$sim_summary$select, 2)
  expect_true(abs(results$sim_summary$select[1] - 0.4) < 0.1)
  expect_true(abs(results$sim_summary$select[2] - 0.4) < 0.1)

  expect_is(      results$sim_summary$futility, "numeric")
  expect_length(  results$sim_summary$futility, 2)
  expect_true(abs(results$sim_summary$futility[1] - 0.2) < 0.1)
  expect_true(abs(results$sim_summary$futility[2] - 0.2) < 0.1)

  expect_is(      results$sim_summary$trad_power, "numeric")
  expect_length(  results$sim_summary$trad_power, 2)
  expect_true(abs(results$sim_summary$trad_power[1] - 0.5) < 0.2)
  expect_true(abs(results$sim_summary$trad_power[2] - 0.5) < 0.2)

  # Check for report generation
  ADTreatSelReportDoc(results)
})

test_that("Success run ADTreatSel with Binary case", {

  # Success run
  results = ADTreatSel(
    list(
      endpoint_type = binaryCase$endpoint_type,
      sample_size = binaryCase$sample_size,
      dropout_rate = binaryCase$dropout_rate,
      control_rate = binaryCase$control_rate,
      treatment_rate = binaryCase$treatment_rate,
      info_frac = binaryCase$info_frac,
      futility_threshold = binaryCase$futility_threshold,
      treatment_count = binaryCase$treatment_count,
      mult_test = binaryCase$mult_test,
      # use default
      #alpha = binaryCase$alpha,
      nsims = binaryCase$nsims
    )
  )
  expect_is(results, "ADTreatSelResults")
  expect_type(  results$sim_results, "double")
  expect_length(results$sim_results, 10000)
  
  expect_type(    results$sim_summary, "list")
  expect_true(abs(results$sim_summary$ad_power - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$overall_futility - 0.1) < 0.1)
  
  expect_is(      results$sim_summary$select, "numeric")
  expect_length(  results$sim_summary$select, 3)
  expect_true(abs(results$sim_summary$select[1] - 0.35) < 0.1)
  expect_true(abs(results$sim_summary$select[2] - 0.33) < 0.1)
  expect_true(abs(results$sim_summary$select[3] - 0.31) < 0.1)
  
  expect_is(      results$sim_summary$futility, "numeric")
  expect_length(  results$sim_summary$futility, 3)
  expect_true(abs(results$sim_summary$futility[1] - 0.1) < 0.1)
  expect_true(abs(results$sim_summary$futility[2] - 0.1) < 0.1)
  expect_true(abs(results$sim_summary$futility[3] - 0.1) < 0.1)

  expect_is(      results$sim_summary$trad_power, "numeric")
  expect_length(  results$sim_summary$trad_power, 3)
  expect_true(abs(results$sim_summary$trad_power[1] - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$trad_power[2] - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$trad_power[3] - 0.8) < 0.2)

  # Check for report generation
  ADTreatSelReportDoc(results)
})

test_that("Success run ADTreatSel with Time-to-event case", {

  # Success run
  results = ADTreatSel(timeToEventCase)
  expect_is(results, "ADTreatSelResults")
  expect_type(  results$sim_results, "double")
  expect_length(results$sim_results, 7000)
  
  expect_type(    results$sim_summary, "list")
  expect_true(abs(results$sim_summary$ad_power - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$overall_futility - 0.1) < 0.1)
  
  expect_is(      results$sim_summary$select, "numeric")
  expect_length(  results$sim_summary$select, 2)
  expect_true(abs(results$sim_summary$select[1] - 0.49) < 0.1)
  expect_true(abs(results$sim_summary$select[2] - 0.47) < 0.1)

  expect_is(      results$sim_summary$futility, "numeric")
  expect_length(  results$sim_summary$futility, 2)
  expect_true(abs(results$sim_summary$futility[1] - 0.1) < 0.1)
  expect_true(abs(results$sim_summary$futility[2] - 0.1) < 0.1)

  expect_is(      results$sim_summary$trad_power, "numeric")
  expect_length(  results$sim_summary$trad_power, 2)
  expect_true(abs(results$sim_summary$trad_power[1] - 0.8) < 0.2)
  expect_true(abs(results$sim_summary$trad_power[2] - 0.8) < 0.2)

  # Check for report generation
  ADTreatSelReportDoc(results)
})

test_that("Success run ADTreatSel with Normal case and short sample_size vector", {

  # Success run
  results = ADTreatSel(
    list(
      narms = normalCase$narms,
      sample_size = c(120, 120),  #normalCase$sample_size = c(120, 120, 120)
      endpoint_type = normalCase$endpoint_type,
      control_mean = normalCase$control_mean,
      control_sd = normalCase$control_sd,
      treatment_mean = c(0.3), # treatment_mean = c(0.3, 0.3),
      treatment_sd = c(1),     # treatment_sd = c(1, 1),
      info_frac = normalCase$info_frac,
      futility_threshold = normalCase$futility_threshold,
      treatment_count = normalCase$treatment_count,
      mult_test = normalCase$mult_test,
      dropout_rate = normalCase$dropout_rate,
      alpha = normalCase$alpha,
      nsims = normalCase$nsims
    )
  )

  expect_is(results, "ADTreatSelResults")
  expect_type(  results$sim_results, "double")
  expect_length(results$sim_results, 4000)
  
  expect_type(    results$sim_summary, "list")
  expect_true(abs(results$sim_summary$ad_power - 0.5) < 0.3)
  expect_true(abs(results$sim_summary$overall_futility - 0.2) < 0.1)
  
  expect_is(      results$sim_summary$select, "numeric")
  expect_length(  results$sim_summary$select, 1)
  expect_true(abs(results$sim_summary$select[1] - 0.7) < 0.2)

  expect_is(      results$sim_summary$futility, "numeric")
  expect_length(  results$sim_summary$futility, 1)
  expect_true(abs(results$sim_summary$futility[1] - 0.2) < 0.1)

  expect_is(      results$sim_summary$trad_power, "numeric")
  expect_length(  results$sim_summary$trad_power, 1)
  expect_true(abs(results$sim_summary$trad_power[1] - 0.5) < 0.2)

  # Check for report generation
  ADTreatSelReportDoc(results)
  GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Success run ADTreatSel with Lower direction", {

  # Success run
  results = ADTreatSel(
    list(
      narms = normalCase$narms,
      sample_size = normalCase$sample_size,
      direction = "Lower",
      endpoint_type = normalCase$endpoint_type,
      control_mean = normalCase$control_mean,
      control_sd = normalCase$control_sd,
      treatment_mean = normalCase$treatment_mean,
      treatment_sd = normalCase$treatment_sd,
      info_frac = normalCase$info_frac,
      futility_threshold = normalCase$futility_threshold,
      treatment_count = normalCase$treatment_count,
      mult_test = normalCase$mult_test,
      dropout_rate = normalCase$dropout_rate,
      alpha = normalCase$alpha,
      nsims = normalCase$nsims
    )
  )

  expect_is(results, "ADTreatSelResults")
  expect_type(  results$sim_results, "double")
  expect_length(results$sim_results, 7000)
  
  expect_type(    results$sim_summary, "list")
  expect_is(      results$sim_summary$select, "numeric")
  expect_is(      results$sim_summary$futility, "numeric")
  expect_is(      results$sim_summary$trad_power, "numeric")

  # Check for report generation
  ADTreatSelReportDoc(results)
})

context("ADTreatSel - Error checks")

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
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
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
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
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
        endpoint_type = normalCase$endpoint_type,
        # wrong
        direction = "Wrong",
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong direction"
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
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
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
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
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
        info_frac = c(0.4, 1.0, 0.9),  # normalCase$info_frac == c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
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
        info_frac = c(0.4, 0.6, 0.9),  # normalCase$info_frac == c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
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
        # wrong
        info_frac = c(0.4, 0.6), # c(0.4, 0.6, 1.0) = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong information fractions at IA1, IA2, FA (incorrect length)"
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
        # wrong
        info_frac = c(0, 0.6, 1.0), # c(0.4, 0.6, 1.0) = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong information fractions at IA1, IA2, FA (value <= 0)"
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
        # wrong
        info_frac = c(0.4, 0.6, 1.1), # c(0.4, 0.6, 1.0) = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong information fractions at IA1, IA2, FA (value > 1)"
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
        # wrong
        info_frac = c(0.4, 0.6, "1.0"), # c(0.4, 0.6, 1.0) = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        treatment_count = normalCase$treatment_count,
        mult_test = normalCase$mult_test,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong information fractions at IA1, IA2, FA (value not a number)"
  )

  # Universal parameter check func
  testParameterErrors = function(params, paramName, paramDesc, 
    checkMissing = TRUE, checkWrong = NA, checkSize = TRUE, checkMin = NA, checkMax = NA) {

    func = ADTreatSel

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
    # Wrong
    if (!is.null(checkWrong) && !is.na(checkWrong)) {
      testParams = params
      testParams[paramName] <- checkWrong
      expect_error(func(testParams), 
        info = paste0("Checking for wrong ", paramDesc))
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
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 1001)

  testParameterErrors(normalCase, 
    'futility_threshold', 
    'Futility threshold at IA1',
    checkMissing = TRUE,
    checkSize = FALSE,
    checkMin = -0.001,
    checkMax = 0.999)

  testParameterErrors(normalCase, 
    'dropout_rate', 
    'Patient dropout rate',
    checkMissing = FALSE,
    checkSize = FALSE,
    checkMin = -0.001,
    checkMax = 1)

  testParameterErrors(normalCase, 
    'nsims', 
    'Number of simulations',
    checkMissing = FALSE,
    checkSize = FALSE,
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
    'Mean effects in the treatment arms',
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
    'Standard deviations in the treatment arms',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(normalCase, 
    'treatment_count', 
    'Number of selected treatments',
    checkMissing = TRUE,
    checkWrong = "Not a number",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = length(normalCase$sample_size))

  testParameterErrors(normalCase, 
    'mult_test', 
    'Multiple testing procedure',
    checkMissing = TRUE,
    checkWrong = "Wrong value",
    checkSize = FALSE,
    checkMin = NA,
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
    'Responses rate in the treatment arms',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = 1)

  testParameterErrors(binaryCase, 
    'treatment_rate', 
    'Responses rate in the treatment arms',
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
    'Median times in the treatment arms',
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

test_that("Input parameters errors check ADTreatSelReportDoc", {

  expect_error(
    ADTreatSelReportDoc(""),
    info = "Checking for wrong parameter type for report generator"
  )

})
