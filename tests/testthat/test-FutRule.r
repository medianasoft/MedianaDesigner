# Test the input of the FutRule function

# Normal case parameters
normalCase = list(
  # Endpoint type
  endpoint_type = "Normal",

  # Direction of favorable outcome (Higher or Lower)
  # Default: Higher
  #direction = "Higher",

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

# Binary case parameters
binaryCase = list(

  # Endpoint type
  endpoint_type = "Binary",

  # Direction of favorable outcome (Higher or Lower)
  direction = "Higher",

  # Number of enrolled patients (control, three treatments)
  sample_size = c(75, 75, 75, 75),

  # Dropout rate
  dropout_rate = 0.05,

  # Response rate in the control arm 
  control_rate = 0.35,

  # Response rates in the treatment arms 
  treatment_rate = c(0.45, 0.5, 0.55),

  # Information fraction
  info_frac = 0.3,

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

  # Number of enrolled patients (control, treatment)
  sample_size = c(125, 125),

  # Target event count at the final analysis
  event_count = 175,

  # Annual patient dropout rate
  dropout_rate = 0.05,

  # Median time in the control arm 
  control_time = 7.5,

  # Median time in the treatment arm
  treatment_time = 10.5,

  # Information fraction
  info_frac = 0.6,

  # Enrollment period
  enrollment_period = 12,

  # Median enrollment time
  enrollment_parameter = 9,

  # One-sided alpha level
  alpha = 0.025,

  # Number of simulations
  nsims = 1000
)

context("FutRule - Success runs")

test_that("Success run FutRule with Normal case", {

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

  expect_equal(optimal_lower, 6, 5, 
    info = paste("optimal_lower(",optimal_lower,") is out of range 6±5"))
  expect_true(
    abs(optimal_upper - 65) < 30,
    info = paste("optimal_upper(",optimal_upper,") is out of range 65±30"))

  # Check for report generation
  FutRuleReportDoc(results)
  GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Success run ModuleD with Binary case", {

  # Success run 
  results = FutRule(binaryCase)

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
    abs(optimal_lower - 19) < 15, 
    info = paste("optimal_lower(",optimal_lower,") is out of range 19±15"))
  expect_true(
    abs(optimal_upper - 95) < 15,
    info = paste("optimal_upper(",optimal_upper,") is out of range 95±15"))

  # Check for report generation
  FutRuleReportDoc(results)
})

test_that("Success run FutRule with Time-to-event case", {

  # Success run 
  results = FutRule(timeToEventCase)

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
    abs(optimal_lower - 5) < 5, 
    info = paste("optimal_lower(",optimal_lower,") is out of range 5±5"))
  expect_true(
    abs(optimal_upper - 50) < 20,
    info = paste("optimal_upper(",optimal_upper,") is out of range 50±20"))

  # Check for report generation
  FutRuleReportDoc(results)
})

test_that("Success run FutRule with default parameters", {

  # Success run 
  results = FutRule(
    list(
      endpoint_type = normalCase$endpoint_type,
      sample_size = normalCase$sample_size,
      #dropout_rate = normalCase$dropout_rate,
      control_mean = normalCase$control_mean,
      control_sd = normalCase$control_sd,
      treatment_mean = normalCase$treatment_mean,
      treatment_sd = normalCase$treatment_sd,
      info_frac = normalCase$info_frac #,
      # alpha = normalCase$alpha,
      # nsims = normalCase$nsims
    )
  )
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
    abs(optimal_upper - 65) < 30,
    info = paste("optimal_upper(",optimal_upper,") is out of range 65±30"))
})

test_that("Success run FutRule with Lower direction", {

  # Success run 
  results = FutRule(
    list(
      endpoint_type = normalCase$endpoint_type,
      sample_size = normalCase$sample_size,
      direction = "Lower",
      #dropout_rate = normalCase$dropout_rate,
      control_mean = normalCase$control_mean,
      control_sd = normalCase$control_sd,
      treatment_mean = normalCase$treatment_mean,
      treatment_sd = normalCase$treatment_sd,
      info_frac = normalCase$info_frac #,
      # alpha = normalCase$alpha,
      # nsims = normalCase$nsims
    )
  )
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

  FutRuleReportDoc(results)
})

context("FutRule - Error checks")

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
        # wrong
        direction = "Wrong",
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
    info = "Checking for wrong direction"
  )

  testParameterErrors = function(params, paramName, paramDesc, 
    checkMissing = TRUE, checkSize = TRUE, checkMin = NA, checkMax = NA) {

    func = FutRule

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
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 1001)

  testParameterErrors(normalCase, 
    'info_frac', 
    'Information fraction',
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

  testParameterErrors(binaryCase, 
    'control_rate', 
    'Response rate in the control arm',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = 1)

  testParameterErrors(binaryCase, 
    'treatment_rate', 
    'Response rates in the treatment arms',
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
    'Target number of events',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = sum(timeToEventCase$sample_size))

})

test_that("Input parameters errors check FutRuleReportDoc", {

  expect_error(
    FutRuleReportDoc(""),
    info = "Checking for wrong parameter type for report generator"
  )

})
