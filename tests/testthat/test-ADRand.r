commonNSim = 50
isTestMultiCore = FALSE

#####################################

# List of all parameters.

normalCase <- list(

  # --[ Endpoint Parameters ]--------------------------------------------------
  # Endpoint distribution type (Normal)
  endpoint_type = "Normal",
  # Direction of favorable outcome (Lower, Higher)
  direction = "Higher",

  # --[ Design Assumptions ]---------------------------------------------------
  # Dose levels in the trial (with 0 corresponding to the placebo arm)
  dose_levels = c(0, 150, 300, 600),
  # Total number of enrolled patients in each trial stage
  stage_sample_size = c(100, 100, 100, 100),

  # Length of the patient enrollment period
  enrollment_period = 52,
  # Median enrollment time, i.e., the time point by which
  # 50% of the patients will be enrolled into the trial
  enrollment_parameter = 36,

  # Length of the treatment period
  treatment_period = 6,
  # Patient dropout rate
  dropout_rate = 0.1,

  #--[ Treatment Assumptions ]-------------------------------------------------
  # Mean and SD in the control arm
  control_mean = 10,
  control_sd = 25,
  # Mean and SD in the treatment arms
  treatment_mean = c(15, 20, 22),
  treatment_sd = c(25, 25, 25),

  # Fixed randomization ratio in the placebo arm
  ratio_placebo = 0.25,

  # Balance parameter for adaptive randomization
  balance = 1,

  # Non-linear parameters of the candidate dose-response models used
  # in the MCPMod method
  exponential_model_parameter = 600,
  emax_model_parameter = 600,
  logistic_model_parameters = c(200, 200),

  # --[ General Parameters ]---------------------------------------------------
  # Clinically relevant improvement over placebo
  delta = 10,
  # One-sided Type I error rate
  alpha = 0.025,
  # Number of simulations
  nsims = commonNSim

)

#####################################

context("ADRand - Success runs")

checkExpectationsForNormalCase = function(res) {
  expect_s3_class(res, "ADRandResults")

  expect_length(res, 2)
  expect_length(res$sim_results, 5)
  
  # Calculate summary
  parameters = res$parameters
  simulations = res$sim_results

  # print(simulations$n)

  # - ComparisonOfTraditionalAndAdaptiveDesigns -----------------------------------------
  rowMax = function(x) {
    row_max = rep(0, nrow(x))
    for (i in 1:nrow(x)) row_max[i] = max(x[i, ])
    return(row_max)
  }
  traditional = mean(rowMax(simulations$traditional) >= simulations$mcpmod_value)
  adaptive = mean(rowMax(simulations$adaptive) >= simulations$mcpmod_value)
  expect_equal(traditional, 0.9, tolerance=0.1)
  expect_equal(adaptive,    0.9, tolerance=0.1)
  # -------------------------------------------------------------------------------------
  # - SampleSizeByTrialArm --------------------------------------------------------------
  col3 = NULL
  n_doses = length(parameters$dose_levels)
  for (i in 1:n_doses) {
    col3 = c(col3, round(summary(simulations$n[, i])[c(1, 3, 4, 6)], 1))
  }
  col3 = unname(col3)
  expect_equal(
    col3,
    # Min. Median   Mean   Max.   Min. Median   Mean   Max.   Min. Median   Mean   Max.   Min. Median   Mean   Max.
    c(89.0,  90.0,  90.2,  92.0,  79.0,  85.5,  87.2, 144.0,  64.0,  90.0,  89.8,  96.0,  63.0,  94.0,  93.6, 102.0),
    tolerance=10
  )
  # -------------------------------------------------------------------------------------
}

test_that("Success run ADRand with Normal case (single core)", {  

  parameters = normalCase
  # Run once with random_seed parameter
  parameters$random_seed = 49283
  # Skip chart generation in tests
  parameters$withoutCharts = TRUE

  # Run simulations
  results = ADRand(parameters)
  checkExpectationsForNormalCase(results)
  expect_length(results$parameters, length(normalCase)+12)

  # Create a simulation report
  temp_file = tempfile("Simulation report.docx", fileext=".docx")
  GenerateReport(results, temp_file)
  expect_true(file.exists(temp_file))
  expect_true(file.size(temp_file) > 35000)
  
    # Check for report generation
  ADRandReportDoc(results)
  GenerateReport(results, tempfile(fileext = ".docx"))
})

if (isTestMultiCore) {
  test_that("Success run ADRand with Normal case (two cores)", {  
    parameters = normalCase
    parameters$ncores = 2
    # Skip chart generation in tests
    parameters$withoutCharts = TRUE

    # Run simulations
    results = ADRand(parameters)
    checkExpectationsForNormalCase(results)
  })
}

test_that("Success run ADRand with Normal case with insignificant changes", {  
  # Remove parameters that accepts default value
  changedNormalCase <- normalCase
  changedNormalCase['direction'] = NULL
  changedNormalCase['dropout_rate'] = NULL
  changedNormalCase['alpha'] = NULL
  changedNormalCase['balance'] = NULL
  changedNormalCase['nsims'] = 25

  # Run simulations
  results = ADRand(changedNormalCase)
  expect_is(results, "ADRandResults")
  expect_equal(length(results), 2)
  expect_equal(length(results$sim_results), 5)
})

test_that("Success run ADRand with Normal case with changed direction", {  
  # Remove parameters that accept default value
  changedNormalCase <- normalCase
  changedNormalCase['direction'] = "Lower"
  changedNormalCase['delta'] = -10
  changedNormalCase['nsims'] = 25

  # Run simulations
  results = ADRand(changedNormalCase)
  expect_is(results, "ADRandResults")
  expect_equal(length(results), 2)
  expect_equal(length(results$sim_results), 5)

  # TODO: Check result
})

# TODO: More success runs

context("ADRand - Error checks")

test_that("Input parameters errors check ADRand", {

  # Errors check
  expect_error(
    ADRand(
      c("Not a list")
    ),
    info = "Checking for wrong parameters collection type"
  )

  # Universal parameter check func
  testParameterErrors = function(params, paramName, paramDesc, 
    checkMissing = TRUE, checkWrong = NA, checkSize = TRUE, checkMin = NA, checkMax = NA) {

    func = ADRand

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
    'endpoint_type', 
    'Endpoint Type',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = NA,
    checkMax = NA)

  # direction = "Higher",
  testParameterErrors(normalCase, 
    'direction', 
    'Direction',
    checkMissing = FALSE,
    checkWrong = "Back",
    checkSize = FALSE,
    checkMin = NA,
    checkMax = NA)

  # dose_levels = c(0, 150, 300, 600),
  testParameterErrors(normalCase, 
    'dose_levels', 
    'Dose levels',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = -1,
    checkMax = 1001)

  # stage_sample_size = c(100, 100, 100, 100),
  testParameterErrors(normalCase, 
    'stage_sample_size', 
    'Stage sample size',
    checkMissing = TRUE,
    checkWrong = FALSE,
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  # enrollment_period = 52,
  testParameterErrors(normalCase, 
    'enrollment_period', 
    'Enrollment period',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  # enrollment_parameter = 36,
  testParameterErrors(normalCase, 
    'enrollment_parameter', 
    'Enrollment Parameter',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = normalCase$enrollment_period)

  # treatment_period = 6,
  testParameterErrors(normalCase, 
    'treatment_period', 
    'Treatment Period',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  # dropout_rate = 0.1,
  testParameterErrors(normalCase, 
    'dropout_rate', 
    'Dropout Rate',
    checkMissing = FALSE,   # default = 0
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = -0.01,
    checkMax = 1)

  # control_mean = 10,
  testParameterErrors(normalCase, 
    'control_mean', 
    'Control Mean',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = NA,
    checkMax = NA)

  # control_sd = 25,
  testParameterErrors(normalCase, 
    'control_sd', 
    'Control SD',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  # treatment_mean = c(15, 20, 22),
  testParameterErrors(normalCase,
    'treatment_mean', 
    'Treatment Mean',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = TRUE,
    checkMin = NA,
    checkMax = NA)

  # treatment_sd = c(25, 25, 25),
  testParameterErrors(normalCase,
    'treatment_sd', 
    'Treatment SD',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  # ratio_placebo = 0.25,
  testParameterErrors(normalCase,
    'ratio_placebo', 
    'Ratio Placebo',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 1)

  # exponential_model_parameter = 600,
  testParameterErrors(normalCase,
    'exponential_model_parameter', 
    'Exponential model parameter',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  # emax_model_parameter = 600,
  testParameterErrors(normalCase,
    'emax_model_parameter', 
    'Emax model parameter',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  # logistic_model_parameters = c(200, 200),
  testParameterErrors(normalCase,
    'logistic_model_parameters', 
    'Logistic model parameters',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = TRUE,
    checkMin = 0,
    checkMax = NA)

  # delta = 10,
  testParameterErrors(normalCase,
    'delta', 
    'Delta',
    checkMissing = TRUE,
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = -0.01,
    checkMax = NA)

  # alpha = 0.025,
  testParameterErrors(normalCase,
    'alpha', 
    'Alpha',
    checkMissing = FALSE,   # default = 0.025
    checkWrong = "String",
    checkSize = FALSE,
    checkMin = 0.001,
    checkMax = 0.5)

  # nsims = 100
  testParameterErrors(normalCase,
    'nsims', 
    'Number of simulations',
    checkMissing = FALSE,   # default = 1000
    checkWrong = 99.5,
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 10001)

})

test_that("Input parameters errors check ADRandReportDoc", {

  expect_error(
    ADRandReportDoc(""),
    info = "Checking for wrong parameter type for report generator"
  )

})
