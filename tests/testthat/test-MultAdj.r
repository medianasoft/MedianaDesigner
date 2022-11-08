commonNSim = 50
isTestMultiCore = FALSE

################################################

# Case 1A 
# One endpoint (n_endpoints = 1), several comparisons (n_comparisons >= 2)
# Normal endpoint (endpoint_type = "Normal")

parametersCase1A = list(
    # Endpoint type
    endpoint_type = "Normal",

    # Direction of beneficial effect
    direction = "Higher",

    # Number of dose-control comparisons
    n_comparisons = 2,

    # Number of endpoints
    n_endpoints = 1,

    # Number of enrolled patients (control, multiple treatments) 
    sample_size = c(120, 120, 120),

    # Patient dropout rate
    dropout_rate = 0.05,

    # Endpoint information 
    control_mean = 0,
    treatment_mean = c(0.25, 0.30),
    control_sd = 1,
    treatment_sd = c(1, 1),

    # Multiple testing procedure
    mult_test = "Chain",

    # Hypothesis weights
    weights = c(0.3, 0.7),

    # Hypothesis transition matrix 
    transition = matrix(c(0, 1, 0, 0), 2, 2, byrow = TRUE),

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = commonNSim
)

################################################

# Case 1B
# One endpoint (n_endpoints = 1), several comparisons (n_comparisons >= 2)
# Binary endpoint (endpoint_type = "Binary")

parametersCase1B = list(
    # Endpoint type
    endpoint_type = "Binary",

    # Direction of beneficial effect
    direction = "Lower",

    # Number of dose-control comparisons
    n_comparisons = 3,

    # Number of endpoints
    n_endpoints = 1,

    # Number of enrolled patients (control, multiple treatments) 
    sample_size = c(100, 150, 150, 150),

    # Patient dropout rate
    dropout_rate = 0.05,

    # Endpoint information 
    control_rate = 0.6,
    treatment_rate = c(0.5, 0.4, 0.3),

    # Multiple testing procedure
    mult_test = "Fixed-sequence",

    # Hypothesis testing sequence
    sequence = 3:1,

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = commonNSim
)

################################################

# Case 2A
# Several endpoints (n_endpoints >= 2), one comparison (n_comparisons = 1)
# Normal endpoint (endpoint_type = "Normal")

# List of all parameters

parametersCase2A = list(
    # Endpoint type
    endpoint_type = "Normal",

    # Direction of beneficial effect
    direction = "Higher",

    # Number of dose-control comparisons
    n_comparisons = 1,

    # Number of endpoints
    n_endpoints = 2,

    # Number of enrolled patients (control, treatment) 
    sample_size = c(120, 120),

    # Patient dropout rate
    dropout_rate = 0.05,

    # Endpoint information 
    control_mean = c(0, 0),
    treatment_mean = c(0.25, 0.30),
    control_sd = c(1, 1),
    treatment_sd = c(1, 1),

    # Multiple testing procedure
    mult_test = "O'Brien",

    # Endpoint correlation matrix
    endpoint_correlation = matrix(c(1, 0.3, 0.3, 1), 2, 2),

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = commonNSim
)

################################################

# Case 2B
# Several endpoints (n_endpoints >= 2), one comparison (n_comparisons = 1)
# Binary endpoint (endpoint_type = "Binary")

parametersCase2B = list(
    # Endpoint type
    endpoint_type = "Binary",

    # Direction of beneficial effect
    direction = "Higher",

    # Number of dose-control comparisons
    n_comparisons = 1,

    # Number of endpoints
    n_endpoints = 3,

    # Number of enrolled patients (control, treatment) 
    sample_size = c(80, 80),

    # Patient dropout rate
    dropout_rate = 0.05,

    # Endpoint information 
    control_rate = c(0.1, 0.1, 0.1),
    treatment_rate = c(0.25, 0.30, 0.35),

    # Multiple testing procedure
    mult_test = "Holm",

    # Endpoint correlation matrix
    endpoint_correlation = matrix(c(1, 0.3, 0.3, 
                                            0.3, 1, 0.3,
                                            0.3, 0.3, 1), 3, 3),

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = commonNSim
)

################################################

# Case 3A
# Several endpoints (n_endpoints >= 2), several comparisons (n_comparisons >= 2)
# Normal endpoint (endpoint_type = "Normal")

parametersCase3A = list(
    # Endpoint type
    endpoint_type = "Normal",

    # Direction of beneficial effect
    direction = "Higher",

    # Number of dose-control comparisons
    n_comparisons = 2,

    # Number of endpoints
    n_endpoints = 3,

    # Number of enrolled patients (control, multiple treatments) 
    sample_size = c(120, 120, 120),

    # Patient dropout rate
    dropout_rate = 0.05,

    # Endpoint information (rows corresponds to endpoints and columns corresponds to treatment-control comparisons)
    control_mean = c(0, 0, 0),
    treatment_mean = matrix(c(0.1, 0.25, 
                              0.20, 0.30, 
                              0.35, 0.40), 3, 2, byrow = TRUE),

    control_sd = c(1, 1, 1),
    treatment_sd = matrix(c(1, 1, 
                            1, 1, 
                            1, 1), 3, 2, byrow = TRUE),

    # Component procedure to be used in the gatekeeping procedure
    mult_test = "Hochberg",

    # Mixture method used in the gatekeeping procedure
    mult_method = "Standard",

    # Truncation parameters in the gatekeeping procedure
    mult_test_gamma = c(0.8, 0.8, 1),

    # Endpoint correlation matrix
    endpoint_correlation = matrix(c(1, 0.3, 0.3, 
                                    0.3, 1, 0.3,
                                    0.3, 0.3, 1), 3, 3),

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = commonNSim
)

################################################

# Case 3B
# Several endpoints (n_endpoints >= 2), several comparisons (n_comparisons >= 2)
# Binary endpoint (endpoint_type = "Binary")

parametersCase3B = list(
    # Endpoint type
    endpoint_type = "Binary",

    # Direction of beneficial effect
    direction = "Higher",

    # Number of dose-control comparisons
    n_comparisons = 3,

    # Number of endpoints
    n_endpoints = 2,

    # Number of enrolled patients (control, multiple treatments) 
    sample_size = c(120, 120, 120, 120),

    # Patient dropout rate
    dropout_rate = 0.05,

    # Endpoint information (rows corresponds to endpoints and columns corresponds to treatment-control comparisons)
    control_rate = c(0.1, 0.1),
    treatment_rate = matrix(c(0.2, 0.3, 0.4, 
                              0.3, 0.4, 0.5), 2, 3, byrow = TRUE),

    # Multiple testing procedure used in the gatekeeping procedure
    mult_test = "Hommel",

    # Mixture method used in the gatekeeping procedure
    mult_method = "Modified",

    # Truncation parameters in the gatekeeping procedure
    mult_test_gamma = c(0.8, 1),

    # Endpoint correlation matrix
    endpoint_correlation = matrix(c(1, 0.3, 
                                    0.3, 1), 2, 2),

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = commonNSim
)

#####################################

context("MultAdj - Success runs")

checkExpectationsForCase1A = function(results) {
  expect_s3_class(results, "MultAdjResults")
  expect_named(results, c("parameters", "sim_results", "sim_summary"))

  expect_equal(nrow(results$sim_results), results$parameters$nsims)

  sim_summary = results$sim_summary
  # print(sim_summary)
  expect_named(sim_summary, c('power', 'adj_power', 'disj_power', 'conj_power'))
  expect_length(sim_summary, 4)
  expect_length(sim_summary$power, 2)
  expect_length(sim_summary$adj_power, 2)

  expect_equal(unname(sim_summary$power[1]), 0.43, tolerance = 0.1)
  expect_equal(unname(sim_summary$power[2]), 0.66, tolerance = 0.1)
  expect_equal(unname(sim_summary$adj_power[1]), 0.29, tolerance = 0.1)
  expect_equal(unname(sim_summary$adj_power[2]), 0.55, tolerance = 0.1)
  expect_equal(unname(sim_summary$disj_power), 0.6, tolerance = 0.1)
  expect_equal(unname(sim_summary$conj_power), 0.25, tolerance = 0.1)
}

test_that("Success run MultAdj with Case 1A (single core)", {  

  # Run simulations
  params = parametersCase1A
  # First run with random seed
  params$random_seed = 49283

  results = MultAdj(params)
  checkExpectationsForCase1A(results)

  # Check for report generation
  GenerateReport(results, tempfile(fileext = ".docx"))
})

if (isTestMultiCore) {
  test_that("Success run MultAdj with Case 1A (two cores)", {  
    # Run simulations
    params = parametersCase1A
    params$ncores = 2
    results = MultAdj(params)
    checkExpectationsForCase1A(results)
  })
}

test_that("Success run MultAdj with Case 1B", {  

  # cat("\nTest Case 1B: begin...\n")

  # Run simulations
  results = MultAdj(parametersCase1B)
  expect_is(results, "MultAdjResults")

  sim_summary = results$sim_summary
  expect_length(sim_summary, 4)
  expect_length(sim_summary$power, 3)
  expect_length(sim_summary$adj_power, 3)
  expect_equal(unname(sim_summary$power[1]), 0.4, tolerance = 0.4)
  expect_equal(unname(sim_summary$power[2]), 0.8, tolerance = 0.3)
  expect_equal(unname(sim_summary$power[3]), 1.0, tolerance = 0.3)
  expect_equal(unname(sim_summary$adj_power[1]), 0.4, tolerance = 0.3)
  expect_equal(unname(sim_summary$adj_power[2]), 0.8, tolerance = 0.3)
  expect_equal(unname(sim_summary$adj_power[3]), 1.0, tolerance = 0.3)
  expect_equal(unname(sim_summary$disj_power), 1.0, tolerance = 0.3)
  expect_equal(unname(sim_summary$conj_power), 0.4, tolerance = 0.3)

  # cat("\nTest Case 1B: done.\n")

  # Check for report generation
  GenerateReport(results, tempfile(fileext = ".docx"))
})

checkExpectationsForCase2A = function(results) {
  expect_s3_class(results, "MultAdjResults")
  expect_named(results, c("parameters", "sim_results", "sim_summary"))

  expect_equal(nrow(results$sim_results), results$parameters$nsims)

  sim_summary = results$sim_summary
  # print(sim_summary)
  expect_named(sim_summary, c('power', 'adj_power'))
  expect_length(sim_summary, 2)
  expect_length(sim_summary$power, 2)
  expect_equal(unname(sim_summary$power[1]), 0.43, tolerance = 0.1)
  expect_equal(unname(sim_summary$power[2]), 0.62, tolerance = 0.1)
  expect_equal(unname(sim_summary$adj_power), 0.67, tolerance = 0.1)
}

test_that("Success run MultAdj with Case 2A (single core)", {  

  # cat("\nTest Case 2A: begin...\n")

  # Run simulations
  results = MultAdj(parametersCase2A)
  checkExpectationsForCase2A(results)

  # cat("\nTest Case 2A: done.\n")

  # Check for report generation
  GenerateReport(results, tempfile(fileext = ".docx"))
})

if (isTestMultiCore) {
  test_that("Success run MultAdj with Case 2A (two cores)", {  

    # cat("\nTest Case 2A-two: begin...\n")

    parameters = parametersCase2A
    parameters$ncores = 2
    # Run simulations
    results = MultAdj(parameters)
    checkExpectationsForCase2A(results)

    # cat("\nTest Case 2A-two: done.\n")
  })
}

test_that("Success run MultAdj with Case 2B", {  
  # Run simulations
  results = MultAdj(parametersCase2B)
  expect_is(results, "MultAdjResults")

  sim_summary = results$sim_summary
  expect_length(sim_summary, 4)
  expect_length(sim_summary$power, 3)
  expect_length(sim_summary$adj_power, 3)
  expect_equal(unname(sim_summary$power[1]), 0.7, tolerance = 0.4)
  expect_equal(unname(sim_summary$power[2]), 0.9, tolerance = 0.3)
  expect_equal(unname(sim_summary$power[3]), 0.9, tolerance = 0.3)
  expect_equal(unname(sim_summary$adj_power[1]), 0.7, tolerance = 0.3)
  expect_equal(unname(sim_summary$adj_power[2]), 0.9, tolerance = 0.3)
  expect_equal(unname(sim_summary$adj_power[3]), 0.9, tolerance = 0.3)
  expect_equal(unname(sim_summary$disj_power), 1, tolerance = 0.3)
  expect_equal(unname(sim_summary$conj_power), 0.6, tolerance = 0.4)

  # Check for report generation
  GenerateReport(results, tempfile(fileext = ".docx"))
})

checkExpectationsForCase3A = function(results) {
  expect_s3_class(results, "MultAdjResults")
  expect_named(results, c("parameters", "sim_results", "sim_summary"))

  expect_equal(nrow(results$sim_results), results$parameters$nsims)

  sim_summary = results$sim_summary
  # print(sim_summary)
  expect_named(sim_summary, c('power', 'adj_power', 'disj_power', 'conj_power'))
  expect_length(sim_summary, 4)
  expect_length(sim_summary$power, 6)
  expect_length(sim_summary$adj_power, 6)
  expect_length(sim_summary$disj_power, 3)
  expect_length(sim_summary$conj_power, 3)
  expect_equal(sum(sim_summary$power), 3.3, tolerance = 0.1)
  expect_equal(sum(sim_summary$adj_power), 0.8, tolerance = 0.2)
  expect_equal(sum(sim_summary$disj_power), 0.7, tolerance = 0.2)
  expect_equal(sum(sim_summary$conj_power)/3, 0.1, tolerance = 0.1)
}

test_that("Success run MultAdj with Case 3A (single core)", {

  # cat("\nTest Case 3A: begin...\n")

  # Run simulations
  results = MultAdj(parametersCase3A)
  checkExpectationsForCase3A(results)

  # cat("\nTest Case 3A: done.\n")

  # Check for report generation
  GenerateReport(results, tempfile(fileext = ".docx"))
})

if (isTestMultiCore) {
  test_that("Success run MultAdj with Case 3A (two cores)", {

    # cat("\nTest Case 3A-two: begin...\n")

    # Run simulations
    parameters = parametersCase3A
    parameters$ncores = 2
    results = MultAdj(parameters)
    checkExpectationsForCase3A(results)

    # cat("\nTest Case 3A-two: done.\n")
  })
}

test_that("Success run MultAdj with Case 3B", {  
  # Run simulations

  # Defaults check (remove parameter)
  params <- parametersCase3B
  params$direction <- NULL
  params$nsims <- NULL
  params$alpha <- NULL

  results = MultAdj(params)
  expect_is(results, "MultAdjResults")

  sim_summary = results$sim_summary
  expect_length(sim_summary, 4)
  expect_length(sim_summary$power, 6)
  expect_length(sim_summary$adj_power, 6)
  expect_length(sim_summary$disj_power, 2)
  expect_length(sim_summary$conj_power, 2)
  expect_equal(sum(sim_summary$power)/6, 0.9, tolerance = 0.3)
  expect_equal(sum(sim_summary$adj_power)/6, 0.9, tolerance = 0.3)
  expect_equal(sum(sim_summary$disj_power)/2, 1.0, tolerance = 0.3)
  expect_equal(sum(sim_summary$conj_power)/2, 0.6, tolerance = 0.3)

  # Check for report generation
  GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Success run MultAdj with default dropout rate", {
  params = parametersCase1A
  params["dropout_rate"] <- NULL
  results = MultAdj(params)
  expect_is(results, "MultAdjResults")
})

test_that("Success run MultAdj with default weights", {
  params = parametersCase1A
  params["weights"] <- NULL
  results = MultAdj(params)
  expect_is(results, "MultAdjResults")
})

#####################################

context("MultAdj - Error checks")

test_that("Input parameters errors check MultAdj", {

  # Errors check
  expect_error(
    MultAdj(
      c("Not a list")
    ),
    info = "Checking for wrong parameters collection type"
  )

  expect_error(
    {
      params = parametersCase1A
      params$n_comparisons = 1
      params$n_endpoints = 1
      MultAdj(params)
    },
    info = "Checking for wrong parameters n_comparisons and n_endpoints"
  )

  expect_error(
    {
      params = parametersCase3A
      params$treatment_mean = c(0.1, 0.25, 
                            0.20, 0.30, 
                            0.35, 0.40)
      MultAdj(params)
    },
    info = "Checking for wrong parameter type treatment_mean (not matrix)"
  )

  expect_error(
    {
      params = parametersCase3A
      params$treatment_mean = matrix(c(0.1, 0.25, 
                                      0.20, 0.30, 
                                      0.35, 0.40), 2, 3, byrow = TRUE)
      MultAdj(params)
    },
    info = "Checking for wrong parameter type treatment_mean (not correct dimentions)"
  )

  expect_error(
    {
      params = parametersCase3A
      params$treatment_sd = c(1,1,1,1,1,1)
      MultAdj(params)
    },
    info = "Checking for wrong parameter type treatment_sd (not matrix)"
  )

  expect_error(
    {
      params = parametersCase3A
      params$treatment_sd = matrix(c(1, 1, 
                                      1, 1, 
                                      1, 1), 2, 3, byrow = TRUE)
      MultAdj(params)
    },
    info = "Checking for wrong parameter type treatment_sd (not correct dimentions)"
  )

  expect_error(
    {
      params = parametersCase3B
      params$treatment_rate = c(0.2, 0.3, 0.4, 
                              0.3, 0.4, 0.5)
      MultAdj(params)
    },
    info = "Checking for wrong parameter type treatment_rate (not matrix)"
  )

  expect_error(
    {
      params = parametersCase3B
      params$treatment_rate = matrix(c(0.2, 0.3, 0.4, 
                              0.3, 0.4, 0.5), 3, 2, byrow = TRUE)
      MultAdj(params)
    },
    info = "Checking for wrong parameter type treatment_rate (not correct dimentions)"
  )

  expect_error(
    {
      params = parametersCase1A
      params$mult_test <- NULL
      MultAdj(params)
    },
    info = "Checking for without parameter mult_test #1"
  )

  expect_error(
    {
      params = parametersCase2A
      params$mult_test <- NULL
      MultAdj(params)
    },
    info = "Checking for without parameter mult_test #2"
  )

  expect_error(
    {
      params = parametersCase3A
      params$mult_test <- NULL
      MultAdj(params)
    },
    info = "Checking for without parameter mult_test #3"
  )

  expect_error(
    {
      params = parametersCase3A
      params$mult_test <- "NotHochbergOrHommel"
      MultAdj(params)
    },
    info = "Checking for wrong parameter mult_test"
  )

  expect_error(
    {
      params = parametersCase1A
      params$weights = c(0.3, 0.8)
      MultAdj(params)
    },
    info = "Checking for error in parameter weights (sum > 1)"
  )

  expect_error(
    {
      params = parametersCase1A
      params$transition = c(0, 1, 0, 0)
      MultAdj(params)
    },
    info = "Checking for error in parameter transition (not matrix)"
  )

  expect_error(
    {
      params = parametersCase1A
      params$transition = matrix(c(0, 1, 0, 0), 1, 4, byrow = TRUE)
      MultAdj(params)
    },
    info = "Checking for error in parameter transition (wrong dimention)"
  )

  expect_error(
    {
      params = parametersCase1A
      params$transition = matrix(c(0.1, 1, 0, 0), 2, 2, byrow = TRUE)
      MultAdj(params)
    },
    info = "Checking for error in parameter transition (row > 1)"
  )

  expect_error(
    {
      params = parametersCase2B
      params$weights = c(0, 0.2, 0.7)
      MultAdj(params)
    },
    info = "Checking for error in parameter weights (elements must be > 0)"
  )

  expect_error(
    {
      params = parametersCase2B
      params$weights = c(0.1, 0.3, 1.0)
      MultAdj(params)
    },
    info = "Checking for error in parameter weights (elements must be < 1)"
  )

  expect_error(
    {
      params = parametersCase2B
      params$weights = c(0.1, 0.2, 0.8)
      MultAdj(params)
    },
    info = "Checking for error in parameter weights (sum of elements must be <= 1)"
  )

  expect_error(
    {
      params = parametersCase1B
      params$sequence = c(1,2,2)
      MultAdj(params)
    },
    info = "Checking for error in parameter sequence"
  )

  expect_error(
    {
      params = parametersCase3A
      params$mult_method = "Random"
      MultAdj(params)
    },
    info = "Checking for error in parameter mult_method"
  )

  expect_error(
    {
      params = parametersCase3A
      params$mult_test_gamma = c(0.8, 1, 1)
      MultAdj(params)
    },
    info = "Checking for error in parameter mult_test_gamma"
  )

  expect_error(
    {
      params = parametersCase2A
      params$endpoint_correlation = c(1, 0.3, 0.3, 1)
      MultAdj(params)
    },
    info = "Checking for error in parameter endpoint_correlation (not matrix)"
  )

  expect_error(
    {
      params = parametersCase2A
      params$endpoint_correlation = matrix(c(1, 0.3, 0.3, 1), 4, 1)
      MultAdj(params)
    },
    info = "Checking for error in parameter endpoint_correlation (wrong size)"
  )

  expect_error(
    {
      params = parametersCase2A
      params$endpoint_correlation = matrix(c(1, 1, 1, 1), 2, 2)
      MultAdj(params)
    },
    info = "Checking for error in parameter endpoint_correlation (det <=0)"
  )

  # Universal parameter check func
  testParameterErrors = function(params, paramName, paramDesc, 
    checkMissing = TRUE, checkWrong = NA, checkSize = TRUE, checkMin = NA, checkMax = NA) {

    func = MultAdj

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

  testParameterErrors(parametersCase1A, 
    'endpoint_type', 
    'Endpoint Type',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = NA,
    checkMax = NA)

  testParameterErrors(parametersCase1A, 
    'direction', 
    'Direction',
    checkMissing = FALSE,
    checkWrong = "Back",
    checkSize = FALSE,
    checkMin = NA,
    checkMax = NA)

  testParameterErrors(parametersCase1A, 
    'n_comparisons', 
    'n_comparisons',
    checkMissing = TRUE,
    checkWrong = 1.5,
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 11)

  testParameterErrors(parametersCase1A, 
    'n_endpoints', 
    'n_endpoints',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 11)

  testParameterErrors(parametersCase1A, 
    'sample_size', 
    'Sample Size',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 10001)

  testParameterErrors(parametersCase1A,
    'dropout_rate', 
    'Dropout Rate',
    checkMissing = FALSE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = -0.1,
    checkMax = 1)

  testParameterErrors(parametersCase1A,
    'nsims', 
    'Number of simulations',
    checkMissing = FALSE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 9,
    checkMax = 10001)

  testParameterErrors(parametersCase1A,
    'alpha', 
    'Alpha',
    checkMissing = FALSE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 0.001,
    checkMax = 0.5)

  testParameterErrors(parametersCase1A,
    'control_mean', 
    'Control Mean',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = NA,
    checkMax = NA)

  testParameterErrors(parametersCase1A,
    'treatment_mean', 
    'Treatment Mean',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = NA,
    checkMax = NA)

  testParameterErrors(parametersCase1A,
    'control_sd', 
    'Control SD',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(parametersCase1A,
    'treatment_sd', 
    'Treatment SD',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = NA)

  testParameterErrors(parametersCase1B,
    'control_rate', 
    'Control Rate',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 1)

  testParameterErrors(parametersCase1B,
    'treatment_rate', 
    'Treatment Rate',
    checkMissing = TRUE,
    checkWrong = "wrong",
    checkSize = FALSE,
    checkMin = 0,
    checkMax = 1)
})

#####################################

test_that("Input parameters errors check MultAdjReportDoc", {

  expect_error(
    MultAdjReportDoc(""),
    info = "Checking for wrong parameter type for report generator"
  )

})
