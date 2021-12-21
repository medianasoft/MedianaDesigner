commonNSim = 100
isTestMultiCore = FALSE

# Base case
baseCase = list(
  # Load the built-in data set with the patient enrollment, event and dropout 
  # information (EventPredData)
  data_set = EventPredData,

  # Future time points for computing event predictions
  time_points = seq(from = 12, to = 24, by = 1),

  # Prior distribution for the event hazard rate based on the 
  # expected median time of 15 months and the uncertainty parameter of 0.3 
  event_prior_distribution = 
      EventPredPriorDistribution(expected = log(2) / 15, uncertainty = 0.3),

  # Prior distribution for the patient dropout hazard rate based on the 
  # expected median time of 80 months and the uncertainty parameter of 0.3
  dropout_prior_distribution = 
      EventPredPriorDistribution(expected = log(2) / 80, uncertainty = 0.3),

  # Prior distribution for the patient enrollment rate based on the
  # expected enrollment rate of 35 patients per month and the uncertainty 
  # parameter of 0.3 
  enrollment_prior_distribution = 
      EventPredPriorDistribution(expected = 35, uncertainty = 0.3),

  # Number of simulations
  nsims = commonNSim
)

context("EventPred - Success runs")

checkExpectationsForBaseCase1 = function(res) {
  expect_s3_class(res, "EventPredResults")
  expect_named(res, c('parameters', 'interim_analysis', 'predictions'))

  # check expected result size
  expect_s3_class(res$interim_analysis, "data.frame")
  expect_equal(dim(res$interim_analysis), c(457,3))
  expect_is(res$predictions, "matrix")
  expect_equal(dim(res$predictions), c(13, 3))

  # check predictions
  predictions = res$predictions
  sum_all_predictions = Reduce('+', predictions)
  expect_equal(sum_all_predictions, 8660.471, tolerance=0.01)
}

test_that("Success run EventPred (one core)", {

  # Success run
  results = EventPred(
    list(
        data_set = baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,
        # missing, use default 1000
        #nsims = baseCase$nsims
        # Skip chart generation in tests
        withoutCharts = TRUE,
        # First run with initial random seed
        random_seed = 49283
    )
  )
  checkExpectationsForBaseCase1(results)

  # Check for report generation
  EventPredReportDoc(results)
  GenerateReport(results, tempfile(fileext = ".docx"))

})

if (isTestMultiCore) {
  test_that("Success run EventPred (two cores)", {
    # Success run
    params = baseCase
    params$ncores = 2
    results = EventPred(params)
    checkExpectationsForBaseCase1(results)
  })
}

context("EventPred - Error checks")

test_that("Input parameters errors check EventPred", {
  # Errors check
  expect_error(
    EventPred(
      c("Not a list")
    ),
    info = "Checking for wrong parameters collection type"
  )

  expect_error(
    EventPred(
      list(
        # missing
        #data_set = baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    info = "Checking for missing dataset"
  )

  # Drop valiable from dataset
  brokenData = subset(EventPredData, select=-c(enrollment))
  expect_error(
    EventPred(
      list(
        # broken dataset
        data_set = brokenData, # baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    "Enrollment variable must be specified.",
    info = "Checking for missing enrollment variable in dataset"
  )

  # Drop valiable from dataset
  brokenData = subset(EventPredData, select=-c(time))
  expect_error(
    EventPred(
      list(
        # broken dataset
        data_set = brokenData, # baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    "Time variable must be specified.",
    info = "Checking for missing time variable in dataset"
  )

  # Drop valiable from dataset
  brokenData = subset(EventPredData, select=-c(event))
  expect_error(
    EventPred(
      list(
        # broken dataset
        data_set = brokenData, # baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    "Event variable must be specified.",
    info = "Checking for missing event variable in dataset"
  )

  # Drop valiable from dataset
  brokenData = subset(EventPredData, select=-c(dropout))
  expect_error(
    EventPred(
      list(
        # broken dataset
        data_set = brokenData, # baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    "Dropout variable must be specified.",
    info = "Checking for missing dropout variable in dataset"
  )

  # Universal parameter check func
  testParameterErrors = function(params, paramName, paramDesc, 
    checkMissing = TRUE, checkWrong = NA, checkSize = TRUE, checkMin = NA, checkMax = NA) {

    func = EventPred

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
    # if (!is.null(checkMin) && !is.na(checkMin)) {
    if (!is.null(checkMin) && !anyNA(checkMin, recursive = FALSE)) {
      testParams = params
      testParams[[paramName]] <- checkMin
      expect_error(func(testParams), 
        info = paste0("Checking for wrong ", paramDesc, " (incorrect value < min)"))
    }
    # Check under max value
    if (!is.null(checkMax) && !is.na(checkMax)) {
      testParams = params
      testParams[[paramName]] <- checkMax
      expect_error(func(testParams), 
        info = paste0("Checking for wrong ", paramDesc, " (incorrect value > max)"))
    }
  }

  # Prepare test values
  data_set = baseCase$data_set
  time_points = baseCase$time_points
  interim_analysis = max(data_set$enrollment + data_set$time)
  interim_analysis_lower_min = interim_analysis - 0.01
  time_points[1] = interim_analysis_lower_min

  testParameterErrors(baseCase,
    "time_points",
    "Future time points for computing event predictions",
    checkMissing = TRUE,
    checkWrong = "Not a number",
    checkSize = FALSE,
    checkMin = time_points,
    checkMax = NA)

  testParameterErrors(baseCase, 
    'event_prior_distribution', 
    'Prior distribution for the event hazard rate',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = c(0,0),
    checkMax = NA)

  testParameterErrors(baseCase, 
    'dropout_prior_distribution', 
    'Prior distribution for the patient dropout hazard rate',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = c(0,0),
    checkMax = NA)

  testParameterErrors(baseCase, 
    'enrollment_prior_distribution', 
    'Prior distribution for the enrollment hazard rate',
    checkMissing = TRUE,
    checkSize = TRUE,
    checkMin = c(0,0),
    checkMax = NA)

  testParameterErrors(baseCase, 
    'nsims', 
    'Number of simulations',
    checkMissing = FALSE,
    checkSize = TRUE,
    checkMin = 0,
    checkMax = 10001)

})

test_that("Input parameters errors check EventPredReport", {

  expect_error(
    EventPredReportDoc(""),
    info = "Checking for wrong parameter type for report generator"
  )

})
