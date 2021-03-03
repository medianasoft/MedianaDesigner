# Test the input of the EventPred function

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
  nsims = 1000
)

test_that("Success run EventPred", {

  # Success run
  results = EventPred(baseCase)
  expect_is(results, "EventPredResults")
  expect_named(results, c('parameters', 'interim_analysis', 'predictions'))

  expect_is(results$interim_analysis, "data.frame")
  expect_length(results$interim_analysis, 3)
  expect_is(results$predictions, "matrix")
  expect_length(results$predictions, 39)

  predictions = results$predictions
  predictionTarget = predictions[13,3]
  expect_true(
    abs(predictionTarget - 350) < 10,
    info = paste("predictionTarget(",predictionTarget,") is out of range 350±10"))

  # Check for report generation
  EventPredReport(results, tempfile(fileext = ".docx"))
})

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

  expect_error(
    EventPred(
      list(
        data_set = baseCase$data_set,
        # missing
        #time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    info = "Checking for missing future time points for computing event predictions"
  )

  expect_error(
    EventPred(
      list(
        data_set = baseCase$data_set,
        time_points = baseCase$time_points,
        # missing
        #event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    info = "Checking for missing prior distribution for the event hazard rate"
  )

  expect_error(
    EventPred(
      list(
        data_set = baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        # missing
        #dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    info = "Checking for missing prior distribution for the patient dropout hazard rate"
  )

  expect_error(
    EventPred(
      list(
        data_set = baseCase$data_set,
        time_points = baseCase$time_points,
        event_prior_distribution = baseCase$event_prior_distribution,
        dropout_prior_distribution = baseCase$dropout_prior_distribution,        
        # missing
        #enrollment_prior_distribution = baseCase$enrollment_prior_distribution,        
        nsims = baseCase$nsims
      )
    ),
    info = "Checking for missing prior distribution for the patient enrollment hazard rate"
  )

})
