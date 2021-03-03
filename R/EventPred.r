EventPred = function(parameters) {

  # Error checks

  if (typeof(parameters) != "list") stop("Function parameters must be a list of named values.", call. = FALSE)

  if (is.null(parameters$data_set)) stop("Data set with the patient enrollment, event and dropout information (data_set): Value must be specified.", call. = FALSE)

  data_set = parameters$data_set

  if (is.null(data_set$enrollment)) stop("Data set with the patient enrollment, event and dropout information (data_set): Enrollment variable must be specified.", call. = FALSE)

  if (is.null(data_set$time)) stop("Data set with the patient enrollment, event and dropout information (data_set): Time variable must be specified.", call. = FALSE)

  if (is.null(data_set$event)) stop("Data set with the patient enrollment, event and dropout information (data_set): Event variable must be specified.", call. = FALSE)

  if (is.null(data_set$dropout)) stop("Data set with the patient enrollment, event and dropout information (data_set): Dropout variable must be specified.", call. = FALSE)

  # Timing of the interim analysis 
  interim_analysis = max(data_set$enrollment + data_set$time)
  parameters$interim_analysis = interim_analysis

  if (is.null(parameters$time_points)) stop("Future time points for computing event predictions (time_points): Value must be specified.", call. = FALSE)

  time_points = 
    ContinuousErrorCheck(parameters$time_points, 
                         NA, 
                         lower_values = interim_analysis,
                         lower_values_sign = ">=",
                         upper_values = NA,
                         upper_values_sign = NA,
                         "Future time points for computing event predictions (time_points)",
                         c("Value"),
                         "double",
                         NA) 

  if (is.null(parameters$event_prior_distribution)) stop("Prior distribution for the event hazard rate (event_prior_distribution): Value must be specified.", call. = FALSE)

  event_prior_distribution = 
    ContinuousErrorCheck(parameters$event_prior_distribution, 
                         2, 
                         lower_values = c(0, 0),
                         lower_values_sign = c(">",">"),
                         upper_values = NA,
                         upper_values_sign = NA,
                         "Prior distribution for the event hazard rate (event_prior_distribution)",
                         c("Value"),
                         "double",
                         NA) 

  if (is.null(parameters$dropout_prior_distribution)) stop("Prior distribution for the patient dropout hazard rate (dropout_prior_distribution): Value must be specified.", call. = FALSE)

  dropout_prior_distribution = 
    ContinuousErrorCheck(parameters$dropout_prior_distribution, 
                         2, 
                         lower_values = c(0, 0),
                         lower_values_sign = c(">",">"),
                         upper_values = NA,
                         upper_values_sign = NA,
                         "Prior distribution for the patient dropout hazard rate (dropout_prior_distribution)",
                         c("Value"),
                         "double",
                         NA) 

  if (is.null(parameters$enrollment_prior_distribution)) stop("Prior distribution for the patient enrollment hazard rate (enrollment_prior_distribution): Value must be specified.", call. = FALSE)

  enrollment_prior_distribution = 
    ContinuousErrorCheck(parameters$enrollment_prior_distribution, 
                         2, 
                         lower_values = c(0, 0),
                         lower_values_sign = c(">",">"),
                         upper_values = NA,
                         upper_values_sign = NA,
                         "Prior distribution for the enrollment hazard rate (enrollment_prior_distribution)",
                         c("Value"),
                         "double",
                         NA) 

  if (!is.null(parameters$nsims)) {

    nsims = 
      ContinuousErrorCheck(parameters$nsims, 
                           1, 
                           lower_values = c(1),
                           lower_values_sign = c(">="),
                           upper_values = c(10000),
                           upper_values_sign = c("<="),
                           "Number of simulations (nsims)",
                           c("Value"),
                           "int",
                           NA) 

  } else {
    parameters$nsims = 1000
  }

  ###########################################################

  # Compute predictions 

  results = EventPredR(parameters)

  n_scenarios = length(time_points)

  # Original event data from the interim analysis data set
  time = parameters$data_set$time + parameters$data_set$enrollment
  event = parameters$data_set$event
  original = data.frame(time, event)
  original = original[order(time), ]
  original$cumsum = cumsum(original$event)

  # Predicted number of events
  predictions = matrix(0, n_scenarios, 3)

  for (i in 1:n_scenarios) {
    predictions[i, ] = c(quantile(results$number_events[[i]], probs = c(0.025, 1 - 0.025)), mean(results$number_events[[i]]))
  }

  results = list(parameters = parameters,
                 interim_analysis = original,
                 predictions = predictions)

  class(results) = "EventPredResults"

  return(results)


}

EventPredPriorDistribution = function(expected,
                                    uncertainty) {

  # Error checks

  expected = 
        ContinuousErrorCheck(expected, 
                             1, 
                             lower_values = 0,
                             lower_values_sign = ">",
                             upper_values = NA,
                             upper_values_sign = NA,
                             "Expected parameter value (expected)",
                             c("Value"),
                             "double",
                             NA) 

  uncertainty = 
        ContinuousErrorCheck(uncertainty, 
                             1, 
                             lower_values = 0,
                             lower_values_sign = ">",
                             upper_values = 1,
                             upper_values_sign = "<",
                             "Uncertainty parameter (uncertainty)",
                             c("Value"),
                             "double",
                             NA) 

  alpha = 1 / uncertainty^2
  beta = alpha / expected

  # Alpha and beta of the prior gamma distribution
  res = c(alpha, beta)

  return(res)

}  

PosteriorDistribution = function(data_set,
                                 interim_analysis,
                                 event_hazard_rate,
                                 dropout_hazard_rate,
                                 enrollment_hazard_rate) {
  
  n_events = sum(data_set$event) 
  n_dropouts = sum(data_set$dropout)
  n_patients = nrow(data_set)
  time = sum(data_set$time)

  # Alpha and beta of the posterior gamma distributions for the event and patient dropout hazard rates and intensity rate for the patient enrollment process
  post_event_rate = c(event_hazard_rate[1] + n_events,
                      event_hazard_rate[2] + time)
  post_dropout_rate = c(dropout_hazard_rate[1] + n_dropouts,
                       dropout_hazard_rate[2] + time)
  enrollment_rate = c(enrollment_hazard_rate[1] + n_patients,
                      enrollment_hazard_rate[2] + interim_analysis)

  res = list(event = post_event_rate, 
             dropout = post_dropout_rate,
             enrollment = enrollment_rate)

  return(res)

}

EventPredR = function(parameters) {

  event_prior_distribution = parameters$event_prior_distribution
  dropout_prior_distribution = parameters$dropout_prior_distribution
  enrollment_prior_distribution = parameters$enrollment_prior_distribution
  data_set = parameters$data_set
  interim_analysis = parameters$interim_analysis
  time_points = parameters$time_points
  nsims = parameters$nsims

  n_scenarios = length(time_points)

  # Compute the posterior distributions of the patient enrollment, event and dropout hazard rates 
  posterior_distributions = PosteriorDistribution(data_set,
                                      interim_analysis, 
                                      event_prior_distribution,
                                      dropout_prior_distribution,
                                      enrollment_prior_distribution) 

  # Sample from the posterior distributions 
  event_sample = rgamma(n = nsims, 
                        shape = posterior_distributions$event[1], 
                        rate = posterior_distributions$event[2])
  dropout_sample = rgamma(n = nsims, 
                         shape = posterior_distributions$dropout[1], 
                         rate = posterior_distributions$dropout[2])
  enrollment_sample = rgamma(n = nsims, 
                         shape = posterior_distributions$enrollment[1], 
                         rate = posterior_distributions$enrollment[2])

  number_events = list()

  current_parameters = list()

  current_parameters$event = data_set$event
  current_parameters$dropout = data_set$dropout
  current_parameters$enrollment = data_set$enrollment
  current_parameters$time = data_set$time
  current_parameters$event_sample = event_sample
  current_parameters$dropout_sample = dropout_sample
  current_parameters$enrollment_sample = enrollment_sample
  current_parameters$interim_analysis = interim_analysis

  # Loop over the number of scenarios
  for (k in 1:n_scenarios) {

    # Create a complete data set using posterior parameters and compute the number of events at the current time point
    current_parameters$time_point = time_points[k]
    number_events[[k]] = EventPredEventCount(current_parameters)$sim_results

  }

  # Save the parameters and results

  results = list(number_events = number_events,
                 posterior_distributions = posterior_distributions)

  return(results)

}

EventPredReportDoc = function(results) {

   #############################################################################

   # Error checks

   if (class(results) != "EventPredResults") stop("The object was not created by the EventPred function", call. = FALSE)

  #############################################################################

  # Empty list of tables to be included in the report
  item_list = list()
  item_index = 1
  table_index = 1
  figure_index = 1

  width = 6.5
  height = 5
  pointsize = 12

  parameters = results$parameters
  interim_analysis = results$interim_analysis
  prediction = results$prediction

  #############################################################################

  report_title = "Event prediction in event-driven trials"

  item_list[[item_index]] = list(type = "paragraph", label = "Description", value = "The simulation report presents a summary of event predictions for Phase II or Phase III trials with an event-driven design. Blinded event data at an interim analysis are used to forecast the number of events at pre-defined time points in the future.")
  item_index = item_index + 1

  ##########################################

  # Trial parameters

  column_names = c("Parameter ", "Value")

  col1 = c("Timing of the interim analysis", "Future time points for computing event predictions")
  col2 = c(parameters$interim_analysis, paste0(parameters$time_points, collapse = ", "))

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Trial parameters")

  column_width = c(3, 3.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1  
  table_index = table_index + 1

  ##########################################

  # Prior distributions

  column_names = c("Parameter ", "Value")

  col1 = c("Event hazard rate",
           "Patient dropout hazard rate",
           "Patient enrollment rate")
  col2 = c(paste0("Alpha = ", round(parameters$event_prior_distribution[1], 2), ", Beta = ", round(parameters$event_prior_distribution[2], 2)),
           paste0("Alpha = ", round(parameters$dropout_prior_distribution[1], 2), ", Beta = ", round(parameters$dropout_prior_distribution[2], 2)),
           paste0("Alpha = ", round(parameters$enrollment_prior_distribution[1], 2), ", Beta = ", round(parameters$enrollment_prior_distribution[2], 2)))

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Prior distribution parameters")

  footnote = "Shape (alpha) and rate (beta) parameters of the prior gamma distributions for the event and patient dropout hazard rates as well as the intensity rate of the patient enrollment process."

  column_width = c(3, 3.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1  
  table_index = table_index + 1

  ##########################################

  # Other parameters

  column_names = c("Parameter ", "Value")

  col1 = c("Number of simulations")
  col2 = c(parameters$nsims)

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Simulation parameters")

  column_width = c(3, 3.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
  item_index = item_index + 1
  table_index = table_index + 1

  ##########################################

  # Event prediction

  column_names = c("Time point", "Mean number of events", "95% predictive interval")

  time_points = parameters$time_points

  n_scenarios = length(time_points)

  col1 = time_points
  col2 = round(prediction[, 3], 1)
  col3 = rep("", 3)

  for (i in 1:n_scenarios) col3[i] = paste0("(", round(prediction[i, 1], 1), ", ", round(prediction[i, 2], 2), ")")

  data_frame = cbind(col1, col2, col3)

  footnote = paste0("Number of events at the interim analysis: ", max(interim_analysis$cumsum), ".")

  title = paste0("Table ", table_index, ". Event prediction at pre-defined time points")

  column_width = c(1.5, rep(5 / 2, 2))
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE, footnote)
  item_index = item_index + 1  
  table_index = table_index + 1

  ##########################################

  filename = "prediction.emf" 

  xlabels = seq(from = 0, to = max(time_points) + 10, by = 10)
  ylabels = seq(from = 0, to = max(prediction) + 10, by = 10)

  emf(file = filename, width = width, height = height, pointsize = pointsize)

  plot(x = interim_analysis$time, y = interim_analysis$cumsum, xlab="",ylab="", xlim=c(min(xlabels), max(xlabels)), ylim = c(min(ylabels), max(ylabels)), col="black", lwd = 2, type="l", axes=FALSE)
  polygon(c(rev(time_points), time_points), c(rev(prediction[, 2]), prediction[, 1]), col = "grey80", border = NA)
  lines(x = time_points, y = prediction[, 3], col="red", lwd = 2)
  axis(1, at = xlabels, labels = as.character(xlabels), tck=0.02, mgp=c(0, 0.4, 0))
  axis(2, at = ylabels, labels = as.character(ylabels), tck=0.02, mgp=c(0, 0.4, 0))
  mtext("Number of events", side=2, line=1.5)
  mtext("Time", side=1, line=1.5)
  box() 
  dev.off()

  footnote = "Black curve: Observed events. Red curve: Predicted mean number of events. Gray band: 95% predictive interval."

  item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Event prediction at pre-defined time points"), 
                           filename = filename,
                           dim = c(width, height),
                           type = "emf_plot",
                           footnote = footnote,
                           page_break = FALSE)

  ##########################################

  report = item_list

  doc = SaveReport(report, report_title)

  return(doc)

}
# End of EventPredReportDoc

