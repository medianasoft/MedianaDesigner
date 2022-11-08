FutRule = function(parameters) {

  # Error checks

  if (typeof(parameters) != "list") stop("Function parameters must be a list of named values.", call. = FALSE)

  if (is.null(parameters$random_seed)) {
    
    random_seed = 49283

  } else {

    random_seed = ContinuousErrorCheck(parameters$random_seed, 
                                 1, 
                                 lower_values = 1,
                                 lower_values_sign = c(">="),
                                 upper_values = 100000,
                                 upper_values_sign = c("<="),
                                 "Seed for the random number generator (random_seed)",
                                 c("Value"),
                                 "int",
                                 NA) 

  }

  parameters$random_seed = random_seed

  # Set the seed of R's random number generator.
  # It also takes effect to Rcpp random generation functions.
  # https://stackoverflow.com/questions/60119621/get-the-same-sample-of-integers-from-rcpp-as-base-r
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  set.seed(random_seed)

  if (is.null(parameters$endpoint_type)) stop("Endpoint type (endpoint_type): Value must be specified.", call. = FALSE)

  if (!tolower(parameters$endpoint_type) %in% tolower(endpoint_list)) stop("Endpoint type (endpoint_type): Value must be Normal, Binary or Time-to-event.", call. = FALSE)

  if (is.null(parameters$direction)) {
        parameters$direction_index = 1
  } else {
        if (!tolower(parameters$direction) %in% c("higher", "lower")) stop("Direction of favorable outcome (direction): Value must be specified.", call. = FALSE)
  }

  if (tolower(parameters$direction) == "higher") parameters$direction_index = 1    
  if (tolower(parameters$direction) == "lower") parameters$direction_index = 2    

  parameters$endpoint_index = 1  

  for (i in 1:length(endpoint_list)) {
      if (tolower(endpoint_list[i]) == tolower(parameters$endpoint_type)) parameters$endpoint_index = i
  }   

  endpoint_index = parameters$endpoint_index

  if (is.null(parameters$sample_size)) stop("Number of enrolled patients (sample_size): Value must be specified.", call. = FALSE)

  sample_size = ContinuousErrorCheck(parameters$sample_size, 
                                     NA, 
                                     lower_values = 0,
                                     lower_values_sign = ">",
                                     upper_values = 1000,
                                     upper_values_sign = "<=",
                                     "Number of enrolled patients (sample_size)",
                                     NA,
                                     "int",
                                     NA) 


  narms = length(sample_size)

  if (is.null(parameters$info_frac)) stop("Information fraction (info_frac): Value must be specified.", call. = FALSE)

  info_frac = 
        ContinuousErrorCheck(parameters$info_frac, 
                             1, 
                             lower_values = c(0.001),
                             lower_values_sign = c(">"),
                             upper_values = c(0.999),
                             upper_values_sign = c("<"),
                             "Information fraction (info_frac)",
                             c("Value"),
                             "double",
                             NA) 

  if (!is.null(parameters$dropout_rate)) {

    dropout_rate = 
      ContinuousErrorCheck(parameters$dropout_rate, 
                           1, 
                           lower_values = c(0),
                           lower_values_sign = c(">="),
                           upper_values = c(1),
                           upper_values_sign = c("<"),
                           "Patient dropout rate (dropout_rate)",
                           c("Value"),
                           "double",
                           NA) 

  } else {
    parameters$dropout_rate = 0
  }

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

  if (!is.null(parameters$ncores)) {
    # nocov start
    # Maximum number of cores
    max_ncores = parallel::detectCores()

    ncores = 
      ContinuousErrorCheck(parameters$ncores, 
                           1, 
                           lower_values = c(1),
                           lower_values_sign = c(">="),
                           upper_values = c(max_ncores),
                           upper_values_sign = c("<="),
                           "Number of cores for parallel calculations (ncores)",
                           c("Value"),
                           "int",
                           NA) 
    # nocov end
  } else {
    parameters$ncores = 1
  }

  # Number of simulations per core
  parameters$nsims_per_core = ceiling(parameters$nsims / parameters$ncores)   

  if (!is.null(parameters$alpha)) {

    alpha = 
          ContinuousErrorCheck(parameters$alpha, 
                               1, 
                               lower_values = c(0.001),
                               lower_values_sign = c(">"),
                               upper_values = c(0.5),
                               upper_values_sign = c("<"),
                               "One-sided Type I error rate (alpha)",
                               c("Value"),
                               "double",
                               NA) 
  } else {
    parameters$alpha = 0.025
  }

  # Treatment effect assumptions 

  if (endpoint_index == 1) {

    if (is.null(parameters$control_mean)) stop("Mean effect in the control arm (control_mean): Value must be specified.", call. = FALSE)

    control_mean = 
          ContinuousErrorCheck(parameters$control_mean, 
                               1, 
                               lower_values = c(NA),
                               lower_values_sign = c(NA),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Mean effect in the control arm (control_mean)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_mean)) stop("Mean effects in the treatment arms (treatment_mean): Value must be specified.", call. = FALSE)

    treatment_mean = 
          ContinuousErrorCheck(parameters$treatment_mean, 
                               narms - 1, 
                               lower_values = c(NA),
                               lower_values_sign = c(NA),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Mean effects in the treatment arms (treatment_mean)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$control_sd)) stop("Standard deviation in the control arm (control_sd): Value must be specified.", call. = FALSE)

    control_sd = 
          ContinuousErrorCheck(parameters$control_sd, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Standard deviation in the control arm (control_sd)",
                               c("Value"),
                               "double",
                               NA) 


    if (is.null(parameters$treatment_sd)) stop("Standard deviations in the treatment arms (treatment_sd): Value must be specified.", call. = FALSE)

    treatment_sd = 
          ContinuousErrorCheck(parameters$treatment_sd, 
                               narms - 1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Standard deviations in the treatment arms (treatment_sd)",
                               c("Value"),
                               "double",
                               NA) 

  }

  if (endpoint_index == 2) {

    if (is.null(parameters$control_rate)) stop("Response rate in the control arm (control_rate): Value must be specified.", call. = FALSE)

    control_rate = 
          ContinuousErrorCheck(parameters$control_rate, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(1),
                               upper_values_sign = c("<"),
                               "Response rate in the control arm (control_rate)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_rate)) stop("Response rates in the treatment arms (treatment_rate): Value must be specified.", call. = FALSE)

    treatment_rate = 
          ContinuousErrorCheck(parameters$treatment_rate, 
                               narms - 1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(1),
                               upper_values_sign = c("<"),
                               "Response rates in the treatment arms (treatment_rate)",
                               c("Value"),
                               "double",
                               NA) 

  }

  if (endpoint_index == 3) {

    if (is.null(parameters$control_time)) stop("Median time in the control arm (control_time): Value must be specified.", call. = FALSE)

    control_time = 
          ContinuousErrorCheck(parameters$control_time, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Median time in the control arm (control_time)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_time)) stop("Median times in the treatment arms (treatment_time): Value must be specified.", call. = FALSE)

    treatment_time = 
          ContinuousErrorCheck(parameters$treatment_time, 
                               narms - 1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Median times in the treatment arms (treatment_time)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$enrollment_period)) stop("Patient enrollment period (enrollment_period): Value must be specified.", call. = FALSE)

    enrollment_period = 
          ContinuousErrorCheck(parameters$enrollment_period, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Patient enrollment period (enrollment_period)",
                               c("Value"),
                               "double",
                               NA)           

    if (is.null(parameters$enrollment_parameter)) stop("Median enrollment time (enrollment_parameter): Value must be specified.", call. = FALSE)

    enrollment_parameter = 
          ContinuousErrorCheck(parameters$enrollment_parameter, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(parameters$enrollment_period),
                               upper_values_sign = c("<"),
                               "Median enrollment time (enrollment_parameter)",
                               c("Value"),
                               "double",
                               NA)   

    event_count = 
          ContinuousErrorCheck(parameters$event_count, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(sum(parameters$sample_size)),
                               upper_values_sign = c("<"),
                               "Target number of events (event_count)",
                               c("Value"),
                               "int",
                               NA)   

  }

  #############################################

  parameters$means = 0
  parameters$sds = 0
  parameters$rates = 0
  parameters$hazard_rates = 0
  parameters$dropout_parameter = 0
  parameters$enrollment_distribution = 2

  # All means and SDs
  if (endpoint_index == 1) {
    parameters$means = c(parameters$control_mean, parameters$treatment_mean)
    parameters$sds = c(parameters$control_sd, parameters$treatment_sd)
    parameters$enrollment_period = 0
    parameters$enrollment_parameter = 0

  }

  # All rates
  if (endpoint_index == 2) {
    parameters$rates = c(parameters$control_rate, parameters$treatment_rate)
    parameters$enrollment_period = 0
    parameters$enrollment_parameter = 0
  }

  # All hazard rates
  if (endpoint_index == 3) {
    parameters$times = c(parameters$control_time, parameters$treatment_time)
    parameters$hazard_rates = log(2) / c(parameters$control_time, parameters$treatment_time)
    parameters$interim = floor(parameters$info_frac * parameters$event_count)
    parameters$final = parameters$event_count
    parameters$dropout_parameter = c(parameters$dropout_rate, dropout_period = 12)

  }

  # Account for dropout rates
  if (endpoint_index != 3) {
    parameters$interim = floor(parameters$info_frac * parameters$sample_size * (1 - parameters$dropout_rate))
    parameters$final = floor(parameters$sample_size * (1 - parameters$dropout_rate))
  }

  ###########################################################

  cp_threshold = seq(from = 0, to = 1, by = 0.01)
  sensitivity = rep(0, length(cp_threshold))
  specificity = rep(0, length(cp_threshold))

  # Compute sensitivity

  simulations = FutRuleNCores(parameters)

  cp = simulations$sim_results  

  for (j in 1:length(cp_threshold)) sensitivity[j] = mean(rowMax(cp) > cp_threshold[j])    

  # Compute specificity

  # Null hypothesis of no effect
  if (endpoint_index == 1) parameters$means = rep(parameters$control_mean, narms)
  if (endpoint_index == 2) parameters$rates = rep(parameters$control_rate, narms)
  if (endpoint_index == 3) parameters$hazard_rates = log(2) / rep(parameters$control_time, narms)

  simulations = FutRuleNCores(parameters)

  cp = simulations$sim_results  

  for (j in 1:length(cp_threshold)) specificity[j] = mean(rowMax(cp) <= cp_threshold[j]) 

  accuracy = (sensitivity + specificity) / 2

  ###########################################################

  # Save the results

  sim_summary = list()

  sim_summary$sensitivity = sensitivity
  sim_summary$specificity = specificity
  sim_summary$accuracy = accuracy
  sim_summary$cp_threshold = cp_threshold

  results = list(parameters = parameters,
                 sim_summary = sim_summary)

  class(results) = "FutRuleResults"

  return(results)

}    

FutRuleReportDoc = function(results) {

   #############################################################################

   # Error checks

   if (!is(results, "FutRuleResults")) stop("The object was not created by the FutRule function.", call. = FALSE)

  #############################################################################

  # Empty list of tables to be included in the report
  item_list = list()
  item_index = 1
  table_index = 1
  figure_index = 1

  width = 6.5
  height = 6
  pointsize = 12

  parameters = results$parameters
  sensitivity = results$sim_summary$sensitivity
  specificity = results$sim_summary$specificity
  accuracy = results$sim_summary$accuracy
  cp_threshold = results$sim_summary$cp_threshold
  narms = length(parameters$sample_size)
  endpoint_index = parameters$endpoint_index

  #############################################################################

  report_title = "Optimal futility stopping rule"

  item_list[[item_index]] = list(type = "paragraph", label = "Description", value = "The simulation report presents operating characteristics of a multi-arm trial design with a single interim analysis. A futility stopping rule will be applied at this interim look and the trial will be stopped early for futility if the predicted probability of success (conditional power) is less than a pre-defined futility threshold in all treatment arms. An optimal value of the futility threshold is computed by maximizing the sensitivity and specificity rates.")
  item_index = item_index + 1

  #############################################################################

  column_names = c("Parameter", "Value")

  if (parameters$direction_index == 1) label = "A higher value of the endpoint indicates a more favorable outcome"

  if (parameters$direction_index == 2) label = "A lower value of the endpoint indicates a more favorable outcome"

  col1 = c("Endpoint type", "Direction of favorable outcome") 
  col2 = c(endpoint_list[endpoint_index], label)

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Primary efficacy endpoint")

  column_width = c(3, 3.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1
  
  #############################################################################

  # Trial arms  
  trial_arms = "Control"
  if (narms >= 3) {
    for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
  } else {
    trial_arms = c(trial_arms, "Treatment")
  }

  # All means and SDs
  if (endpoint_index == 1) {
    means = c(parameters$control_mean, parameters$treatment_mean)
    sds = c(parameters$control_sd, parameters$treatment_sd)
  }

  # All rates
  if (endpoint_index == 2) {
    rates = c(parameters$control_rate, parameters$treatment_rate)
  }

  # All hazard rates
  if (endpoint_index == 3) {
    times = c(parameters$control_time, parameters$treatment_time)
  }

  #############################################################################

  column_names = c("Trial arm", "Sample size")

  col1 = trial_arms
  col2 = parameters$sample_size

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Number of enrolled patients")

  column_width = c(5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Trial arm", "Parameter", "Value")

  col1 = NULL
  col2 = NULL
  col3 = NULL

  for (i in 1:narms) {

    if (endpoint_index == 1) {
      col1 = c(col1, trial_arms[i], "")
      col2 = c(col2, "Mean", "SD")
      col3 = c(col3, means[i], sds[i])
    }

    if (endpoint_index == 2) {
      col1 = c(col1, trial_arms[i])
      col2 = c(col2, "Rate (%)")
      col3 = c(col3, 100 * rates[i])
    }

    if (endpoint_index == 3) {
      col1 = c(col1, trial_arms[i])
      col2 = c(col2, "Median time")
      col3 = c(col3, times[i])
    }

  }

  data_frame = data.frame(col1, col2, col3)
  title = paste0("Table ", table_index, ". Treatment effect assumptions")

  column_width = c(2, 2, 2.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  if (endpoint_index %in% c(1, 2)) column_names = c("Decision point", "Total number of patients", "Information fraction (%)")

  if (endpoint_index %in% c(3)) column_names = c("Decision point", "Total number of events", "Information fraction (%)")

  col1 = c("Interim analysis", 
           "Final analysis")
  if (endpoint_index %in% c(1, 2)) col2 = c(sum(parameters$interim), sum(parameters$final))
  if (endpoint_index %in% c(3)) col2 = c(round(parameters$event_count * c(parameters$info_frac, 1)))
  col3 = c(100 * parameters$info_frac, 100)
  
  data_frame = data.frame(col1, col2, col3)
  if (endpoint_index %in% c(1, 2)) title = paste0("Table ", table_index, ". Number of patients at the interim and final analyses")
  if (endpoint_index %in% c(3)) title = paste0("Table ", table_index, ". Number of events at the interim and final analyses")

  if (endpoint_index %in% c(1, 2)) footnote = "The number of patients at the interim and final analyses may be reduced due to patient dropout." else footnote = NULL 

  column_width = c(2, 2.5, 2)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Parameter", "Value")

  if (endpoint_index %in% c(1, 2)) {
    col1 = c("Dropout rate at the end of the treatment period (%)")
    col2 = c(100 * parameters$dropout_rate)
  }
  if (endpoint_index %in% c(3)) {
    col1 = c("Patient enrollment period", "Median enrollment time", "Annual dropout rate (%)")
    col2 = c(parameters$enrollment_period, parameters$enrollment_parameter, 100 * parameters$dropout_rate)
  }

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Other design parameters")

  column_width = c(5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Parameter", "Value")

  col1 = c("One-sided Type I error rate", "Number of simulations")
  col2 = c(sprintf("%0.3f", sum(parameters$alpha)),
           sprintf("%d", parameters$nsims))

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Simulation parameters")

  column_width = c(5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  if (is.null(parameters$withoutCharts)) {   # skip chart generation on tests
    # nocov start
    filename = "sensspec.emf"

    emf(file = filename, width = width, height = height, pointsize = pointsize)
    plot(x = 100 * cp_threshold, y = 100 * sensitivity, xlab="",ylab="", xlim=c(0, 100), ylim=c(0, 100), axes=FALSE, col="red", lwd = 2, type="l")
    # 1=bottom, 2=left, 3=top, 4=right
    lines(x = 100 * cp_threshold, y = 100 * specificity, col="blue", lwd = 2, type="l")
    labels = seq(from = 0, to = 100, by = 20)
    axis(1,at=labels, labels=labels, tck = 0.02, mgp = c(0, 0.4, 0))
    axis(2,at=labels, labels=labels, tck = 0.02, mgp = c(0, 0.4, 0))
    mtext("Futility threshold (%)", side=1, line=1.5)
    mtext("Sensitivity/specificity rate (%)", side=2, line=1.75)
    box()
    dev.off()

    footnote = "Red curve: Sensitivity rate (probability of correctly retaining at least one treatment arm at the interim analysis, evaluated under the alternative hypothesis of beneficial effect, i.e., all treatments are effective). Blue curve: Specificity rate (probability of correctly stopping all treatment arms at the interim analysis due to futility, evaluated under the null hypothesis of no effect, i.e., all treatments are ineffective)."

    item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Sensitivity and specificity rates as functions of the futility threshold"), 
                              filename = filename,
                              dim = c(width, height),
                              type = "emf_plot",
                              page_break = TRUE,
                              footnote = footnote)

    item_index = item_index + 1
    figure_index = figure_index + 1
    # nocov end
  }

  #############################################################################

  if (is.null(parameters$withoutCharts)) {   # skip chart generation on tests
    # nocov start
    filename = "accuracy.emf"

    index = which.max(accuracy)
    optimal_point = round(100 * cp_threshold[index], 1)
    level = 0.95 * max(accuracy)
    zone = cp_threshold[accuracy >= level]
    optimal_lower = round(100 * min(zone), 1)
    optimal_upper = round(100 * max(zone), 1)

    emf(file = filename, width = width, height = height)
    plot(x = 100 * cp_threshold, y = 100 * accuracy, xlab="",ylab="", xlim=c(0, 100), ylim=c(50, 100), axes=FALSE, col="red", lwd = 2, type="l")
    box()
    # 1=bottom, 2=left, 3=top, 4=right
    xlabels = seq(from = 0, to = 100, by = 20)
    ylabels = seq(from = 50, to = 100, by = 10)
    axis(1,at=xlabels, labels=xlabels, tck = 0.02, mgp = c(0, 0.4, 0))
    axis(2,at=ylabels, labels=ylabels, tck = 0.02, mgp = c(0, 0.4, 0))
    abline(v = optimal_point, lty = "solid")
    abline(v = c(optimal_lower, optimal_upper), lty = "dashed")
    mtext("Futility threshold (%)", side=1, line=1.5)
    mtext("Accuracy rate (%)", side=2, line=1.75)
    box()
    dev.off()

    footnote = paste0("The accuracy rate is defined as the average of the sensitivity and specificity rates and an optimal futility threshold is defined as the threshold that maximizes the accuracy rate. Optimal futility threshold: ", optimal_point, "%. 95% optimal interval: (", optimal_lower, "%, ", optimal_upper, "%).")

    item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Accuracy rate as a function of the futility threshold"), 
                              filename = filename,
                              dim = c(width, height),
                              type = "emf_plot",
                              page_break = FALSE,
                              footnote = footnote)

    item_index = item_index + 1
    figure_index = figure_index + 1
    # nocov end
  }

  #############################################################################

  report = item_list

  doc = SaveReport(report, report_title)

  return(doc)

}
# End of FutRuleReportDoc

FutRuleNCores = function(parameters) {

  # Run simulations on multiple cores to compute key characteristics

  ncores = parameters$ncores

  if (ncores > 1) {
    # nocov start
    cl = parallel::makeCluster(ncores)

    # Export all functions in the global environment to each node
    parallel::clusterExport(cl, ls(envir = .GlobalEnv))

    doParallel::registerDoParallel(cl)
    simulation_list = foreach(counter=(1:ncores), .packages = c("MedianaDesigner")) %dorng% { 
      FutRuleSingleCore(parameters)
    }
    stopCluster(cl)

    # Combine the simulation results across the cores 

    sim_results = NULL

    for (i in 1:ncores) {

      sim_results = rbind(sim_results, simulation_list[[i]]$sim_results)

    }

    simulations = list(sim_results = sim_results)
    # nocov end
  } else {

    simulations = FutRuleSingleCore(parameters)

  }

  return(simulations)
}

FutRuleSingleCore = function(parameters) {
  params_for_core = parameters
  params_for_core$nsims = params_for_core$nsims_per_core
  simulations = FutRuleC(params_for_core)
  return(simulations)
}