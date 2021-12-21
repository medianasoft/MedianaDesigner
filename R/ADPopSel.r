ADPopSel = function(parameters) {

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

  parameters$endpoint_index = 1  

  for (i in 1:length(endpoint_list)) {
      if (tolower(endpoint_list[i]) == tolower(parameters$endpoint_type)) parameters$endpoint_index = i
  }   

  endpoint_index = parameters$endpoint_index

  if (is.null(parameters$direction)) {
        parameters$direction_index = 1
  } else {
        if (!tolower(parameters$direction) %in% c("higher", "lower")) stop("Direction of favorable outcome (direction): Value must be specified.", call. = FALSE)
  }

  if (tolower(parameters$direction) == "higher") parameters$direction_index = 1    
  if (tolower(parameters$direction) == "lower") parameters$direction_index = 2    

  if (is.null(parameters$sample_size)) stop("Number of enrolled patients (sample_size): Value must be specified.", call. = FALSE)

  sample_size = ContinuousErrorCheck(parameters$sample_size, 
                                     2, 
                                     lower_values = 0,
                                     lower_values_sign = ">",
                                     upper_values = 1000,
                                     upper_values_sign = "<=",
                                     "Number of enrolled patients (sample_size)",
                                     NA,
                                     "int",
                                     NA) 

  narms = length(sample_size)

  if (is.null(parameters$prevalence)) stop("Prevalence of biomarker-positive patients (prevalence): Value must be specified.", call. = FALSE)

  prevalence = 
        ContinuousErrorCheck(parameters$prevalence, 
                             1, 
                             lower_values = c(0),
                             lower_values_sign = c(">"),
                             upper_values = c(0.999),
                             upper_values_sign = c("<"),
                             "Prevalence of biomarker-positive patients (prevalence)",
                             c("Value"),
                             "double",
                             NA) 


  if (is.null(parameters$info_frac)) stop("Information fractions at IA1, IA2, FA (info_frac): Value must be specified.", call. = FALSE)

  info_frac = ContinuousErrorCheck(parameters$info_frac, 
                                     3, 
                                     lower_values = 0,
                                     lower_values_sign = ">",
                                     upper_values = 1,
                                     upper_values_sign = "<=",
                                     "Information fractions at IA1, IA2, FA (info_frac)",
                                     NA,
                                     "double",
                                     NA) 

  if (info_frac[1] >= info_frac[2]) stop("Information fractions at IA1, IA2, FA (info_frac): Information fraction at IA1 must be less than Information fraction at IA2.", call. = FALSE)

  if (info_frac[2] >= info_frac[3]) stop("Information fractions at IA1, IA2, FA (info_frac): Information fraction at IA2 must be less than Information fraction at FA", call. = FALSE)

  if (info_frac[3] != 1) stop("Information fractions at IA1, IA2, FA (info_frac): Information fraction at FA must be equal to 1.", call. = FALSE)
  
  if (is.null(parameters$futility_threshold)) stop("Futility threshold at IA1 (futility_threshold): Value must be specified.", call. = FALSE)

  futility_threshold = 
        ContinuousErrorCheck(parameters$futility_threshold, 
                             1, 
                             lower_values = c(0),
                             lower_values_sign = c(">="),
                             upper_values = c(0.999),
                             upper_values_sign = c("<"),
                             "Futility threshold at IA1 (futility_threshold)",
                             c("Value"),
                             "double",
                             NA) 

  if (is.null(parameters$influence)) stop("Influence threshold at IA2 (influence): Value must be specified.", call. = FALSE)

  influence = 
        ContinuousErrorCheck(parameters$influence, 
                             1, 
                             lower_values = c(0),
                             lower_values_sign = c(">"),
                             upper_values = c(0.999),
                             upper_values_sign = c("<"),
                             "Influence threshold at IA2 (influence)",
                             c("Value"),
                             "double",
                             NA) 

  if (is.null(parameters$interaction)) stop("Interaction threshold at IA2 (interaction): Value must be specified.", call. = FALSE)

  interaction = 
        ContinuousErrorCheck(parameters$interaction, 
                             1, 
                             lower_values = c(1),
                             lower_values_sign = c(">"),
                             upper_values = c(2),
                             upper_values_sign = c("<"),
                             "Interaction threshold at IA2 (interaction)",
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

    if (is.null(parameters$control_mean)) stop("Mean effects in the control arm (control_mean): Value must be specified.", call. = FALSE)

    control_mean = 
          ContinuousErrorCheck(parameters$control_mean, 
                               2, 
                               lower_values = c(NA),
                               lower_values_sign = c(NA),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Mean effects in the control arm (control_mean)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_mean)) stop("Mean effects in the treatment arm (treatment_mean): Value must be specified.", call. = FALSE)

    treatment_mean = 
          ContinuousErrorCheck(parameters$treatment_mean, 
                               2, 
                               lower_values = c(NA),
                               lower_values_sign = c(NA),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Mean effects in the treatment arm (treatment_mean)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$control_sd)) stop("Standard deviations in the control arm (control_sd): Value must be specified.", call. = FALSE)

    control_sd = 
          ContinuousErrorCheck(parameters$control_sd, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Standard deviations in the control arm (control_sd)",
                               c("Value"),
                               "double",
                               NA) 


    if (is.null(parameters$treatment_sd)) stop("Standard deviations in the treatment arm (treatment_sd): Value must be specified.", call. = FALSE)

    treatment_sd = 
          ContinuousErrorCheck(parameters$treatment_sd, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Standard deviations in the treatment arm (treatment_sd)",
                               c("Value"),
                               "double",
                               NA) 

  }

  if (endpoint_index == 2) {

    if (is.null(parameters$control_rate)) stop("Response rates in the control arm (control_rate): Value must be specified.", call. = FALSE)

    control_rate = 
          ContinuousErrorCheck(parameters$control_rate, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(1),
                               upper_values_sign = c("<"),
                               "Response rates in the control arm (control_rate)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_rate)) stop("Response rates in the treatment arm (treatment_rate): Value must be specified.", call. = FALSE)

    treatment_rate = 
          ContinuousErrorCheck(parameters$treatment_rate, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(1),
                               upper_values_sign = c("<"),
                               "Responses rate in the treatment arm (treatment_rate)",
                               c("Value"),
                               "double",
                               NA) 

  }

  if (endpoint_index == 3) {

    if (is.null(parameters$control_time)) stop("Median times in the control arm (control_time): Value must be specified.", call. = FALSE)

    control_time = 
          ContinuousErrorCheck(parameters$control_time, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Median times in the control arm (control_time)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_time)) stop("Median times in the treatment arm (treatment_time): Value must be specified.", call. = FALSE)

    treatment_time = 
          ContinuousErrorCheck(parameters$treatment_time, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Median times in the treatment arm (treatment_time)",
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

    if (is.null(parameters$event_count)) stop("Target number of events at FA (event_count): Value must be specified.", call. = FALSE)

    event_count = 
          ContinuousErrorCheck(parameters$event_count, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(sum(parameters$sample_size)),
                               upper_values_sign = c("<"),
                               "Target number of events at FA (event_count)",
                               c("Value"),
                               "int",
                               NA)   

  }

  #############################################

  parameters$means = c(0, 0)
  parameters$sds = c(0, 0)
  parameters$rates = c(0, 0)
  parameters$hazard_rates = c(0, 0)
  parameters$dropout_parameter = 0
  parameters$enrollment_distribution = 2
  parameters$sample_size_ia1 = 0
  parameters$sample_size_ia2 = 0
  parameters$sample_size_fa = 1
  parameters$event_count_ia1 = 0
  parameters$event_count_ia2 = 0
  parameters$event_count_fa = 0
  parameters$max_sample_size = 0

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
    parameters$dropout_parameter = c(parameters$dropout_rate, dropout_period = 12)
    parameters$event_count_ia1 = floor(parameters$event_count[1] * info_frac[1])
    parameters$event_count_ia2 = floor(parameters$event_count[1] * info_frac[2])
    parameters$event_count_fa = floor(parameters$event_count * info_frac[3])
  }

  # Sample size in each subpopulation (biomarker-negative patients in the control arm, biomarker-positive patients in the control arm, biomarker-negative patients in the treatment arm, biomarker-positive patients in the treatment arm)
  sample_size_pop = rep(0, 4)
  sample_size_pop[1] = round(sample_size[1] * (1 - prevalence))
  sample_size_pop[2] = sample_size[1] - sample_size_pop[1]
  sample_size_pop[3] = round(sample_size[2] * (1 - prevalence))
  sample_size_pop[4] = sample_size[2] - sample_size_pop[3]

  parameters$sample_size_pop = sample_size_pop

  # Total sample size after accounting for dropout rates
  if (endpoint_index != 3) {
    parameters$sample_size_ia1 = floor(parameters$sample_size_pop * (1 - parameters$dropout_rate) * info_frac[1])
    parameters$sample_size_ia2 = floor(parameters$sample_size_pop * (1 - parameters$dropout_rate) * info_frac[2])
    parameters$sample_size_fa = floor(parameters$sample_size_pop * (1 - parameters$dropout_rate) * info_frac[3])
    parameters$max_sample_size = max(parameters$sample_size_fa)
  }



  ###########################################################

  # Run simulations on multiple cores to compute key characteristics

  if (parameters$ncores > 1) {
    # nocov start
    cl = parallel::makeCluster(parameters$ncores)

    # Export all functions in the global environment to each node
    parallel::clusterExport(cl, ls(envir = .GlobalEnv))

    doParallel::registerDoParallel(cl)
    simulation_list = foreach(counter=(1:parameters$ncores), .packages = c("MedianaDesigner")) %dorng% { 
      ADPopSelSingleCore(parameters)
    }
    stopCluster(cl)

    # Combine the simulation results across the cores 

    sim_results = NULL

    for (i in 1:parameters$ncores) {

      sim_results = rbind(sim_results, simulation_list[[i]]$sim_results)

    }
    # nocov end
  } else {

    simulations = ADPopSelSingleCore(parameters)
    sim_results = simulations$sim_results

  }

  # Add column names
  column_names = c("trad_sign_outcome_op", "trad_sign_outcome_bplus", "adapt_sign_outcome_op", "adapt_sign_outcome_bplus", "ia2_cp", "futility_flag", "effect_size_bminus", "effect_size_bplus", "selection_flag_op_only", "selection_flag_bplus_only", "selection_flag_both", "ia1_time", "ia2_time", "fa_time_op_only", "fa_time_bplus_only") 
  colnames(sim_results) = column_names

  sim_summary = list()

  nsims = nrow(sim_results)

  sim_summary$trad_power = colMeans(sim_results[, 1:2])
  ad_power = rep(0, 3)
  ad_power[1:2] = colMeans(sim_results[, 3:4])
  ad_power[3] = mean(sim_results[, 3] + sim_results[, 4] >= 1)  
  sim_summary$ad_power = ad_power
  sim_summary$futility = mean(sim_results[, 6])

  sim_summary$hypothesis_selection = colMeans(sim_results[, 9:11])
  look_time = rep(0, 4)
  look_time[1:2] = colMeans(sim_results[, 12:13])
  look_time[3] = mean(sim_results[sim_results[, 14] != 0, 14])
  look_time[4] = mean(sim_results[sim_results[, 15] != 0, 15])
  sim_summary$look_time = look_time

  results = list(parameters = parameters,
                 sim_results = sim_results,
                 sim_summary = sim_summary)

  class(results) = "ADPopSelResults"

  return(results)
}    

ADPopSelReportDoc = function(results) {

   #############################################################################

   # Error checks

   if (class(results) != "ADPopSelResults") stop("The object was not created by the ADPopSel function", call. = FALSE)

  #############################################################################

  statistics = c("Lower quartile", "Median", "Mean", "Upper quartile")

  #############################################################################

  # Empty list of tables to be included in the report
  item_list = list()
  item_index = 1
  table_index = 1
  figure_index = 1

  width = 6.5
  height = 5

  parameters = results$parameters
  sim_results = results$sim_results
  sim_summary = results$sim_summary
  endpoint_index = parameters$endpoint_index
  narms = 2

  # Trial arms  
  trial_arms = c("Control", "Treatment")

  # Populations
  populations = c("Biomarker-negative", "Biomarker-positive")

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

  report_title = "Adaptive designs with data-driven population selection"

  item_list[[item_index]] = list(type = "paragraph", label = "Description", value = "The simulation report presents key operating characteristics of an adaptive design for a two-arm Phase III clinical trial with two interim analyses and two pre-defined patient populations (overall population and subpopulation of biomarker-positive patients). The first interim analysis supports early stopping for futility in the overall population and the second interim analysis identifies the most promising patient population or populations.")

  item_index = item_index + 1

  #############################################################################

  column_names = c("Trial arm", "Population", "Sample size")

  col1 = c(trial_arms[1], "", trial_arms[2], "")
  col2 = c(populations, populations)
  col3 = parameters$sample_size_pop

  data_frame = data.frame(col1, col2, col3)
  title = paste0("Table ", table_index, ". Number of enrolled patients")

  column_width = c(2, 3, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1

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

  column_names = c("Trial arm", "Population", "Parameter", "Value")

  if (endpoint_index == 1) {
    col1 = c(trial_arms[1], rep("", 3), trial_arms[2], rep("", 3))
    col2 = c(populations[1], "", populations[2], "", populations[1], "", populations[2], "")
    col3 = rep(c("Mean", "SD"), 2)
    col4 = c(means[1], sds[1], means[2], sds[2], means[3], sds[3], means[4], sds[4])
  }

  if (endpoint_index == 2) {
    col1 = c(trial_arms[1], "", trial_arms[2], "")
    col2 = c(populations, populations)
    col3 = rep("Rate (%)", 4)
    col4 = 100 * rates
  }

  if (endpoint_index == 3) {
    col1 = c(trial_arms[1], "", trial_arms[2], "")
    col2 = c(populations, populations)
    col3 = rep("Median time", 4)
    col4 = times
  }

  data_frame = data.frame(col1, col2, col3, col4)
  title = paste0("Table ", table_index, ". Treatment effect assumptions")

  column_width = c(1, 2.5, rep(3 / 2, 2))
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  if (endpoint_index %in% c(1, 2)) column_names = c("Decision point", "Total number of patients", "Information fraction (%)")

  if (endpoint_index %in% c(3)) column_names = c("Decision point", "Total number of events", "Information fraction (%)")

  if (endpoint_index %in% c(1, 2)) {
    col1 = c("Interim analysis 1", 
           "Interim analysis 2",
           "Final analysis")
    col2 = c(sum(parameters$sample_size_ia1), sum(parameters$sample_size_ia2), sum(parameters$sample_size_fa))
    col3 = round(100 * parameters$info_frac[1:3], 1)    
  }
  if (endpoint_index %in% c(3)) {
    col1 = c("Interim analysis 1", 
           "Interim analysis 2",
           "Final analysis (OP)",
           "Final analysis (BPP)")
    col2 = c(parameters$event_count_ia1, parameters$event_count_ia2, parameters$event_count_fa)
    col3 = round(c(100 * parameters$info_frac[1:3], 100), 1)    
  }

  data_frame = data.frame(col1, col2, col3)
  if (endpoint_index %in% c(1, 2)) title = paste0("Table ", table_index, ". Number of patients at the interim and final analyses")
  if (endpoint_index %in% c(3)) title = paste0("Table ", table_index, ". Number of events at the interim and final analyses")

  if (endpoint_index %in% c(1, 2)) footnote = "The number of patients at the interim and final analyses may be reduced due to patient dropout." else footnote = "OP: Overall population. BPP: Biomarker-positive population."

  column_width = c(2.5, 2, 2)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Parameter", "PPS (%)")

  col1 = c("Futility threshold")
  col2 = 100 * c(parameters$futility_threshold)
  
  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Decision rule at Interim analysis 1")

  footnote = "PPS: Predicted probability of success at Interim analysis 1. The trial will be stopped for futility at Interim analysis 1 if the predicted probability of success in the overall population is less than the futility threshold."

  column_width = c(5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Parameter", "Value")

  col1 = c("Influence threshold", "Interaction threshold")
  col2 = c(parameters$influence, parameters$interaction)
  
  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Decision rules  at Interim analysis 2")

  footnote = paste0("Influence condition: Only OP is selected if the effect size in BNP is less than the influence threshold. Interaction condition: Both OP and BPP are selected if the ratio of the effect sizes in BNP and BPP is greater than the interaction threshold and the influence condition is met. OP: Overall population. BPP: Biomarker-positive population. BNP: Biomarker-negative population")

  column_width = c(5, 1.5)
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

  column_names = c("Outcome", "Probability (%)")

  col1 = c("Futility stopping at Interim analysis 1",
           "Only OP is selected at Interim analysis 2",
           "Only BPP is selected at Interim analysis 2",
           "Both populations are selected at Interim analysis 2")
  col2 = round(100 * c(sim_summary$futility, sim_summary$hypothesis_selection), 1)

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Simulation results: Futility stopping and population selection")

  footnote = "OP: Overall population. BPP: Biomarker-positive population."

  column_width = c(4.5, 2)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Population", "Power (%)")

  col1 = c("Overall population")
  col2 = round(100 * sim_summary$trad_power[1], 1)

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Simulation results: Traditional design")

  footnote = "Power for a traditional design conducted in the overall population."

  column_width = c(4.5, 2)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Population", "Power (%)")

  col1 = c("Overall population", "Biomarker-positive population", "Either population")
  col2 = round(100 * sim_summary$ad_power, 1)

  footnote = "Probability of establishing a significant treatment effect in the overall population, biomarker-positive population or at least one of the two patient populations at the final analysis."

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Simulation results: Adaptive design")

  column_width = c(4.5, 2)

  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  report = item_list

  doc = SaveReport(report, report_title)

  return(doc)

}
# End of ADPopSelReportDoc

ADPopSelSingleCore = function(parameters) {

  params_for_run = parameters
  params_for_run$nsims = params_for_run$nsims_per_core
  simulations = ADPopSelC(params_for_run)

  return(simulations)

}
