MultAdj = function(parameters) {

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

  if (!tolower(parameters$endpoint_type) %in% tolower(c("normal", "binary"))) stop("Endpoint type (endpoint_type): Value must be Normal or Binary.", call. = FALSE)

  if (tolower(parameters$endpoint_type) == "normal") endpoint_index = 1  
  if (tolower(parameters$endpoint_type) == "binary") endpoint_index = 2  

  parameters$endpoint_index = endpoint_index

  # Default value of direction
  if (is.null(parameters$direction)) {

    parameters$direction = "Higher"
    parameters$direction_index = 1  

  } else {

    direction_list = c("Higher", "Lower")

    if (!tolower(parameters$direction) %in% tolower(direction_list)) stop("Direction of beneficial effect (direction): Value must be Higher or Lower.", call. = FALSE)

    for (i in 1:length(direction_list)) {
        if (tolower(direction_list[i]) == tolower(parameters$direction)) parameters$direction_index = i
    }   

  }

  if (is.null(parameters$n_comparisons)) stop("Number of dose-control comparisons (n_comparisons): Value must be specified.", call. = FALSE)

  n_comparisons = ContinuousErrorCheck(parameters$n_comparisons, 
                                     1, 
                                     lower_values = 1,
                                     lower_values_sign = ">=",
                                     upper_values = 10,
                                     upper_values_sign = "<=",
                                     "Number of dose-control comparisons (n_comparisons)",
                                     NA,
                                     "int",
                                     NA) 

  if (is.null(parameters$n_endpoints)) stop("Number of endpoints (n_endpoints): Value must be specified.", call. = FALSE)

  n_endpoints = ContinuousErrorCheck(parameters$n_endpoints, 
                                     1, 
                                     lower_values = 1,
                                     lower_values_sign = ">=",
                                     upper_values = 10,
                                     upper_values_sign = "<=",
                                     "Number of endpoints (n_endpoints)",
                                     NA,
                                     "int",
                                     NA) 

  if (n_comparisons == 1 & n_endpoints == 1) stop("Number of dose-control comparisons (n_comparisons) must be greater than 1 or Number of endpoints (n_endpoints) must be greater than 1.", call. = FALSE)

  if (is.null(parameters$sample_size)) stop("Number of enrolled patients (sample_size): Value must be specified.", call. = FALSE)

  sample_size = ContinuousErrorCheck(parameters$sample_size, 
                                     n_comparisons + 1, 
                                     lower_values = 0,
                                     lower_values_sign = ">",
                                     upper_values = 10000,
                                     upper_values_sign = "<=",
                                     "Number of enrolled patients (sample_size)",
                                     NA,
                                     "int",
                                     NA) 

  if (is.null(parameters$dropout_rate)) {
    parameters$dropout_rate = 0
  }

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

  if (is.null(parameters$nsims)) {
    parameters$nsims = 1000
  }

  nsims = 
    ContinuousErrorCheck(parameters$nsims, 
                          1, 
                          lower_values = c(10),
                          lower_values_sign = c(">="),
                          upper_values = c(10000),
                          upper_values_sign = c("<="),
                          "Number of simulations (nsims)",
                          c("Value"),
                          "int",
                          NA) 

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

  if (is.null(parameters$alpha)) {
    parameters$alpha = 0.025
  }

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

  # Treatment effect assumptions 

  if (endpoint_index == 1) {

    if (is.null(parameters$control_mean)) stop("Mean effect in the control arm (control_mean): Value must be specified.", call. = FALSE)

    control_mean = 
          ContinuousErrorCheck(parameters$control_mean, 
                               n_endpoints, 
                               lower_values = c(NA),
                               lower_values_sign = c(NA),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Mean effect in the control arm (control_mean)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_mean)) stop("Mean effect in the treatment arm or arms (treatment_mean): Value must be specified.", call. = FALSE)

    treatment_mean = 
          ContinuousErrorCheck(parameters$treatment_mean, 
                               n_comparisons * n_endpoints, 
                               lower_values = c(NA),
                               lower_values_sign = c(NA),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Mean effect in the treatment arm (treatment_mean)",
                               c("Value"),
                               "double",
                               NA) 

    if (n_comparisons >= 2 & n_endpoints >= 2) {

        if (!is.matrix(parameters$treatment_mean)) stop("Mean effect in the treatment arm (treatment_mean) must be a matrix.", call. = FALSE)

        if (any(dim(parameters$treatment_mean) != c(n_endpoints, n_comparisons))) stop("Mean effect in the treatment arm (treatment_mean) does not have the correction dimensions.", call. = FALSE) 

    } 

    if (is.null(parameters$control_sd)) stop("Standard deviation in the control arm (control_sd): Value must be specified.", call. = FALSE)

    control_sd = 
          ContinuousErrorCheck(parameters$control_sd, 
                               n_endpoints, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Standard deviation in the control arm (control_sd)",
                               c("Value"),
                               "double",
                               NA) 


    if (is.null(parameters$treatment_sd)) stop("Standard deviation in the treatment arm or arms (treatment_sd): Value must be specified.", call. = FALSE)

    treatment_sd = 
          ContinuousErrorCheck(parameters$treatment_sd, 
                               n_comparisons * n_endpoints, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Standard deviation in the treatment arm (treatment_sd)",
                               c("Value"),
                               "double",
                               NA) 

    if (n_comparisons >= 2 & n_endpoints >= 2) {

        if (!is.matrix(parameters$treatment_sd)) stop("Standard deviation in the treatment arm (treatment_sd) must be a matrix.", call. = FALSE)

        if (any(dim(parameters$treatment_sd) != c(n_endpoints, n_comparisons))) stop("Standard deviation in the treatment arm (treatment_sd) does not have the correction dimensions.", call. = FALSE) 

    } 

  }

  if (endpoint_index == 2) {

    if (is.null(parameters$control_rate)) stop("Response rate in the control arm (control_rate): Value must be specified.", call. = FALSE)

    control_rate = 
          ContinuousErrorCheck(parameters$control_rate, 
                               n_endpoints, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(1),
                               upper_values_sign = c("<"),
                               "Response rate in the control arm (control_rate)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_rate)) stop("Response rate in the treatment arm (treatment_rate): Value must be specified.", call. = FALSE)

    treatment_rate = 
          ContinuousErrorCheck(parameters$treatment_rate, 
                               n_comparisons * n_endpoints, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(1),
                               upper_values_sign = c("<"),
                               "Response rate in the treatment arm (treatment_rate)",
                               c("Value"),
                               "double",
                               NA) 

    if (n_comparisons >= 2 & n_endpoints >= 2) {

        if (!is.matrix(parameters$treatment_rate)) stop("Response rate in the treatment arm (treatment_rate) must be a matrix.", call. = FALSE)

        if (any(dim(parameters$treatment_rate) != c(n_endpoints, n_comparisons))) stop("Response rate in the treatment arm (treatment_rate) does not have the correction dimensions.", call. = FALSE) 

    }       

  }

  if (n_comparisons >= 2 & n_endpoints == 1) {

    # Analysis of multiple treatment-control comparisons
    parameters$mult_test_type = 1 

    n_hypotheses = n_comparisons

    mult_test_list = c("Bonferroni", "Holm", "Hochberg", "Hommel", "Fixed-sequence", "Chain")

    if (!tolower(parameters$mult_test) %in% tolower(mult_test_list)) stop("Multiple testing procedure (mult_test): Value must be Bonferroni, Holm, Hochberg, Hommel, Fixed-sequence or Chain.", call. = FALSE)

    for (i in 1:length(mult_test_list)) {
        if (tolower(mult_test_list[i]) == tolower(parameters$mult_test)) mult_test_index = i
    }   

    parameters$mult_test_index = mult_test_index

  }

  if (n_comparisons == 1 & n_endpoints >= 2) {

    # Analysis of multiple endpoints
    parameters$mult_test_type = 2 

    n_hypotheses = n_endpoints

    mult_test_list = c("Bonferroni", "Holm", "Hochberg", "Hommel", "Fixed-sequence", "Chain", "O'Brien")

    if (!tolower(parameters$mult_test) %in% tolower(mult_test_list)) stop("Multiple testing procedure (mult_test): Value must be Bonferroni, Holm, Hochberg, Hommel, Fixed-sequence, Chain or O'Brien.", call. = FALSE)

    for (i in 1:length(mult_test_list)) {
        if (tolower(mult_test_list[i]) == tolower(parameters$mult_test)) mult_test_index = i
    }   

    parameters$mult_test_index = mult_test_index

  }

  if (n_comparisons >= 2 & n_endpoints >= 2) {

    # Analysis of multiple endpoints and multiple treatment-control comparisons
    parameters$mult_test_type = 3 

    n_hypotheses = n_comparisons * n_endpoints

    mult_test_list = c("Holm", "Hochberg", "Hommel")

    if (!tolower(parameters$mult_test) %in% tolower(mult_test_list)) stop("Multiple testing procedure (mult_test): Value must be Holm, Hochberg or Hommel.", call. = FALSE)

    for (i in 1:length(mult_test_list)) {
        if (tolower(mult_test_list[i]) == tolower(parameters$mult_test)) mult_test_index = i
    }   

    parameters$mult_test_index = mult_test_index

  }

  if (parameters$mult_test_type %in% c(1, 2)) {

    # Chain procedure
    if (mult_test_index == 6) {

        if (is.null(parameters$weights)) {

          weights = rep(1 / n_hypotheses, n_hypotheses)

        } else {

          weights = 
                ContinuousErrorCheck(parameters$weights, 
                                     n_hypotheses, 
                                     lower_values = c(0),
                                     lower_values_sign = c(">="),
                                     upper_values = c(1),
                                     upper_values_sign = c("<="),
                                     "Hypothesis weights (weights)",
                                     c("Value"),
                                     "double",
                                     NA)           
        }

        if (sum(weights) > 1) stop("Hypothesis weights (weights): Sum of weights must be less than or equal to 1.", call. = FALSE)

        parameters$weights = weights  

        transition = 
            ContinuousErrorCheck(as.vector(parameters$transition), 
                                 n_hypotheses * n_hypotheses, 
                                 lower_values = c(0),
                                 lower_values_sign = c(">="),
                                 upper_values = c(1),
                                 upper_values_sign = c("<="),
                                 "Hypothesis transition matrix (transition)",
                                 c("Value"),
                                 "double",
                                 NA) 

        if (!is.matrix(parameters$transition)) stop("Hypothesis transition matrix (transition) must be a matrix.", call. = FALSE)

        if (any(dim(parameters$transition) != c(n_hypotheses, n_hypotheses))) stop("Hypothesis transition matrix (transition) does not have the correction dimensions.", call. = FALSE) 

        for (i in 1:n_hypotheses) if (sum(parameters$transition[i, ]) > 1) stop("Hypothesis transition matrix (transition): The sum of elements in each row must be less than or equal to 1.", call. = FALSE)   

        parameters$ctransition = t(parameters$transition)

    } 

    # All other procedures but fixed-sequence or chain
    if (mult_test_index %in% c(1, 2, 3, 4)) {

      if (is.null(parameters$weights)) {

        weights = rep(1 / n_hypotheses, n_hypotheses)

      } else {

        weights = 
              ContinuousErrorCheck(parameters$weights, 
                                   n_hypotheses, 
                                   lower_values = c(0),
                                   lower_values_sign = c(">"),
                                   upper_values = c(1),
                                   upper_values_sign = c("<"),
                                   "Hypothesis weights (weights)",
                                   c("Value"),
                                   "double",
                                   NA)           
      }

      if (sum(weights) > 1) stop("Sum of hypothesis weights (weights) must be less than or equal to 1.", call. = FALSE)

      parameters$weights = weights  

      parameters$transition = 0 
      parameters$ctransition = 0 

    }

    # Fixed-sequence procedure
    if (mult_test_index == 5) {

      sequence = 
            ContinuousErrorCheck(parameters$sequence, 
                                 n_hypotheses, 
                                 lower_values = c(1),
                                 lower_values_sign = c(">="),
                                 upper_values = c(n_hypotheses),
                                 upper_values_sign = c("<="),
                                 "Hypothesis testing sequence (sequence)",
                                 c("Value"),
                                 "int",
                                 NA)           

      if (any(sort(parameters$sequence) != 1:n_hypotheses)) stop("Hypothesis testing sequence for the fixed-sequence procedure (sequence) is not correctly specified.", call. = FALSE)

      parameters$weights = parameters$sequence
      parameters$transition = 0 
      parameters$ctransition = 0 

    }

    # Global test
    if (mult_test_index == 7) {

      parameters$weights = 0
      parameters$transition = 0 
      parameters$ctransition = 0 

    }

  }

  if (parameters$mult_test_type == 3) {

    mult_method_list = c("Standard", "Modified", "Enhanced")

    if (!tolower(parameters$mult_method) %in% tolower(mult_method_list)) stop("Mixture method (mult_method): Value must be Standard, Modified or Enhanced.", call. = FALSE)

    for (i in 1:length(mult_method_list)) {
        if (tolower(mult_method_list[i]) == tolower(parameters$mult_method)) mult_method_index = i
    }   

    parameters$mult_method_index = mult_method_index

    # Truncation parameters
    mult_test_gamma = 
        ContinuousErrorCheck(as.vector(parameters$mult_test_gamma), 
                             n_endpoints, 
                             lower_values = c(0),
                             lower_values_sign = c(">"),
                             upper_values = c(1),
                             upper_values_sign = c("<="),
                             "Truncation parameters (mult_test_gamma)",
                             c("Value"),
                             "double",
                             NA) 

    for (i in 1:(n_endpoints - 1)) {
        if (mult_test_gamma[i] == 1) stop("Truncation parameters (mult_test_gamma): Value must be less than 1 in all families except for the last one.", call. = FALSE)
    }   

    parameters$mult_test_gamma[n_endpoints] = 1

  }

  if (parameters$mult_test_type %in% c(2, 3)) {

    # Correlation matrix
    endpoint_correlation = 
        ContinuousErrorCheck(as.vector(parameters$endpoint_correlation), 
                             n_endpoints * n_endpoints, 
                             lower_values = c(0),
                             lower_values_sign = c(">="),
                             upper_values = c(1),
                             upper_values_sign = c("<="),
                             "Endpoint correlation matrix (endpoint_correlation)",
                             c("Value"),
                             "double",
                             NA) 

    parameters$corr_sum = sum(parameters$endpoint_correlation)    

    if (!is.matrix(parameters$endpoint_correlation)) stop("Endpoint correlation matrix (endpoint_correlation) must be a matrix.", call. = FALSE)

    if (any(dim(parameters$endpoint_correlation) != c(n_endpoints, n_endpoints))) stop("Endpoint correlation matrix (endpoint_correlation) does not have the correction dimensions.", call. = FALSE)

    if (det(parameters$endpoint_correlation) <= 0) stop("Endpoint correlation matrix (endpoint_correlation) must be positive definite.", call. = FALSE)

  }

  #############################################

  parameters$n_hypotheses = n_hypotheses

  parameters$means = 0
  parameters$sds = 0
  parameters$rates = 0
  parameters$hazard_rates = 0
  parameters$dropout_parameter = 0
  parameters$enrollment_distribution = 2
  parameters$max_sample_size = max(parameters$sample_size)

  # All means and SDs
  if (endpoint_index == 1) {
    parameters$means = c(parameters$control_mean, parameters$treatment_mean)
    parameters$sds = c(parameters$control_sd, parameters$treatment_sd)
    parameters$enrollment_period = 0
    parameters$enrollment_parameter = 0
    parameters$event_count = 0
    parameters$control_rate = 0
    if (parameters$mult_test_type != 3) parameters$treatment_rate = 0 else parameters$treatment_rate = matrix(0, 1, 1)
    parameters$control_hazard_rate = 0
    if (parameters$mult_test_type != 3) parameters$treatment_hazard_rate = 0 else parameters$treatment_hazard_rate = matrix(0, 1, 1)
  }

  # All rates
  if (endpoint_index == 2) {
    parameters$rates = c(parameters$control_rate, parameters$treatment_rate)
    parameters$enrollment_period = 0
    parameters$enrollment_parameter = 0
    parameters$event_count = 0
    parameters$control_mean = rep(0, n_endpoints)
    if (parameters$mult_test_type != 3) parameters$treatment_mean = rep(0, n_endpoints) else parameters$treatment_mean = matrix(rep(0, n_endpoints), 1, n_endpoints)
    parameters$control_sd = rep(1, n_endpoints)
    if (parameters$mult_test_type != 3) parameters$treatment_sd = rep(1, n_endpoints) else parameters$treatment_sd = matrix(rep(1, n_endpoints), 1, n_endpoints)
    parameters$control_hazard_rate = 0
    if (parameters$mult_test_type != 3) parameters$treatment_hazard_rate = 0 else parameters$treatment_hazard_rate = matrix(0, 1, 1)

  }

  # Total sample size after accounting for dropout rates
  if (endpoint_index != 3) {
    parameters$sample_size_adj = floor(parameters$sample_size * (1 - parameters$dropout_rate))
  }

  ###########################################################

  # Run simulations to compute key characteristics

  # Analysis of multiple treatment-control comparisons
  if (parameters$mult_test_type == 1) {

    simulations = MultAdj1NCores(parameters)

    sim_results = simulations$sim_results

    # Add column names
    column_names = NULL
    for (i in 1:n_comparisons) column_names = c(column_names, paste0("pvalue", i))
    for (i in 1:n_comparisons) column_names = c(column_names, paste0("adjpvalue", i))
    colnames(sim_results) = column_names

    sim_summary = list()

    # Unadjusted power
    pvalue = sim_results[, 1:n_comparisons]
    sim_summary$power = colMeans(pvalue <= alpha)

    # Adjusted overall power
    adj_pvalue = sim_results[, (n_comparisons + 1):(2 * n_comparisons)]
    sim_summary$adj_power = colMeans(adj_pvalue <= alpha)

    # Disjunctive and conjunctive power
    disj_power = 0
    conj_power = 0
    for (i in 1:nsims) {
      disj_power = disj_power + any(adj_pvalue[i, ] <= alpha)
      conj_power = conj_power + all(adj_pvalue[i, ] <= alpha)
    }

    sim_summary$disj_power = disj_power / nsims
    sim_summary$conj_power = conj_power / nsims

  }

  # Analysis of multiple endpoints
  if (parameters$mult_test_type == 2) {

    simulations = MultAdj2NCores(parameters)

    sim_results = simulations$sim_results

    # Add column names
    column_names = NULL
    for (i in 1:n_endpoints) column_names = c(column_names, paste0("pvalue", i))
    for (i in 1:n_endpoints) column_names = c(column_names, paste0("adjpvalue", i))
    colnames(sim_results) = column_names

    sim_summary = list()

    # Unadjusted power
    pvalue = sim_results[, 1:n_endpoints]
    sim_summary$power = colMeans(pvalue <= alpha)

    # Adjusted marginal power
    if (parameters$mult_test_index <= 6) {

      adj_pvalue = sim_results[, (n_endpoints + 1):(2 * n_endpoints)]
      sim_summary$adj_power = colMeans(adj_pvalue <= alpha)

      # Disjunctive and conjunctive power
      disj_power = 0
      conj_power = 0
      for (i in 1:nsims) {
        disj_power = disj_power + any(adj_pvalue[i, ] <= alpha)
        conj_power = conj_power + all(adj_pvalue[i, ] <= alpha)
      }

      sim_summary$disj_power = disj_power / nsims
      sim_summary$conj_power = conj_power / nsims

    }

    # Adjusted overall power
    if (parameters$mult_test_index == 7) {

      adj_pvalue = sim_results[, n_endpoints + 1]
      sim_summary$adj_power = mean(adj_pvalue <= alpha)

    }

  }

  # Analysis of multiple treatment-control comparisons and multiple endpoints
  if (parameters$mult_test_type == 3) {

    simulations = MultAdj3NCores(parameters)

    sim_results = simulations$sim_results

    # Add column names
    column_names = NULL
    for (i in 1:n_hypotheses) column_names = c(column_names, paste0("pvalue", i))
    for (i in 1:n_hypotheses) column_names = c(column_names, paste0("adjpvalue", i))
    colnames(sim_results) = column_names

    sim_summary = list()

    # Unadjusted power
    pvalue = sim_results[, 1:n_hypotheses]
    sim_summary$power = colMeans(pvalue <= alpha)

    # Adjusted overall power
    adj_pvalue = sim_results[, (n_hypotheses + 1):(2 * n_hypotheses)]
    sim_summary$adj_power = colMeans(adj_pvalue <= alpha)

    # Disjunctive and conjunctive power for each endpoint
    disj_power = rep(0, n_endpoints)
    conj_power = rep(0, n_endpoints)
    for (i in 1:nsims) {
      k = 0
      for (j in 1:n_endpoints) {
        disj_power[j] = disj_power[j] + any(adj_pvalue[i, (k + 1):(k + n_comparisons)] <= alpha)
        conj_power[j] = conj_power[j] + all(adj_pvalue[i, (k + 1):(k + n_comparisons)] <= alpha)
        k = k + n_comparisons
      }
    }

    sim_summary$disj_power = disj_power / nsims
    sim_summary$conj_power = conj_power / nsims

  }

  results = list(parameters = parameters,
                 sim_results = sim_results,
                 sim_summary = sim_summary)

  class(results) = "MultAdjResults"

  return(results)

}    

MultAdjReportDoc = function(results) {

   #############################################################################

   # Error checks

   if (class(results) != "MultAdjResults") stop("The object was not created by the MultAdj function", call. = FALSE)

  #############################################################################

  # Empty list of tables to be included in the report
  item_list = list()
  item_index = 1
  table_index = 1
  figure_index = 1

  width = 6.5
  height = 5

  parameters = results$parameters
  sim_summary = results$sim_summary
  endpoint_index = parameters$endpoint_index
  ncomparisons = parameters$n_comparisons
  nendpoints = parameters$n_endpoints
  nhypotheses = parameters$n_hypotheses
  mult_test_type = parameters$mult_test_type
  mult_test_index = parameters$mult_test_index

  # Trial arms  
  trial_arms = "Control"
  if (ncomparisons >= 2) {
    for (i in 1:ncomparisons) trial_arms = c(trial_arms, paste0("Treatment ", i))
  } else {
    trial_arms = c(trial_arms, "Treatment")
  }

  # Endpoints  
  endpoints = rep("", nendpoints)
  if (nendpoints >= 2) {
    for (i in 1:nendpoints) endpoints[i] = paste0("Endpoint ", i)
  } 

  # Hypotheses
  short_hypotheses = rep("", nhypotheses)
  for (i in 1:nhypotheses) short_hypotheses[i] = paste0("H", i)

  hypotheses = rep("", nhypotheses)
  for (i in 1:nhypotheses) hypotheses[i] = paste0("Hypothesis ", i)

  # Treatment effect assumptions
  if (mult_test_type %in% c(1, 2)) {

    # All means and SDs
    if (endpoint_index == 1) {
      means = c(parameters$control_mean, parameters$treatment_mean)
      sds = c(parameters$control_sd, parameters$treatment_sd)
    }

    # All rates
    if (endpoint_index == 2) {
      rates = c(parameters$control_rate, parameters$treatment_rate)
    }

  }

  # Treatment effect assumptions
  if (mult_test_type == 3) {

    # All means and SDs
    if (endpoint_index == 1) {
      means = cbind(parameters$control_mean, parameters$treatment_mean)
      sds = cbind(parameters$control_sd, parameters$treatment_sd)
    }

    # All rates
    if (endpoint_index == 2) {
      rates = cbind(parameters$control_rate, parameters$treatment_rate)
    }

  }

  #############################################################################

  report_title = "Traditional design with multiple objectives"

  endpoint_labels = c("normally distributed", "binary")

  if (mult_test_type == 1) item_list[[item_index]] = list(type = "paragraph", label = "Description", value = paste0("The simulation report presents key operating characteristics of a multiplicity adjustment for multiple treatment-control comparisons in a trial with a ", tolower(endpoint_labels[endpoint_index]), " endpoint."))

  if (mult_test_type == 2) item_list[[item_index]] = list(type = "paragraph", label = "Description", value = paste0("The simulation report presents key operating characteristics of a multiplicity adjustment or a global testing procedure for multiple endpoints in a trial with a ", tolower(endpoint_labels[endpoint_index]), " endpoint."))

  if (mult_test_type == 3) item_list[[item_index]] = list(type = "paragraph", label = "Description", value = paste0("The simulation report presents key operating characteristics of a multiple testing procedure (gatekeeping procedure) for multiple endpoints and multiple treatment-control comparisons in a trial with a ", tolower(endpoint_labels[endpoint_index]), " endpoint."))

  item_index = item_index + 1

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

  if (mult_test_type == 1) {

    column_names = c("Trial arm", "Parameter", "Value")

    col1 = NULL
    col2 = NULL
    col3 = NULL

    for (i in 1:(ncomparisons + 1)) {

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

    }

    data_frame = data.frame(col1, col2, col3)

    column_width = c(2.5, 2.5, 1.5)

  }

  if (mult_test_type == 2) {

    column_names = c("Endpoint", "Trial arm", "Parameter", "Value")

    col1 = NULL
    col2 = NULL
    col3 = NULL
    col4 = NULL

    for (i in 1:nendpoints) {

      if (endpoint_index == 1) {
        col1 = c(col1, endpoints[i], rep("", 3))
        col2 = c(col2, "Control", "", "Treatment", "")       
        col3 = c(col3, "Mean", "SD", "Mean", "SD")
        col4 = c(col4, parameters$control_mean[i], parameters$control_sd[i], parameters$treatment_mean[i], parameters$treatment_sd[i])
      }

      if (endpoint_index == 2) {
        col1 = c(col1, endpoints[i], "")
        col2 = c(col2, "Control", "Treatment")       
        col3 = c(col3, "Rate (%)", "Rate (%)")
        col4 = c(col4, 100 * parameters$control_rate[i], 100 * parameters$treatment_rate[i])
      }

    }

    data_frame = data.frame(col1, col2, col3, col4)

    column_width = c(1.5, 1.5, 2, 1.5)

  }

  if (mult_test_type == 3) {

    column_names = c("Endpoint", "Trial arm", "Parameter", "Value")

    col1 = NULL
    col2 = NULL
    col3 = NULL
    col4 = NULL

    for (i in 1:nendpoints) {

      for (j in 1:(ncomparisons + 1)) {

        if (endpoint_index == 1) {
          col1 = c(col1, endpoints[i], "")
          col2 = c(col2, trial_arms[j], "")       
          col3 = c(col3, "Mean", "SD")
          col4 = c(col4, means[i, j], sds[i, j])

        }

        if (endpoint_index == 2) {
          col1 = c(col1, endpoints[i])
          col2 = c(col2, trial_arms[j])       
          col3 = c(col3, "Rate (%)")
          col4 = c(col4, 100 * rates[i, j])
        }

      }

    }

    data_frame = data.frame(col1, col2, col3, col4)

    column_width = c(1.5, 1.5, 2, 1.5)

  }

  if (parameters$direction_index == 1) footnote = "A higher value of the endpoint indicates a beneficial effect." 

  if (parameters$direction_index == 2) footnote = "A lower value of the endpoint indicates a beneficial effect." 

  title = paste0("Table ", table_index, ". Treatment effect assumptions")
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  if (mult_test_type == 1) {

    column_names = c("Hypothesis", "Definition")

    col1 = hypotheses
    col2 = rep("", ncomparisons)
    for (i in 1:ncomparisons) col2[i] = paste0("Null hypothesis of no difference between Treatment ", i, " and control")

    data_frame = data.frame(col1, col2)

    column_width = c(1.5, 5)

    title = paste0("Table ", table_index, ". Hypothesis definitions")
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  if (mult_test_type == 2) {

    column_names = c("Hypothesis", "Definition")

    col1 = hypotheses
    col2 = rep("", nendpoints)
    for (i in 1:nendpoints) col2[i] = paste0("Null hypothesis of no difference between the treatment and control with respect to Endpoint ", i)

    data_frame = data.frame(col1, col2)

    column_width = c(1.5, 5)

    title = paste0("Table ", table_index, ". Hypothesis definitions")
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  if (mult_test_type == 3) {

    column_names = c("Hypothesis", "Definition")

    col1 = hypotheses
    col2 = rep("", ncomparisons * nendpoints)
    k = 1
    for (i in 1:nendpoints) {
      for (j in 1:ncomparisons) {
        col2[k] = paste0("Null hypothesis of no difference between Treatment ", j, " and control with respect to Endpoint ", i)
        k = k + 1
      }
    }

    data_frame = data.frame(col1, col2)

    column_width = c(1.5, 5)

    title = paste0("Table ", table_index, ". Hypothesis definitions")
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  #############################################################################

  if (mult_test_type %in% c(2, 3)) {

    column_names = c("Endpoint", endpoints)

    data_frame = cbind(endpoints, parameters$endpoint_correlation)

    column_width = c(1, rep(5.5 / nendpoints, nendpoints))

    title = paste0("Table ", table_index, ". Endpoint correlation matrix")
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  #############################################################################

  if (mult_test_type %in% c(1, 2)) {

    if (mult_test_index <= 6) {

      column_names = c("Parameter", "Value")

      col1 = c("Multiple testing procedure")
      col2 = c(parameters$mult_test)
      
      data_frame = data.frame(col1, col2)

      column_width = c(5, 1.5)

      title = paste0("Table ", table_index, ". Multiple testing procedure")
      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
      item_index = item_index + 1
      table_index = table_index + 1

      # Fixed-sequence procedure
      if (mult_test_index == 5) {

        column_names = short_hypotheses

        data_frame = as.data.frame(matrix(0, 1, nhypotheses))

        data_frame[1, ] = round(parameters$sequence)

        title = paste0("Table ", table_index, ". Multiple testing procedure: Hypothesis testing sequence")

        column_width = c(rep(6.5 / nhypotheses, nhypotheses))
        item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
        item_index = item_index + 1
        table_index = table_index + 1

      }

      # All other procedures but fixed-sequence
      if (mult_test_index %in% c(1, 2, 3, 4, 6)) {

        column_names = short_hypotheses

        data_frame = as.data.frame(matrix(0, 1, nhypotheses))

        data_frame[1, ] = round(parameters$weights, 4)

        title = paste0("Table ", table_index, ". Multiple testing procedure: Initial hypothesis weights")

        column_width = c(rep(6.5 / nhypotheses, nhypotheses))
        item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
        item_index = item_index + 1
        table_index = table_index + 1

      }

      # Chain procedure
      if (mult_test_index == 6) {

        data_frame = as.data.frame(matrix(0, nhypotheses, nhypotheses + 1))
        data_frame[, 1] = hypotheses
        temp = round(parameters$transition, 4)
        temp = matrix(temp, nhypotheses, nhypotheses)
        data_frame[, 2:(nhypotheses + 1)] = temp

        column_names = c("Hypothesis", short_hypotheses)

        title = paste0("Table ", table_index, ". Multiple testing procedure: Transition parameters")

        column_width = c(1.25, rep(5.25 / nhypotheses, nhypotheses))
        item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
        item_index = item_index + 1
        table_index = table_index + 1

      }
    }

    if (mult_test_index == 7) {

      column_names = c("Parameter", "Value")

      col1 = c("Global testing procedure")
      col2 = c(parameters$mult_test)
      
      data_frame = data.frame(col1, col2)

      column_width = c(5, 1.5)

      title = paste0("Table ", table_index, ". Global testing procedure")
      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
      item_index = item_index + 1
      table_index = table_index + 1

    }


  }

  if (mult_test_type == 3) {

    column_names = c("Parameter", "Value")

    col1 = c("Component procedure", "Mixture method", "Truncation parameters")
    col2 = c(parameters$mult_test, parameters$mult_method, paste0(parameters$mult_test_gamma, collapse = ", "))
    
    data_frame = data.frame(col1, col2)

    column_width = c(5, 1.5)

    title = paste0("Table ", table_index, ". Multiple testing procedure (Gatekeeping procedure)")
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  #############################################################################

  column_names = c("Parameter", "Value")

  if (endpoint_index %in% c(1, 2)) {
    col1 = c("Dropout rate (%)")
    col2 = c(100 * parameters$dropout_rate)
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

  # Traditional multiplicity adjustments or gatekeeping procedures
  if (mult_test_type == 1 | (mult_test_type == 2 & mult_test_index <= 6) | mult_test_type == 3) {

    column_names = c("Hypothesis", "Power (%)", "Adjusted power (%)")

    col1 = hypotheses
    col2 = round(100 * sim_summary$power, 1)
    col3 = round(100 * sim_summary$adj_power, 1)
    
    data_frame = data.frame(col1, col2, col3)
    title = paste0("Table ", table_index, ". Simulation results: Hypothesis-specific power")

    column_width = c(2.5, 2, 2)

    footnote = "Power: Probability of rejecting each hypothesis of no effect without a multiplicity adjustment. Adjusted power: Probability of rejecting each hypothesis of no effect using a multiplicity adjustment based on the specified multiple testing procedure."

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  # Global testing procedure
  if (mult_test_type == 2 & mult_test_index == 7) {

    column_names = c("Hypothesis", "Power (%)")

    col1 = hypotheses
    col2 = round(100 * sim_summary$power, 1)
    
    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Simulation results: Hypothesis-specific power")

    column_width = c(2.5, 4)

    footnote = "Power: Probability of rejecting each hypothesis of no effect without a multiplicity adjustment."

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

    column_names = c("Overall power", "Power (%)")

    col1 = c("Global power")
    col2 = round(100 * sim_summary$adj_power, 1)
    
    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Simulation results: Overall power")

    column_width = c(2.5, 4)

    footnote = "Global power: Probability of rejecting at least one hypothesis of no effect using the specified global testing approach."

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  # Traditional multiplicity adjustments
  if (mult_test_type == 1 | (mult_test_type == 2 & mult_test_index <= 6)) {

    column_names = c("Overall power", "Power (%)")

    col1 = c("Disjunctive power", "Conjunctive power")
    col2 = c(round(100 * sim_summary$disj_power, 1), round(100 * sim_summary$conj_power, 1))
    
    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Simulation results: Overall power")

    column_width = c(2.5, 4)

    footnote = "Disjunctive power: Probability of rejecting at least one hypothesis of no effect using a multiplicity adjustment based on the specified multiple testing procedure. Conjunctive power: Probability of rejecting all hypotheses of no effect using a multiplicity adjustment based on the specified multiple testing procedure."

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  # Gatekeeping procedures
  if (mult_test_type == 3) {

    column_names = c("Endpoint family", "Overall power", "Power (%)")

    col1 = NULL
    col2 = NULL
    col3 = NULL

    for (i in 1:nendpoints) {

      col1 = c(col1, paste0("Endpoint ", i), "")
      col2 = c(col2, "Disjunctive power", "Conjunctive power")
      col3 = c(col3, round(100 * sim_summary$disj_power[i], 1), round(100 * sim_summary$conj_power[i], 1))
    
    }

    data_frame = data.frame(col1, col2, col3)
    title = paste0("Table ", table_index, ". Simulation results: Overall power")

    column_width = c(2.5, 2, 2)

    footnote = "Disjunctive power: Probability of rejecting at least one hypothesis of no effect within each endpoint family using a multiplicity adjustment based on the specified multiple testing procedure. Conjunctive power: Probability of rejecting all hypotheses of no effect within each endpoint family using a multiplicity adjustment based on the specified multiple testing procedure."

    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

  }

  #############################################################################

  report = item_list

  doc = SaveReport(report, report_title)

  return(doc)

}
# End of MultAdjReportDoc

# --==[ MultAdj1 ]==--
MultAdj1NCores = function(parameters) {
  ncores = parameters$ncores

  # Run simulations on multiple cores to compute key characteristics

  if (ncores > 1) {
    # nocov start
    cl = parallel::makeCluster(ncores)

    # Export all functions in the global environment to each node
    parallel::clusterExport(cl, ls(envir = .GlobalEnv))

    doParallel::registerDoParallel(cl)
    simulation_list = foreach(counter=(1:ncores), .packages = c("MedianaDesigner")) %dorng% { 
      MultAdj1SingleCore(parameters)
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

    simulations = MultAdj1SingleCore(parameters)

  }

  return(simulations)
}

MultAdj1SingleCore = function(parameters) {
  params_for_core = parameters
  params_for_core$nsims = parameters$nsims_per_core

  simulations = MultAdjC1(params_for_core)
  return(simulations)
}

# --==[ MultAdj2 ]==--
MultAdj2NCores = function(parameters) {
  ncores = parameters$ncores

  # Run simulations on multiple cores to compute key characteristics

  if (ncores > 1) {
    # nocov start
    cl = parallel::makeCluster(ncores)

    # Export all functions in the global environment to each node
    parallel::clusterExport(cl, ls(envir = .GlobalEnv))

    doParallel::registerDoParallel(cl)
    simulation_list = foreach(counter=(1:ncores), .packages = c("MedianaDesigner")) %dorng% { 
      MultAdj2SingleCore(parameters)
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

    simulations = MultAdj2SingleCore(parameters)

  }

  return(simulations)
}

MultAdj2SingleCore = function(parameters) {
  params_for_core = parameters
  params_for_core$nsims = parameters$nsims_per_core

  simulations = MultAdjC2(params_for_core)
  return(simulations)
}

# --==[ MultAdj3 ]==--
MultAdj3NCores = function(parameters) {
  ncores = parameters$ncores

  # Run simulations on multiple cores to compute key characteristics

  if (ncores > 1) {
    # nocov start
    cl = parallel::makeCluster(ncores)

    # Export all functions in the global environment to each node
    parallel::clusterExport(cl, ls(envir = .GlobalEnv))

    doParallel::registerDoParallel(cl)
    simulation_list = foreach(counter=(1:ncores), .packages = c("MedianaDesigner")) %dorng% { 
      MultAdj3SingleCore(parameters)
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

    simulations = MultAdj3SingleCore(parameters)

  }

  return(simulations)
}

MultAdj3SingleCore = function(parameters) {
  params_for_core = parameters
  params_for_core$nsims = parameters$nsims_per_core

  simulations = MultAdjC3(params_for_core)
  return(simulations)
}
