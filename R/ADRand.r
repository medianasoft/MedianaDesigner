rowMax = function(x) {

  row_max = rep(0, nrow(x))

  for (i in 1:nrow(x)) row_max[i] = max(x[i, ])

  return(row_max)

}

# nocov start

# Truncated exponential distribution
TruncExpDist = function(x, rate, max) {

    x = x / max
    if (abs(rate) < 0.0001) y = x else y = (1 - exp(- rate * x)) / (1 - exp(- rate))

    return(y)

}

# Find the parameter for a given median enrollment time 
TruncatedExponentialMedian = function(median, max) {

        lower_bound = -9.0
        upper_bound = 10.0
        midpoint = (lower_bound + upper_bound) / 2

        while (upper_bound - lower_bound >= 0.001) {

            if (TruncExpDist(median, midpoint, max) < 0.5) 
                { 
                    lower_bound = midpoint
                } else {   
                    upper_bound = midpoint        
                }
            midpoint = (lower_bound + upper_bound) / 2.0

        }

        return(midpoint)

}

# nocov end

ADRand = function(parameters) {
  
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

    set.seed(random_seed)
        
    if (is.null(parameters$endpoint_type)) stop("Endpoint type (endpoint_type): Value must be specified.", call. = FALSE)

    if (!tolower(parameters$endpoint_type) %in% c("normal")) stop("Endpoint type (endpoint_type): Value must be Normal.", call. = FALSE)

    parameters$endpoint_index = 1  

    if (is.null(parameters$direction)) {
          parameters$direction_index = 1
    } else {
          if (!tolower(parameters$direction) %in% c("higher", "lower")) stop("Direction of favorable outcome (direction): Value must be specified.", call. = FALSE)
    }

    if (tolower(parameters$direction) == "higher") parameters$direction_index = 1    
    if (tolower(parameters$direction) == "lower") parameters$direction_index = 2    

    if (is.null(parameters$dose_levels)) stop("Dose levels in the trial (dose_levels): Value must be specified.", call. = FALSE)

    dose_levels = 
          ContinuousErrorCheck(parameters$dose_levels, 
                               NA, 
                               lower_values = 0,
                               lower_values_sign = ">=",
                               upper_values = 1000,
                               upper_values_sign = "<=",
                               "Dose levels in the trial (dose_levels)",
                               c("Value"),
                               "double",
                               NA) 

    n_doses = length(dose_levels)
    
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
                               n_doses - 1, 
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
                               n_doses - 1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = c(NA),
                               upper_values_sign = c(NA),
                               "Standard deviations in the treatment arms (treatment_sd)",
                               c("Value"),
                               "double",
                               NA) 

    parameters$mean = c(parameters$control_mean, parameters$treatment_mean)      
    parameters$sd = c(parameters$control_sd, parameters$treatment_sd)      

    if (is.null(parameters$stage_sample_size)) stop("Total number of enrolled patients in each stage (stage_sample_size): Value must be specified.", call. = FALSE)

    stage_sample_size = 
          ContinuousErrorCheck(parameters$stage_sample_size, 
                               NA, 
                               lower_values = 0,
                               lower_values_sign = ">",
                               upper_values = NA,
                               upper_values_sign = NA,
                               "Total number of enrolled patients in each stage (stage_sample_size)",
                               c("Value"),
                               "int",
                               NA) 

    if (is.null(parameters$ratio_placebo)) stop("Fixed randomization ratio in the placebo arm (ratio_placebo): Value must be specified.", call. = FALSE)

    ratio_placebo = 
          ContinuousErrorCheck(parameters$ratio_placebo, 
                               1, 
                               lower_values = 0,
                               lower_values_sign = ">",
                               upper_values = 1,
                               upper_values_sign = "<",
                               "Fixed randomization ratio in the placebo arm (ratio_placebo)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_period)) stop("Treatment period (treatment_period): Value must be specified.", call. = FALSE)

    treatment_period = 
          ContinuousErrorCheck(parameters$treatment_period, 
                               1, 
                               lower_values = 0,
                               lower_values_sign = ">",
                               upper_values = NA,
                               upper_values_sign = NA,
                               "Treatment period (treatment_period)",
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
                               upper_values = c(enrollment_period),
                               upper_values_sign = c("<"),
                               "Median enrollment time (enrollment_parameter)",
                               c("Value"),
                               "double",
                               NA)   

    if (is.null(parameters$delta)) stop("Clinically relevant improvement over placebo (delta): Value must be specified.", call. = FALSE)

    if (parameters$direction_index == 1) { 

      delta = 
            ContinuousErrorCheck(parameters$delta, 
                                 1, 
                                 lower_values = c(0),
                                 lower_values_sign = c(">="),
                                 upper_values = NA,
                                 upper_values_sign = NA,
                                 "Clinically relevant improvement over placebo (delta)",
                                 c("Value"),
                                 "double",
                                 NA) 

    }

    if (parameters$direction_index == 2) { 

      delta = 
            ContinuousErrorCheck(parameters$delta, 
                                 1, 
                                 lower_values = NA,
                                 lower_values_sign = NA,
                                 upper_values = c(0),
                                 upper_values_sign = c("<="),
                                 "Clinically relevant improvement over placebo (delta)",
                                 c("Value"),
                                 "double",
                                 NA) 

    }

    parameters$linear_model_parameter = 0
    
    if (is.null(parameters$exponential_model_parameter)) stop("Non-linear parameter for the exponential model (exponential_model_parameter): Value must be specified.", call. = FALSE)

    exponential_model_parameter = 
          ContinuousErrorCheck(parameters$exponential_model_parameter, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = NA,
                               upper_values_sign = NA,
                               "Non-linear parameter for the exponential model (exponential_model_parameter)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$emax_model_parameter)) stop("Non-linear parameter for the Emax model (emax_model_parameter): Value must be specified.", call. = FALSE)

    emax_model_parameter = 
          ContinuousErrorCheck(parameters$emax_model_parameter, 
                               1, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = NA,
                               upper_values_sign = NA,
                               "Non-linear parameter for the Emax model (emax_model_parameter)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$logistic_model_parameters)) stop("Non-linear parameters for the logistic model (logistic_model_parameters): Value must be specified.", call. = FALSE)

    logistic_model_parameters = 
          ContinuousErrorCheck(parameters$logistic_model_parameters, 
                               2, 
                               lower_values = c(0),
                               lower_values_sign = c(">"),
                               upper_values = NA,
                               upper_values_sign = NA,
                               "Non-linear parameters for the logistic model (logistic_model_parameters)",
                               c("Value"),
                               "double",
                               NA) 

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
      parameters$nsims = 1000     # nocov
    }

    if (!is.null(parameters$balance)) {

      balance = 
            ContinuousErrorCheck(parameters$balance, 
                                 1, 
                                 lower_values = c(0),
                                 lower_values_sign = c(">="),
                                 upper_values = c(3),
                                 upper_values_sign = c("<="),
                                 "Balance parameter for adaptive randomization (balance)",
                                 c("Value"),
                                 "double",
                                 NA) 
    } else {
      parameters$balance = 1
    }

    parameters$n_per_arm = round(sum(parameters$stage_sample_size) / n_doses)

    ###########################################################

    parameters$model_index = 1:4

    parameters$non_linear_vector = c(parameters$linear_model_parameter, 
                          parameters$exponential_model_parameter,
                          parameters$emax_model_parameter,
                          parameters$logistic_model_parameters)

    ###########################################################

    # Run simulations

    withCallingHandlers({        

        sim_results = ADRandC(parameters)

        },

        warning = function(c) {
          # nocov start
          msg <- conditionMessage(c)
          if ( grepl("the line search step became smaller than the minimum value allowed", msg, fixed = TRUE) ) {
            invokeRestart("muffleWarning")
          }
          # nocov end
        }

    )

    # Number of models
    n_models = 4

    # Number of stages in the trial
    n_stages = length(parameters$stage_sample_size)

    # Add column names
    column_names = c("placebo")
    for (i in 1:(n_doses - 1)) column_names = c(column_names, paste0("dose", i))
    colnames(sim_results$n) = column_names

    column_names = NULL
    for (i in 1:n_models) column_names = c(column_names, paste0("model", i))
    colnames(sim_results$traditional) = column_names

    column_names = NULL
    for (i in 1:n_models) column_names = c(column_names, paste0("model", i))
    colnames(sim_results$adaptive) = column_names

    column_names = NULL
    for (i in 1:n_stages) column_names = c(column_names, paste0("stage", i))
    colnames(sim_results$stage_sample_size) = column_names

    # Compute the critical values for the MCPMod method
    mcpmod_value = rep(0, parameters$nsims)

    for (i in 1:parameters$nsims) mcpmod_value[i] = MCPModCriticalValue(parameters, sim_results$n[i, ])
    sim_results$mcpmod_value = mcpmod_value

    ###########################################################

    # Return the results

    results = list(parameters = parameters,
                   sim_results = sim_results)

    class(results) = "ADRandResults"

    return(results)

}    


# Standardize a vector
Standardize = function(vec) {

  results = (vec - mean(vec)) / as.numeric(sqrt(t(vec - mean(vec)) %*% (vec - mean(vec))))

  return(results)

}

# Dose-response functions
DRFunction = function(model_index, coef, x) {

    # Linear model
    if (model_index == 1) {
        y = coef[1] + coef[2] * x
    }

    # Exponential model
    if (model_index == 2) {
        y = coef[1] + coef[2] * (exp(x / coef[3]) - 1)
    }

    # Emax model
    if (model_index == 3) {
        y = coef[1] + coef[2] * x / (coef[3] + x)
    }

    # Logistic model
    if (model_index == 4) {
            den = 1.0 + exp((coef[3] - x) / coef[4])
            y = coef[1] + coef[2] / den
    }

    return(y)

}

# Compute the model parameters to match the placebo and maximum effects
ComputeDRFunctionParameters = function(model_index, placebo_effect, max_effect, max_dose, parameters) {

    # Linear model
    if (model_index == 1) {
        coef = rep(0, 2)
        coef[1] = placebo_effect
        coef[2] = max_effect / max_dose
    }

    # Exponential model
    if (model_index == 2) {
        coef = rep(0, 3)
        coef[1] = placebo_effect
        coef[2] = max_effect / (exp(max_dose / parameters[1]) - 1)
        coef[3] = parameters[1]
    }

    # Emax model
    if (model_index == 3) {
        coef = rep(0, 3)
        coef[1] = placebo_effect
        coef[2] = max_effect * (parameters[1] + max_dose) / max_dose
        coef[3] = parameters[1]
    }

    # Logistic model
    if (model_index == 4) {
        coef = rep(0, 4)
        temp_coef = c(0, 1, parameters[1], parameters[2])
        temp =  max_effect / (DRFunction(4, temp_coef, max_dose) - DRFunction(4, temp_coef, 0))
        coef[1] = placebo_effect- temp * DRFunction(4, temp_coef, 0)
        coef[2] = temp
        coef[3] = parameters[1]
        coef[4] = parameters[2]
    }

    return(coef)

}

MCPModCriticalValue = function(parameters, n_groups) {

  model_index = parameters$model_index
  doses = parameters$dose_levels
  non_linear_vector = parameters$non_linear_vector
  alpha = parameters$alpha

  # Total number of models
  n_models = length(model_index)

  n_doses = length(doses)
  n_patients = sum(n_groups)

  max_dose = max(doses)
  diag_vec = rep(0, n_doses)

  dr_model = rep(0, n_doses)

  one_vec = rep(1, n_doses)

  Sinv = diag(n_groups)
  S = diag(1 / n_groups)

  #####################################################

  # Compute the model-specific optimal contrasts

  opt_contrast = matrix(0, n_doses, n_models)

  # Set up dose-response models based on the initial values
  for (i in 1:n_models) {

      if (model_index[i] == 1) non_linear_parameters = non_linear_vector[1]
      if (model_index[i] == 2) non_linear_parameters = non_linear_vector[2]
      if (model_index[i] == 3) non_linear_parameters = non_linear_vector[3]
      if (model_index[i] == 4) non_linear_parameters = non_linear_vector[4:5]

      # Parameters of a standardized model
      parameter_values = ComputeDRFunctionParameters(model_index[i], 0, 1, max_dose, non_linear_parameters)

      for (j in 1:n_doses) {

        dr_model[j] = DRFunction(model_index[i], parameter_values, doses[j]) 

      }

      contrast = Sinv %*% (dr_model - as.numeric(t(dr_model) %*% Sinv %*% one_vec) / as.numeric(t(one_vec) %*% Sinv %*% one_vec)) 

      opt_contrast[, i] = Standardize(contrast)

  }

  cov_mat = t(opt_contrast) %*% S %*% opt_contrast
  diag_mat = diag(sqrt(diag(cov_mat)))
  corr_matrix = solve(diag_mat) %*% cov_mat %*% solve(diag_mat)

  if (det(corr_matrix) > 1E-10) {

    # nocov start
    crit_value = qmvt(p = 1 - alpha, tail = "lower.tail", df = Inf, corr = corr_matrix, maxpts = 30000, abseps = 0.001, releps = 0, algorithm = GenzBretz())$quantile
    # nocov end

  } else crit_value = qnorm(1 - alpha)

  return(crit_value)

}

ADRandReportDoc = function(results) {

    #############################################################################

    # Error checks

    if (class(results) != "ADRandResults") stop("The object was not created by the ADRand function.", call. = FALSE)

    #############################################################################

    parameters = results$parameters
    simulations = results$sim_results

    dose_levels = parameters$dose_levels  
    n_doses = length(dose_levels)
    n_stages = length(parameters$stage_sample_size)

    # Trial arms
    trial_arms = rep("", n_doses)
    for (i in 1:n_doses) {
        if (i == 1) trial_arms[i] = "Placebo" else trial_arms[i] = paste0("Dose ", i - 1, " (", dose_levels[i], ")")
    }

    #############################################################################

    # Empty list of tables to be included in the report
    item_list = list()
    item_index = 1
    table_index = 1
    figure_index = 1

    #############################################################################

    report_title = "Adaptive designs with adaptive randomization"

    item_list[[item_index]] = list(type = "paragraph", label = "Description", value = "The simulation report presents key operating characteristics of an adaptive design for a dose-finding Phase II clinical trial with multiple interim analyses aimed at updating the randomization scheme based on the accumulating efficacy data.")

    item_index = item_index + 1

    #############################################################################

    column_names = c("Parameter", "Value")

    # nocov start
    if (parameters$direction_index == 1) label = "A higher value of the endpoint indicates a more favorable outcome"

    if (parameters$direction_index == 2) label = "A lower value of the endpoint indicates a more favorable outcome"
    # nocov end

    col1 = c("Endpoint type", "Direction of favorable outcome") 
    col2 = c(endpoint_list[parameters$endpoint_index], label)

    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Primary efficacy endpoint")

    column_width = c(3, 3.5)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Stage", "Planned number of enrolled patients")

    stage_sample_size = parameters$stage_sample_size
    n_stages = length(stage_sample_size)

    col1 = rep("", n_stages)
    col2 = stage_sample_size
    for (i in 1:n_stages) col1[i] = paste0("Stage ", i)

    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Trial stages")

    column_width = c(2.5, 4)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Trial arm", "Parameter", "Value") 

    col1 = NULL
    col2 = NULL
    col3 = NULL

    for (i in 1:n_doses) {

      col1 = c(col1, trial_arms[i], "")
      col2 = c(col2, "Mean", "SD")
      col3 = c(col3, parameters$mean[i], parameters$sd[i])

    }

    data_frame = data.frame(col1, col2, col3)
    title = paste0("Table ", table_index, ". Treatment effect assumptions")

    column_width = c(2.5, 2, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = c("Fixed randomization ratio in the placebo arm (%)",
             "Balance parameter for adaptive randomization",  
             "Clinically relevant improvement over placebo")
    col2 = c(100 * parameters$ratio_placebo,
             parameters$balance,
             parameters$delta)

    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Decision rule parameters")

    column_width = c(4.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Model", "Parameters")

    col1 = c("Exponential", "Emax", "Logistic")
    col2 = c(paste0("Delta = ", parameters$exponential_model_parameter),   
             paste0("ED50 = ", parameters$emax_model_parameter),
             paste0("ED50 = ", parameters$logistic_model_parameters[1], ", Delta = ", parameters$logistic_model_parameters[2]))


    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Parameters of candidate dose-response models used in the MCPMod method")

    footnote = "This table defines non-linear parameters of the candidate dose-response models and therefore no parameters are specified for the linear model."

    column_width = c(4.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = c("Length of the patient enrollment period", 
             "Median enrollment time",
             "Patient dropout rate (%)",
             "Length of the treatment period")
    col2 = c(parameters$enrollment_period,
             parameters$enrollment_parameter, 
             100 * parameters$dropout_rate,
             parameters$treatment_period)

    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Trial design parameters")

    footnote = "Median enrollment time: Time point by which 50% of the patients will be enrolled into the trial."

    column_width = c(4.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Parameter", "Value")

    col1 = c("One-sided Type I error rate", "Number of simulations")
    col2 = c(sprintf("%0.3f", parameters$alpha),
           sprintf("%d", parameters$nsims))

    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Simulation parameters")

    column_width = c(4.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    if (is.null(parameters$withoutCharts)) {   # skip chart generation on tests

      # Plot candidate dose-response models

      width = 6.5
      height = 5
      filename = "models.emf"

      # Linear model
      coef1 = ComputeDRFunctionParameters(1, 0, 1, max(dose_levels), NA)

      # Exponential model
      coef2 = ComputeDRFunctionParameters(2, 0, 1, max(dose_levels), parameters$exponential_model_parameter)

      # Emax model
      coef3 = ComputeDRFunctionParameters(3, 0, 1, max(dose_levels), parameters$
        emax_model_parameter)

      # Logistic model
      coef4 = ComputeDRFunctionParameters(4, 0, 1, max(dose_levels), parameters$logistic_model_parameters)

      x = seq(from = 0, to = max(dose_levels), length = 100)
      y1 = numeric(length(x))
      y2 = numeric(length(x))
      y3 = numeric(length(x))
      y4 = numeric(length(x))

      for (i in 1:length(x)) {

        y1[i] = DRFunction(1, coef1, x[i])
        y2[i] = DRFunction(2, coef2, x[i])
        y3[i] = DRFunction(3, coef3, x[i])
        y4[i] = DRFunction(4, coef4, x[i])

      }
   
      emf(file = filename, width = width, height = height)
      plot(x = x, y = y1, xlab="Dose", ylab="Response", xlim = c(0, max(dose_levels)), ylim = c(0, 1), type="l", lwd = 2, col = "black") 
      lines(x = x, y = y2, col = "blue", lwd = 2)
      lines(x = x, y = y3, col = "red", lwd = 2)
      lines(x = x, y = y4, col = "darkgreen", lwd = 2)
      dev.off() 

      item_list[[item_index]] =  list(label = paste0("Figure ", figure_index, ". Candidate dose-response models used in the MCPMod method."), 
                              filename = filename,
                              dim = c(width, height),
                              type = "emf_plot",
                              footnote = "Black curve: Linear model, Blue curve: Exponential model, Red curve: Emax model, Green curve: Logistic model.",
                              page_break = TRUE)

      item_index = item_index + 1
      figure_index = figure_index + 1

    }

    #############################################################################

    column_names = c("Stage", "Statistic", "Sample size")

    col1 = NULL
    col2 = NULL
    col3 = NULL

    for (i in 1:n_stages) {

      col1 = c(col1, paste0("Stage ", i), rep("", 3))
      col2 = c(col2, "Min", "Median", "Mean", "Max")
      col3 = c(col3, round(summary(simulations$stage_sample_size[, i])[c(1, 3, 4, 6)], 1))

    }

    data_frame = data.frame(col1, col2, col3)
    title = paste0("Table ", table_index, ". Simulation results: Sample size by stage")

    column_width = c(2, 2.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Dose", "Statistic", "Sample size")

    col1 = NULL
    col2 = NULL
    col3 = NULL

    for (i in 1:n_doses) {

      col1 = c(col1, trial_arms[i], rep("", 3))
      col2 = c(col2, "Min", "Median", "Mean", "Max")
      col3 = c(col3, round(summary(simulations$n[, i])[c(1, 3, 4, 6)], 1))

    }

    data_frame = data.frame(col1, col2, col3)
    title = paste0("Table ", table_index, ". Simulation results: Sample size by trial arm")

    column_width = c(2, 2.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    column_names = c("Design", "Power (%)") 

    traditional = mean(rowMax(simulations$traditional) >= simulations$mcpmod_value, na.rm = TRUE) 
    adaptive = mean(rowMax(simulations$adaptive) >= simulations$mcpmod_value, na.rm = TRUE) 

    col1 = c("Traditional", "Adaptive")
    col2 = c(round(100 * traditional, 1), 
             round(100 * adaptive, 1))

    data_frame = data.frame(col1, col2)
    title = paste0("Table ", table_index, ". Simulation results: Comparison of traditional and adaptive designs")

    footnote = "Probability of a statistically significant dose-response relationship based on the MCPMod method."

    column_width = c(4.5, 2)
    item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
    item_index = item_index + 1
    table_index = table_index + 1

    #############################################################################

    report = item_list

    doc = SaveReport(report, report_title)

    return(doc)

}
# End of ADRandReportDoc

