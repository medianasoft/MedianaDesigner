ComputePower = function(data_frame, convergence, alpha) {   

    ncomps = ncol(data_frame)
    nsims = nrow(data_frame)

    power = rep(0, ncomps + 1)
    outcomes = matrix(0, nsims, ncomps + 1)
    for (i in 1:ncomps) {
      outcomes[, i] = as.numeric(data_frame[, i] <= alpha)
      outcomes[convergence[, i] == 0, i] = NA
      power[i] = mean(outcomes[, i], na.rm = TRUE)
    }
    outcomes[, ncomps + 1] = as.numeric(rowSums(outcomes[, 1:ncomps]) >= 1)
    power[ncomps + 1] = mean(outcomes[, ncomps + 1], na.rm = TRUE)

    return(power)

}

AntiLogit = function(x) {
    return(1 / (1 + exp(-x)))
}

ClustRand = function(parameters) {

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

    if (is.null(parameters$method_type)) stop("Data analysis method (method_type): Value must be specified.", call. = FALSE)

    if (!tolower(parameters$method_type) %in% tolower(c("gee", "glmem"))) stop("Data analysis method (method_type): Value must be GEE or GLMEM", call. = FALSE)

    if (tolower(parameters$method_type) == "gee") method_index = 1  
    if (tolower(parameters$method_type) == "glmem") method_index = 2  

    parameters$method_index = method_index

    if (is.null(parameters$direction)) {
        parameters$direction_index = 1
    } else {
        if (!tolower(parameters$direction) %in% c("higher", "lower")) stop("Direction of favorable outcome (direction): Value must be specified.", call. = FALSE)
    }

    if (tolower(parameters$direction) == "higher") parameters$direction_index = 1    
    if (tolower(parameters$direction) == "lower") parameters$direction_index = 2    

    if (is.null(parameters$cluster_scheme)) stop("Cluster scheme (cluster_scheme): Value must be specified.", call. = FALSE)

    if (!tolower(parameters$cluster_scheme) %in% tolower(c("fixed", "random"))) stop("Cluster scheme (cluster_scheme): Value must be Fixed or Random.", call. = FALSE)

    if (tolower(parameters$cluster_scheme) == "fixed") cluster_index = 1  
    if (tolower(parameters$cluster_scheme) == "random") cluster_index = 2  

    parameters$cluster_index = cluster_index

    if (is.null(parameters$sample_size)) stop("Number of completers in each trial arm (sample_size): Value must be specified.", call. = FALSE)

    sample_size = ContinuousErrorCheck(parameters$sample_size, 
                                     NA, 
                                     lower_values = 0,
                                     lower_values_sign = ">",
                                     upper_values = 1000,
                                     upper_values_sign = "<=",
                                     "Number of completers in each trial arm (sample_size)",
                                     NA,
                                     "int",
                                     NA) 

    narms = length(sample_size)
    parameters$narms = narms

    mult_test_index = 0

    if (narms >= 3) {

      mult_test_list = c("Bonferroni", "Holm", "Hochberg")

      if (is.null(parameters$mult_test)) stop("Multiple testing procedure (mult_test): Value must be specified.", call. = FALSE)

      if (!tolower(parameters$mult_test) %in% tolower(mult_test_list)) stop("Multiple testing procedure (mult_test): Value must be Bonferroni, Holm or Hochberg.", call. = FALSE)

      for (i in 1:length(mult_test_list)) {
          if (tolower(mult_test_list[i]) == tolower(parameters$mult_test)) mult_test_index = i
      }   

    } 

    parameters$mult_test_index = mult_test_index

    parameters$treatment_cluster_size_vector = rep(0, 2)
    parameters$treatment_cluster_size_matrix = matrix(0, 2, 2)

    parameters$treatment_cluster_cum_vector = rep(0, 2)
    parameters$treatment_cluster_cum_matrix = matrix(0, 2, 2)

    # Fixed cluster sizes
    if (cluster_index == 1) {

      # Vectors of cluster sizes in trial arms

      if (is.null(parameters$control_cluster_size)) stop("Vector of cluster sizes in the control arm (control_cluster_size): Value must be specified.", call. = FALSE)

      control_cluster_size = ContinuousErrorCheck(parameters$control_cluster_size, 
                                       NA, 
                                       lower_values = 5,
                                       lower_values_sign = ">=",
                                       upper_values = 100,
                                       upper_values_sign = "<=",
                                       "Vector of cluster sizes in the control arm (control_cluster_size)",
                                       NA,
                                       "int",
                                       NA) 

      if (sum(parameters$control_cluster_size) != parameters$sample_size[1]) stop(
        "Cluster sizes in the control arm (control_cluster_size) are incorrectly specified.", call. = FALSE)

      if (is.null(parameters$treatment_cluster_size)) stop("Vector of cluster sizes in the treatment arm (treatment_cluster_size): Value must be specified.", call. = FALSE)

      if (narms == 2) {

        treatment_cluster_size = ContinuousErrorCheck(parameters$treatment_cluster_size, 
                                         NA, 
                                         lower_values = 5,
                                         lower_values_sign = ">=",
                                         upper_values = 100,
                                         upper_values_sign = "<=",
                                         "Vector of cluster sizes in the treatment arm (treatment_cluster_size)",
                                         NA,
                                         "int",
                                         NA) 

        if (sum(parameters$treatment_cluster_size) != parameters$sample_size[2]) stop(
          "Cluster sizes in the treatment arm (treatment_cluster_size) are incorrectly specified.", call. = FALSE)

        parameters$treatment_cluster_size_vector = parameters$treatment_cluster_size

      } else {

        if (!is(parameters$treatment_cluster_size, "matrix")) stop(
          "A matrix of cluster sizes in the treatment arms (treatment_cluster_size) must be specified.", call. = FALSE)

        if (nrow(parameters$treatment_cluster_size) != narms - 1) stop(
          "A matrix of cluster sizes in the treatment arms (treatment_cluster_size) must be specified.", call. = FALSE)

        treatment_cluster_size = ContinuousErrorCheck(parameters$treatment_cluster_size, 
                                         NA, 
                                         lower_values = 5,
                                         lower_values_sign = ">=",
                                         upper_values = 100,
                                         upper_values_sign = "<=",
                                         "Matrix of cluster sizes in the treatment arms (treatment_cluster_size)",
                                         NA,
                                         "int",
                                         NA) 

        for (i in 1:nrow(parameters$treatment_cluster_size)) if (sum(parameters$treatment_cluster_size[i, ]) != parameters$sample_size[i + 1]) stop(
          "Cluster sizes in the treatment arms (treatment_cluster_size) are incorrectly specified.", call. = FALSE)

        parameters$treatment_cluster_size_matrix = parameters$treatment_cluster_size

      }  

      parameters$control_cluster_cum = 0

    }

    # Random cluster sizes
    if (cluster_index == 2) {
      
      if (is.null(parameters$control_cluster_proportion)) stop("Vector of relative cluster sizes in the control arm (control_cluster_proportion): Value must be specified.", call. = FALSE)

      control_cluster_proportion = ContinuousErrorCheck(parameters$control_cluster_proportion, 
                                       NA, 
                                       lower_values = 0,
                                       lower_values_sign = ">",
                                       upper_values = 1,
                                       upper_values_sign = "<",
                                       "Vector of relative cluster sizes in the control arm (control_cluster_proportion)",
                                       NA,
                                       "double",
                                       NA) 

      if (abs(sum(parameters$control_cluster_proportion)-1) > 0.0001) stop(
        "Relative cluster sizes in the control arm (control_cluster_proportion) are incorrectly specified.", call. = FALSE)

      if (narms == 2) {

        treatment_cluster_proportion = ContinuousErrorCheck(parameters$treatment_cluster_proportion, 
                                         NA, 
                                         lower_values = 0,
                                         lower_values_sign = ">",
                                         upper_values = 1,
                                         upper_values_sign = "<",
                                         "Vector of relative cluster sizes in the treatment arm (treatment_cluster_proportion)",
                                         NA,
                                         "double",
                                         NA) 

        if (sum(parameters$treatment_cluster_proportion) != 1) stop(
          "Relative cluster sizes in the control arm (treatment_cluster_proportion) are incorrectly specified.", call. = FALSE)

        parameters$treatment_cluster_proportion_vector = parameters$treatment_cluster_proportion

      } else {

        if (!is(parameters$treatment_cluster_proportion, "matrix")) stop(
          "A matrix of relative cluster sizes in the treatment arms (treatment_cluster_proportion) must be specified.", call. = FALSE)

        if (nrow(parameters$treatment_cluster_proportion) != narms - 1) stop(
          "A matrix of relative cluster sizes in the treatment arms (treatment_cluster_proportion) must be specified.", call. = FALSE)

        treatment_cluster_proportion = ContinuousErrorCheck(parameters$treatment_cluster_proportion, 
                                         NA, 
                                         lower_values = 0,
                                         lower_values_sign = ">",
                                         upper_values = 1,
                                         upper_values_sign = "<",
                                         "Matrix of relative cluster sizes in the treatment arms (treatment_cluster_proportion)",
                                         NA,
                                         "double",
                                         NA) 

        for (i in 1:(narms - 1)) {
          if (sum(parameters$treatment_cluster_proportion[i, ]) != 1) stop(
          "Relative cluster sizes in the treatment arms (treatment_cluster_proportion) are incorrectly specified.", call. = FALSE)
        }

        parameters$treatment_cluster_proportion_matrix = parameters$treatment_cluster_proportion

      }  

      # Create vectors of cimulative proportions
      parameters$control_cluster_cum = cumsum(parameters$control_cluster_proportion)

      if (narms == 2) {
        parameters$treatment_cluster_cum_vector = cumsum(parameters$treatment_cluster_proportion)
      } else {
        parameters$treatment_cluster_cum_matrix = matrix(0, narms - 1, ncol(parameters$treatment_cluster_proportion))
        for (i in 1:(narms - 1)) {
          parameters$treatment_cluster_cum_matrix[i, ] = cumsum(parameters$treatment_cluster_proportion[i, ])        
        }
      }

      parameters$control_cluster_size = 0

    }

    if (!is.null(parameters$nsims)) {

      nsims = ContinuousErrorCheck(parameters$nsims, 
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
      parameters$nsims = 1000       # nocov
    }

    if (!is.null(parameters$ncores)) {
      # nocov start
      # Maximum number of cores
      max_ncores = parallel::detectCores()

      ncores = ContinuousErrorCheck(parameters$ncores, 
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

    # Compute descriptive statistics (arm-specific effects, cluster sizes) for each simulation run

    if (!is.null(parameters$descriptive_statistics)) {

      if (!is.logical(parameters$descriptive_statistics)) stop(
          "Descriptive statistics flag (descriptive_statistics): Value must be TRUE or FALSE.", call. = FALSE)
    
    } else {
      parameters$descriptive_statistics = FALSE # nocov
    }

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

    if (is.null(parameters$control_icc)) stop("Intracluster correlation coefficient in the control arm (control_icc): Value must be specified.", call. = FALSE)

    control_icc = 
          ContinuousErrorCheck(parameters$control_icc, 
                               1, 
                               lower_values = 0,
                               lower_values_sign = ">",
                               upper_values = 1,
                               upper_values_sign = "<",
                               "Intracluster correlation coefficient in the control arm (control_icc)",
                               c("Value"),
                               "double",
                               NA) 

    if (is.null(parameters$treatment_icc)) stop("Intracluster correlation coefficient in the treatment arms (treatment_icc): Value must be specified.", call. = FALSE)

    treatment_icc = 
          ContinuousErrorCheck(parameters$treatment_icc, 
                               narms - 1, 
                               lower_values = 0,
                               lower_values_sign = ">",
                               upper_values = 1,
                               upper_values_sign = "<",
                               "Intracluster correlation coefficient in the treatment arm (treatment_icc)",
                               c("Value"),
                               "double",
                               NA) 


    # Normal endpoint     
    if (endpoint_index == 1) {

      if (is.null(parameters$control_between_cluster_sd)) stop("Between-cluster standard deviation in the control arm (control_sd): Value must be specified.", call. = FALSE)

      control_between_cluster_sd = 
            ContinuousErrorCheck(parameters$control_between_cluster_sd, 
                                 1, 
                                 lower_values = c(0),
                                 lower_values_sign = c(">"),
                                 upper_values = c(NA),
                                 upper_values_sign = c(NA),
                                 "Between-cluster standard deviation in the control arm (control_between_cluster_sd)",
                                 c("Value"),
                                 "double",
                                 NA) 

      if (is.null(parameters$treatment_between_cluster_sd)) stop("Between-cluster standard deviation in the treatment arms (treatment_sd): Value must be specified.", call. = FALSE)

      treatment_between_cluster_sd = 
            ContinuousErrorCheck(parameters$treatment_between_cluster_sd, 
                                 narms - 1, 
                                 lower_values = c(0),
                                 lower_values_sign = c(">"),
                                 upper_values = c(NA),
                                 upper_values_sign = c(NA),
                                 "Between-cluster standard deviation in the treatment arm (treatment_between_cluster_sd)",
                                 c("Value"),
                                 "double",
                                 NA) 

      # Within-cluster standard deviation in the control arm 
      parameters$control_within_cluster_sd = control_between_cluster_sd * (1 - control_icc) / control_icc

      # Within-cluster standard deviation in the treatment arm 
      parameters$treatment_within_cluster_sd = treatment_between_cluster_sd * (1 - treatment_icc) / treatment_icc

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

      if (is.null(parameters$treatment_mean)) stop("Mean effect in the treatment arms (treatment_mean): Value must be specified.", call. = FALSE)

      treatment_mean = 
            ContinuousErrorCheck(parameters$treatment_mean, 
                                 narms - 1, 
                                 lower_values = c(NA),
                                 lower_values_sign = c(NA),
                                 upper_values = c(NA),
                                 upper_values_sign = c(NA),
                                 "Mean effect in the treatment arm (treatment_mean)",
                                 c("Value"),
                                 "double",
                                 NA) 

        parameters$control_alpha = 0
        parameters$control_beta = 0
        parameters$treatment_alpha = 0
        parameters$treatment_beta = 0
        parameters$control_rate = 0
        parameters$treatment_rate = 0

    }

    # Binary endpoint
    if (endpoint_index == 2) {

      if (is.null(parameters$control_rate)) stop("Response rate in the control arm (control_rate): Value must be specified.", call. = FALSE)

      control_rate = 
            ContinuousErrorCheck(parameters$control_rate, 
                                 1, 
                                 lower_values = 0,
                                 lower_values_sign = ">",
                                 upper_values = 1,
                                 upper_values_sign = "<",
                                 "Response rate in the control arm (control_rate)",
                                 c("Value"),
                                 "double",
                                 NA) 

      if (is.null(parameters$treatment_rate)) stop("Response rate in the treatment arms (treatment_rate): Value must be specified.", call. = FALSE)

      treatment_rate = 
            ContinuousErrorCheck(parameters$treatment_rate, 
                                 narms - 1, 
                                 lower_values = 0,
                                 lower_values_sign = ">",
                                 upper_values = 1,
                                 upper_values_sign = "<",
                                 "Response rate in the treatment arm (treatment_rate)",
                                 c("Value"),
                                 "double",
                                 NA) 

      # Alpha and beta parameters of the beta distribution in the control arm
      parameters$control_alpha = (1 / control_icc - 1) * control_rate 
      parameters$control_beta = (1 / control_icc - 1) * (1 - control_rate) 

      # Alpha and beta parameters of the beta distribution in the treatment arm
      parameters$treatment_alpha = (1 / treatment_icc - 1) * treatment_rate 
      parameters$treatment_beta = (1 / treatment_icc - 1) * (1 - treatment_rate) 

      parameters$control_mean = 0
      parameters$treatment_mean = 0

      parameters$control_between_cluster_sd = 0
      parameters$treatment_between_cluster_sd = 0

      parameters$control_within_cluster_sd = 0
      parameters$treatment_within_cluster_sd = 0

    }

    #############################################

    parameters$means = 0
    parameters$within_cluster_sds = 0
    parameters$between_cluster_sds = 0
    parameters$rates = 0

    # All means and SDs
    if (endpoint_index == 1) {
        parameters$means = c(parameters$control_mean, parameters$treatment_mean)
        parameters$within_cluster_sds = c(parameters$control_within_cluster_sd, parameters$treatment_within_cluster_sd)
        parameters$between_cluster_sds = c(parameters$control_between_cluster_sd, parameters$treatment_between_cluster_sd)
    }

    # All rates
    if (endpoint_index == 2) {
        parameters$rates = c(parameters$control_rate, parameters$treatment_rate)
    }

    # Maximum number of iterations when fitting models
    parameters$max_iterations = 100

    # Minimum increment for stopping the model fitting process
    parameters$min_increment = 0.001

    # Parameters of multiple testing procedures
    parameters$weight = rep(1 / (narms - 1), narms - 1)
    parameters$transition = 0

    #########################################################

    # Run simulations on multple cores to compute operating characteristics

    if (parameters$ncores > 1) {
        # nocov start
        cl = parallel::makeCluster(parameters$ncores)

        # Export all functions in the global environment to each node
        parallel::clusterExport(cl, ls(envir = .GlobalEnv))

        doParallel::registerDoParallel(cl)
        simulation_list = foreach(counter=(1:parameters$ncores), .packages = c("ClustRand")) %dorng% { 
          ClustRandSingleCore(parameters)
        }
        stopCluster(cl)

        # Combine the simulation results across the cores 

        pval_results = NULL
        coef_results = NULL
        cluster_size_results = NULL

        for (i in 1:parameters$ncores) {
          pval_results = rbind(pval_results, simulation_list[[i]]$pval_results)
          coef_results = rbind(coef_results, simulation_list[[i]]$coef_results)
          cluster_size_results = c(cluster_size_results, simulation_list[[i]]$cluster_size_results)
        }
        # nocov end
    } else {

        simulations = ClustRandSingleCore(parameters)
        pval_results = simulations$pval_results
        coef_results = simulations$coef_results
        cluster_size_results = simulations$cluster_size_results

    }

    pval_results = as.data.frame(pval_results)

    if (parameters$descriptive_statistics) {

      coef_results = as.data.frame(coef_results)
      cluster_size_results = as.numeric(cluster_size_results)

    }

    sim_summary = list()

    # GEE
    if (parameters$method_index == 1) {

      # Extract p-values for each treatment-control comparison
      if (narms == 2) {

        pval_convergence = pval_results[, 1]
        sandwich = pval_results[, 2]
        kc = pval_results[, 3]
        md = pval_results[, 4]

        # Compute overall power
        sim_summary$power_sandwich = mean(sandwich[pval_convergence == 1] <= parameters$alpha)
        sim_summary$power_kc = mean(kc[pval_convergence == 1] <= parameters$alpha)
        sim_summary$power_md = mean(md[pval_convergence == 1] <= parameters$alpha)

        if (parameters$descriptive_statistics) {

            # Summary of random cluster size estimates
            if (parameters$cluster_scheme == "Random") {

              cluster_size = cluster_size_results

              cluster_size_summary = rep(0, 5)
              cluster_size_summary[1] = mean(cluster_size, na.rm = TRUE)
              cluster_size_summary[2:5] = quantile(cluster_size, probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)

              sim_summary$cluster_size_summary = cluster_size_summary

            }          

            coef_convergence = coef_results[, 1]
            intercept = coef_results[, 2]
            slope = coef_results[, 3]
            if (parameters$direction_index == 2) {
              intercept = -intercept
              slope = -slope
            }

            # Compute effect estimates 
            nsims = nrow(coef_results)
            effect = matrix(NA, nsims, narms) 

            for (i in 1:nsims) {

              if (coef_convergence[i] == 1) {

                if (endpoint_index == 1) { 
                  effect[i, 1] = intercept[i]
                  effect[i, 2] = intercept[i] + slope[i]
                }

                if (endpoint_index == 2) {                  # nocov start
                  effect[i, 1] = AntiLogit(intercept[i]) 
                  effect[i, 2] = AntiLogit(intercept[i] + slope[i]) 
                }                                           # nocov end

              }

            }

            # Compute the mean and 95% CI
            eff_summary = matrix(0, 2, 5)
            for (i in 1:2) {
              eff_summary[i, 1] = mean(effect[, i], na.rm = TRUE)
              eff_summary[i, 2:5] = quantile(effect[, i], probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)
            }

            sim_summary$eff_summary = eff_summary

        }

      } else {

        pval_convergence = pval_results[, seq(from = 1, to = 1 + (narms - 2) * 4, by = 4)]
        sandwich = pval_results[, seq(from = 2, to = 2 + (narms - 2) * 4, by = 4)]
        kc = pval_results[, seq(from = 3, to = 3 + (narms - 2) * 4, by = 4)]
        md = pval_results[, seq(from = 4, to = 4 + (narms - 2) * 4, by = 4)]

        # Compute treatment-specific and overall power
        sim_summary$power_sandwich = ComputePower(sandwich, pval_convergence, parameters$alpha)  
        sim_summary$power_kc = ComputePower(kc, pval_convergence, parameters$alpha)    
        sim_summary$power_md = ComputePower(md, pval_convergence, parameters$alpha)    

        if (parameters$descriptive_statistics) {

          # Summary of random cluster size estimates
          if (parameters$cluster_scheme == "Random") {

            cluster_size = cluster_size_results

            cluster_size_summary = rep(0, 5)
            cluster_size_summary[1] = mean(cluster_size, na.rm = TRUE)
            cluster_size_summary[2:5] = quantile(cluster_size, probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)

            sim_summary$cluster_size_summary = cluster_size_summary

          }          

          coef_convergence = coef_results[, seq(from = 1, to = 1 + (narms - 2) * 3, by = 3)]
          intercept = coef_results[, seq(from = 2, to = 2 + (narms - 2) * 3, by = 3)]
          slope = coef_results[, seq(from = 3, to = 3 + (narms - 2) * 3, by = 3)]

          if (parameters$direction_index == 2) {
            intercept = -intercept
            slope = -slope
          }

          # Compute effect estimates 
          nsims = nrow(coef_results)
          effect = matrix(NA, nsims, narms) 

          for (i in 1:nsims) {

            if (coef_convergence[i, 1] == 1) {

              if (endpoint_index == 1) effect[i, 1] = intercept[i, 1]
              if (endpoint_index == 2) effect[i, 1] = AntiLogit(intercept[i, 1]) 

            }

            for (j in 2:narms) {

              if (coef_convergence[i, j - 1] == 1) {

                if (endpoint_index == 1) effect[i, j] = intercept[i, j - 1] + slope[i, j - 1]
                if (endpoint_index == 2) effect[i, j] = AntiLogit(intercept[i, j - 1] + slope[i, j - 1])

            }
          
          }

        }  

        # Compute the mean and 95% CI
        eff_summary = matrix(0, narms, 5)
        for (i in 1:narms) {
          eff_summary[i, 1] = mean(effect[, i], na.rm = TRUE)
          eff_summary[i, 2:5] = quantile(effect[, i], probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)
        }

        sim_summary$eff_summary = eff_summary

        }
      }  

    }

    # GLMEM
    if (parameters$method_index == 2) {

      # Extract p-values for each treatment-control comparison
      if (narms == 2) {

        pval_convergence = pval_results[, 1]
        pval = pval_results[, 2]

        # Compute overall power
        sim_summary$power = mean(pval[pval_convergence == 1] <= parameters$alpha)

        if (parameters$descriptive_statistics) {

          # Summary of random cluster size estimates
          if (parameters$cluster_scheme == "Random") {

            cluster_size = cluster_size_results

            cluster_size_summary = rep(0, 5)
            cluster_size_summary[1] = mean(cluster_size, na.rm = TRUE)
            cluster_size_summary[2:5] = quantile(cluster_size, probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)

            sim_summary$cluster_size_summary = cluster_size_summary

          }          

          # Summary of effect estimates
          coef_convergence = coef_results[, 1]
          intercept = coef_results[, 2]
          slope = coef_results[, 3]
          if (parameters$direction_index == 2) {
            intercept = -intercept
            slope = -slope
          }

          nsims = nrow(coef_results)
          effect = matrix(NA, nsims, narms) 

          for (i in 1:nsims) {

            if (coef_convergence[i] == 1) {

              if (endpoint_index == 1) {      # nocov start
                effect[i, 1] = intercept[i]
                effect[i, 2] = intercept[i] + slope[i]
              }                               # nocov end

              if (endpoint_index == 2) { 
                effect[i, 1] = AntiLogit(intercept[i]) 
                effect[i, 2] = AntiLogit(intercept[i] + slope[i]) 
              }

            }

          }

          eff_summary = matrix(0, narms, 5)
          for (i in 1:narms) {
            eff_summary[i, 1] = mean(effect[, i], na.rm = TRUE)
            eff_summary[i, 2:5] = quantile(effect[, i], probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)
          }

          sim_summary$eff_summary = eff_summary

        }

      } else {

        pval_convergence = pval_results[, seq(from = 1, to = 1 + (narms - 2) * 2, by = 2)]
        pval = pval_results[, seq(from = 2, to = 2 + (narms - 2) * 2, by = 2)]

        # Compute treatment-specific and overall power
        sim_summary$power = ComputePower(pval, pval_convergence, parameters$alpha)  

        if (parameters$descriptive_statistics) {

          # Summary of random cluster size estimates
          if (parameters$cluster_scheme == "Random") {

            cluster_size = cluster_size_results

            cluster_size_summary = rep(0, 5)
            cluster_size_summary[1] = mean(cluster_size, na.rm = TRUE)
            cluster_size_summary[2:5] = quantile(cluster_size, probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)

            sim_summary$cluster_size_summary = cluster_size_summary

          }          

          # Summary of effect estimates
          coef_convergence = coef_results[, seq(from = 1, to = 1 + (narms - 2) * 3, by = 3)]
          intercept = coef_results[, seq(from = 2, to = 2 + (narms - 2) * 3, by = 3)]
          slope = coef_results[, seq(from = 3, to = 3 + (narms - 2) * 3, by = 3)]

          if (parameters$direction_index == 2) {  # nocov start
            intercept = -intercept
            slope = -slope
          }                                       # nocov end

          nsims = nrow(coef_results)
          effect = matrix(NA, nsims, narms) 

          for (i in 1:nsims) {

            if (coef_convergence[i, 1] == 1) {

              if (endpoint_index == 1) effect[i, 1] = intercept[i, 1]
              if (endpoint_index == 2) effect[i, 1] = AntiLogit(intercept[i, 1])  # nocov

            }

            for (j in 2:narms) {

              if (coef_convergence[i, j - 1] == 1) {

                if (endpoint_index == 1) effect[i, j] = intercept[i, j - 1] + slope[i, j - 1]
                if (endpoint_index == 2) effect[i, j] = AntiLogit(intercept[i, j - 1] + slope[i, j - 1])    # nocov

            }
          
          }

        }  

        eff_summary = matrix(0, narms, 5)
        for (i in 1:narms) {
          eff_summary[i, 1] = mean(effect[, i], na.rm = TRUE)
          eff_summary[i, 2:5] = quantile(effect[, i], probs = c(0.025, 0.25, 0.75, 0.975), na.rm = TRUE)
        }

        sim_summary$eff_summary = eff_summary

        }
      }  
    }

    results = list(parameters = parameters,
                   sim_summary = sim_summary,
                   pval_results = pval_results,
                   coef_results = coef_results,
                   cluster_size_results = cluster_size_results)

    class(results) = "ClustRandResults"

    return(results)
}    

ClustRandReportDoc = function(results) {

   #############################################################################

   # Error checks

   if (!is(results, "ClustRandResults")) stop("The object was not created by the ClustRand function", call. = FALSE)

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
  narms = parameters$narms

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
  }

  # All rates
  if (endpoint_index == 2) {
    rates = c(parameters$control_rate, parameters$treatment_rate)
  }

  # Other parameters
  between_cluster_sds = c(parameters$control_between_cluster_sd, parameters$treatment_between_cluster_sd)
  iccs = c(parameters$control_icc, parameters$treatment_icc)    

  #############################################################################

  report_title = "Cluster-randomized design"

  item_list[[item_index]] = list(type = "paragraph", label = "Description", value = "The simulation report presents key operating characteristics of a cluster-randomized design for a Phase III clinical trial.")

  item_index = item_index + 1

  #############################################################################

  column_names = c("Trial arm", "Sample size")

  col1 = trial_arms
  col2 = parameters$sample_size

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Number of completers")

  footnote = "Completers are defined as patients who complete the trial and are included in the final analysis."

  column_width = c(5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
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

  if (parameters$cluster_index == 1) {

      column_names = c("Trial arm", "Number of clusters", "Cluster sizes")

      col1 = trial_arms
      if (narms == 2) {
        col2 = c(length(parameters$control_cluster_size), length(parameters$treatment_cluster_size)) 
        col3 = c(paste0(parameters$control_cluster_size, collapse = ", "), paste0(parameters$treatment_cluster_size, collapse = ", "))
      } else {
        col2 = c(length(parameters$control_cluster_size), rep(ncol(parameters$treatment_cluster_size), narms - 1))         
        col3 = c(paste0(parameters$control_cluster_size, collapse = ", "))
        for (i in 1:(narms - 1)) col3 = c(col3, paste0(parameters$treatment_cluster_size[i, ], collapse = ", "))
      }   

      data_frame = data.frame(col1, col2, col3)
      title = paste0("Table ", table_index, ". Cluster scheme parameters")

      footnote = "A design with pre-defined cluster sizes is assumed."

      column_width = c(1.5, 1.5, 3.5)
      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
      item_index = item_index + 1
      table_index = table_index + 1

  }

  if (parameters$cluster_index == 2) {

      column_names = c("Trial arm", "Number of clusters", "Relative cluster sizes")

      col1 = trial_arms
      if (narms == 2) {
        col2 = c(length(parameters$control_cluster_proportion), length(parameters$treatment_cluster_proportion)) 
        col3 = c(paste0(parameters$control_cluster_proportion, collapse = ", "), paste0(parameters$treatment_cluster_proportion, collapse = ", "))
      } else {
        col2 = c(length(parameters$control_cluster_proportion), rep(ncol(parameters$treatment_cluster_proportion), narms - 1))         
        col3 = c(paste0(parameters$control_cluster_proportion, collapse = ", "))
        for (i in 1:(narms - 1)) col3 = c(col3, paste0(parameters$treatment_cluster_proportion[i, ], collapse = ", "))
      }   

      data_frame = data.frame(col1, col2, col3)
      title = paste0("Table ", table_index, ". Cluster scheme parameters")

      footnote = paste0("A design with random cluster sizes is assumed. The cluster sizes are generated using a generalized Bernoulli distribution based on the relative cluster sizes.")

      column_width = c(1.5, 1.5, 3.5)
      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
      item_index = item_index + 1
      table_index = table_index + 1

  }

  #############################################################################

  column_names = c("Trial arm", "Parameter", "Value")

  col1 = NULL
  col2 = NULL
  col3 = NULL

  for (i in 1:narms) {

    if (endpoint_index == 1) {
      col1 = c(col1, trial_arms[i])
      col2 = c(col2, "Mean")
      col3 = c(col3, means[i])
    }

    if (endpoint_index == 2) {
      col1 = c(col1, trial_arms[i])
      col2 = c(col2, "Rate (%)")
      col3 = c(col3, 100 * rates[i])
    }

  }

  data_frame = data.frame(col1, col2, col3)
  title = paste0("Table ", table_index, ". Treatment effect assumptions")

  column_width = c(1.5, 3.5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Trial arm", "Parameter", "Value")

  col1 = NULL
  col2 = NULL
  col3 = NULL

  for (i in 1:narms) {

    if (parameters$endpoint_index == 1) {

      col1 = c(col1, trial_arms[i], "")
      col2 = c(col2, "Between-cluster standard deviation", "Intra-cluster correlation coefficient")
      col3 = c(col3, between_cluster_sds[i], iccs[i])

    }

    if (parameters$endpoint_index == 2) {

      col1 = c(col1, trial_arms[i])
      col2 = c(col2, "Intra-cluster correlation coefficient")
      col3 = c(col3, iccs[i])

    }


  }

  data_frame = data.frame(col1, col2, col3)
  title = paste0("Table ", table_index, ". Variability parameters")

  column_width = c(1.5, 3.5, 1.5)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  column_names = c("Parameter", "Value")

  if (parameters$method_index == 1) method_type = "Generalized estimating equations (GEE)" 
  if (parameters$method_index == 2) method_type = "Generalized linear mixed effects model (GLMEM)" 

  col1 = c("Data analysis method")
  col2 = c(method_type)

  if (parameters$mult_test_index >= 1) {
    col1 = c(col1, "Multiple testing procedure")
    col2 = c(col2, parameters$mult_test)    
  }

  col1 = c(col1, "One-sided Type I error rate", "Number of simulations")
  col2 = c(col2, sprintf("%0.3f", sum(parameters$alpha)), sprintf("%d", parameters$nsims))

  data_frame = data.frame(col1, col2)
  title = paste0("Table ", table_index, ". Analysis and simulation parameters")

  column_width = c(2.5, 4)
  item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, TRUE)
  item_index = item_index + 1
  table_index = table_index + 1

  #############################################################################

  if (parameters$descriptive_statistics) {

      column_names = c("Parameter", "Trial arm", "Statistic", "Value")

      eff_summary = sim_summary$eff_summary
      cluster_size_summary = sim_summary$cluster_size_summary

      if (endpoint_index == 1) {

        col1 = c("Mean", rep("", 3 * narms - 1))
        col2 = NULL
        col3 = NULL
        col4 = NULL
        for (i in 1:narms) {

          col2 = c(col2, trial_arms[i], rep("", 2))
          col3 = c(col3, "Mean", "25%P, 75%P", "2.5%P, 97.5%P")
          col4 = c(col4, round(eff_summary[i, 1], 2), paste0(round(eff_summary[i, 3], 2), ", ", round(eff_summary[i, 4], 2)), paste0(round(eff_summary[i, 2], 2), ", ", round(eff_summary[i, 5], 2)))

        }

      }

      if (endpoint_index == 2) {

        col1 = c("Rate (%)", rep("", 3 * narms - 1))
        col2 = NULL
        col3 = NULL
        col4 = NULL
        for (i in 1:narms) {

          col2 = c(col2, trial_arms[i], rep("", 2))
          col3 = c(col3, "Mean", "25%P, 75%P", "2.5%P, 97.5%P")
          col4 = c(col4, round(100 * eff_summary[i, 1], 1), paste0(round(100 * eff_summary[i, 3], 1), ", ", round(100 * eff_summary[i, 4], 1)), paste0(round(100 * eff_summary[i, 2], 1), ", ", round(100 * eff_summary[i, 5], 1)))

        }

      }

      if (parameters$cluster_scheme == "Random") {

        col1 = c(col1, "Cluster size", rep("", 2))
        col2 = c(col2, "All trial arms", rep("", 2))
        col3 = c(col3, "Mean", "25%P, 75%P", "2.5%P, 97.5%P")
        col4 = c(col4, round(cluster_size_summary[1], 1), paste0(round(cluster_size_summary[3], 1), ", ", round(cluster_size_summary[4], 1)), paste0(round(cluster_size_summary[2], 1), ", ", round(cluster_size_summary[5], 1)))

      }

      data_frame = data.frame(col1, col2, col3, col4)
      title = paste0("Table ", table_index, ". Simulation results: Descriptive statistics")

      column_width = c(2, 1.5, 2, 1)

      footnote = "Descriptive statistics computed from each simulation run: 2.5%P (2.5th percentile), 25%P (25th percentile), 75%P (75th percentile) and 97.5%P (97.5th percentile)."

      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE, footnote)
      item_index = item_index + 1
      table_index = table_index + 1

  }

  #############################################################################

  # GEE
  if (parameters$method_index == 1) {

    if (parameters$narms == 2) {

      column_names = c("Analysis method", "Power (%)")

      col1 = c("Standard analysis", 
               "Bias-corrected analysis (Kauermann-Carroll correction)",
               "Bias-corrected analysis (Mancl-DeRouen correction)") 
               
      col2 = c(round(100 * sim_summary$power_sandwich, 1),
               round(100 * sim_summary$power_kc, 1),
               round(100 * sim_summary$power_md, 1))
      
      data_frame = data.frame(col1, col2)
      title = paste0("Table ", table_index, ". Simulation results: Power calculations")

      column_width = c(5.5, 1)

      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
      item_index = item_index + 1
      table_index = table_index + 1
    
    } else {

      column_names = c("Comparison", "Analysis method", "Power (%)")

      col1 = NULL
      col2 = NULL
      col3 = NULL

      for (i in 1:narms) {

        if (i < narms) comp = paste0(trial_arms[i + 1], " vs control") else comp = "Overall effect"
        col1 = c(col1, comp, "", "")

        col2 = c(col2, "Standard analysis", 
                 "Bias-corrected analysis (Kauermann-Carroll correction)",
                 "Bias-corrected analysis (Mancl-DeRouen correction)") 
                 
        col3 = c(col3, round(100 * sim_summary$power_sandwich[i], 1),
                 round(100 * sim_summary$power_kc[i], 1),
                 round(100 * sim_summary$power_md[i], 1))

      }

      data_frame = data.frame(col1, col2, col3)
      title = paste0("Table ", table_index, ". Simulation results: Power calculations")

      column_width = c(1.5, 4, 1)

      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
      item_index = item_index + 1
      table_index = table_index + 1

    }

  }

  # GLMEM
  if (parameters$method_index == 2) {

    if (parameters$narms == 2) {

      column_names = c("Comparison", "Power (%)")

      col1 = "Treatment vs control"
      col2 = round(100 * sim_summary$power, 1)
      
      data_frame = data.frame(col1, col2)
      title = paste0("Table ", table_index, ". Simulation results: Power calculations")

      column_width = c(5.5, 1)

      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
      item_index = item_index + 1
      table_index = table_index + 1
    
    } else {

      column_names = c("Comparison", "Power (%)")

      col1 = NULL
      col2 = NULL

      for (i in 1:narms) {

        if (i < narms) comp = paste0(trial_arms[i + 1], " vs control") else comp = "Overall effect"
        col1 = c(col1, comp)
                 
        col2 = c(col2, round(100 * sim_summary$power[i], 1))

      }

      data_frame = data.frame(col1, col2)
      title = paste0("Table ", table_index, ". Simulation results: Power calculations")

      column_width = c(5.5, 1)

      item_list[[item_index]] = CreateTable(data_frame, column_names, column_width, title, FALSE)
      item_index = item_index + 1
      table_index = table_index + 1

    }

  }

  #############################################################################

  report = item_list

  doc = SaveReport(report, report_title)

  return(doc)

}
# End of ClustRandReportDoc

ClustRandSingleCore = function(parameters) {

  params_for_run = parameters
  params_for_run$nsims = params_for_run$nsims_per_core

  # GEE
  if (params_for_run$method_index == 1) simulations = ClustRandGEEC(params_for_run)

  # GLMEM
  if (params_for_run$method_index == 2) simulations = ClustRandGLMEMR(params_for_run)

  return(simulations)

}

ClustRandGLMEMR = function(parameters) { 

  nsims = parameters$nsims
  narms = parameters$narms

  control_cluster_size = parameters$control_cluster_size

  if (narms == 2) treatment_size = length(parameters$treatment_cluster_size) else treatment_size = ncol(parameters$treatment_cluster_size)

  pval_results = NULL
  coef_results = NULL
  cluster_size_results = NULL

  for (sim in 1:nsims) {

    # Random cluster sizes
    if (parameters$cluster_index == 2) {

        # Generate a vector of random cluster sizes
        control_cluster_size = ExportRandomClusterSize(parameters$sample_size[1], parameters$control_cluster_cum)

    }

    control_size = length(control_cluster_size)

    # Generate data for the common control arm
    x_control = NULL
    y_control = NULL
    id_control = NULL

    for (i in 1:control_size) {

      # Normal endpoint
      if (parameters$endpoint_index == 1) {

          cluster_term = rnorm(n = 1, mean = 0, sd = parameters$between_cluster_sds[1])
          control_sample = rnorm(n = control_cluster_size[i], mean = parameters$means[1] + cluster_term, sd = parameters$within_cluster_sds[1])

      }

      # Binary endpoint
      if (parameters$endpoint_index == 2) {

        cluster_term = rbeta(n = 1, shape1 = parameters$control_alpha, shape2 = parameters$control_beta)
        control_sample = rbinom(n = control_cluster_size[i], size = 1, prob = cluster_term)

      }

      x_control = c(x_control, rep(0, control_cluster_size[i]))
      y_control = c(y_control, control_sample)
      id_control = c(id_control, rep(i, control_cluster_size[i]))

    }

    p = rep(NA, narms - 1)
    convergence = rep(1, narms - 1)
    intercept = rep(NA, narms - 1)
    slope = rep(NA, narms - 1)
    cluster_size = control_cluster_size

    # Generate data for each treatment arm
    for (k in 1:(narms - 1)) {

      x = x_control
      y = y_control
      id = id_control

      # Fixed cluster sizes
      if (parameters$cluster_index == 1) {

        if (narms == 2) treatment_cluster_size = parameters$treatment_cluster_size else treatment_cluster_size = parameters$treatment_cluster_size[k, ] 

      }

      # Random cluster sizes
      if (parameters$cluster_index == 2) {

          if (narms == 2) proportions = parameters$treatment_cluster_cum_vector else proportions = parameters$treatment_cluster_cum_matrix[k, ]     # nocov

          # Generate a vector of random cluster sizes
          treatment_cluster_size = ExportRandomClusterSize(parameters$sample_size[k + 1], proportions)

      }

      cluster_size = c(cluster_size, treatment_cluster_size)

      treatment_size = length(treatment_cluster_size)

      for (i in 1:treatment_size) {

        # Normal endpoint
        if (parameters$endpoint_index == 1) {

            cluster_term = rnorm(n = 1, mean = 0, sd = parameters$between_cluster_sds[k + 1])
            current_treatment_sample = rnorm(n = treatment_cluster_size[i], mean = parameters$means[k + 1] + cluster_term, sd = parameters$within_cluster_sds[k + 1])

        }

        # Binary endpoint
        if (parameters$endpoint_index == 2) {

            cluster_term = rbeta(n = 1, shape1 = parameters$treatment_alpha[k], shape2 = parameters$treatment_beta[k])
            current_treatment_sample = rbinom(n = treatment_cluster_size[i], size = 1, prob = cluster_term)

        }

        x = c(x, rep(1, treatment_cluster_size[i]))
        y = c(y, current_treatment_sample)
        id = c(id, rep(control_size + i, treatment_cluster_size[i]))

      }

      # Direction of favorable outcome
      if (parameters$direction_index == 2) {
          if (parameters$endpoint_index == 1) y = -y
          if (parameters$endpoint_index == 2) y = 1 - y
      }

      data_set = data.frame(x, y, id)

      # Analyze the data for the current treatment-control comparison using generalized linear mixed-effects models

      # Normal endpoint
      if (parameters$endpoint_index == 1) {

          result = tryCatch(
            {
              lmerTest::lmer(y ~ x + (1|id), data = data_set, REML = TRUE)
            },
            warning = function(cond) {     # nocov start
              w = -1
              class(w) = "try-error"
              return(w)
            })                             # nocov end

          if (!inherits(result, "try-error")) {

              model_fit  = result
              summary_fit = summary(model_fit, ddf = "Kenward-Roger")

              # Estimated model parameters
              if (nrow(summary_fit$coefficients) == 2) {
                intercept[k] = summary_fit$coefficients[1, 1]
                slope[k] = summary_fit$coefficients[2, 1]
                delta = slope[k]
                # One-sided p-value
                if (delta > 0) p[k] = summary_fit$coefficients[2, 5] / 2 else p[k] = 1 - summary_fit$coefficients[2, 5] / 2
              }

              if (nrow(summary_fit$coefficients) == 1) {                      # nocov start
                if (rownames(summary_fit$coefficients) == "(Intercept)") {
                   intercept[k] = summary_fit$coefficients[1, 1] 
                }
                if (rownames(summary_fit$coefficients) == "x") {
                   slope[k] = summary_fit$coefficients[2, 1]
                }
              }                                                               # nocov end

          } else {
            convergence[k] = 0    # nocov
          }

      }

      # Binary endpoint
      if (parameters$endpoint_index == 2) {

          result = tryCatch({

              glmer(y ~ x + (1|id), data = data_set, family = binomial(link = 'logit'))

            },
            warning = function(cond) {
              w = -1
              class(w) = "try-error"
              return(w)
            })

          if (!inherits(result, "try-error")) {

              model_fit  = result
              summary_fit = summary(model_fit)

              intercept[k] = summary_fit$coefficients[1, 1]
              slope[k] = summary_fit$coefficients[2, 1]
              delta = slope[k]

              # One-sided p-value
              if (delta > 0) p[k] = summary_fit$coefficients[2, 4] / 2 else p[k] = 1 - summary_fit$coefficients[2, 4] / 2

          } else {
            convergence[k] = 0
          }

      }

    }

    # Apply a multiplicity adjustment
    if (parameters$mult_test_index >= 1) { 

        p = ExportTradMultAdj(parameters$mult_test_index, p, parameters$weight, parameters$transition)

    }

    # One-sided p-values
    pval_vec = NULL
    for (k in 1:(narms - 1)) pval_vec = c(pval_vec, convergence[k], p[k])
    pval_results = rbind(pval_results, pval_vec)

    # Descriptive statistics for each simulation run
    if (parameters$descriptive_statistics) {
      coef_vec = NULL
      cluster_size_vec = NULL
      for (k in 1:(narms - 1)) {
          
        coef_vec = c(coef_vec, convergence[k], intercept[k], slope[k])          
        cluster_size_vec = c(cluster_size_vec, cluster_size)          

      }                

      coef_results = rbind(coef_results, coef_vec)
      cluster_size_results = c(cluster_size_results, cluster_size_vec)
    }

  }
  # End of simulations

  return(list(pval_results = pval_results,
              coef_results = coef_results,
              cluster_size_results = cluster_size_results))

}