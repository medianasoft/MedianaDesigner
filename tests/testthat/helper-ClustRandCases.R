globalTestNCores = 1
globalTestNSims  = 50

######################################
# Normal endpoint with 2 arms and fixed cluster sizes
parametersNormal2ArmsFixedClusterSize = list(
    # Endpoint type
    endpoint_type = "Normal",

    # Direction of favorable outcome
    direction = "Lower",

    # Number of completers in the trial arms (control, multiple treatments)
    sample_size = c(100, 100),

    # Cluster scheme
    cluster_scheme = "Fixed",

    # Vector of cluster sizes in the control arm
    control_cluster_size = rep(10, 10),

    # Vector of cluster sizes in the treatment arm
    treatment_cluster_size = rep(10, 10),

    # Mean in the control arm 
    control_mean = 1,

    # Intra-cluster correlation coefficient in the control arm 
    control_icc = 0.6,

    # Between-cluster standard deviation in the control arm 
    control_between_cluster_sd = 1.2,

    # Mean in the treatment arm
    treatment_mean = 0.3,

    # Intra-cluster correlation coefficient in the treatment arm
    treatment_icc = 0.6,

    # Between-cluster standard deviation in the treatment arm 
    treatment_between_cluster_sd = 1.2,

    # Data analysis method (generalized estimating equations (GEE) or generalized linear mixed effects model (GLMEM))
    method_type = "GEE",

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = globalTestNSims,

    # Number of cores for parallel calculations
    ncores = globalTestNCores,

    # Compute descriptive statistics (arm-specific effects, ICC, cluster sizes) for each simulation run
    descriptive_statistics = TRUE
)

parametersNormal2ArmsRandomClusterSize = parametersNormal2ArmsFixedClusterSize
parametersNormal2ArmsRandomClusterSize$control_cluster_proportion = c(0.5, 0.5)
parametersNormal2ArmsRandomClusterSize$treatment_cluster_proportion = c(0.5, 0.5)
parametersNormal2ArmsRandomClusterSize$cluster_scheme = "Random"

parametersNormal2ArmsRandomClusterSizeGLMEM = parametersNormal2ArmsRandomClusterSize
parametersNormal2ArmsRandomClusterSizeGLMEM$method_type = "GLMEM"

######################################
# Normal endpoint with 3 arms and fixed cluster sizes
parametersNormal3ArmsFixedClusterSize = list(
    # Endpoint type
    endpoint_type = "Normal",

    # Direction of favorable outcome
    direction = "Higher",

    # Number of completers in the trial arms (control, multiple treatments)
    sample_size = c(100, 100, 100),

    # Cluster scheme
    cluster_scheme = "Fixed",

    # Vector of cluster sizes in the control arm
    control_cluster_size = rep(10, 10),

    # Matrix of cluster sizes in the treatment arms 
    treatment_cluster_size = rbind(rep(10, 10),
                                   rep(10, 10)),

    # Mean in the control arm 
    control_mean = 0,

    # Intra-cluster correlation coefficient in the control arm 
    control_icc = 0.5,

    # Between-cluster standard deviation in the control arm 
    control_between_cluster_sd = 1.2,

    # Mean in the treatment arms 
    treatment_mean = c(1, 1.1),

    # Intra-cluster correlation coefficient in the treatment arms
    treatment_icc = c(0.5, 0.5),

    # Between-cluster standard deviation in the treatment arms 
    treatment_between_cluster_sd = c(1.2, 1.2),

    # Data analysis method (generalized estimating equations (GEE) or generalized linear mixed effects model (GLMEM))
    method_type = "GLMEM",

    # Multiple testing procedure
    mult_test = "Bonferroni",

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = globalTestNSims,

    # Number of cores for parallel calculations
    ncores = globalTestNCores,

    # Compute descriptive statistics (arm-specific effects, ICC, cluster sizes) for each simulation run
    descriptive_statistics = TRUE
)

######################################
# Normal endpoint with 3 arms and random cluster sizes
parametersNormal3ArmsRandomClusterSize = list(
    # Endpoint type
    endpoint_type = "Normal",

    # Direction of favorable outcome
    direction = "Higher",

    # Number of completers in the trial arms (control, multiple treatments)
    sample_size = c(100, 100, 100),

    # Cluster scheme 
    cluster_scheme = "Random",

    # Vector of relative cluster sizes in the control arm
    control_cluster_proportion = rep(0.1, 10),

    # Matrix of cluster sizes in the treatment arms 
    treatment_cluster_proportion = rbind(rep(0.1, 10),
                                         rep(0.1, 10)),

    # Mean in the control arm 
    control_mean = 0,

    # Intra-cluster correlation coefficient in the control arm 
    control_icc = 0.5,

    # Between-cluster standard deviation in the control arm 
    control_between_cluster_sd = 1.2,

    # Mean in the treatment arms 
    treatment_mean = c(1, 1.1),

    # Intra-cluster correlation coefficient in the treatment arms
    treatment_icc = c(0.5, 0.5),

    # Between-cluster standard deviation in the treatment arms 
    treatment_between_cluster_sd = c(1.2, 1.2),

    # Data analysis method (generalized estimating equations (GEE) or generalized linear mixed effects model (GLMEM))
    method_type = "GEE",

    # Multiple testing procedure
    mult_test = "Hochberg",

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = globalTestNSims,

    # Number of cores for parallel calculations
    ncores = globalTestNCores,

    # Compute descriptive statistics (arm-specific effects, ICC, cluster sizes) for each simulation run
    descriptive_statistics = TRUE
)

parametersNormal3ArmsRandomClusterSizeGLMEM = parametersNormal3ArmsRandomClusterSize
parametersNormal3ArmsRandomClusterSizeGLMEM$method_type = "GLMEM"

############################
# Binary endpoint with 2 arms and fixed cluster sizes
parametersBinary2ArmsFixedClusterSize = list(
    # Endpoint type
    endpoint_type = "Binary",

    # Direction of favorable outcome
    direction = "Lower",

    # Number of completers in the trial arms (control, multiple treatments)
    sample_size = c(100, 100),

    # Cluster scheme (fixed or random cluster sizes)
    cluster_scheme = "Fixed",

    # Vector of cluster sizes in the control arm 
    control_cluster_size = rep(10, 10),

    # Vector of cluster sizes in the treatment arm 
    treatment_cluster_size = rep(10, 10),

    # Response rate in the control arm 
    control_rate = 0.6,

    # Intracluster correlation coefficient in the control arm 
    control_icc = 0.3,

    # Response rate in the treatment arms 
    treatment_rate = 0.3,

    # Intracluster correlation coefficient in the treatment arms
    treatment_icc = 0.3,

    # Data analysis method (generalized estimating equations (GEE) or generalized linear mixed effects model (GLMEM))
    method_type = "GLMEM",

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = globalTestNSims,

    # Number of cores for parallel calculations
    ncores = globalTestNCores,

    # Compute descriptive statistics (arm-specific effects, ICC, cluster sizes) for each simulation run
    descriptive_statistics = TRUE
)

############################
# Binary endpoint with 4 arms and random cluster sizes
parametersBinary4ArmsRandomClusterSize = list(
    # Endpoint type
    endpoint_type = "Binary",

    # Direction of favorable outcome
    direction = "Lower",

    # Number of completers in the trial arms (control, multiple treatments)
    sample_size = c(100, 100, 100, 100),

    # Cluster scheme 
    cluster_scheme = "Random",

    # Vector of relative cluster sizes in the control arm
    control_cluster_proportion = rep(0.1, 10),

    # Matrix of cluster sizes in the treatment arms 
    treatment_cluster_proportion = rbind(rep(0.1, 10),
                                         rep(0.1, 10),
                                         rep(0.1, 10)),

    # Response rate in the control arm 
    control_rate = 0.6,

    # Intracluster correlation coefficient in the control arm 
    control_icc = 0.3,

    # Response rate in the treatment arms 
    treatment_rate = c(0.3, 0.25, 0.2),

    # Intracluster correlation coefficient in the treatment arms
    treatment_icc = c(0.3, 0.3, 0.3),

    # Data analysis method (generalized estimating equations (GEE) or generalized linear mixed effects model (GLMEM))
    method_type = "GEE",

    # Multiple testing procedure
    mult_test = "Hochberg",

    # One-sided alpha level
    alpha = 0.025,

    # Number of simulations
    nsims = globalTestNSims,

    # Number of cores for parallel calculations
    #ncores = globalTestNCores      # Use default value - 1 core

    # Compute descriptive statistics (arm-specific effects, ICC, cluster sizes) for each simulation run
    descriptive_statistics = TRUE
)
