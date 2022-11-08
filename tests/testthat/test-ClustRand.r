test_that("Success run ClustRand with normal endpoint", {
    # Run simulations to compute key characteristics of the cluster-randomized design
    results = ClustRand(parametersNormal2ArmsFixedClusterSize)

    expect_s3_class(results, "ClustRandResults")
    expect_equal(ncol(results$pval_results), 4)
    expect_equal(nrow(results$pval_results), results$parameters$nsims)
    expect_equal(ncol(results$coef_results), 3)
    expect_equal(nrow(results$coef_results), results$parameters$nsims)

    sim_summary = results$sim_summary
    expect_equivalent(sim_summary$power_sandwich, 0.3, tolerance=0.01)
    expect_equivalent(sim_summary$power_kc, 0.24, tolerance=0.01)
    expect_equivalent(sim_summary$power_md, 0.14, tolerance=0.01)

    # Check for report generation
    GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Success run ClustRand with binary endpoint", {
    # Run simulations to compute key characteristics of the cluster-randomized design
    results = ClustRand(parametersBinary2ArmsFixedClusterSize)

    expect_s3_class(results, "ClustRandResults")

    sim_summary = results$sim_summary
    expect_equivalent(sim_summary$power, 0.66, tolerance=0.02)

    # Check for report generation
    GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Success run ClustRand with binary endpoint and some parameter modifications", {
    # Run simulations to compute key characteristics of the cluster-randomized design
    parameters = parametersBinary2ArmsFixedClusterSize
    # parameters$method_type = "glmem"
    parameters$direction <- NULL
    parameters$alpha <- NULL
    results = ClustRand(parameters)

    expect_s3_class(results, "ClustRandResults")

    sim_summary = results$sim_summary
    expect_equal(sim_summary$power, 0.0, tolerance=0.01)

    # Check for report generation
    GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Success runs ClustRand for additional cases", {
    results = ClustRand(parametersNormal2ArmsRandomClusterSize)
    expect_s3_class(results, "ClustRandResults")
    GenerateReport(results, tempfile(fileext = ".docx"))

    results = ClustRand(parametersNormal2ArmsRandomClusterSizeGLMEM)
    expect_s3_class(results, "ClustRandResults")
    GenerateReport(results, tempfile(fileext = ".docx"))

    results = ClustRand(parametersNormal3ArmsFixedClusterSize)
    expect_s3_class(results, "ClustRandResults")
    GenerateReport(results, tempfile(fileext = ".docx"))

    results = ClustRand(parametersNormal3ArmsRandomClusterSize)
    expect_s3_class(results, "ClustRandResults")
    GenerateReport(results, tempfile(fileext = ".docx"))

    results = ClustRand(parametersNormal3ArmsRandomClusterSizeGLMEM)
    expect_s3_class(results, "ClustRandResults")
    GenerateReport(results, tempfile(fileext = ".docx"))

    results = ClustRand(parametersBinary4ArmsRandomClusterSize)
    expect_s3_class(results, "ClustRandResults")
    GenerateReport(results, tempfile(fileext = ".docx"))
})

test_that("Parameters range checks for ClustRand", {

    expect_error({
        ClustRand(c("Not a list"))
    })

    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'random_seed',
        checkMissing = FALSE, 
        checkWrong = "Not int", 
        checkSize = FALSE, 
        checkMin = 0, 
        checkMax = 100001
    )

    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'endpoint_type',
        checkMissing = TRUE, 
        checkWrong = "Wrong",
        checkSize = FALSE, 
        checkMin = NA, 
        checkMax = NA
    )

    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'method_type',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = NA,
        checkMax = NA
    )

    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'direction',
        checkMissing = FALSE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = NA,
        checkMax = NA
    )

    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'cluster_scheme',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = NA,
        checkMax = NA
    )

    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'control_cluster_size',
        checkMissing = TRUE,
        checkWrong = rep(11, 10),
        checkSize = FALSE,
        checkMin = 4,
        checkMax = 101
    )

    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'treatment_cluster_size',
        checkMissing = TRUE,
        checkWrong = rep(11, 10),
        checkSize = FALSE,
        checkMin = 4,
        checkMax = 101
    )

    testParameterErrors(parametersNormal3ArmsFixedClusterSize, 
        'treatment_cluster_size',
        checkMissing = TRUE,
        checkWrong = matrix(c(50,50,50,50), nrow=4, ncol=1),
        checkSize = TRUE,
        checkMin = matrix(c(4,50,50,50), nrow=2, ncol=2),
        checkMax = matrix(c(50,50,50,100), nrow=2, ncol=2)
    )

    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'sample_size',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = TRUE,
        checkMin = c(0, 0),
        checkMax = c(1001, 1001)
    )

    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'control_cluster_proportion',
        checkMissing = TRUE,
        checkWrong = c(0.5, 0.6),   # Sum is not equal 1
        checkSize = TRUE,
        checkMin = c(0, 0),
        checkMax = c(1, 1)
    )

    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'treatment_cluster_proportion',
        checkMissing = TRUE,
        checkWrong = c(0.5, 0.6),   # Sum is not equal 1
        checkSize = TRUE,
        checkMin = c(0, 0),
        checkMax = c(1, 1)
    )

    testParameterErrors(parametersNormal3ArmsRandomClusterSize, 
        'treatment_cluster_proportion',
        checkWrong = c(0.5, 0.5, 0.5, 0.5)
    )
    testParameterErrors(parametersNormal3ArmsRandomClusterSize, 
        'treatment_cluster_proportion',
        checkWrong = matrix(c(0.5, 0.5, 0.5, 0.5), nrow=4, ncol=1)
    )
    testParameterErrors(parametersNormal3ArmsRandomClusterSize, 
        'treatment_cluster_proportion',
        checkMissing = TRUE,
        checkWrong = matrix(c(0.5, 0.5, 0.5, 0.6), nrow=2, ncol=2),   # Sum is not equal 1
        checkSize = TRUE,
        checkMin = matrix(c(0, 0.5, 0.5, 0.5), nrow=2, ncol=2),
        checkMax = matrix(c(0.5, 0.5, 0.5, 1), nrow=2, ncol=2)
    )

    testParameterErrors(parametersNormal3ArmsRandomClusterSize, 
        'mult_test',
        checkMissing = TRUE,
        checkWrong = "Wrong"
    )

    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'nsims',
        checkMissing = FALSE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0,
        checkMax = 100000
    )

    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'alpha',
        checkMissing = FALSE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0.001,
        checkMax = 0.5
    )

    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'control_icc',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0,
        checkMax = 1
    )

    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'treatment_icc',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0,
        checkMax = 1
    )

    # control_between_cluster_sd > 0
    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'control_between_cluster_sd',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0,
        checkMax = NA
    )

    # treatment_between_cluster_sd > 0
    testParameterErrors(parametersNormal2ArmsRandomClusterSize, 
        'treatment_between_cluster_sd',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0,
        checkMax = NA
    )

    # Normal : 
    # control_mean not null, type=double
    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'control_mean',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = NA,
        checkMax = NA
    )

    # treatment_mean same checks
    testParameterErrors(parametersNormal2ArmsFixedClusterSize, 
        'treatment_mean',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = NA,
        checkMax = NA
    )

    # Binary :
    # control_rate (0, 1)
    testParameterErrors(parametersBinary2ArmsFixedClusterSize, 
        'control_rate',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0,
        checkMax = 1
    )
    
    # treatment_rate (0, 1)
    testParameterErrors(parametersBinary2ArmsFixedClusterSize, 
        'treatment_rate',
        checkMissing = TRUE,
        checkWrong = "Wrong",
        checkSize = FALSE,
        checkMin = 0,
        checkMax = 1
    )

    # descriptive_statistics
    testParameterErrors(parametersBinary2ArmsFixedClusterSize, 
        'descriptive_statistics',
        checkWrong = "Wrong"
    )

    # Report generation input error check
    expect_error(
        ClustRandReportDoc(list())
    )

})