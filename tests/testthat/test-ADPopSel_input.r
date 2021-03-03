# Test the input of the ADPopSel function

# Normal case
normalCase = list(

  # Number of enrolled patients in each trial arm
  sample_size = c(120, 120),

  # Prevalence
  prevalence = 0.4,

  # Primary endpoint's type
  endpoint_type = "Normal",

  # Mean and SD in the control arm 
  control_mean = c(0, 0),
  control_sd = c(1, 1),
  
  # Mean and SD in the treatment arm 
  treatment_mean = c(0.2, 0.4),
  treatment_sd = c(1, 1),
  
  # Information fractions at IA1, IA2, FA
  info_frac = c(0.4, 0.6, 1.0),

  # Futility threshold
  futility_threshold = 0.1,

  # Influence threshold for population selection
  influence = 0.1,

  # Interaction threshold for population selection
  interaction = 1.3,

  # Dropout rate at the end of the treatment period
  dropout_rate = 0.05,

  # One-sided Type I error rate
  alpha = 0.025,

  # Number of simulations
  nsims = 1000
)

# Run simulations to compute key characteristics

test_that("Success run ADPopSel", {

  # Success run
  results = ADPopSel(normalCase)
  expect_is(results, "ADPopSelResults")

  expect_type(  results$sim_results, "double")
  expect_length(results$sim_results, 19000)
  
  expect_type(    results$sim_summary, "list")
  expect_true(abs(results$sim_summary$futility - 0.2) < 0.2)
  
  expect_is(      results$sim_summary$trad_power, "numeric")
  expect_length(  results$sim_summary$trad_power, 2)
  expect_true(abs(results$sim_summary$trad_power[1] - 0.5) < 0.1)
  expect_true(abs(results$sim_summary$trad_power[2] - 0.4) < 0.1)

  expect_is(      results$sim_summary$ad_power, "numeric")
  expect_length(  results$sim_summary$ad_power, 3)
  expect_true(abs(results$sim_summary$ad_power[1] - 0.4) < 0.1)
  expect_true(abs(results$sim_summary$ad_power[2] - 0.3) < 0.1)
  expect_true(abs(results$sim_summary$ad_power[3] - 0.5) < 0.1)

  expect_is(      results$sim_summary$hypothesis_selection, "numeric")
  expect_length(  results$sim_summary$hypothesis_selection, 3)
  expect_true(abs(results$sim_summary$hypothesis_selection[1] - 0.3) < 0.1)
  expect_true(abs(results$sim_summary$hypothesis_selection[2] - 0.3) < 0.1)
  expect_true(abs(results$sim_summary$hypothesis_selection[3] - 0.3) < 0.1)

  expect_is(      results$sim_summary$look_time, "numeric")
  expect_length(  results$sim_summary$look_time, 4)
  expect_true(abs(results$sim_summary$look_time[1] - 0) < 0.1)
  expect_true(abs(results$sim_summary$look_time[2] - 0) < 0.1)
  expect_true(is.nan(results$sim_summary$look_time[3]))
  expect_true(is.nan(results$sim_summary$look_time[4]))

  # Check for report generation
  ADPopSelReport(results, tempfile(fileext = ".docx"))
})

test_that("Input parameters errors check ADPopSel", {
  # Errors check
  expect_error(
    ADPopSel(
      c("Not a list")
    ),
    info = "Checking for wrong parameters collection type"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        # missing
        #endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing endpoint type"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        # wrong
        endpoint_type = "WrongType", # normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for wrong endpoint type"
  )

  expect_error(
    ADPopSel(
      list(
        # missing
        #sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing number of enrolled patients"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        # missing
        #prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing prevalence of biomarker-positive patients"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # missing
        #info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing information fractions at IA1, IA2, FA"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong
        info_frac = c(0.8, 0.6, 1.0), # normalCase$info_frac = c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for error information fractions (IA1 must be < IA2)"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong
        info_frac = c(0.4, 1.6, 1.0), # normalCase$info_frac = c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for error information fractions (IA2 must be < FA)"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        # wrong
        info_frac = c(0.4, 0.6, 1.1), # normalCase$info_frac = c(0.4, 0.6, 1.0),
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for error information fractions (FA must be 1)"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        # missing
        #futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing futility threshold at IA1"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        # missing
        #influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing influence threshold at IA2"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        # missing
        #interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing interaction threshold at IA2"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        # missing
        #control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effects in the control arm"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        # missing
        #treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing mean effects in the treatment arm"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        # missing
        #control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviations in the control arm"
  )

  expect_error(
    ADPopSel(
      list(
        sample_size = normalCase$sample_size,
        prevalence = normalCase$prevalence,
        endpoint_type = normalCase$endpoint_type,
        control_mean = normalCase$control_mean,
        control_sd = normalCase$control_sd,
        treatment_mean = normalCase$treatment_mean,
        # missing
        #treatment_sd = normalCase$treatment_sd,
        info_frac = normalCase$info_frac,
        futility_threshold = normalCase$futility_threshold,
        influence = normalCase$influence,
        interaction = normalCase$interaction,
        dropout_rate = normalCase$dropout_rate,
        alpha = normalCase$alpha,
        nsims = normalCase$nsims
      )
    ),
    info = "Checking for missing standard deviations in the treatment arm"
  )

})
