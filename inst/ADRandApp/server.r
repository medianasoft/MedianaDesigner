
shinyServer(function(input, output, session) {

  Simulation = reactive({

    #########################################################

    # Design parameters

    parameters = list()

    endpoint_index = as.numeric(input$endpoint_index)

    direction_index = as.numeric(input$direction_index)
    if (direction_index == 1) parameters$direction = "Higher"
    if (direction_index == 2) parameters$direction = "Lower"

    parameters$dose_levels = as.numeric(input$dose_levels[1,])
    parameters$stage_sample_size = as.numeric(input$stage_sample_size[1,])

    parameters$enrollment_period = as.numeric(input$enrollment_period)
    parameters$enrollment_parameter = as.numeric(input$enrollment_parameter)

    parameters$treatment_period = as.numeric(input$treatment_period)
    parameters$dropout_rate = as.numeric(input$dropout_rate) / 100

    input_mean = as.numeric(input$treatment_assumptions[1,])
    input_sd = as.numeric(input$treatment_assumptions[2,])

    parameters$control_mean = input_mean[1]
    parameters$control_sd = input_sd[1]

    parameters$treatment_mean = input_mean[2:length(input_mean)]
    parameters$treatment_sd = input_sd[2:length(input_sd)]

    parameters$ratio_placebo = as.numeric(input$ratio_placebo_percent) / 100

    parameters$balance = as.numeric(input$balance)

    parameters$exponential_model_parameter = as.numeric(input$exponential_model_parameter)
    parameters$emax_model_parameter = as.numeric(input$emax_model_parameter)
    parameters$logistic_model_parameters = c(
      as.numeric(input$logistic_model_ed50), 
      as.numeric(input$logistic_model_delta)
    )

    parameters$delta = as.numeric(input$delta)
    parameters$alpha = as.numeric(input$alpha)
    parameters$nsims = as.numeric(input$nsims)

    if (endpoint_index == 1) {
      parameters$endpoint_type = "Normal"
    }

    withProgress(message = "Running simulations", value = 1, {

      # Run simulations

      results = ADRand(parameters) 

    })

    # Return the list of results
    results

  })

  output$DownloadResults = downloadHandler(

    filename = function() {
      "Report.docx"
    },

    content = function(file) {

      # Run simulations  
      results = Simulation()

      doc = ADRandReportDoc(results)

      # Save the report
      xfile = paste0(file, ".docx")
      print(doc, target = xfile)          
      file.rename(xfile, file)
    }
  )

  # Create a matrix for entering dose levels 
  output$DoseLevels = renderUI({

      # Trial arms
      narms = as.numeric(input$n_arms)

      trial_arms = rep("Placebo", narms)
      for (i in 2:narms) trial_arms[i] = paste0("Dose ", i - 1)

      value = matrix(0, 1, narms)
      value[1, ] = 100 * (0:(narms - 1))
      rownames(value) = c("Dose")
      colnames(value) = trial_arms

      matrixInput("dose_levels",
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering sample sizes in each trial stage
  output$StageSampleSize = renderUI({

      # Trial stages
      nstages = as.numeric(input$n_stages)

      stages = rep("", nstages)
      for (i in 1:nstages) stages[i] = paste0("Stage ", i)

      value = matrix(0, 1, nstages)
      value[1, ] = rep(50, nstages)
      rownames(value) = c("Sample size")
      colnames(value) = stages

      matrixInput("stage_sample_size",
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering treatment effect assumptions
  output$TreatmentEffectAssumptions = renderUI({

      # Trial arms  
      narms = as.numeric(input$n_arms)

      trial_arms = rep("Placebo", narms)
      for (i in 2:narms) trial_arms[i] = paste0("Dose ", i - 1)

      endpoint_index = as.numeric(input$endpoint_index)  

      if (endpoint_index == 1) {
          value = matrix(0, 2, narms)
          value[1, ] = 5 * (0:(narms - 1))
          value[2, ] = rep(25, narms)
          rownames(value) = c("Mean", "Standard deviation")
      }

      colnames(value) = trial_arms

      matrixInput("treatment_assumptions", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Plots of candidate dose-response models used in the MCPMod method
  output$DoseResponseModels = renderImage({

    filename = tempfile(fileext='.png')

    max_dose_levels = max(input$dose_levels[1, ])

    # Linear model
    coef1 = ComputeDRFunctionParameters(1, 0, 1, max_dose_levels, NA)

    # Exponential model
    coef2 = ComputeDRFunctionParameters(2, 0, 1, max_dose_levels, as.numeric(input$exponential_model_parameter))

    # Emax model
    coef3 = ComputeDRFunctionParameters(3, 0, 1, max_dose_levels, as.numeric(input$emax_model_parameter))

    # Logistic model
    coef4 = ComputeDRFunctionParameters(4, 0, 1, max_dose_levels, as.numeric(c(input$logistic_model_ed50, input$logistic_model_delta)))

    x = seq(from = 0, to = max_dose_levels, length = 100)
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

    width = 800
    height = 400

    xlabels = as.numeric(input$dose_levels[1, ])
    ylabels = seq(from = 0, to = 100, by = 20)
    png(file = filename, width = width, height = height, bg = "transparent")
    plot(x = x, y = y1, xlab="Dose", ylab="Response", xlim = c(0, max_dose_levels), ylim = c(0, 1), type="l", lwd = 2, col = "black") 
    lines(x = x, y = y2, col = "blue", lwd = 2)
    lines(x = x, y = y3, col = "red", lwd = 2)
    lines(x = x, y = y4, col = "darkgreen", lwd = 2)
    dev.off()

    # Return a list containing the filename
    list(src = filename,
         contentType = 'image/png',
         width = width,
         height = height,
         alt = "Candidate dose-response models")
  }, deleteFile = TRUE)   


  output$SampleSizeByStage = renderTable({

    results = Simulation()

    parameters = results$parameters
    simulations = results$simulations

    column_names = c("Stage", "Statistic", "Sample size")

    col1 = NULL
    col2 = NULL
    col3 = NULL

    n_stages = dim(simulations$stage_sample_size)[2]
    for (i in 1:n_stages) {

      col1 = c(col1, paste0("Stage ", i), rep("", 3))
      col2 = c(col2, "Min", "Median", "Mean", "Max")
      col3 = c(col3, round(summary(simulations$stage_sample_size[, i])[c(1, 3, 4, 6)], 1))

    }

    data_frame = data.frame(col1, col2, col3)
    
    colnames(data_frame) = column_names

    data_frame

  })


  output$SampleSizeByTrialArm = renderTable({

    results = Simulation()

    parameters = results$parameters
    simulations = results$simulations

    column_names = c("Dose", "Statistic", "Sample size")

    col1 = NULL
    col2 = NULL
    col3 = NULL

    n_doses = length(parameters$dose_levels)
    # Trial arms
    trial_arms = rep("", n_doses)
    for (i in 1:n_doses) {
        if (i == 1) trial_arms[i] = "Placebo" else trial_arms[i] = paste0("Dose ", i - 1, " (", parameters$dose_levels[i], ")")
    }

    for (i in 1:n_doses) {

      col1 = c(col1, trial_arms[i], rep("", 3))
      col2 = c(col2, "Min", "Median", "Mean", "Max")
      col3 = c(col3, round(summary(simulations$n[, i])[c(1, 3, 4, 6)], 1))

    }

    data_frame = data.frame(col1, col2, col3)

    colnames(data_frame) = column_names

    data_frame

  })


  output$ComparisonOfTraditionalAndAdaptiveDesigns = renderTable({

    results = Simulation()

    parameters = results$parameters
    simulations = results$simulations

    rowMax = function(x) {
      row_max = rep(0, nrow(x))
      for (i in 1:nrow(x)) row_max[i] = max(x[i, ])
      return(row_max)
    }

    column_names = c("Design", "Power (%)") 

    traditional = mean(rowMax(simulations$traditional) >= simulations$mcpmod_value) 
    adaptive = mean(rowMax(simulations$adaptive) >= simulations$mcpmod_value) 

    col1 = c("Traditional", "Adaptive")
    col2 = c(round(100 * traditional, 1), 
             round(100 * adaptive, 1))

    data_frame = data.frame(col1, col2)

    colnames(data_frame) = column_names

    data_frame

  })


  observeEvent(input$jump_to_panel2, {
        updateTabItems(session, "sidebar",
                          selected = "endpoint_parameters")
  })

  observeEvent(input$jump_to_panel3, {
        updateTabItems(session, "sidebar",
                          selected = "interim_parameters")
  })

  observeEvent(input$jump_to_panel4, {
        updateTabItems(session, "sidebar",
                          selected = "general_parameters")
  })

  observeEvent(input$jump_to_panel5, {
        updateTabItems(session, "sidebar",
                          selected = "simulation")
  })

  observeEvent(input$jump_to_panel6, {
        updateTabItems(session, "sidebar",
                          selected = "report")
  })
})
