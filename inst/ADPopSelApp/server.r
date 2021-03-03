
shinyServer(function(input, output, session) {

  Simulation = reactive({

    #########################################################

    # Design parameters

    parameters = list()

    endpoint_index = as.numeric(input$endpoint_index)

    parameters$sample_size = as.numeric(input$n)   

    assumptions = input$assumptions

    parameters$info_frac = c(as.numeric(input$info_frac1 / 100), as.numeric(input$info_frac2 / 100), 1)   

    assumptions = input$assumptions

    if (endpoint_index == 1) {
      parameters$endpoint_type = "Normal"
      parameters$control_mean = as.numeric(assumptions[1, 1:2])
      parameters$treatment_mean = as.numeric(assumptions[1, 3:4])
      parameters$control_sd = as.numeric(assumptions[2, 1:2])
      parameters$treatment_sd = as.numeric(assumptions[2, 3:4])
    }

    if (endpoint_index == 2) {
      parameters$endpoint_type = "Binary"
      parameters$control_rate = as.numeric(assumptions[1, 1:2] / 100)
      parameters$treatment_rate = as.numeric(assumptions[1, 3:4]  / 100)
    }

    if (endpoint_index == 3) {
      parameters$endpoint_type = "Time-to-event"  
      parameters$control_time = as.numeric(assumptions[1, 1:2])
      parameters$treatment_time = as.numeric(assumptions[1, 3:4])          
      parameters$event_count = c(as.numeric(input$event_count1), as.numeric(input$event_count2))   
      parameters$enrollment_period = as.numeric(input$enrollment_period)   
      parameters$enrollment_parameter = as.numeric(input$enrollment_parameter)   
    }
 
    parameters$dropout_rate = as.numeric(input$dropout_rate / 100)

    parameters$futility_threshold = as.numeric(input$futility_threshold / 100)
    parameters$prevalence = as.numeric(input$prevalence / 100)

    parameters$influence = as.numeric(input$influence)
    parameters$interaction = as.numeric(input$interaction)

    parameters$alpha = as.numeric(input$alpha)
    parameters$nsims = as.numeric(input$nsims)

    #########################################################

    withProgress(message = "Running simulations", value = 1, {

      # Run simulations

      results = ADPopSel(parameters) 

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

      doc = ReportDoc(results)

      # Save the report
      xfile = paste0(file, ".docx")
      print(doc, target = xfile)          
      file.rename(xfile, file)
    }
  )

  # Create a matrix for entering sample sizes
  output$SampleSize = renderUI({

      # Trial arms  
      trial_arms = c("Control", "Treatment")

      narms = 2

      value = matrix(0, 1, narms)
      value[1, ] = rep(120, narms)

      rownames(value) = "Sample size"
      colnames(value) = trial_arms

      matrixInput("n", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering treatment effect assumptions
  output$TreatmentEffectAssumptions = renderUI({

      # Trial arms  
      trial_arms = c("Control arm (Biomarker-negative)", 
                     "Control arm (Biomarker-positive)",
                     "Treatment arm (Biomarker-negative)",
                     "Treatment arm (Biomarker-positive)")

      narms = 4

      endpoint_index = as.numeric(input$endpoint_index)  

      if (endpoint_index == 1) {
          value = matrix(0, 2, narms)
          value[1, ] = c(0, 0, 0.2, 0.4)
          value[2, ] = rep(1, narms)
          rownames(value) = c("Mean", "Standard deviation")
      }

      if (endpoint_index == 2) {
          value = matrix(0, 1, narms)
          value[1, ] = c(10, 10, 20, 40)
          rownames(value) = c("Rate (%)")
      }

      if (endpoint_index == 3) {
          value = matrix(0, 1, narms)
          value[1, ] = c(10, 10, 12, 14)
          rownames(value) = c("Median time")
      }

      colnames(value) = trial_arms

      matrixInput("assumptions", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })


  output$FutilityStopping = renderTable({

    results = Simulation()

    parameters = results$parameters
    sim_summary = results$sim_summary

    column_names = c("Outcome", "Probability (%)")

    col1 = c("Futility stopping at Interim analysis 1",
             "Only OP is selected at Interim analysis 2",
             "Only BPP is selected at Interim analysis 2",
             "Both populations are selected at Interim analysis 2")
    col2 = round(100 * c(sim_summary$futility, sim_summary$hypothesis_selection), 1)

    data_frame = data.frame(col1, col2)
    
    colnames(data_frame) = column_names

    data_frame

  })


  output$TraditionalDesign = renderTable({

    results = Simulation()

    parameters = results$parameters
    sim_summary = results$sim_summary

    column_names = c("Population", "Power (%)")

    col1 = c("Overall population")
    col2 = round(100 * sim_summary$trad_power[1], 1)

    data_frame = data.frame(col1, col2)

    colnames(data_frame) = column_names

    data_frame

  })

  output$AdaptiveDesign = renderTable({

    results = Simulation()

    parameters = results$parameters
    sim_summary = results$sim_summary

    column_names = c("Population", "Power (%)")

    col1 = c("Overall population", "Biomarker-positive population", "Either population")
    col2 = round(100 * sim_summary$ad_power, 1)

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
