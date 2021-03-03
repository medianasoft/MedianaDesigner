shinyServer(function(input, output, session) {

  Simulation = reactive({

    #########################################################

    # Design parameters

    parameters = list()

    endpoint_index = as.numeric(input$endpoint_index)

    parameters$sample_size = as.numeric(input$n)   

    assumptions = input$assumptions

    parameters$info_frac = c(as.numeric(input$info_frac1 / 100), as.numeric(input$info_frac2 / 100), 1, as.numeric(input$info_frac3 / 100))   

    if (endpoint_index == 1) {
      parameters$endpoint_type = "Normal"
      parameters$control_mean = assumptions[1, 1]
      parameters$treatment_mean = as.numeric(assumptions[1, 2])
      parameters$control_sd = assumptions[2, 1]
      parameters$treatment_sd = as.numeric(assumptions[2, 2])
    }

    if (endpoint_index == 2) {
      parameters$endpoint_type = "Binary"
      parameters$control_rate = assumptions[1, 1] / 100
      parameters$treatment_rate = as.numeric(assumptions[1, 2]  / 100)
    }

    if (endpoint_index == 3) {
      parameters$endpoint_type = "Time-to-event"  
      parameters$control_time = assumptions[1, 1]
      parameters$treatment_time = as.numeric(assumptions[1, 2])          
      parameters$event_count = as.numeric(input$event_count)
      parameters$enrollment_period = as.numeric(input$enrollment_period)   
      parameters$enrollment_parameter = as.numeric(input$enrollment_parameter)   
    }
 
    parameters$dropout_rate = as.numeric(input$dropout_rate / 100)

    parameters$futility_threshold = as.numeric(input$futility_threshold / 100)
    parameters$promising_interval = c(as.numeric(input$lower_limit / 100), as.numeric(input$upper_limit / 100))
    parameters$target_power = as.numeric(input$target_power / 100)

    parameters$alpha = as.numeric(input$alpha)
    parameters$nsims = as.numeric(input$nsims)

    #########################################################

    withProgress(message = "Running simulations", value = 1, {

      # Run simulations

      results = ADSSMod(parameters) 

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

      narms = 2  

      # Trial arms  
      trial_arms = c("Control", "Treatment")

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

      narms = 2 
      endpoint_index = as.numeric(input$endpoint_index)  

      # Trial arms  
      trial_arms = "Control"
      if (narms >= 3) {
        for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
      } else {
        trial_arms = c(trial_arms, "Treatment")
      }

      if (endpoint_index == 1) {
          value = matrix(0, 2, narms)
          value[1, ] = c(0, rep(0.3, narms - 1))
          value[2, ] = rep(1, narms)
          rownames(value) = c("Mean", "Standard deviation")
      }

      if (endpoint_index == 2) {
          value = matrix(0, 1, narms)
          value[1, ] = c(10, rep(30, narms - 1))
          rownames(value) = c("Rate (%)")
      }

      if (endpoint_index == 3) {
          value = matrix(0, 1, narms)
          value[1, ] = c(10, rep(14, narms - 1))
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


  output$OutcomeProbabilities = renderTable({

    results = Simulation()

    parameters = results$parameters
    sim_results = results$sim_results
    sim_summary = results$sim_summary

    column_names = c("Parameter", "Value")

    col1 = c("Probability of stopping for futility at Interim analysis 1 (%)",
             "Probability of increasing the number of events at Interim analysis 2 (%)",
             "Traditional design: Power (%)",
             "Adaptive design: Power (%)") 
    col2 = c(sprintf("%0.1f", 100 * sim_summary$futility),
             sprintf("%0.1f", 100 * sim_summary$increase),
             sprintf("%0.1f", 100 * sim_summary$trad_power),
             sprintf("%0.1f", 100 * sim_summary$ad_power))
    
    data_frame = data.frame(col1, col2)

    colnames(data_frame) = column_names

    data_frame

  })

  output$Comparison = renderTable({

    results = Simulation()

    parameters = results$parameters
    sim_results = results$sim_results
    sim_summary = results$sim_summary

    column_names = c("Interval", "Design", "Power (%)")

    col1 = c("Unfavorable interval", "", "Promising interval", "", "Favorable interval", "")
    col2 = rep(c("Traditional design", "Adaptive design"), 3)
    col3 = c(sprintf("%0.1f", 100 * sim_summary$trad_under),
             sprintf("%0.1f", 100 * sim_summary$ad_under),
             sprintf("%0.1f", 100 * sim_summary$trad_prom),
             sprintf("%0.1f", 100 * sim_summary$ad_prom),
             sprintf("%0.1f", 100 * sim_summary$trad_over),
             sprintf("%0.1f", 100 * sim_summary$ad_over))
    
    data_frame = data.frame(col1, col2, col3)
    
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
