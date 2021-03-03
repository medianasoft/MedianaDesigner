shinyServer(function(input, output, session) {

  Simulation = reactive({

    #########################################################

    # Design parameters

    parameters = list()

    endpoint_index = as.numeric(input$endpoint_index)
    narms = as.numeric(input$narms)
    parameters$narms = narms

    parameters$sample_size = as.numeric(input$sample_size)   

    assumptions = input$assumptions

    if (endpoint_index == 1) {
      parameters$endpoint_type = "Normal"
      parameters$control_mean = assumptions[1, 1]
      parameters$treatment_mean = as.numeric(assumptions[1, 2:narms])
      parameters$control_sd = assumptions[2, 1]
      parameters$treatment_sd = as.numeric(assumptions[2, 2:narms])
    }

    if (endpoint_index == 2) {
      parameters$endpoint_type = "Binary"
      parameters$control_rate = assumptions[1, 1] / 100
      parameters$treatment_rate = as.numeric(assumptions[1, 2:narms]  / 100)
    }

    if (endpoint_index == 3) {
      parameters$endpoint_type = "Time-to-event"  
      parameters$control_time = assumptions[1, 1]
      parameters$treatment_time = as.numeric(assumptions[1, 2:narms])          
      parameters$event_count = as.numeric(input$event_count)   
      parameters$enrollment_period = as.numeric(input$enrollment_period)   
      parameters$enrollment_parameter = as.numeric(input$enrollment_parameter)   
    }
 
    parameters$dropout_rate = as.numeric(input$dropout_rate / 100)
    parameters$info_frac = as.numeric(input$info_frac / 100)

    parameters$alpha = as.numeric(input$alpha)
    parameters$nsims = as.numeric(input$nsims)

    #########################################################

    withProgress(message = "Running simulations", value = 1, {

      # Run simulations

      results = FutRule(parameters) 

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

    narms = as.numeric(input$narms)  

    # Trial arms  
    trial_arms = "Control"
    if (narms >= 3) {
      for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
    } else {
      trial_arms = c(trial_arms, "Treatment")
    }

    value = matrix(0, 1, narms)
    value[1, ] = rep(50, narms)

    rownames(value) = "Sample size"
    colnames(value) = trial_arms

    matrixInput("sample_size", 
      class = "numeric",
      rows = list(names = TRUE),
      cols = list(names = TRUE),
      value = value
    )

  })

  # Create a matrix for entering treatment effect assumptions
  output$TreatmentEffectAssumptions = renderUI({

    narms = as.numeric(input$narms)  
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

  output$Sensitivity = renderImage({

    results = Simulation()

    sim_summary = results$sim_summary

    sensitivity = sim_summary$sensitivity
    cp_threshold = sim_summary$cp_threshold

    width = 600
    height = 500

    filename = tempfile(fileext='.png')
    png(file = filename, width = width, height = height, bg = "transparent")

    plot(x = 100 * cp_threshold, y = 100 * sensitivity, xlab="",ylab="", xlim=c(0, 100), ylim=c(0, 100), axes=FALSE, col="red", lwd = 2, type="l")
    box()
    xlabels = seq(from = 0, to = 100, by = 20)
    ylabels = seq(from = 0, to = 100, by = 20)
    axis(1,at=xlabels, labels=xlabels, tck = 0.02, mgp = c(0, 0.4, 0))
    axis(2,at=ylabels, labels=ylabels, tck = 0.02, mgp = c(0, 0.4, 0))
    mtext("Futility threshold (%)", side=1, line=1.5)
    mtext("Sensitivity rate (%)", side=2, line=1.75)
    dev.off()

    # Return a list containing the filename
    list(src = filename,
         contentType = 'image/png',
         width = width,
         height = height,
         alt = "Sensitivity rate")

  }, deleteFile = TRUE) 

  output$Specificity = renderImage({

    results = Simulation()

    sim_summary = results$sim_summary

    specificity = sim_summary$specificity
    cp_threshold = sim_summary$cp_threshold

    width = 600
    height = 500

    filename = tempfile(fileext='.png')
    png(file = filename, width = width, height = height, bg = "transparent")

    plot(x = 100 * cp_threshold, y = 100 * specificity, xlab="",ylab="", xlim=c(0, 100), ylim=c(0, 100), axes=FALSE, col="red", lwd = 2, type="l")
    box()
    xlabels = seq(from = 0, to = 100, by = 20)
    ylabels = seq(from = 0, to = 100, by = 20)
    axis(1,at=xlabels, labels=xlabels, tck = 0.02, mgp = c(0, 0.4, 0))
    axis(2,at=ylabels, labels=ylabels, tck = 0.02, mgp = c(0, 0.4, 0))
    mtext("Futility threshold (%)", side=1, line=1.5)
    mtext("Specificity rate (%)", side=2, line=1.75)
    dev.off()

    # Return a list containing the filename
    list(src = filename,
         contentType = 'image/png',
         width = width,
         height = height,
         alt = "Specificity rate")

  }, deleteFile = TRUE) 

  output$Accuracy = renderImage({

    results = Simulation()

    sim_summary = results$sim_summary

    accuracy = sim_summary$accuracy
    cp_threshold = sim_summary$cp_threshold

    index = which.max(accuracy)
    optimal_point = round(100 * cp_threshold[index], 1)
    level = 0.95 * max(accuracy)
    zone = cp_threshold[accuracy >= level]
    optimal_lower = round(100 * min(zone), 1)
    optimal_upper = round(100 * max(zone), 1)

    width = 600
    height = 500

    filename = tempfile(fileext='.png')
    png(file = filename, width = width, height = height, bg = "transparent")

    plot(x = 100 * cp_threshold, y = 100 * accuracy, xlab="",ylab="", xlim=c(0, 100), ylim=c(50, 100), axes=FALSE, col="red", lwd = 2, type="l")
    box()
    xlabels = seq(from = 0, to = 100, by = 20)
    ylabels = seq(from = 50, to = 100, by = 10)
    axis(1,at=xlabels, labels=xlabels, tck = 0.02, mgp = c(0, 0.4, 0))
    axis(2,at=ylabels, labels=ylabels, tck = 0.02, mgp = c(0, 0.4, 0))
    abline(v = optimal_point, lty = "solid")
    abline(v = c(optimal_lower, optimal_upper), lty = "dashed")
    mtext("Futility threshold (%)", side=1, line=1.5)
    mtext("Accuracy rate (%)", side=2, line=1.75)
    dev.off()

    # Return a list containing the filename
    list(src = filename,
         contentType = 'image/png',
         width = width,
         height = height,
         alt = "Accuracy rate")

  }, 
  deleteFile = TRUE) 

  output$OptimalThreshold = renderTable({

    results = Simulation()

    sim_summary = results$sim_summary

    accuracy = sim_summary$accuracy
    cp_threshold = sim_summary$cp_threshold

    index = which.max(accuracy)
    optimal_point = round(100 * cp_threshold[index], 1)
    level = 0.95 * max(accuracy)
    zone = cp_threshold[accuracy >= level]
    optimal_lower = round(100 * min(zone), 1)
    optimal_upper = round(100 * max(zone), 1)

    column_names = c("Parameter", "Value")

    col1 = c("Optimal futility threshold (%)", "95% optimal interval (%)")
    col2 = c(optimal_point, paste0("(", optimal_lower, ", ", optimal_upper, ")"))
    
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
                          selected = "general_parameters")
  })

  observeEvent(input$jump_to_panel4, {
        updateTabItems(session, "sidebar",
                          selected = "simulation")
  })

  observeEvent(input$jump_to_panel5, {
        updateTabItems(session, "sidebar",
                          selected = "report")
  })

})
