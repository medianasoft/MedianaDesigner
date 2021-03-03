shinyServer(function(input, output, session) {

  GetDataSet = reactive({

      validate(
        need(input$data_set, message = "Error: Please load the data set.")
      )

      data_set_object = input$data_set

      if (is.null(data_set_object)) return(NULL)

      tryCatch({
          data_set = read.csv(data_set_object$datapath, header = TRUE, na = ".", stringsAsFactors = FALSE)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )
              
      return(data_set)
  })

  Simulation = reactive({

    #########################################################

    # Design parameters

    parameters = list()

    # Open the trial data set  

    data_set = GetDataSet() 
    parameters$data_set = data_set 

    parameters$time_points = as.numeric(input$time_points) 

    parameters$event_prior_distribution = EventPredPriorDistribution(expected = log(2) / as.numeric(input$event_median_time), uncertainty = as.numeric(input$event_uncertainty))

    parameters$dropout_prior_distribution = EventPredPriorDistribution(expected = log(2) / as.numeric(input$dropout_median_time), uncertainty = as.numeric(input$dropout_uncertainty))

    parameters$enrollment_prior_distribution = EventPredPriorDistribution(expected = as.numeric(input$enrollment_median_time), uncertainty = as.numeric(input$enrollment_uncertainty))

    parameters$nsims = as.numeric(input$nsims)

    #########################################################

    withProgress(message = "Running simulations", value = 1, {

      # Run simulations

      results = EventPred(parameters) 

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

  # Create a matrix for entering time points
  output$FinalAnalysis = renderUI({

      n_time_points = as.numeric(input$n_time_points)  

      data_set = GetDataSet() 

      interim_analysis = max(data_set$enrollment + data_set$time) 

      labels = rep("", n_time_points)
      for (i in 1:n_time_points) labels[i] = paste0("Time point ", i)

      value = matrix(0, n_time_points, 1)
      value[, 1] = interim_analysis + ((1:n_time_points) - 1) * 3

      colnames(value) = "Time"
      rownames(value) = labels

      matrixInput("time_points", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )


  })

output$EventPredictionPlot = renderImage({

    results = Simulation()

    interim_analysis = results$interim_analysis
    prediction = results$prediction
    parameters = results$parameters
    time_points = parameters$time_points

    xlabels = seq(from = 0, to = max(time_points) + 10, by = 10)
    ylabels = seq(from = 0, to = max(prediction) + 10, by = 10)

    width = 600
    height = 500

    filename = tempfile(fileext='.png')
    png(file = filename, width = width, height = height, bg = "transparent")

    plot(x = interim_analysis$time, y = interim_analysis$cumsum, xlab="",ylab="", xlim=c(min(xlabels), max(xlabels)), ylim = c(min(ylabels), max(ylabels)), col="black", lwd = 2, type="l", axes=FALSE)
    polygon(c(rev(time_points), time_points), c(rev(prediction[, 2]), prediction[, 1]), col = "grey80", border = NA)
    lines(x = time_points, y = prediction[, 3], col="red", lwd = 2)
    axis(1, at = xlabels, labels = as.character(xlabels), tck=0.02, mgp=c(0, 0.4, 0))
    axis(2, at = ylabels, labels = as.character(ylabels), tck=0.02, mgp=c(0, 0.4, 0))
    mtext("Number of events", side=2, line=1.5)
    mtext("Time", side=1, line=1.5)
    box() 
    dev.off()

    # Return a list containing the filename
    list(src = filename,
         contentType = 'image/png',
         width = width,
         height = height,
         alt = "Event prediction")

  }, 
  deleteFile = TRUE) 

  output$EventPrediction = renderTable({

    results = Simulation()

    parameters = results$parameters
    original = results$original
    prediction = results$prediction

    column_names = c("Time point", "Mean number of events", "95% predictive interval")

    time_points = parameters$time_points

    n_scenarios = length(time_points)

    col1 = time_points
    col2 = round(prediction[, 3], 1)
    col3 = rep("", 3)

    for (i in 1:n_scenarios) col3[i] = paste0("(", round(prediction[i, 1], 1), ", ", round(prediction[i, 2], 2), ")")

    data_frame = cbind(col1, col2, col3)

    colnames(data_frame) = column_names

    data_frame

  })


  observeEvent(input$jump_to_panel2, {
        updateTabItems(session, "sidebar",
                          selected = "general_parameters")
  })

  observeEvent(input$jump_to_panel3, {
        updateTabItems(session, "sidebar",
                          selected = "simulation")
  })

  observeEvent(input$jump_to_panel4, {
        updateTabItems(session, "sidebar",
                          selected = "report")
  })

})
