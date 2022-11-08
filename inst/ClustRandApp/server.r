shinyServer(function(input, output, session) {

  Simulation = reactive({

    #########################################################

    # Design parameters

    parameters = list()

    narms = as.numeric(input$narms)
    parameters$narms = narms

    endpoint_index = as.numeric(input$endpoint_index)

    direction_index = as.numeric(input$direction_index)
    if (direction_index == 1) parameters$direction = "Higher"
    if (direction_index == 2) parameters$direction = "Lower"

    parameters$sample_size = as.numeric(input$n)   

    cluster_index = as.numeric(input$cluster_index)
    if (cluster_index == 1) parameters$cluster_scheme = "Fixed"
    if (cluster_index == 2) parameters$cluster_scheme = "Random"

    # Fixed cluster sizes
    if (cluster_index == 1) {

      parameters$control_cluster_size = as.numeric(input$cluster_size_fixed1)   
      parameters$treatment_cluster_size = t(input$cluster_size_fixed2)

    }

    # Random cluster sizes
    if (cluster_index == 2) {

      parameters$control_cluster_proportion = as.numeric(input$cluster_size_random1)

      # Matrix of cluster sizes in the treatment arms 
      parameters$treatment_cluster_proportion = t(input$cluster_size_random2)

    }

    assumptions = input$assumptions

    if (endpoint_index == 1) {
      parameters$endpoint_type = "Normal"
      parameters$control_mean = assumptions[1, 1]
      parameters$treatment_mean = as.numeric(assumptions[1, 2:narms])
      parameters$control_between_cluster_sd = as.numeric(input$variability[2, 1]) 
      parameters$treatment_between_cluster_sd = as.numeric(input$variability[2, 2:narms]) 
    }

    if (endpoint_index == 2) {
      parameters$endpoint_type = "Binary"
      parameters$control_rate = assumptions[1, 1] / 100
      parameters$treatment_rate = as.numeric(assumptions[1, 2:narms]  / 100)
    } 

    parameters$control_icc = as.numeric(input$variability[1, 1]) 

    parameters$treatment_icc = as.numeric(input$variability[1, 2:narms]) 

    method_index = as.numeric(input$method_index)
    if (method_index == 1) parameters$method_type = "GEE"
    if (method_index == 2) parameters$method_type = "GLMEM"

    mult_test_list = c("Bonferroni", "Holm", "Hochberg")

    parameters$mult_test = mult_test_list[as.numeric(input$mult_test)]

    parameters$alpha = as.numeric(input$alpha)
    parameters$nsims = as.numeric(input$nsims)

    parameters$descriptive_statistics = as.logical(input$descriptive_statistics)

    #########################################################

    withProgress(message = "Running simulations", value = 1, {

      # Run simulations

      results = ClustRand(parameters) 

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
      value[1, ] = rep(100, narms)

      rownames(value) = "Sample size"
      colnames(value) = trial_arms

      matrixInput("n", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering cluster sizes
  output$ClusterSizeFixed1 = renderUI({

      n_clusters = as.numeric(input$n_clusters_fixed1)

      # Cluster labels  
      labels = rep("", n_clusters)

      for (i in 1:n_clusters) labels[i] = paste0("Cluster ", i)

      value = matrix(0, n_clusters, 1)
      value[, 1] = rep(10, n_clusters)

      rownames(value) = labels
      colnames(value) = "Cluster size"

      matrixInput("cluster_size_fixed1", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering cluster sizes
  output$ClusterSizeFixed2 = renderUI({

      n_clusters = as.numeric(input$n_clusters_fixed2)

      narms = as.numeric(input$narms)

      # Cluster labels  
      labels = rep("", n_clusters)

      for (i in 1:n_clusters) labels[i] = paste0("Cluster ", i)

      value = matrix(10, n_clusters, narms - 1)

      rownames(value) = labels

      if (narms >= 3) {
        col_names = rep("", narms - 1)
        for (i in 1:(narms - 1)) col_names[i] = paste0("Cluster size (Treatment ", i, ")")
      } else {
        col_names = "Cluster size"
      }

      colnames(value) = col_names

      matrixInput("cluster_size_fixed2", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering relative cluster sizes
  output$ClusterSizeRandom1 = renderUI({

      n_clusters = as.numeric(input$n_clusters_random1)

      # Cluster labels  
      labels = rep("", n_clusters)

      for (i in 1:n_clusters) labels[i] = paste0("Cluster ", i)

      value = matrix(0, n_clusters, 1)
      value[, 1] = rep(0.1, n_clusters)

      rownames(value) = labels
      colnames(value) = "Relative cluster size"

      matrixInput("cluster_size_random1", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering relative cluster sizes
  output$ClusterSizeRandom2 = renderUI({

      n_clusters = as.numeric(input$n_clusters_random2)

      narms = as.numeric(input$narms)

      # Cluster labels  
      labels = rep("", n_clusters)

      for (i in 1:n_clusters) labels[i] = paste0("Cluster ", i)

      value = matrix(0.1, n_clusters, narms - 1)

      rownames(value) = labels

      if (narms >= 3) {
        col_names = rep("", narms - 1)
        for (i in 1:(narms - 1)) col_names[i] = paste0("Relative cluster size (Treatment ", i, ")")
      } else {
        col_names = "Relative cluster size"
      }

      colnames(value) = col_names

      matrixInput("cluster_size_random2", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )
  })

  # Create a matrix for entering the variability parameters
  output$VariabilityParameters = renderUI({

      endpoint_index = as.numeric(input$endpoint_index)  

      narms = as.numeric(input$narms)

      # Trial arms  
      trial_arms = "Control"
      if (narms >= 3) {
        for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
      } else {
        trial_arms = c(trial_arms, "Treatment")
      }

      if (endpoint_index == 1) {

        value = matrix(0, 2, narms)
        value[1, ] = rep(0.5, narms)
        value[2, ] = rep(1, narms)

        rownames(value) = c("Intra-cluster correlation coefficient", "Between-cluster standard deviation")

      }

      if (endpoint_index == 2) {

        value = matrix(0, 1, narms)
        value[1, ] = rep(0.5, narms)

        rownames(value) = c("Intra-cluster correlation coefficient")

      }

      colnames(value) = trial_arms

      matrixInput("variability", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  })

  # Create a matrix for entering treatment effect assumptions
  output$TreatmentEffectAssumptions = renderUI({

      endpoint_index = as.numeric(input$endpoint_index)  

      narms = as.numeric(input$narms)

      # Trial arms  
      trial_arms = "Control"
      if (narms >= 3) {
        for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
      } else {
        trial_arms = c(trial_arms, "Treatment")
      }

      value = matrix(0, 1, narms)

      if (endpoint_index == 1) {
          value[1, ] = c(0, rep(0.5, narms - 1))
          rownames(value) = c("Mean")
      }

      if (endpoint_index == 2) {
          value[1, ] = c(10, rep(30, narms - 1))
          rownames(value) = c("Rate (%)")
      }

      colnames(value) = trial_arms

      matrixInput("assumptions", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )


  })

  output$DescriptiveStatistics = renderTable({

    results = Simulation()

    column_names = c("Parameter", "Trial arm", "Statistic", "Value")

    eff_summary = results$sim_summary$eff_summary
    cluster_size_summary = results$sim_summary$cluster_size_summary
    parameters = results$parameters

    narms = parameters$narms

    # Trial arms  
    trial_arms = "Control"
    if (narms >= 3) {
      for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
    } else {
      trial_arms = c(trial_arms, "Treatment")
    }

    if (parameters$endpoint_index == 1) {

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

    if (parameters$endpoint_index == 2) {

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

    colnames(data_frame) = column_names

    data_frame

  })

  output$Power = renderTable({

    results = Simulation()

    parameters = results$parameters
    sim_summary = results$sim_summary

    narms = parameters$narms

    # Trial arms  
    trial_arms = "Control"
    if (narms >= 3) {
      for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
    } else {
      trial_arms = c(trial_arms, "Treatment")
    }

    # GEE
    if (parameters$method_index == 1) {

      if (narms == 2) {

        column_names = c("Analysis method", "Power (%)")

        col1 = c("Standard analysis", 
                 "Bias-corrected analysis (Kauermann-Carroll correction)",
                 "Bias-corrected analysis (Mancl-DeRouen correction)") 
                 
        col2 = c(round(100 * sim_summary$power_sandwich, 1),
                 round(100 * sim_summary$power_kc, 1),
                 round(100 * sim_summary$power_md, 1))
        
        data_frame = data.frame(col1, col2)
      
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

      }

    }  

    # GLMEM
    if (parameters$method_index == 2) {

      if (narms == 2) {

        col1 = "Treatment vs control"
        col2 = round(100 * sim_summary$power, 1)
        
      } else {

        col1 = NULL
        col2 = NULL

        for (i in 1:narms) {

          if (i < narms) comp = paste0(trial_arms[i + 1], " vs control") else comp = "Overall effect"
          col1 = c(col1, comp)
                   
          col2 = c(col2, round(100 * sim_summary$power[i], 1))

        }

      }  

      data_frame = data.frame(col1, col2)
      column_names = c("Comparison", "Power (%)")

    }

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
