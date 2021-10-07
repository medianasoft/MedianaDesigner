
shinyServer(function(input, output, session) {

  Simulation = reactive({

    #########################################################

    # Design parameters

    parameters = list()

    endpoint_index = as.numeric(input$endpoint_index)
    n_comparisons = as.numeric(input$n_comparisons)
    n_endpoints = as.numeric(input$n_endpoints)

    parameters$n_comparisons = n_comparisons
    parameters$n_endpoints = n_endpoints

    parameters$sample_size = as.numeric(input$n)   

    assumptions1 = input$assumptions1

    if (endpoint_index == 1) {
      parameters$endpoint_type = "Normal"
      if (n_endpoints == 1 & n_comparisons >= 2) {
        parameters$control_mean = as.numeric(assumptions1[1, 1])
        parameters$treatment_mean = as.numeric(assumptions1[1, 2:(n_comparisons + 1)])
        parameters$control_sd = as.numeric(assumptions1[2, 1])
        parameters$treatment_sd = as.numeric(assumptions1[2, 2:(n_comparisons + 1)])
      }
      if (n_endpoints >= 2 & n_comparisons == 1) {
        control_mean = rep(0, n_endpoints)
        control_sd = rep(0, n_endpoints)
        treatment_mean = rep(0, n_endpoints)
        treatment_sd = rep(0, n_endpoints)
        for (i in 1:n_endpoints) {
          assumption = input[[paste0("assumptions", i)]]
          control_mean[i] = assumption[1, 1]
          control_sd[i] = assumption[2, 1]
          treatment_mean[i] = assumption[1, 2]
          treatment_sd[i] = assumption[2, 2]
        }
        parameters$control_mean = control_mean
        parameters$control_sd = control_sd
        parameters$treatment_mean = treatment_mean
        parameters$treatment_sd = treatment_sd
      }
      if (n_endpoints >= 2 & n_comparisons >= 2) {
        control_mean = rep(0, n_endpoints)
        control_sd = rep(0, n_endpoints)
        treatment_mean = matrix(0, n_endpoints, n_comparisons)
        treatment_sd = matrix(0, n_endpoints, n_comparisons)
        for (i in 1:n_endpoints) {
          assumption = input[[paste0("assumptions", i)]]
          control_mean[i] = assumption[1, 1]
          control_sd[i] = assumption[2, 1]
          for (j in 2:(n_comparisons + 1)) {
            treatment_mean[i, j - 1] = assumption[1, j]
            treatment_sd[i, j - 1] = assumption[2, j]
          }
        }
        parameters$control_mean = control_mean
        parameters$control_sd = control_sd
        parameters$treatment_mean = treatment_mean
        parameters$treatment_sd = treatment_sd
      }
    }

    if (endpoint_index == 2) {
      parameters$endpoint_type = "Binary"
      if (n_endpoints == 1 & n_comparisons >= 2) {
        parameters$control_rate = as.numeric(assumptions[1, 1] / 100)
        parameters$treatment_rate = as.numeric(assumptions[1, 2:(n_comparisons + 1)] / 100)
      }
      if (n_endpoints >= 2 & n_comparisons == 1) {
        control_rate = rep(0, n_endpoints)
        treatment_rate = rep(0, n_endpoints)
        for (i in 1:n_endpoints) {
          assumption = input[[paste0("assumptions", i)]]
          control_rate[i] = assumption[1, 1]
          treatment_rate[i] = assumption[1, 2]
        }
        parameters$control_rate = control_rate
        parameters$treatment_rate = treatment_rate        
      }
      if (n_endpoints >= 2 & n_comparisons >= 2) {
        control_rate = rep(0, n_endpoints)
        treatment_rate = matrix(0, n_endpoints, n_comparisons)
        for (i in 1:n_endpoints) {
          assumption = get(paste0("input$assumptions", i))
          control_rate[i] = assumption[1, 1]
          for (j in 2:(n_comparisons + 1)) {
            treatment_rate[i, j - 1] = assumption[1, j]
          }
        }
        parameters$control_rate = control_rate
        parameters$treatment_rate = treatment_rate        
      }
    }

    direction_index = as.numeric(input$direction_index)
    if (direction_index == 1) parameters$direction = "Higher"
    if (direction_index == 2) parameters$direction = "Lower"

    if (n_endpoints == 1 & n_comparisons >= 2) {

      mult_test_list = c("Bonferroni", "Holm", "Hochberg", "Hommel", "Fixed-sequence", "Chain")
      comp_mult_test = as.numeric(input$comp_mult_test)
      parameters$mult_test = mult_test_list[comp_mult_test]

      if (comp_mult_test <= 4 | comp_mult_test == 6) {
        parameters$weights = as.numeric(input$comp_weights)        
      }

      if (comp_mult_test == 5) {
        parameters$sequence = as.numeric(input$comp_sequence)        
      }

      if (comp_mult_test == 6) {
        parameters$transition = input$comp_transition
      }

    }

    if (n_endpoints >= 2 & n_comparisons == 1) {

      mult_test_list = c("Bonferroni", "Holm", "Hochberg", "Hommel", "Fixed-sequence", "Chain", "O'Brien")
      end_mult_test = as.numeric(input$end_mult_test)
      parameters$mult_test = mult_test_list[end_mult_test]

      if (end_mult_test <= 4 | end_mult_test == 6) {
        parameters$weights = as.numeric(input$end_weights)        
      }

      if (end_mult_test == 5) {
        parameters$sequence = as.numeric(input$end_sequence)        
      }

      if (end_mult_test == 6) {
        parameters$transition = input$end_transition
      }

    }

    if (n_endpoints >= 2 & n_comparisons >= 2) {

      mult_test_list = c("Holm", "Hochberg", "Hommel")
      gate_mult_test = as.numeric(input$gate_mult_test)
      parameters$mult_test = mult_test_list[gate_mult_test]

      mult_method_list = c("Standard", "Modified", "Enhanced")
      mult_method = as.numeric(input$mult_method)
      parameters$mult_method = mult_method_list[mult_method]

      parameters$mult_test_gamma = as.numeric(input$mult_test_gamma)        

    }

    if (n_endpoints >= 2) parameters$endpoint_correlation = input$correlation

    parameters$dropout_rate = as.numeric(input$dropout_rate / 100)

    parameters$alpha = as.numeric(input$alpha)
    parameters$nsims = as.numeric(input$nsims)

    #########################################################

    withProgress(message = "Running simulations", value = 1, {

      # Run simulations

      results = MultAdj(parameters) 

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

      narms = as.numeric(input$n_comparisons) + 1 

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

  # Create a matrix for entering treatment effect assumptions
  output$TreatmentEffectAssumptions = renderUI({

      narms = as.numeric(input$n_comparisons) + 1 
      n_endpoints = as.numeric(input$n_endpoints) 
      endpoint_index = as.numeric(input$endpoint_index)
      direction_index = as.numeric(input$direction_index)

      # Trial arms  
      trial_arms = "Control"
      if (narms >= 3) {
        for (i in 2:narms) trial_arms = c(trial_arms, paste0("Treatment ", i - 1))
      } else {
        trial_arms = c(trial_arms, "Treatment")
      }

      endpoint_label = rep(", n_endpoints")
      for (i in 1:n_endpoints) endpoint_label[i] = paste0("Endpoint ", i)

      lapply(c(1:n_endpoints), function(i) {

        if (endpoint_index == 1) {
          value = matrix(0, 2, narms)
          for (j in 1:narms) {
            if (direction_index == 1) {
              if (j == 1) value[1, j] = 0 else value[1, j] = 0.3
            }
            if (direction_index == 2) {
              if (j == 1) value[1, j] = 0.3 else value[1, j] = 0
            }
            value[2, j] = 1
          }
          rownames(value) = c("Mean", "SD")
        }

        if (endpoint_index == 2) {
          value = matrix(0, 1, narms)
          for (j in 1:narms) {
            if (direction_index == 1) {
            if (j == 1) value[1, j] = 10 else value[1, j] = 30
            }
            if (direction_index == 2) {
            if (j == 1) value[1, j] = 30 else value[1, j] = 10
            }
          }
          rownames(value) = c("Rate (%)")
        }

        colnames(value) = trial_arms

        list(
          tags$h4(endpoint_label[i]),   
          matrixInput(paste0("assumptions", i), 
            class = "numeric",
            rows = list(names = TRUE),
            cols = list(names = TRUE),
            value = value
          )
        )

      })

  })  

  # Create a matrix for entering pairwise endpoint correlations
  output$EndpointCorrelation = renderUI({

      n_endpoints = as.numeric(input$n_endpoints) 

      value = matrix(1, n_endpoints, n_endpoints)
      for (i in 1:n_endpoints) { 
        for (j in 1:n_endpoints) { 
          if (i != j) value[i, j] = 0.3
        }
      }

      endpoint_label = rep(", n_endpoints")
      for (i in 1:n_endpoints) endpoint_label[i] = paste0("Endpoint ", i)

      rownames(value) = endpoint_label
      colnames(value) = endpoint_label

      matrixInput("correlation", 
        class = "numeric",
        rows = list(names = TRUE),
        cols = list(names = TRUE),
        value = value
      )

  }) 

  # Definitions of all null hypotheses
  output$HypothesisList = renderUI({

      n_comparisons = as.numeric(input$n_comparisons) 
      n_endpoints = as.numeric(input$n_endpoints) 
      n_hypotheses = n_comparisons * n_endpoints

      # Assign hypothesis definitions
      hypothesis_list = rep("", n_hypotheses)

      if (n_comparisons >= 2 & n_endpoints == 1) {

          for (i in 1:n_comparisons) hypothesis_list[i] = paste0("H", i,": Null hypothesis of no difference between Treatment ", i, " and control")
      }

      if (n_comparisons == 1 & n_endpoints >= 2) {

          for (i in 1:n_endpoints) hypothesis_list[i] = paste0("H", i,": Null hypothesis of no difference between the treatment and control with respect to Endpoint ", i)
      }

      if (n_comparisons >= 2 & n_endpoints >= 2) {

        k = 1
        for (i in 1:n_endpoints) {
          for (j in 1:n_comparisons) {
            hypothesis_list[k] = paste0("H", i,": Null hypothesis of no difference between Treatment ", j, " and control with respect to Endpoint ", i)
            k = k + 1
          }
        }

      }

      lapply(1:n_hypotheses, function(i) {

        tags$p(hypothesis_list[i])

      })

  })

  # Create a matrix for entering initial hypothesis weights in trials with a single endpoint and several dose-control comparisons
  output$CompHypothesisWeights = renderUI({

      n_hypotheses = as.numeric(input$n_comparisons) 

      value = matrix(0, n_hypotheses, 1)
      hypothesis_label = rep("", n_hypotheses)
      for (i in 1:n_hypotheses) {
        hypothesis_label[i] = paste0("H", i)
        value[i, 1] = round(1 / n_hypotheses, 3)
      }
      colnames(value) = c("Weight")
      rownames(value) = hypothesis_label

      list(
        matrixInput("comp_weights", 
          class = "numeric",
          rows = list(names = TRUE),
          cols = list(names = TRUE),
          value = value
        )
      )

  })  

  # Create a matrix for entering the hypothesis testing sequence in trials with a single endpoint and several dose-control comparisons
  output$CompSequence = renderUI({

      n_hypotheses = as.numeric(input$n_comparisons) 

      value = matrix(0, n_hypotheses, 1)
      hypothesis_label = rep("", n_hypotheses)
      for (i in 1:n_hypotheses) {
        hypothesis_label[i] = paste0("H", i)
        value[i, 1] = i
      }
      colnames(value) = c("Sequence")
      rownames(value) = hypothesis_label

      list(
        matrixInput("comp_sequence", 
          class = "numeric",
          rows = list(names = TRUE),
          cols = list(names = TRUE),
          value = value
        )
      )

  }) 

  # Create a matrix for entering initial hypothesis weights in trials with several endpoints and a single dose-control comparison
  output$EndHypothesisWeights = renderUI({

      n_hypotheses = as.numeric(input$n_endpoints) 

      value = matrix(0, n_hypotheses, 1)
      hypothesis_label = rep("", n_hypotheses)
      for (i in 1:n_hypotheses) {
        hypothesis_label[i] = paste0("H", i)
        value[i, 1] = round(1 / n_hypotheses, 3)
      }
      colnames(value) = c("Weight")
      rownames(value) = hypothesis_label

      list(
        matrixInput("end_weights", 
          class = "numeric",
          rows = list(names = TRUE),
          cols = list(names = TRUE),
          value = value
        )
      )

  }) 

  # Create a matrix for entering the hypothesis testing sequence in trials with several endpoints and a single dose-control comparison
  output$EndSequence = renderUI({

      n_hypotheses = as.numeric(input$n_endpoints) 

      value = matrix(0, n_hypotheses, 1)
      hypothesis_label = rep("", n_hypotheses)
      for (i in 1:n_hypotheses) {
        hypothesis_label[i] = paste0("H", i)
        value[i, 1] = i
      }
      colnames(value) = c("Sequence")
      rownames(value) = hypothesis_label

      list(
        matrixInput("end_sequence", 
          class = "numeric",
          rows = list(names = TRUE),
          cols = list(names = TRUE),
          value = value
        )
      )

  }) 

  # Create a matrix for entering transition parameters in trials with a single endpoint and several dose-control comparisons
  output$CompTransition = renderUI({

      n_hypotheses = as.numeric(input$n_comparisons) 

      value = matrix(0, n_hypotheses, n_hypotheses, byrow = TRUE)

      hypothesis_label = rep("", n_hypotheses)
      for (i in 1:n_hypotheses) hypothesis_label[i] = paste0("H", i)

      for (i in 1:(n_hypotheses - 1)) value[i, i + 1] = 1

      colnames(value) = hypothesis_label
      rownames(value) = hypothesis_label

      list(
        matrixInput("comp_transition", 
          class = "numeric",
          rows = list(names = TRUE),
          cols = list(names = TRUE),
          value = value
        )
      )

  })  

  # Create a matrix for entering transition parameters in trials with several endpoints and a single dose-control comparison
  output$EndTransition = renderUI({

      n_hypotheses = as.numeric(input$n_endpoints) 

      value = matrix(0, n_hypotheses, n_hypotheses, byrow = TRUE)

      hypothesis_label = rep("", n_hypotheses)
      for (i in 1:n_hypotheses) hypothesis_label[i] = paste0("H", i)

      for (i in 1:(n_hypotheses - 1)) value[i, i + 1] = 1

      colnames(value) = hypothesis_label
      rownames(value) = hypothesis_label

      list(
        matrixInput("end_transition", 
          class = "numeric",
          rows = list(names = TRUE),
          cols = list(names = TRUE),
          value = value
        )
      )

  })  

  # Create a matrix for entering the family-specific truncation parameters in trials with several endpoints and several dose-control comparisons
  output$TruncationParameters = renderUI({

      n_endpoints = as.numeric(input$n_endpoints) 

      value = matrix(0, n_endpoints, 1)
      endpoint_label = rep("", n_endpoints)
      for (i in 1:n_endpoints) {
        endpoint_label[i] = paste0("Family ", i)
        value[i, 1] = 0.8
      }
      value[n_endpoints, 1] = 1
      colnames(value) = c("Truncation parameter")
      rownames(value) = endpoint_label

      list(
        matrixInput("mult_test_gamma", 
          class = "numeric",
          rows = list(names = TRUE),
          cols = list(names = TRUE),
          value = value
        )
      )

  }) 

  # Hypothesis-specific power for traditional multiplicity adjustments or gatekeeping procedures
  output$HypothesisPower1 = renderTable({

    results = Simulation()

    n_comparisons = as.numeric(input$n_comparisons) 
    n_endpoints = as.numeric(input$n_endpoints) 
    n_hypotheses = n_comparisons * n_endpoints

    hypothesis_list = rep("", n_hypotheses)
    for (i in 1:n_hypotheses) hypothesis_list[i] = paste0("H", i)

    sim_summary = results$sim_summary

print(sim_summary)

    column_names = c("Hypothesis", "Power (%)", "Adjusted power (%)")

    col1 = hypothesis_list
    col2 = round(100 * sim_summary$power, 1)
    col3 = round(100 * sim_summary$adj_power, 1)
    
    data_frame = data.frame(col1, col2, col3)

    colnames(data_frame) = column_names

    data_frame

  })

  # Hypothesis-specific power for global testing procedures
  output$HypothesisPower2 = renderTable({

    results = Simulation()

    n_hypotheses = as.numeric(input$n_endpoints) 

    hypothesis_list = rep("", n_hypotheses)
    for (i in 1:n_hypotheses) hypothesis_list[i] = paste0("H", i)

    sim_summary = results$sim_summary

    column_names = c("Hypothesis", "Power (%)")

    col1 = hypothesis_list
    col2 = round(100 * sim_summary$power, 1)
    
    data_frame = data.frame(col1, col2)

    colnames(data_frame) = column_names

    data_frame

  })

  # Overall power for global testing procedures
  output$OverallPower1 = renderTable({

    results = Simulation()

    sim_summary = results$sim_summary

    column_names = c("Overall power", "Power (%)")

    col1 = c("Disjunctive power", "Conjunctive power")
    col2 = c(round(100 * sim_summary$disj_power, 1), round(100 * sim_summary$conj_power, 1))
    
    data_frame = data.frame(col1, col2)

    colnames(data_frame) = column_names

    data_frame

  })

  # Global power for traditional multiplicity adjustments
  output$OverallPower2 = renderTable({

    results = Simulation()

    sim_summary = results$sim_summary

    column_names = c("Overall power", "Power (%)")

    col1 = c("Global power")
    col2 = round(100 * sim_summary$adj_power, 1)
    
    data_frame = data.frame(col1, col2)

    colnames(data_frame) = column_names

    data_frame

  })

  # Global power for gatekeeping procedures
  output$OverallPower3 = renderTable({

    results = Simulation()

    sim_summary = results$sim_summary
    n_endpoints = as.numeric(input$n_endpoints) 

    column_names = c("Endpoint family", "Overall power", "Power (%)")

    col1 = NULL
    col2 = NULL
    col3 = NULL

    for (i in 1:n_endpoints) {

      col1 = c(col1, paste0("Endpoint ", i), "")
      col2 = c(col2, "Disjunctive power", "Conjunctive power")
      col3 = c(col3, round(100 * sim_summary$disj_power[i], 1), round(100 * sim_summary$conj_power[i], 1))
    
    }

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
                          selected = "adjustment_parameters")
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
