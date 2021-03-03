library(shiny)
suppressPackageStartupMessages(library(shinydashboard))
library(MedianaDesigner)
library(shinyMatrix)

shinyUI(
  dashboardPage(
    dashboardHeader(title = "MedianaDesigner"),
    dashboardSidebar(
      sidebarMenu(id = "sidebar",
        menuItem("Design parameters", tabName = "design_parameters", icon = icon("sliders-h")),
        menuItem("Endpoint parameters", tabName = "endpoint_parameters", icon = icon("sliders-h")),
        menuItem("Interim parameters", tabName = "interim_parameters", icon = icon("sliders-h")),
        menuItem("General  parameters", tabName = "general_parameters", icon = icon("sliders-h")),
        menuItem("Simulation results", tabName = "simulation", icon = icon("table")), 
        menuItem("Simulation report", tabName = "report", icon = icon("file"))
      )
    ),
    dashboardBody(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "main.css?2")
      ),

      tabItems(

        tabItem(tabName = "design_parameters",

          fluidRow(
            box(
              title = "Adaptive design with sample size or event count re-estimation",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12,

              "The web application computes key operating characteristics of an adaptive design for a two-arm Phase III clinical trial with two interim analyses. The first interim analysis supports early stopping for futility and the second interim analysis supports an option to increase the sample size or number of events."
            )
          ),

          fluidRow(

          column(12, class = "col-md-12",
              box(
                title = "Sample size",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                uiOutput("SampleSize"),

                tags$p(class = "help-block",
                  "Number of enrolled patients in each trial arm.")

              )

            )
          ),

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Endpoint parameters."),
                actionButton("jump_to_panel2", "Next tab")
              )
          )          

        ),

        tabItem(tabName = "endpoint_parameters",

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Endpoint type",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                selectInput("endpoint_index", label = "Primary endpoint's type",
                          c("Normal" = 1, "Binary" = 2, "Time-to-event" = 3))


              )
            ),

            conditionalPanel(
              condition = "input.endpoint_index == 3",

              column(6, class = "col-md-6",

                box(
                  title = "Number of events",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,
                  
                  numericInput(inputId = "event_count", label = "Target number of events", value = 100),

                  tags$p(class = "help-block",
                  "Target number of events at the final analysis before event count re-estimation.")

                )
              )


            )


          ),


          fluidRow(

          column(12, class = "col-md-12",
              box(
                title = "Treatment effect assumptions",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                uiOutput("TreatmentEffectAssumptions"),

                tags$p(class = "help-block",
                  "Treatment effect assumptions for the selected primary endpoint under the alternative hypothesis of beneficial effect, i.e., the experimental treatment is effective.")

              )

            )
          ),


          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Interim parameters."),
                actionButton("jump_to_panel3", "Next tab")
              )
          )          

        ),


        tabItem(tabName = "interim_parameters",

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Interim analysis 1",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "info_frac1", label = "Interim fraction (%)", value = 40, min = 10, max = 90),

                conditionalPanel(
                  condition = "input.endpoint_index == 1 || input.endpoint_index == 2",
                  tags$p(class = "help-block",
                     "Fraction of the total number of patients at Interim analysis 1.")

                ),

                conditionalPanel(
                  condition = "input.endpoint_index == 3",
                  tags$p(class = "help-block",
                     "Fraction of the total number of events at Interim analysis 1.")

                ),

                numericInput(inputId = "futility_threshold", label = "Futility threshold (%)", value = 10),

                tags$p(class = "help-block",
                  "The trial will be stopped for futility at Interim analysis 1 if the predicted probability of success is less than the futility threshold.")



              )

            ),

            column(6, class = "col-md-6",

              box(
                title = "Interim analysis 2",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "info_frac2", label = "Interim fraction (%)", value = 60, min = 10, max = 90),

                conditionalPanel(
                  condition = "input.endpoint_index == 1 || input.endpoint_index == 2",
                  tags$p(class = "help-block",
                     "Fraction of the total number of patients at Interim analysis 2.")

                ),

                conditionalPanel(
                  condition = "input.endpoint_index == 3",
                  tags$p(class = "help-block",
                     "Fraction of the total number of events at Interim analysis 2.")

                ),

                numericInput(inputId = "lower_limit", label = "Promising interval (lower limit) (%)", value = 50),

                numericInput(inputId = "upper_limit", label = "Promising interval (upper limit) (%)", value = 90),

                numericInput(inputId = "target_power", label = "Target probability of success at Final analysis (%)", value = 90),

                conditionalPanel(
                  condition = "input.endpoint_index == 1 || input.endpoint_index == 2",
                tags$p(class = "help-block",
                  "The sample size will be increased after Interim analysis 2 to achieve the target probability of success at Final analysis if the predicted probability of success lies within the promising interval.")
                ),

                conditionalPanel(
                  condition = "input.endpoint_index == 3",
                tags$p(class = "help-block",
                  "The target number of events will be increased after Interim analysis 2 to achieve the target probability of success at Final analysis if the predicted probability of success lies within the promising interval.")
                ),

                numericInput(inputId = "info_frac3", label = "Maximum information fraction (%)", value = 120),

                conditionalPanel(
                  condition = "input.endpoint_index == 1 || input.endpoint_index == 2",
                  tags$p(class = "help-block",
                     "Fraction of the total number of patients at Final analysis after sample size re-estimation.")

                ),

                conditionalPanel(
                  condition = "input.endpoint_index == 3",
                  tags$p(class = "help-block",
                     "Fraction of the total number of events at Final analysis after event count re-estimation.")

                )


              )

            )

          ),

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to General parameters."),
                actionButton("jump_to_panel4", "Next tab")
              )
          )          

        ),


        tabItem(tabName = "general_parameters",

          fluidRow(

            column(6, class = "col-md-6",
              box(
                title = "Patient dropout rate",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "dropout_rate", label = "Patient dropout rate (%)", value = 5, min = 0, max = 50),

                conditionalPanel(
                  condition = "input.endpoint_index == 1 || input.endpoint_index == 2",
                  tags$p(class = "help-block",
                     "Dropout rate at the end of the treatment period.")

                ),

                conditionalPanel(
                  condition = "input.endpoint_index == 3",
                  tags$p(class = "help-block",
                     "Annual dropout rate.")

                )


              )


            )            

          ),

          conditionalPanel(
            condition = "input.endpoint_index == 3",

            fluidRow(
              column(6, class = "col-md-6",

                box(
                  title = "Patient enrollment period",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,
                  
                  numericInput(inputId = "enrollment_period", label = "Length of the patient enrollment period", value = 12)

                )

              ),

              column(6, class = "col-md-6",
                box(
                  title = "Patient enrollment pattern",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,

                  numericInput(inputId = "enrollment_parameter", label = "Median enrollment time", value = 9),

                  tags$p(class = "help-block",
                    "Median enrollment time is defined as the time point by which 50% of the patients will be enrolled.")


                )


              )            

            )

          ),

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Alpha",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "alpha", label = "Type I error rate", value = 0.025),

                tags$p(class = "help-block",
                  "One-sided Type I error rate.")

              )

            ),

            column(6, class = "col-md-6",
              box(
                title = "Number of simulations",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "nsims", label = "Number of simulation runs", value = 10000)

              )


            )            

          ),

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Simulation results."),
                actionButton("jump_to_panel5", "Next tab")
              )
          )          

        ),

        tabItem(tabName = "simulation", 
          fluidRow(
              box(
                title = "Summary of simulation results",
                status = "primary",
                solidHeader = FALSE,
                collapsible = FALSE,
                width = 12
                ),  
            

            column(12, class="col-lg-12",
              box(
                title = "Outcome probabilities",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("OutcomeProbabilities")

              )
            ),

            column(12, class="col-lg-12",
              box(
                title = "Comparison of traditional and adaptive designs",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("Comparison")

              )
            )



          ),

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Simulation report."),
                actionButton("jump_to_panel6", "Next tab")
              )
          )                     

        ),

        tabItem(tabName = "report", 
          fluidRow(
            box(
              title = "Create a simulation report",
              background = "red",
              width = 12,

                tags$p("Click the Download button to create and save a detailed simulation report in a Microsoft Word format."),
                downloadButton("DownloadResults", "Download")
              )
            )
          )
        )

      )  

)
)