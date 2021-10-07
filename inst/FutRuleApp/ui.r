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
              title = "Optimal selection of a futility stopping rule",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12,

              "The web application computes operating characteristics of a multi-arm trial design with a single interim analysis. A futility stopping rule will be applied at this interim look and the trial will be stopped early for futility if the predicted probability of success (conditional power) is less than a pre-defined futility threshold in all treatment arms. An optimal value of the futility threshold is derived by maximizing the sensitivity and specificity rates."
            )
          ),

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Number of trial arms",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "narms", label = "Number of trial arms", value = 3, min = 2, max = 4),

                tags$p(class = "help-block",
                  "Total number of trial arms, including the control arm.")


              )

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
                          c("Normal" = 1, "Binary" = 2, "Time-to-event" = 3)),

                selectInput("direction_index", label = "Direction of favorable outcome",
                          c("Higher" = 1, "Lower" = 2))

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
                  "Target number of events at the final analysis.")

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
                  "Treatment effect assumptions for the selected primary endpoint under the alternative hypothesis of beneficial effect, i.e., all treatments are effective.")

              )

            )
          ),


          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to General parameters."),
                actionButton("jump_to_panel3", "Next tab")
              )
          )          

        ),


        tabItem(tabName = "general_parameters",

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Interim fraction",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "info_frac", label = "Interim fraction (%)", value = 50, min = 10, max = 90),

                conditionalPanel(
                  condition = "input.endpoint_index == 1 || input.endpoint_index == 2",
                  tags$p(class = "help-block",
                     "Fraction of the total number of patients at the interim analysis.")

                ),

                conditionalPanel(
                  condition = "input.endpoint_index == 3",
                  tags$p(class = "help-block",
                     "Fraction of the total number of events at the interim analysis.")

                )


              )

            ),

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
                actionButton("jump_to_panel4", "Next tab")
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
                title = "Sensitivity and specificity rates",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                imageOutput("SensSpec", height = "auto"),

                tags$p(class = "help-block",
                  "Red curve: Sensitivity rate (probability of correctly retaining at least one treatment arm at the interim analysis, evaluated under the alternative hypothesis of beneficial effect, i.e., all treatments are effective). Blue curve: Specificity rate (probability of correctly stopping all treatment arms at the interim analysis due to futility, evaluated under the null hypothesis of no effect, i.e., all treatments are ineffective).")

              )
            ),            
 
            column(12, class="col-lg-12",

              box(
                title = "Accuracy rate",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                imageOutput("Accuracy", height = "auto"),

                tags$p(class = "help-block",
                  "The accuracy rate is defined as the average of the sensitivity and specificity rates. An optimal futility threshold (solid line) is defined as the threshold that maximizes the accuracy rate. The dashed lines are drawn at the lower and upper limits of the 95% optimal interval.")

              )
            ),

            column(12, class="col-lg-12",
              box(
                title = "Optimal futility threshold",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("OptimalThreshold"),

                tags$p(class = "help-block",
                  "Optimal futility threshold is defined as the threshold that maximizes the accuracy rate.")

              )
            )


          ),

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Simulation report."),
                actionButton("jump_to_panel5", "Next tab")
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