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
              title = "Blinded event prediction in event-driven trials",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12,

              "The web application computes event predictions for Phase II or Phase III trials with an event-driven design. Blinded event data at an interim analysis are used to forecast the number of events at pre-defined time points in the future."
            )
          ),

          fluidRow(
            column(12, class = "col-md-12",

              box(
                title = "Event data",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                tags$p("Select a data set with blinded patient enrollment, event and dropout data at the interim analysis (comma-separated value file):"),
                fileInput('data_set', '', accept = c('.csv') ),

                tags$p(class = "help-block",
                  "The data set is required to include four variables (enrollment, time of patient enrollment; time, time to event or last contact; event, event indicator; dropout, patient dropout indicator) with a single record per patient.")

              )

            )

          ),

          fluidRow(

            column(6, class = "col-md-6",

              box(
                title = "Future time points for computing event predictions",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "n_time_points", label = "Number of time points", 5, min = 2, max = 10),

                uiOutput("FinalAnalysis"),

                tags$p(class = "help-block",
                  "Future time points for computing event predictions in appropriate time units, e.g., weeks or months.")


              )

            )

          ),

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to General parameters."),
                actionButton("jump_to_panel2", "Next tab")
              )
          )          

        ),

        tabItem(tabName = "general_parameters",

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Event rate: Prior distribution",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "event_median_time", label = "Expected median time", value = 15),

                  tags$p(class = "help-block",
                     "Expected median time for the event of interest."),

                numericInput(inputId = "event_uncertainty", label = "Uncertainty parameter", value = 0.3),

                  tags$p(class = "help-block",
                     "Uncertainty parameter for the event hazard rate.")

              )

            ),

            column(6, class = "col-md-6",

              box(
                title = "Patient dropout rate: Prior distribution",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "dropout_median_time", label = "Expected median time", value = 80),

                  tags$p(class = "help-block",
                     "Expected median time for the patient dropout process."),

                numericInput(inputId = "dropout_uncertainty", label = "Uncertainty parameter", value = 0.3),

                  tags$p(class = "help-block",
                     "Coefficient of variation for the patient dropout hazard rate.")

              )

            )     

          ),

          fluidRow(

            column(6, class = "col-md-6",

              box(
                title = "Patient enrollment rate: Prior distribution",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "enrollment_median_time", label = "Enrollment rate", value = 35),

                  tags$p(class = "help-block",
                     "Expected patient enrollment rate."),

                numericInput(inputId = "enrollment_uncertainty", label = "Uncertainty parameter", value = 0.3),

                  tags$p(class = "help-block",
                     "Coefficient of variation for the patient enrollment rate.")

              )

            )            

          ),

          fluidRow(

            column(6, class = "col-md-6",
              box(
                title = "Number of simulations",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "nsims", label = "Number of simulation runs", value = 1000)

              )


            )            

          ),

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Simulation results."),
                actionButton("jump_to_panel3", "Next tab")
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
                title = "Event prediction at pre-defined time points",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                tableOutput("EventPrediction")

              )
            ),            

            column(12, class="col-lg-12",

              box(
                title = "Event prediction at pre-defined time points",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                imageOutput("EventPredictionPlot", height = "auto"),

                tags$p(class = "help-block",
                  "Black curve: Observed events. Red curve: Predicted mean number of events. Gray band: 95% predictive interval.")

              )
            )
          ),            

          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Simulation report."),
                actionButton("jump_to_panel4", "Next tab")
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