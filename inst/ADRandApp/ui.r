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
              title = "Adaptive trials with response-adaptive randomization",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12,

              "The web application computes key operating characteristics 
              of an adaptive design for a dose-finding Phase II clinical trial 
              with multiple interim analyses aimed at updating the randomization 
              scheme based on the accumulating efficacy data."
            )
          ),

          fluidRow(
            column(12, class = "col-md-12",

              box(
                title = "Number of trial arms",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "n_arms", label = "Number of trial arms", value = 3, min = 2, max = 6),

                tags$p(class = "help-block",
                  "Total number of trial arms, including the placebo arm.")
              )

            )
          ),

          fluidRow(
            column(12, class = "col-md-12",

              box(
                title = "Dose levels",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                uiOutput("DoseLevels")
              )

            )
          ),

          fluidRow(
            column(12, class = "col-md-12",

              box(
                title = "Number of trial stages",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "n_stages", label = "Number of trial stages", value = 3, min = 2, max = 5)
              )

            )
          ),

          fluidRow(

            column(12, class = "col-md-12",
              box(
                title = "Sample size in each stage",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                uiOutput("StageSampleSize"),

                tags$p(class = "help-block",
                  "Number of enrolled patients in each trial stage.")

              )
            )

          ),

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Patient enrollment period",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "enrollment_period", label = "Length of the patient enrollment period", value = 36)

              )

            ),

            column(6, class = "col-md-6",
              box(
                title = "Patient enrollment pattern",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "enrollment_parameter", label = "Median enrollment time", value = 24),

                tags$p(class = "help-block",
                  "Median enrollment time is defined as the time point by which 50% of the patients will be enrolled.")

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

            column(width = 6,
              box(
                title = "Endpoint type",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                selectInput("endpoint_index", label = "Primary endpoint's type", c("Normal" = 1)),

                selectInput("direction_index", label = "Direction of favorable outcome", c("Higher" = 1, "Lower" = 2))

              )

            ),

            column(width = 6,
              box(
                title = "Treatment period",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                  numericInput(inputId = "treatment_period", label = "Length of the treatment period", value = 6)

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

                uiOutput("TreatmentEffectAssumptions")

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

            column(width = 6,
              box(
                title = "Randomization ratio",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "ratio_placebo_percent", label = "Fixed randomization ratio in the placebo arm (%)", value = 25, min = 1, max = 99)

              )

            ),

            column(width = 6,
              box(
                title = "Clinically meaningful improvement",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "delta", label = "Threshold for clinically meaningful improvement over placebo", value = 10)

              )

            )


          ),

          fluidRow(

            column(width = 6,
              box(
                title = "Balance parameter",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "balance", label = "Balance parameter for adaptive randomization", value = 1, min = 0, max = 3)

              )

            )
          ),

          fluidRow(

            column(width = 12,
              box(
                title = "Parameters of candidate dose-response models",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(
                  inputId = "exponential_model_parameter", 
                  label = "Exponential model: Delta", 
                  value = 200, min = 0
                ),

                numericInput(
                  inputId = "emax_model_parameter", 
                  label = "Emax model: ED50", 
                  value = 200, min = 0, max = 0
                ),

                numericInput(
                  inputId = "logistic_model_ed50", 
                  label = "Logistic model: ED50", 
                  value = 200, min = 0, max = 0
                ),
                numericInput(
                  inputId = "logistic_model_delta", 
                  label = "Logistic model: Delta", 
                  value = 50, min = 0, max = 0
                ),

                tags$p(class = "help-block",
                  "Non-linear parameters of the candidate dose-response models used in the MCPMod method. No parameters need to be specified for the linear model."
                )

              )

            )

          ),

          fluidRow(
            column(width = 12,

              box(
                title = "Candidate dose-response models used in the MCPMod method",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                imageOutput("DoseResponseModels", height = "auto"),

                tags$p(class = "help-block",
                  "Black curve: Linear model, Blue curve: Exponential model, Red curve: Emax model, Green curve: Logistic model.")

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

                  tags$p(class = "help-block",
                     "Dropout rate at the end of the treatment period.")

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

                numericInput(inputId = "nsims", label = "Number of simulation runs", value = 100)

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
                title = "Sample size by stage",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("SampleSizeByStage")

              )
            ),

            column(12, class="col-lg-12",
              box(
                title = "Sample size by trial arm",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("SampleSizeByTrialArm")

              )
            ),

            column(12, class="col-lg-12",
              box(
                title = "Comparison of traditional and adaptive designs",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("ComparisonOfTraditionalAndAdaptiveDesigns"),

                tags$p(class = "help-block",
                  "Probability of a statistically significant dose-response relationship based on the MCPMod method.")

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