library(shiny)
suppressPackageStartupMessages(library(shinydashboard))
library(MultAdj)
library(shinyMatrix)

shinyUI(
  dashboardPage(
    dashboardHeader(title = "MedianaDesigner"),
    dashboardSidebar(
      sidebarMenu(id = "sidebar",
        menuItem("Design parameters", tabName = "design_parameters", icon = icon("sliders-h")),
        menuItem("Endpoint parameters", tabName = "endpoint_parameters", icon = icon("sliders-h")),
        menuItem("Multiplicity adjustment", tabName = "adjustment_parameters", icon = icon("sliders-h")),
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
              title = "Traditional design with multiple outcomes",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12,

              tags$h5("The web application computes key operating characteristics of a traditional design with multiple objectives, including"),

              tags$h5("Trials with a single endpoint and several dose-control comparisons."),
              tags$h5("Trials with a single dose-control comparison and several endpoints"),
              tags$h5("Trials with several dose-control comparisons and several endpoints")

            )
          ),


          fluidRow(

          column(6, class = "col-md-6",
              box(
                title = "Dose-placebo comparisons",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "n_comparisons", label = "Number of dose-placebo comparisons", value = 3, min = 1, max = 5)

              )

            ),

          column(6, class = "col-md-6",
              box(
                title = "Clinical endpoints",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                numericInput(inputId = "n_endpoints", label = "Number of endpoints", value = 1, min = 1, max = 5)

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
                
                selectInput("endpoint_index", label = "Primary endpoint's type", c("Normal" = 1, "Binary" = 2)),

                selectInput("direction_index", label = "Direction of favorable outcome", c("Higher" = 1, "Lower" = 2))

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

          conditionalPanel(
            condition = "input.n_endpoints >= 2",

            fluidRow(
              column(12, class="col-lg-12",

                box(
                  title = "Endpoint correlation matrix",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,
                  
                  uiOutput("EndpointCorrelation"),

                  tags$p(class = "help-block",
                    "Pairwise correlation coefficients for the endpoint-specified test statistics.")

                )
              )
            
            )          

          ),


          fluidRow(
              box(
                title = "Next step",
                background = "red",
                width = 12,
                tags$p("Proceed to Multiplicity adjustment."),
                actionButton("jump_to_panel3", "Next tab")
              )
          )          

        ),

        tabItem(tabName = "adjustment_parameters",

            fluidRow(

              column(12, class="col-lg-12",

                box(
                  title = "Hypothesis definitions",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,

                  uiOutput("HypothesisList")

                )
              )

          ),    

          conditionalPanel(
            condition = "input.n_comparisons >= 2 && input.n_endpoints == 1",

            fluidRow(

                column(6, class = "col-md-6",
                    box(
                      title = "Multiple testing procedure",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      selectInput("comp_mult_test", label = "Multiple testing procedure", c("Bonferroni" = 1, "Holm" = 2, "Hochberg" = 3, "Hommel" = 4, "Fixed-sequence" = 5, "Chain" = 6))

                    )

                  ),

                conditionalPanel(
                  condition = "input.comp_mult_test <= 4 || input.comp_mult_test == 6",

                  column(6, class = "col-md-6",

                      box(
                        title = "Initial hypothesis weights",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width = NULL,
                        
                        uiOutput("CompHypothesisWeights"),

                        tags$p(class = "help-block",
                          "The hypothesis weights range from 0 to 1 and must add up to 1.")

                      )
                    )

                ),

                conditionalPanel(
                  condition = "input.comp_mult_test == 5",

                  column(6, class = "col-md-6",

                      box(
                        title = "Hypothesis testing sequence",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width = NULL,
                        
                        uiOutput("CompSequence"),

                        tags$p(class = "help-block",
                          "The hypothesis testing sequence for the fixed-sequence procedure.")

                      )
                    )
                )
            ),

            conditionalPanel(
              condition = "input.comp_mult_test == 6",

              fluidRow(
                  column(12, class = "col-lg-12",

                      box(
                        title = "Transition parameters",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width = NULL,
                        
                        uiOutput("CompTransition"),

                        tags$p(class = "help-block",
                          "The transition parameters range from 0 to 1, the diagonal entries must be 0 and the transition parameters in each row must add up to 1.")

                      )
                    )
              )
            )
          
          ),

          conditionalPanel(
            condition = "input.n_comparisons == 1 && input.n_endpoints >= 2",

            fluidRow(

                column(6, class = "col-md-6",
                    box(
                      title = "Multiple or global testing procedure",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      selectInput("end_mult_test", label = "Multiple or global testing procedure", c("Bonferroni" = 1, "Holm" = 2, "Hochberg" = 3, "Hommel" = 4, "Fixed-sequence" = 5, "Chain" = 6, "O'Brien" = 7))

                    )

                  ),

                conditionalPanel(
                  condition = "input.end_mult_test <= 4 || input.end_mult_test == 6",

                  column(6, class = "col-md-6",

                      box(
                        title = "Initial hypothesis weights",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width = NULL,
                        
                        uiOutput("EndHypothesisWeights"),

                        tags$p(class = "help-block",
                          "The hypothesis weights range from 0 to 1 and must add up to 1.")

                      )
                    )

                ),            

                conditionalPanel(
                  condition = "input.end_mult_test == 5",

                  column(6, class = "col-md-6",

                      box(
                        title = "Hypothesis testing sequence",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width = NULL,
                        
                        uiOutput("EndSequence"),

                        tags$p(class = "help-block",
                          "The hypothesis testing sequence for the fixed-sequence procedure.")

                      )
                    )
                )

            ),

            conditionalPanel(
              condition = "input.end_mult_test == 6",

              fluidRow(
                  column(12, class = "col-lg-12",

                      box(
                        title = "Transition parameters",
                        status = "primary",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        width = NULL,
                        
                        uiOutput("EndTransition"),

                        tags$p(class = "help-block",
                          "The transition parameters range from 0 to 1, the diagonal entries must be 0 and the transition parameters in each row must add up to 1.")

                      )
                    )
              )
            )            
          
          ),

          conditionalPanel(
            condition = "input.n_comparisons >= 2 && input.n_endpoints >= 2",

            fluidRow(

                column(6, class = "col-md-6",
                    box(
                      title = "Gatekeeping procedure",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      selectInput("gate_mult_test", label = "Component of gatekeeping procedure", c("Hochberg" = 1, "Hommel" = 2)),

                      selectInput("mult_method", label = "Mixture method used in the gatekeeping procedure", c("Standard" = 1, "Modified" = 2, "Enhanced" = 3))

                    )

                  ),

                column(6, class = "col-md-6",
                    box(
                      title = "Truncation parameters",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      uiOutput("TruncationParameters"),

                      tags$p(class = "help-block",
                        "The family-specific truncation parameters element  range from 0 to 1. The last parameter may be equal to 1 whereas the other parameters must be strictly less than 1.")

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
                
                numericInput(inputId = "alpha", label = "Overall Type I error rate", value = 0.025),

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

              # Traditional multiplicity adjustments or gatekeeping procedures
              conditionalPanel(
                condition = "(input.n_comparisons >= 2 && input.n_endpoints == 1) | (input.n_comparisons == 1 && input.n_endpoints >= 2 && input.end_mult_test <= 6) | (input.n_comparisons >= 2 && input.n_endpoints >= 2)",

                  column(12, class="col-lg-12",
                    box(
                      title = "Hypothesis-specific power",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      tableOutput("HypothesisPower1"),

                      tags$p(class = "help-block",
                        "Power: Probability of rejecting each hypothesis of no effect without a multiplicity adjustment. Adjusted power: Probability of rejecting each hypothesis of no effect using a multiplicity adjustment based on the specified multiple testing procedure.")

                    )
                  )

              ),

              # Global testing procedure
              conditionalPanel(
                condition = "input.n_comparisons == 1 && input.n_endpoints >= 2 && input.end_mult_test == 7",

                  column(12, class="col-lg-12",
                    box(
                      title = "Hypothesis-specific power",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      tableOutput("HypothesisPower2"),

                      tags$p(class = "help-block",
                        "Power: Probability of rejecting each hypothesis of no effect without a multiplicity adjustment.")

                    )
                  ),

                  column(12, class="col-lg-12",
                    box(
                      title = "Overall power",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      tableOutput("OverallPower2"),

                      tags$p(class = "help-block",
                        "Global power: Probability of rejecting at least one hypothesis of no effect using the specified global testing approach.")

                    )
                  )

              ),

              # Traditional multiplicity adjustments 
              conditionalPanel(
                condition = "(input.n_comparisons >= 2 && input.n_endpoints == 1) | (input.n_comparisons == 1 && input.n_endpoints >= 2 && input.end_mult_test <= 6)",

                  column(12, class="col-lg-12",
                    box(
                      title = "Overall power",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      tableOutput("OverallPower1"),

                      tags$p(class = "help-block",
                        "Disjunctive power: Probability of rejecting at least one hypothesis of no effect using a multiplicity adjustment based on the specified multiple testing procedure. Conjunctive power: Probability of rejecting all hypotheses of no effect using a multiplicity adjustment based on the specified multiple testing procedure.")

                    )
                  )

              ),

              # Gatekeeping procedures
              conditionalPanel(
                condition = "input.n_comparisons >= 2 && input.n_endpoints >= 2",

                  column(12, class="col-lg-12",
                    box(
                      title = "Overall power",
                      status = "primary",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      width = NULL,

                      tableOutput("OverallPower3"),

                      tags$p(class = "help-block",
                        "Disjunctive power: Probability of rejecting at least one hypothesis of no effect within each endpoint family using a multiplicity adjustment based on the specified multiple testing procedure. Conjunctive power: Probability of rejecting all hypotheses of no effect within each endpoint family using a multiplicity adjustment based on the specified multiple testing procedure.")

                    )
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