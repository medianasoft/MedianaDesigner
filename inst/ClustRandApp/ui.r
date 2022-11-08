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
              title = "Cluster-randomized design",
              solidHeader = FALSE,
              collapsible = FALSE,
              width = 12,

              "The web application computes key operating characteristics of a cluster-randomized design for a Phase III clinical trial."
            )
          ),

          fluidRow(

            column(3, class = "col-md-3",

              box(
                title = "Number of trial arms",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                numericInput(inputId = "narms", label = "Number of trial arms", value = 3, min = 2, max = 5),

                tags$p(class = "help-block",
                  "Total number of trial arms, including the control arm.")


              )

            ),

            column(9, class = "col-md-9",
                box(
                  title = "Sample size",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,

                  uiOutput("SampleSize"),

                  tags$p(class = "help-block",
                    "Number of completers in each trial arm. Completers are defined as patients who complete the trial and are included in the final analysis.")

                )

              )
          ),

          fluidRow(

            column(6, class = "col-md-6",
                box(
                  title = "Cluster scheme",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,

                selectInput("cluster_index", label = "Cluster scheme",
                          c("Fixed" = 1, "Random" = 2),

                    tags$p(class = "help-block",
                    "A trial design with fixed (pre-defined) or random cluster sizes.")

                )

              )

            )

          ),

          conditionalPanel(
            condition = "input.cluster_index == 1",

            fluidRow(

                column(3, class = "col-md-3",

                  box(
                    title = "Cluster sizes in the control arm",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = NULL,
                    
                    numericInput(inputId = "n_clusters_fixed1", label = "Number of clusters in the control arm", value = 10, min = 5),

                    uiOutput("ClusterSizeFixed1")


                  )
                ),

                column(9, class = "col-md-9",

                  box(
                    title = "Cluster sizes in the treatment arms",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = NULL,
                    
                    numericInput(inputId = "n_clusters_fixed2", label = "Number of clusters in each treatment arm", value = 10, min = 5),

                    uiOutput("ClusterSizeFixed2")


                  )
                )              

              )

           ),   

            conditionalPanel(
              condition = "input.cluster_index == 2",

            fluidRow(

                column(3, class = "col-md-3",

                  box(
                    title = "Relative cluster sizes in the control arm",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = NULL,
                    
                    numericInput(inputId = "n_clusters_random1", label = "Number of clusters in the control arm", value = 10, min = 5),

                    uiOutput("ClusterSizeRandom1"),

                    tags$p(class = "help-block",
                    "The cluster membership in the control arm is determined using a generalized Bernoulli distribution based on the relative cluster sizes.")

                  )
                ),

                column(9, class = "col-md-9",

                  box(
                    title = "Relative cluster sizes in the treatment arms",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = NULL,
                    
                    numericInput(inputId = "n_clusters_random2", label = "Number of clusters in each treatment arm", value = 10, min = 5),

                    uiOutput("ClusterSizeRandom2"),

                    tags$p(class = "help-block",
                    "The cluster membership in each treatment arm is determined using a generalized Bernoulli distribution based on the relative cluster sizes.")


                  )
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
            column(3, class = "col-md-3",

              box(
                title = "Endpoint type",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                selectInput("endpoint_index", label = "Primary endpoint's type",
                          c("Normal" = 1, "Binary" = 2)),

                selectInput("direction_index", label = "Direction of favorable outcome",
                          c("Higher" = 1, "Lower" = 2))

                )

            ),

            column(9, class = "col-md-9",
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

            column(9, class = "col-md-9",
                box(
                  title = "Variability parameters",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = NULL,

                  uiOutput("VariabilityParameters")

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
                title = "Data analysis method",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                selectInput("method_index", label = "Data analysis method",
                          c("GEE" = 1, "GLMEM" = 2)),

                tags$p(class = "help-block",
                  "GEE: Generalized estimating equations, GLMEM: Generalized linear mixed effects model. The current implementation of the GLMEM method is much less efficient than that for the GEE method and results in longer simulation run times.")

              )

            ),

            conditionalPanel(
              condition = "input.narms >= 3",

                column(6, class = "col-md-6",

                  box(
                    title = "Data analysis method",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = NULL,

                    selectInput("mult_test", label = "Multiple testing procedure",
                              c("Bonferroni" = 1, "Holm" = 2, "Hochberg" = 3))

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

                numericInput(inputId = "nsims", label = "Number of simulation runs", value = 1000)

              )


            )            

          ),

          fluidRow(
            column(6, class = "col-md-6",

              box(
                title = "Descriptive statistics",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,
                
                checkboxInput("descriptive_statistics", "Compute descriptive statistics", FALSE),

                tags$p(class = "help-block",
                  "Key descriptive statistics (arm-specific effects and cluster sizes)  will be computed from each simulation run.")

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
            
            conditionalPanel(
             condition = "input.descriptive_statistics == 1", 

                column(12, class="col-lg-12",
                  box(
                    title = "Descriptive statistics",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = NULL,

                    tableOutput("DescriptiveStatistics"),

                    tags$p(class = "help-block",
                    "Descriptive statistics computed from each simulation run: 2.5%P (2.5th percentile), 25%P (25th percentile), 75%P (75th percentile) and 97.5%P (97.5th percentile).")


                  )
                )

            ),

            column(12, class="col-lg-12",
              box(
                title = "Power calculations",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                width = NULL,

                tableOutput("Power")

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