# =============================================================================
# ui.R - Survival Analysis Toolkit
# User interface definition
# =============================================================================

fluidPage(
  theme = bslib::bs_theme(version = 5),
  titlePanel(
    div(
      "Survival Analysis Toolkit",
      div(style = "font-size: 14px; font-weight: normal; color: #666; margin-top: 5px;",
          "Yiwen Xu | Biostatistics II Final Project")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Upload Data"),
      fileInput("csv_file", "Select .csv file", accept = ".csv"),
      verbatimTextOutput("csv_preview"),
      
      hr(),
      h4("2. Column Mapping"),
      uiOutput("col_mapping_ui"),
      
      div(style = "margin-top: 5px;",
          actionLink("select_all_numeric", "Select all numeric", 
                     style = "font-size: 0.85em; margin-right: 15px;"),
          actionLink("clear_covars", "Clear all", 
                     style = "font-size: 0.85em; color: #999;")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel("KM Curve", uiOutput("km_tab_ui")),
        tabPanel("Log-Rank Test", uiOutput("logrank_tab_ui")),
        tabPanel("Cox Model", uiOutput("cox_tab_ui")),
        tabPanel(
          "Advanced",
          tabsetPanel(
            tabPanel("Forest Plot", uiOutput("forest_tab_ui")),
            tabPanel("PH Diagnostics", uiOutput("cox_diag_tab_ui")),
            tabPanel("Cutpoint Analysis", uiOutput("cutpoint_tab_ui")),
            tabPanel("Adjusted Curves", uiOutput("adj_tab_ui")),
            tabPanel("Competing Risks", uiOutput("cr_tab_ui")),
            tabPanel("Multi-State", uiOutput("ms_tab_ui"))
          )
        ),
        tabPanel("Help", uiOutput("help_tab_ui"))
      )
    )
  )
)