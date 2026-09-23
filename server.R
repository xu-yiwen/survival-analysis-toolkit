# =============================================================================
# server.R - Survival Analysis Toolkit
# Server logic and analysis functions
# =============================================================================

function(input, output, session) {
  
  library(cmprsk)
  
  # ---------------------------------------------------------------------------
  # Utility Functions
  # ---------------------------------------------------------------------------
  
  # Safely quote column names for formula construction
  safe_colname <- function(colname) {
    if (grepl("[^a-zA-Z0-9_.]", colname) || grepl("^[0-9]", colname)) {
      return(paste0("`", gsub("`", "\\\\`", colname), "`"))
    }
    return(colname)
  }
  
  # Convert character columns to factors for modeling
  prepare_model_data <- function(d, vars) {
    if (length(vars) > 0) {
      for (v in vars) {
        if (v %in% names(d) && is.character(d[[v]])) {
          d[[v]] <- as.factor(d[[v]])
        }
      }
    }
    d
  }
  
  # Build survival formula string safely
  build_surv_formula_str <- function(time_var, event_var, rhs = "1") {
    time_safe <- safe_colname(time_var)
    event_safe <- safe_colname(event_var)
    
    lhs <- paste0("Surv(", time_safe, ", ", event_safe, ")")
    
    if (length(rhs) == 1 && rhs == "1") {
      return(paste(lhs, "~ 1"))
    }
    
    rhs_safe <- paste(sapply(rhs, safe_colname), collapse = " + ")
    paste(lhs, "~", rhs_safe)
  }
  
  # Build survival formula object safely
  build_surv_formula <- function(time_var, event_var, rhs = "1") {
    as.formula(build_surv_formula_str(time_var, event_var, rhs))
  }
  
  # Format numeric values for display
  fmt_num <- function(x, digits = 3) {
    if (is.null(x) || length(x) == 0 || is.na(x)) return("NA")
    formatC(x, format = "f", digits = digits)
  }
  
  # Format p-values for display
  fmt_p <- function(p, digits = 3) {
    if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
    if (p < 0.001) return("<0.001")
    formatC(p, format = "f", digits = digits)
  }
  
  # Guess if a variable is continuous or categorical
  guess_var_type <- function(x, max_unique_categorical = 5, max_unique_integer = 10) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return("categorical")
    if (is.factor(x) || is.character(x) || is.logical(x)) return("categorical")
    if (!is.numeric(x)) return("categorical")
    
    ux <- unique(x)
    n_ux <- length(ux)
    
    if (n_ux <= max_unique_categorical) return("categorical")
    if (all(abs(ux - round(ux)) < 1e-8) && n_ux <= max_unique_integer) return("categorical")
    
    "continuous"
  }
  
  # Create group labels for cutpoint analysis
  make_group_labels <- function(x, cut_value = NULL) {
    var_type <- guess_var_type(x)
    
    if (var_type == "categorical") {
      levs <- sort(unique(as.character(x[!is.na(x)])))
      return(list(type = "categorical", labels = levs))
    }
    
    if (is.null(cut_value)) {
      stop("cut_value is required for continuous variables.")
    }
    
    low_lab  <- paste0("Not greater than ", signif(cut_value, 4))
    high_lab <- paste0("Greater than ", signif(cut_value, 4))
    
    list(type = "continuous", labels = c(low_lab, high_lab))
  }
  
  # ---------------------------------------------------------------------------
  # Data Input and Preview
  # ---------------------------------------------------------------------------
  
  df_raw <- reactive({
    req(input$csv_file)
    read.csv(input$csv_file$datapath, stringsAsFactors = FALSE)
  })
  
  output$csv_preview <- renderPrint({
    req(df_raw())
    str(df_raw())
  })
  
  # ---------------------------------------------------------------------------
  # Column Mapping UI
  # ---------------------------------------------------------------------------
  
  output$col_mapping_ui <- renderUI({
    req(df_raw())
    cols <- as.character(unlist(colnames(df_raw())))
    
    t_guess <- cols[grepl("time|surv|os|pfs|day", cols, ignore.case = TRUE)][1]
    e_guess <- cols[grepl("event|status|dead", cols, ignore.case = TRUE)][1]
    
    if (is.na(t_guess) && length(cols) >= 1) t_guess <- cols[1]
    if (is.na(e_guess) && length(cols) >= 2) e_guess <- cols[2]
    
    tagList(
      selectInput("time_var", "Time variable", choices = cols, selected = t_guess),
      selectInput("event_var", "Event indicator", choices = cols, selected = e_guess),
      selectizeInput("covars", "Covariates", choices = cols, multiple = TRUE),
      selectInput("strata", "Stratify by", choices = c("None" = " ", cols))
    )
  })
  
  # Covariate selection shortcuts
  observeEvent(input$select_all_numeric, {
    req(df_raw())
    num_cols <- names(df_raw())[sapply(df_raw(), is.numeric)]
    exclude <- c(input$time_var, input$event_var, input$strata)
    num_cols <- setdiff(num_cols, exclude)
    updateSelectizeInput(session, "covars", selected = num_cols)
  })
  
  observeEvent(input$clear_covars, {
    updateSelectizeInput(session, "covars", selected = character(0))
  })
  
  observeEvent(input$strata, {
    req(input$time_var, input$event_var)
    exclude <- c(input$time_var, input$event_var)
    if (!is.null(input$strata) && input$strata != " ") exclude <- c(exclude, input$strata)
    exclude <- exclude[!is.null(exclude) & exclude != " "]
    updateSelectInput(session, "covars", choices = setdiff(names(df_raw()), exclude), selected = input$covars)
  }, ignoreInit = TRUE)
  
  observeEvent(input$covars, {
    req(input$time_var, input$event_var)
    exclude <- c(input$time_var, input$event_var)
    if (!is.null(input$covars) && length(input$covars) > 0) exclude <- c(exclude, input$covars)
    exclude <- exclude[!is.null(exclude) & exclude != " "]
    updateSelectInput(session, "strata", 
                      choices = c("Please select" = " ", setdiff(names(df_raw()), exclude)), 
                      selected = input$strata)
  }, ignoreInit = TRUE)
  
  # ---------------------------------------------------------------------------
  # Data Cleaning and Validation
  # ---------------------------------------------------------------------------
  
  clean_data <- reactive({
    req(df_raw(), input$time_var, input$event_var)
    
    df <- as.data.frame(df_raw())
    
    time_col <- df[[input$time_var]]
    validate(need(is.numeric(time_col), "Time variable must be numeric."))
    
    evt <- standardize_event(df[[input$event_var]])
    df[[input$event_var]] <- evt
    
    validate(need(all(evt %in% c(0,1), na.rm = TRUE), 
                  "Event must be binary (0/1) after conversion."))
    
    sel <- c(input$time_var, input$event_var)
    if (!is.null(input$covars) && length(input$covars) > 0) sel <- c(sel, input$covars)
    if (!is.null(input$strata) && input$strata != " ") sel <- c(sel, input$strata)
    
    sel <- unique(sel[!is.na(sel) & sel %in% names(df)])
    df <- df[, sel, drop = FALSE]
    
    n0 <- nrow(df)
    df <- na.omit(df)
    if (nrow(df) < n0) {
      shiny::showNotification(
        sprintf("Removed %d rows with missing data.", n0 - nrow(df)), 
        type = "warning", duration = 5
      )
    }
    
    n_ev <- sum(df[[input$event_var]] == 1)
    n_cv <- ifelse(is.null(input$covars), 0, length(input$covars))
    validate(
      need(n_cv == 0 || n_ev >= 10 * n_cv, 
           sprintf("Too few events. Need >=10 per covariate. Current: %d events, %d covariates.", 
                   n_ev, n_cv))
    )
    
    df
  })
  
  # ---------------------------------------------------------------------------
  # Download Handlers
  # ---------------------------------------------------------------------------
  
  output$dl_km <- downloadHandler(
    filename = function() paste0("km_curve_", Sys.Date(), ".png"),
    content = function(file) {
      d <- clean_data()
      dat <- d
      names(dat)[names(dat) == input$time_var]  <- ".time"
      names(dat)[names(dat) == input$event_var] <- ".event"
      
      if (!is.null(input$strata) && input$strata != " ") {
        names(dat)[names(dat) == input$strata] <- ".strata"
        dat$.strata <- as.factor(dat$.strata)
        fit <- survival::survfit(survival::Surv(.time, .event) ~ .strata, data = dat)
        clean_labels <- levels(dat$.strata)
        g <- survminer::ggsurvplot(
          fit, data = dat, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, ggtheme = theme_bw(),
          legend.title = input$strata, legend.labs = clean_labels
        )
      } else {
        fit <- survival::survfit(survival::Surv(.time, .event) ~ 1, data = dat)
        g <- survminer::ggsurvplot(
          fit, data = dat, risk.table = TRUE, pval = TRUE, 
          conf.int = TRUE, ggtheme = theme_bw()
        )
      }
      
      # 用 png() 设备，不要用 ggsave
      png(file, width = 8, height = 6, units = "in", res = 300)
      print(g)
      dev.off()
    },
    contentType = "image/png"
  )
  
  output$dl_cox <- downloadHandler(
    filename = function() paste0("cox_results_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- clean_data()
      d <- prepare_model_data(d, input$covars)
      f <- build_surv_formula(input$time_var, input$event_var, input$covars)
      fit <- survival::coxph(f, data = d)
      tb <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
      write.csv(tb, file, row.names = FALSE)
    }
  )
  
  output$dl_forest <- downloadHandler(
    filename = function() paste0("forest_plot_", Sys.Date(), ".png"),
    content = function(file) {
      d <- clean_data()
      d <- prepare_model_data(d, input$covars)
      f <- build_surv_formula(input$time_var, input$event_var, input$covars)
      fit <- survival::coxph(f, data = d)
      g <- survminer::ggforest(fit, data = d)
      
      png(file, width = 10, height = 8, units = "in", res = 300)
      print(g)
      dev.off()
    },
    contentType = "image/png"
  )
  
  output$dl_cutpoint <- downloadHandler(
    filename = function() paste0("cutpoint_km_", Sys.Date(), ".png"),
    content = function(file) {
      d <- clean_data()
      d2 <- d
      d2$.time  <- d2[[input$time_var]]
      d2$.event <- d2[[input$event_var]]
      
      var_type <- guess_var_type(d2[[input$cut_var]])
      
      if (var_type == "categorical") {
        d2$.group <- as.factor(d2[[input$cut_var]])
        fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = d2)
        g <- survminer::ggsurvplot(fit, data = d2, risk.table = TRUE, pval = TRUE, 
                                   conf.int = TRUE, ggtheme = theme_bw())
      } else {
        cp <- survminer::surv_cutpoint(
          d2, time = input$time_var, event = input$event_var,
          variables = input$cut_var, minprop = input$cut_minprop
        )
        cut_value <- as.numeric(cp$cutpoint[1, "cutpoint"])
        labs <- make_group_labels(d2[[input$cut_var]], cut_value)$labels
        d2$.group <- factor(ifelse(d2[[input$cut_var]] <= cut_value, labs[1], labs[2]), levels = labs)
        fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = d2)
        g <- survminer::ggsurvplot(fit, data = d2, risk.table = TRUE, pval = TRUE,
                                   conf.int = TRUE, ggtheme = theme_bw())
      }
      ggsave(file, plot = print(g), width = 8, height = 6, dpi = 300)
    }
  )
  
  output$dl_adj <- downloadHandler(
    filename = function() paste0("adjusted_curves_", Sys.Date(), ".png"),
    content = function(file) {
      d <- clean_data()
      d[[input$adj_var]] <- as.factor(d[[input$adj_var]])
      
      covars <- setdiff(input$covars %||% character(0), input$adj_var)
      for (v in covars) {
        if (v %in% names(d) && is.character(d[[v]])) d[[v]] <- as.factor(d[[v]])
      }
      
      f <- if (length(covars) > 0) {
        as.formula(paste0("Surv(`", input$time_var, "`, `", input$event_var, "`) ~ `",
                          input$adj_var, "` + ", paste(sprintf("`%s`", covars), collapse = " + ")))
      } else {
        as.formula(paste0("Surv(`", input$time_var, "`, `", input$event_var, "`) ~ `", input$adj_var, "`"))
      }
      
      fit <- survival::coxph(f, data = d)
      g <- survminer::ggadjustedcurves(fit, data = d, variable = input$adj_var, method = "average")
      ggsave(file, plot = g, width = 8, height = 6, dpi = 300)
    }
  )
  
  output$dl_cr <- downloadHandler(
    filename = function() paste0("competing_risks_", Sys.Date(), ".png"),
    content = function(file) {
      d <- as.data.frame(df_raw())
      time <- as.numeric(d[[input$cr_time_var]])
      status <- trimws(as.character(d[[input$cr_status_var]]))
      event_code <- trimws(as.character(input$cr_event_code))
      
      fstatus <- ifelse(status == event_code, 1L,
                        ifelse(status %in% c("0", "", "NA") | is.na(status), 0L, 2L))
      
      group <- if (input$cr_group_var != " ") factor(d[[input$cr_group_var]]) else NULL
      ci <- if (is.null(group)) {
        cmprsk::cuminc(ftime = time, fstatus = fstatus)
      } else {
        cmprsk::cuminc(ftime = time, fstatus = fstatus, group = group)
      }
      
      png(file, width = 8, height = 6, units = "in", res = 300)
      plot(ci, xlab = "Time", ylab = "Cumulative incidence", lwd = 2)
      dev.off()
    }
  )
  
  # ---------------------------------------------------------------------------
  # Tab UI Components
  # ---------------------------------------------------------------------------
  
  output$km_tab_ui <- renderUI({
    req(input$csv_file)
    tagList(
      plotOutput("km_plot", height = "500px"),
      uiOutput("km_interp"),
      br(),
      downloadButton("dl_km", "Download PNG", class = "btn-sm")
    )
  })
  
  output$logrank_tab_ui <- renderUI({
    req(input$csv_file)
    tagList(
      verbatimTextOutput("logrank_out"),
      uiOutput("logrank_interp")
    )
  })
  
  output$cox_tab_ui <- renderUI({
    req(input$csv_file)
    tagList(
      tableOutput("cox_table"),
      uiOutput("cox_interp"),
      br(),
      downloadButton("dl_cox", "Download Table (CSV)", class = "btn-sm")
    )
  })
  
  output$cox_diag_tab_ui <- renderUI({
    req(input$csv_file)
    sidebarLayout(
      sidebarPanel(
        helpText("Proportional hazards assumption diagnostics."),
        verbatimTextOutput("cox_diag_text")
      ),
      mainPanel(
        plotOutput("cox_diag_plot", height = "550px"),
        uiOutput("cox_diag_interp")
      )
    )
  })
  
  output$forest_tab_ui <- renderUI({
    req(input$csv_file)
    sidebarLayout(
      sidebarPanel(
        helpText("Forest plot of hazard ratios from the Cox model."),
        verbatimTextOutput("forest_note")
      ),
      mainPanel(
        plotOutput("forest_plot", height = "650px"),
        br(),
        uiOutput("forest_interp"),
        downloadButton("dl_forest", "Download PNG", class = "btn-sm")
      )
    )
  })
  
  output$cutpoint_tab_ui <- renderUI({
    req(input$csv_file)
    sidebarLayout(
      sidebarPanel(
        uiOutput("cut_var_ui"),
        numericInput("cut_minprop", "Minimum proportion in each group",
                     value = 0.10, min = 0.05, max = 0.40, step = 0.05)
      ),
      mainPanel(
        tableOutput("cutpoint_table"),
        plotOutput("cutpoint_km_plot", height = "550px"),
        uiOutput("cutpoint_interp"),
        br(),
        downloadButton("dl_cutpoint", "Download PNG", class = "btn-sm")
      )
    )
  })
  
  output$adj_tab_ui <- renderUI({
    req(input$csv_file)
    d <- as.data.frame(df_raw())
    cols <- names(d)
    
    sidebarLayout(
      sidebarPanel(
        selectInput("adj_var", "Group variable", choices = cols,
                    selected = if ("celltype" %in% cols) "celltype" else cols[1]),
        helpText("Adjusted survival curves using marginal approach.")
      ),
      mainPanel(
        plotOutput("adj_surv_plot", height = "600px"),
        uiOutput("adj_interp"),
        br(),
        downloadButton("dl_adj", "Download PNG", class = "btn-sm")
      )
    )
  })
  
  output$cr_tab_ui <- renderUI({
    req(input$csv_file)
    d <- as.data.frame(df_raw())
    cols <- names(d)
    
    sidebarLayout(
      sidebarPanel(
        selectInput("cr_time_var", "Time variable", choices = cols,
                    selected = if ("time" %in% cols) "time" else cols[1]),
        selectInput("cr_status_var", "Status variable", choices = cols,
                    selected = if ("status" %in% cols) "status" else cols[2]),
        selectInput("cr_group_var", "Group variable", choices = c("None" = " ", cols)),
        selectInput("cr_event_code", "Event of interest", choices = character(0))
      ),
      mainPanel(
        verbatimTextOutput("cr_gray_test"),
        plotOutput("cr_plot", height = "600px"),
        uiOutput("cr_interp"),
        br(),
        downloadButton("dl_cr", "Download PNG", class = "btn-sm")
      )
    )
  })
  
  output$ms_tab_ui <- renderUI({
    req(input$csv_file)
    d <- as.data.frame(df_raw())
    cols <- names(d)
    
    sidebarLayout(
      sidebarPanel(
        selectInput("ms_time_var", "Time variable", choices = cols,
                    selected = if ("time" %in% cols) "time" else cols[1]),
        selectInput("ms_state_var", "State variable", choices = cols,
                    selected = if ("status" %in% cols) "status" else cols[2]),
        selectInput("ms_censor_level", "Censoring level", choices = character(0)),
        selectInput("ms_group_var", "Group variable", choices = c("None" = " ", cols))
      ),
      mainPanel(
        plotOutput("ms_plot", height = "600px"),
        uiOutput("ms_interp")
      )
    )
  })
  
  output$help_tab_ui <- renderUI({
    tagList(
      h4("Quick Guide"),
      hr(),
      h5("Data Format"),
      p("Upload a CSV file with one row per subject. Required columns:"),
      tags$ul(
        tags$li(tags$b("Time:"), " Numeric follow-up time"),
        tags$li(tags$b("Event:"), " Binary indicator (0 = censored, 1 = event)"),
        tags$li(tags$b("Covariates:"), " Predictors (numeric or categorical)")
      ),
      h5("Analysis Modules"),
      tags$ul(
        tags$li(tags$b("KM Curve:"), " Kaplan-Meier survival estimates"),
        tags$li(tags$b("Log-Rank Test:"), " Compare survival between groups"),
        tags$li(tags$b("Cox Model:"), " Multivariable proportional hazards regression"),
        tags$li(tags$b("Forest Plot:"), " Visualize hazard ratios"),
        tags$li(tags$b("PH Diagnostics:"), " Check proportional hazards assumption"),
        tags$li(tags$b("Cutpoint Analysis:"), " Find optimal cutpoint for continuous variables"),
        tags$li(tags$b("Adjusted Curves:"), " Survival curves adjusted for covariates"),
        tags$li(tags$b("Competing Risks:"), " Cumulative incidence with competing events"),
        tags$li(tags$b("Multi-State:"), " State occupation probabilities")
      ),
      h5("Notes"),
      p("Missing data rows are automatically removed."),
      p("Cox model requires at least 10 events per covariate."),
      p("Results include automatic statistical interpretation.")
    )
  })
  
  # ---------------------------------------------------------------------------
  # Analysis Outputs - Kaplan-Meier and Log-Rank
  # ---------------------------------------------------------------------------
  
  output$km_plot <- renderPlot({
    req(clean_data())
    
    d <- clean_data()
    dat <- d
    names(dat)[names(dat) == input$time_var]  <- ".time"
    names(dat)[names(dat) == input$event_var] <- ".event"
    
    if (!is.null(input$strata) && input$strata != " ") {
      names(dat)[names(dat) == input$strata] <- ".strata"
      dat$.strata <- as.factor(dat$.strata)
      
      fit <- survival::survfit(survival::Surv(.time, .event) ~ .strata, data = dat)
      clean_labels <- levels(dat$.strata)
      
      g <- survminer::ggsurvplot(
        fit, data = dat, risk.table = TRUE, pval = TRUE,
        conf.int = TRUE, ggtheme = theme_bw(),
        legend.title = input$strata, legend.labs = clean_labels
      )
    } else {
      fit <- survival::survfit(survival::Surv(.time, .event) ~ 1, data = dat)
      g <- survminer::ggsurvplot(
        fit, data = dat, risk.table = TRUE, pval = TRUE,
        conf.int = TRUE, ggtheme = theme_bw()
      )
    }
    print(g)
  })
  
  output$km_interp <- renderUI({
    req(clean_data())
    
    d <- clean_data()
    dat <- d
    names(dat)[names(dat) == input$time_var]  <- ".time"
    names(dat)[names(dat) == input$event_var] <- ".event"
    
    if (!is.null(input$strata) && input$strata != " ") {
      names(dat)[names(dat) == input$strata] <- ".strata"
      dat$.strata <- as.factor(dat$.strata)
      fit <- survival::survfit(survival::Surv(.time, .event) ~ .strata, data = dat)
      s <- summary(fit)$table
      
      if (is.matrix(s) || is.data.frame(s)) {
        raw_names <- rownames(s)
        clean_names <- gsub("^\\.strata=", "", raw_names)
        
        txt <- character()
        for (i in seq_along(clean_names)) {
          med <- if ("median" %in% colnames(s)) s[i, "median"] else NA
          txt <- c(txt, paste0("<b>", clean_names[i], "</b>: median survival = ", fmt_num(med)))
        }
        
        HTML(paste0(
          "<b>Interpretation:</b> The KM curves describe survival patterns across groups defined by <b>",
          input$strata, "</b>.<br>", paste(txt, collapse = "<br>")
        ))
      } else {
        HTML(paste0("<b>Interpretation:</b> The KM curve summarizes time-to-event patterns across levels of <b>",
                    input$strata, "</b>."))
      }
    } else {
      fit <- survival::survfit(survival::Surv(.time, .event) ~ 1, data = dat)
      med <- summary(fit)$table["median"]
      HTML(paste0(
        "<b>Interpretation:</b> This KM curve summarizes the overall survival experience. ",
        "Estimated median survival: <b>", fmt_num(med), "</b>."
      ))
    }
  })
  
  output$logrank_out <- renderPrint({
    req(clean_data(), !is.null(input$strata) && input$strata != " ")
    
    d <- clean_data()
    f <- as.formula(paste0("survival::Surv(", input$time_var, ", ", input$event_var, ") ~ ", input$strata))
    res <- survival::survdiff(f, data = d)
    pval <- 1 - pchisq(res$chisq, length(res$n) - 1)
    
    cat("Log-Rank Test Results:\n")
    print(res)
    cat("\nP-value:", format.pval(pval, digits = 3), "\n")
  })
  
  output$logrank_interp <- renderUI({
    req(clean_data(), !is.null(input$strata) && input$strata != " ")
    
    d <- clean_data()
    d$strata_group <- as.factor(d[[input$strata]])
    
    f <- as.formula(paste0("Surv(", input$time_var, ", ", input$event_var, ") ~ strata_group"))
    res <- survival::survdiff(f, data = d)
    pval <- 1 - pchisq(res$chisq, length(res$n) - 1)
    
    msg <- if (pval < 0.05) {
      paste0("Survival differs across <b>", input$strata, "</b> groups (p = ", fmt_p(pval), ").")
    } else {
      paste0("No significant survival difference across <b>", input$strata, "</b> groups (p = ", fmt_p(pval), ").")
    }
    HTML(paste0("<b>Interpretation:</b> ", msg))
  })
  
  # ---------------------------------------------------------------------------
  # Analysis Outputs - Cox Model
  # ---------------------------------------------------------------------------
  
  output$cox_table <- renderTable({
    req(clean_data(), !is.null(input$covars) && length(input$covars) > 0)
    
    d <- clean_data()
    d <- prepare_model_data(d, input$covars)
    f <- build_surv_formula(input$time_var, input$event_var, input$covars)
    
    fit <- tryCatch(
      survival::coxph(f, data = d),
      error = function(e) validate(need(FALSE, paste("Cox model failed:", e$message)))
    )
    
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
      dplyr::rename(HR = estimate, `95% CI Low` = conf.low, `95% CI High` = conf.high)
  }, digits = 3)
  
  output$cox_interp <- renderUI({
    req(clean_data(), !is.null(input$covars) && length(input$covars) > 0)
    
    d <- clean_data()
    for (v in input$covars) {
      if (v %in% names(d) && is.character(d[[v]])) d[[v]] <- as.factor(d[[v]])
    }
    
    rhs <- paste(input$covars, collapse = " + ")
    f <- as.formula(paste0("Surv(", input$time_var, ", ", input$event_var, ") ~ ", rhs))
    
    fit <- tryCatch(
      survival::coxph(f, data = d),
      error = function(e) return(HTML(paste0("<b>Interpretation:</b> Cox model failed: ", e$message)))
    )
    
    if (inherits(fit, "html")) return(fit)
    
    tb <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    sig <- tb[tb$p.value < 0.05, , drop = FALSE]
    nsig <- tb[tb$p.value >= 0.05, , drop = FALSE]
    
    if (nrow(sig) == 0) {
      return(HTML("<b>Interpretation:</b> No covariate is statistically significant at the 0.05 level."))
    }
    
    lines <- character()
    for (i in 1:nrow(sig)) {
      term <- sig$term[i]
      hr <- as.numeric(sig$estimate[i])
      low <- as.numeric(sig$conf.low[i])
      high <- as.numeric(sig$conf.high[i])
      p <- as.numeric(sig$p.value[i])
      risk <- if (hr > 1) "higher" else "lower"
      
      lines <- c(lines, paste0(
        "<li><b>", term, "</b>: HR = ", fmt_num(hr),
        " (95% CI: ", fmt_num(low), "-", fmt_num(high),
        "), p = ", fmt_p(p), " \u2014 associated with ", risk, " hazard.</li>"
      ))
    }
    
    ns_lines <- if (nrow(nsig) > 0) {
      paste0("<br><small>Non-significant (p >= 0.05): ", 
             paste(nsig$term, collapse = ", "), ".</small>")
    } else ""
    
    HTML(paste0(
      "<b>Interpretation:</b> Significant covariates in the multivariable Cox model:<ul>",
      paste(lines, collapse = ""), "</ul>", ns_lines
    ))
  })
  
  # ---------------------------------------------------------------------------
  # Analysis Outputs - Advanced Features
  # ---------------------------------------------------------------------------
  
  output$cut_var_ui <- renderUI({
    req(clean_data())
    d <- clean_data()
    num_vars <- names(d)[vapply(d, is.numeric, logical(1))]
    num_vars <- setdiff(num_vars, c(input$time_var, input$event_var))
    
    if (length(num_vars) == 0) return(helpText("No numeric covariate available."))
    selectInput("cut_var", "Continuous variable", choices = num_vars, selected = num_vars[1])
  })
  
  output$cox_diag_text <- renderPrint({
    req(clean_data())
    validate(need(!is.null(input$covars) && length(input$covars) > 0,
                  "Please select at least one covariate first."))
    
    d <- clean_data()
    for (v in input$covars) {
      if (v %in% names(d) && is.character(d[[v]])) d[[v]] <- as.factor(d[[v]])
    }
    
    rhs <- paste(input$covars, collapse = " + ")
    f <- as.formula(paste0("survival::Surv(", input$time_var, ", ", input$event_var, ") ~ ", rhs))
    fit <- survival::coxph(f, data = d)
    print(survival::cox.zph(fit))
  })
  
  output$cox_diag_plot <- renderPlot({
    req(clean_data())
    validate(need(!is.null(input$covars) && length(input$covars) > 0,
                  "Please select at least one covariate first."))
    
    d <- clean_data()
    for (v in input$covars) {
      if (v %in% names(d) && is.character(d[[v]])) d[[v]] <- as.factor(d[[v]])
    }
    
    rhs <- paste(input$covars, collapse = " + ")
    f <- as.formula(paste0("survival::Surv(", input$time_var, ", ", input$event_var, ") ~ ", rhs))
    fit <- survival::coxph(f, data = d)
    print(survminer::ggcoxzph(survival::cox.zph(fit)))
  })
  
  output$cox_diag_interp <- renderUI({
    req(clean_data())
    validate(need(!is.null(input$covars) && length(input$covars) > 0,
                  "Please select at least one covariate first."))
    
    d <- clean_data()
    for (v in input$covars) {
      if (v %in% names(d) && is.character(d[[v]])) d[[v]] <- as.factor(d[[v]])
    }
    
    rhs <- paste(input$covars, collapse = " + ")
    f <- as.formula(paste0("survival::Surv(", input$time_var, ", ", input$event_var, ") ~ ", rhs))
    
    fit <- tryCatch(survival::coxph(f, data = d), error = function(e) NULL)
    if (is.null(fit)) return(HTML("<b>Interpretation:</b> Cox model fitting failed."))
    
    zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
    if (is.null(zph)) return(HTML("<b>Interpretation:</b> PH test could not be computed."))
    
    tb <- as.data.frame(zph$table)
    global_p <- if ("GLOBAL" %in% rownames(tb)) tb["GLOBAL", "p"] else NA
    
    var_tb <- tb[rownames(tb) != "GLOBAL", , drop = FALSE]
    bad_vars <- if (nrow(var_tb) > 0 && "p" %in% colnames(var_tb)) {
      rownames(var_tb)[var_tb$p < 0.05]
    } else character(0)
    
    msg1 <- if (!is.na(global_p)) {
      if (global_p < 0.05) {
        paste0("The <b>global test</b> suggests PH violation (p = ", fmt_p(global_p), "). ")
      } else {
        paste0("The <b>global test</b> does not show strong evidence against PH (p = ", fmt_p(global_p), "). ")
      }
    } else ""
    
    msg2 <- if (length(bad_vars) > 0) {
      paste0("Potential time-varying effects: <b>", paste(bad_vars, collapse = ", "), "</b>.")
    } else if (nrow(var_tb) > 0) {
      "No individual covariate shows significant PH violation."
    } else ""
    
    HTML(paste0("<b>Interpretation:</b> ", msg1, msg2))
  })
  
  output$forest_note <- renderPrint({
    cat("Forest plot from current multivariable Cox model.\n")
  })
  
  output$forest_plot <- renderPlot({
    req(clean_data())
    validate(need(!is.null(input$covars) && length(input$covars) > 0,
                  "Please select at least one covariate first."))
    
    d <- clean_data()
    d <- prepare_model_data(d, input$covars)
    f <- build_surv_formula(input$time_var, input$event_var, input$covars)
    fit <- survival::coxph(f, data = d)
    print(survminer::ggforest(fit, data = d))
  })
  
  output$forest_interp <- renderUI({
    req(clean_data(), !is.null(input$covars) && length(input$covars) > 0)
    HTML(paste0(
      "<b>Interpretation:</b> Forest plot of hazard ratios from the Cox model. ",
      "HR > 1 indicates worse survival; HR < 1 indicates better survival. ",
      "Confidence intervals crossing 1 are not statistically significant."
    ))
  })
  
  output$cutpoint_table <- renderTable({
    req(clean_data(), input$cut_var)
    d <- clean_data()
    validate(need(input$cut_var %in% names(d), "Invalid cutpoint variable."))
    validate(need(is.numeric(d[[input$cut_var]]), "Variable must be numeric."))
    
    cp <- survminer::surv_cutpoint(
      d, time = input$time_var, event = input$event_var,
      variables = input$cut_var, minprop = input$cut_minprop
    )
    as.data.frame(cp$cutpoint)
  }, digits = 3)
  
  output$cutpoint_km_plot <- renderPlot({
    req(clean_data(), input$cut_var)
    d <- clean_data()
    validate(need(input$cut_var %in% names(d), "Invalid cutpoint variable."))
    
    d2 <- d
    d2$.time  <- d2[[input$time_var]]
    d2$.event <- d2[[input$event_var]]
    
    var_type <- guess_var_type(d2[[input$cut_var]])
    
    if (var_type == "categorical") {
      d2$.group <- as.factor(d2[[input$cut_var]])
      fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = d2)
      clean_labels <- levels(d2$.group)
      
      g <- survminer::ggsurvplot(
        fit, data = d2, risk.table = TRUE, pval = TRUE,
        conf.int = TRUE, ggtheme = theme_bw(),
        legend.title = input$cut_var, legend.labs = clean_labels
      )
      print(g)
      return()
    }
    
    validate(need(is.numeric(d2[[input$cut_var]]), "Variable must be numeric."))
    
    cp <- survminer::surv_cutpoint(
      d2, time = input$time_var, event = input$event_var,
      variables = input$cut_var, minprop = input$cut_minprop
    )
    cut_value <- as.numeric(cp$cutpoint[1, "cutpoint"])
    labs <- make_group_labels(d2[[input$cut_var]], cut_value)$labels
    
    d2$.group <- factor(
      ifelse(d2[[input$cut_var]] <= cut_value, labs[1], labs[2]),
      levels = labs
    )
    fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = d2)
    
    g <- survminer::ggsurvplot(
      fit, data = d2, risk.table = TRUE, pval = TRUE,
      conf.int = TRUE, ggtheme = theme_bw(),
      legend.title = input$cut_var, legend.labs = labs
    )
    print(g)
  })
  
  output$cutpoint_interp <- renderUI({
    req(clean_data(), input$cut_var)
    d <- clean_data()
    validate(need(input$cut_var %in% names(d), "Invalid cutpoint variable."))
    
    var_type <- guess_var_type(d[[input$cut_var]])
    
    if (var_type == "categorical") {
      return(HTML(paste0("<b>Interpretation:</b> <b>", input$cut_var,
                         "</b> is categorical \u2014 groups shown are original categories.")))
    }
    
    cp <- survminer::surv_cutpoint(
      d, time = input$time_var, event = input$event_var,
      variables = input$cut_var, minprop = input$cut_minprop
    )
    cut_value <- as.numeric(cp$cutpoint[1, "cutpoint"])
    
    x <- d[[input$cut_var]]
    n_low <- sum(x <= cut_value, na.rm = TRUE)
    n_high <- sum(x > cut_value, na.rm = TRUE)
    pct_low <- round(100 * n_low / length(x), 1)
    pct_high <- round(100 * n_high / length(x), 1)
    
    HTML(paste0(
      "<b>Interpretation:</b> Optimal cutpoint for <b>", input$cut_var, "</b> is <b>", fmt_num(cut_value), "</b>. ",
      "Splits into <b>Not greater than ", fmt_num(cut_value), "</b> (n = ", n_low, ", ", pct_low, "%) ",
      "and <b>Greater than ", fmt_num(cut_value), "</b> (n = ", n_high, ", ", pct_high, "%)."
    ))
  })
  
  output$adj_surv_plot <- renderPlot({
    req(clean_data(), input$adj_var)
    d <- clean_data()
    validate(need(input$adj_var %in% names(d), "Invalid group variable."))
    validate(need(guess_var_type(d[[input$adj_var]]) == "categorical",
                  "Choose a categorical variable."))
    
    d[[input$adj_var]] <- as.factor(d[[input$adj_var]])
    covars <- setdiff(input$covars %||% character(0), input$adj_var)
    for (v in covars) {
      if (v %in% names(d) && is.character(d[[v]])) d[[v]] <- as.factor(d[[v]])
    }
    
    f <- if (length(covars) > 0) {
      as.formula(paste0("Surv(`", input$time_var, "`, `", input$event_var, "`) ~ `",
                        input$adj_var, "` + ", paste(sprintf("`%s`", covars), collapse = " + ")))
    } else {
      as.formula(paste0("Surv(`", input$time_var, "`, `", input$event_var, "`) ~ `", input$adj_var, "`"))
    }
    
    fit <- survival::coxph(f, data = d)
    print(survminer::ggadjustedcurves(fit, data = d, variable = input$adj_var, method = "average"))
  })
  
  output$adj_interp <- renderUI({
    req(clean_data(), input$adj_var)
    d <- clean_data()
    validate(need(input$adj_var %in% names(d), "Invalid group variable."))
    
    covars <- setdiff(input$covars %||% character(0), input$adj_var)
    
    if (length(covars) == 0) {
      return(HTML(paste0("<b>Interpretation:</b> Unadjusted survival curves for <b>", input$adj_var, "</b>.")))
    }
    
    HTML(paste0(
      "<b>Interpretation:</b> Adjusted survival curves for <b>", input$adj_var, "</b>, ",
      "marginalizing over: <b>", paste(covars, collapse = ", "), "</b>."
    ))
  })
  
  # ---------------------------------------------------------------------------
  # Competing Risks and Multi-State
  # ---------------------------------------------------------------------------
  
  observeEvent(input$cr_status_var, {
    req(df_raw(), input$cr_status_var)
    d <- as.data.frame(df_raw())
    
    status_vec <- d[[input$cr_status_var]]
    status_clean <- trimws(as.character(status_vec))
    unique_codes <- sort(unique(status_clean[!is.na(status_clean) & status_clean != ""]))
    
    numeric_codes <- suppressWarnings(as.numeric(unique_codes))
    
    if (!anyNA(numeric_codes) && length(numeric_codes) >= 2) {
      sorted_codes <- unique_codes[order(numeric_codes)]
      event_choices <- sorted_codes[-1]
    } else {
      event_choices <- unique_codes
    }
    
    if (length(event_choices) > 0) {
      updateSelectInput(session, "cr_event_code", 
                        choices = event_choices, selected = event_choices[1])
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$ms_state_var, {
    req(df_raw(), input$ms_state_var)
    d <- as.data.frame(df_raw())
    levs <- sort(unique(trimws(as.character(d[[input$ms_state_var]]))))
    levs <- levs[!is.na(levs) & levs != ""]
    if (length(levs) > 0) {
      updateSelectInput(session, "ms_censor_level", choices = levs, selected = levs[1])
    }
  }, ignoreInit = TRUE)
  
  output$cr_gray_test <- renderPrint({
    req(df_raw(), input$cr_time_var, input$cr_status_var, input$cr_group_var, input$cr_event_code)
    d <- as.data.frame(df_raw())
    
    time <- suppressWarnings(as.numeric(d[[input$cr_time_var]]))
    validate(need(!anyNA(time), "Time must be numeric."))
    
    status <- trimws(as.character(d[[input$cr_status_var]]))
    event_code <- trimws(as.character(input$cr_event_code))
    
    fstatus <- ifelse(status == event_code, 1L,
                      ifelse(status %in% c("0", "", "NA") | is.na(status), 0L, 2L))
    
    group <- if (input$cr_group_var != " ") factor(d[[input$cr_group_var]]) else NULL
    
    ci <- if (is.null(group)) {
      cmprsk::cuminc(ftime = time, fstatus = fstatus)
    } else {
      cmprsk::cuminc(ftime = time, fstatus = fstatus, group = group)
    }
    
    if (!is.null(ci$Tests)) {
      print(ci$Tests)
    } else {
      cat("No group comparison was run.\n")
    }
  })
  
  output$cr_plot <- renderPlot({
    req(df_raw(), input$cr_time_var, input$cr_status_var, input$cr_group_var, input$cr_event_code)
    d <- as.data.frame(df_raw())
    
    time <- suppressWarnings(as.numeric(d[[input$cr_time_var]]))
    validate(need(!anyNA(time), "Time must be numeric."))
    
    status <- trimws(as.character(d[[input$cr_status_var]]))
    event_code <- trimws(as.character(input$cr_event_code))
    
    unique_codes <- sort(unique(status[!is.na(status) & status != ""]))
    numeric_codes <- suppressWarnings(as.numeric(unique_codes))
    
    if (!anyNA(numeric_codes) && length(numeric_codes) >= 2) {
      sorted_codes <- unique_codes[order(numeric_codes)]
      censor_codes <- sorted_codes[1]
    } else {
      censor_codes <- c("0", "censor", "censored", "alive")
    }
    
    fstatus <- ifelse(
      status == event_code, 1L,
      ifelse(status %in% censor_codes | status == "" | is.na(status), 0L, 2L)
    )
    
    group <- if (input$cr_group_var != " ") factor(d[[input$cr_group_var]]) else NULL
    
    ci <- if (is.null(group)) {
      cmprsk::cuminc(ftime = time, fstatus = fstatus)
    } else {
      cmprsk::cuminc(ftime = time, fstatus = fstatus, group = group)
    }
    
    plot(ci, xlab = "Time", ylab = "Cumulative incidence", lwd = 2)
  })
  
  output$cr_interp <- renderUI({
    req(df_raw(), input$cr_time_var, input$cr_status_var, input$cr_event_code)
    d <- as.data.frame(df_raw())
    
    time <- suppressWarnings(as.numeric(d[[input$cr_time_var]]))
    status <- trimws(as.character(d[[input$cr_status_var]]))
    event_code <- trimws(as.character(input$cr_event_code))
    
    n_event <- sum(status == event_code, na.rm = TRUE)
    n_other <- sum(!status %in% c(event_code, "0", "", NA) & !is.na(status), na.rm = TRUE)
    n_censor <- sum(status %in% c("0", "") | is.na(status), na.rm = TRUE)
    
    base_msg <- paste0(
      "<b>Interpretation:</b> Cumulative incidence of event '<b>", event_code, "</b>' (n = ", n_event, "). "
    )
    
    if (input$cr_group_var == " ") {
      msg <- paste0(base_msg, "Competing events: n = ", n_other, "; Censored: n = ", n_censor, ".")
      return(HTML(msg))
    }
    
    group <- factor(d[[input$cr_group_var]])
    fstatus <- ifelse(status == event_code, 1L,
                      ifelse(status %in% c("0", "", "NA") | is.na(status), 0L, 2L))
    
    ci <- tryCatch(
      cmprsk::cuminc(ftime = time, fstatus = fstatus, group = group),
      error = function(e) NULL
    )
    
    if (!is.null(ci) && !is.null(ci$Tests)) {
      pval <- ci$Tests[1, "pv"]
      test_msg <- if (!is.na(pval) && pval < 0.05) {
        paste0("Gray's test: significant difference across <b>", input$cr_group_var, "</b> (p = ", fmt_p(pval), ").")
      } else {
        paste0("Gray's test: no significant difference (p = ", fmt_p(pval), ").")
      }
      msg <- paste0(base_msg, test_msg)
    } else {
      msg <- paste0(base_msg, "Cumulative incidence curves shown for each group.")
    }
    
    HTML(msg)
  })
  
  output$ms_plot <- renderPlot({
    req(df_raw(), input$ms_time_var, input$ms_state_var, input$ms_censor_level)
    d <- as.data.frame(df_raw())
    
    time <- suppressWarnings(as.numeric(d[[input$ms_time_var]]))
    validate(need(!anyNA(time), "Time must be numeric."))
    
    state <- trimws(as.character(d[[input$ms_state_var]]))
    censor_level <- trimws(as.character(input$ms_censor_level))
    
    validate(need(censor_level %in% unique(state), "Invalid censoring level."))
    
    state <- factor(state)
    state <- stats::relevel(state, ref = censor_level)
    validate(need(nlevels(state) >= 3, "Multi-state needs at least 3 levels."))
    
    d$.time <- time
    d$.state <- state
    
    if (input$ms_group_var != " ") {
      validate(need(input$ms_group_var %in% names(d), "Invalid group variable."))
      d$.group <- factor(d[[input$ms_group_var]])
      f <- survival::Surv(.time, .state) ~ .group
    } else {
      f <- survival::Surv(.time, .state) ~ 1
    }
    
    fit <- survival::survfit(f, data = d)
    
    plot(fit, xlab = "Time", ylab = "State occupation probability",
         mark.time = TRUE, lty = 1, col = 1:ncol(summary(fit)$pstate))
    
    if (!is.null(fit$states)) {
      legend("topright", legend = fit$states, lty = 1, bty = "n")
    }
  })
  
  output$ms_interp <- renderUI({
    req(df_raw(), input$ms_time_var, input$ms_state_var, input$ms_censor_level)
    HTML(paste0("<b>Interpretation:</b> State occupation probabilities over time",
                if (input$ms_group_var != " ") paste0(" across levels of <b>", input$ms_group_var, "</b>") else "",
                "."))
  })
}