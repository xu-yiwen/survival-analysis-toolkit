# =============================================================================
# global.R - Survival Analysis Toolkit
# =============================================================================

library(shiny)
library(readr)
library(survival)
library(survminer)
library(broom)
library(dplyr)
library(ggplot2)

options(shiny.maxRequestSize = 50 * 1024^2)

# -----------------------------------------------------------------------------
# Utility Functions
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Standardize event variable to binary 0/1 (1 = event, 0 = censored)
standardize_event <- function(x) {
  x <- trimws(tolower(as.character(x)))
  
  if (all(na.omit(x) %in% c("0", "1"))) {
    return(as.integer(x == "1"))
  }
  
  if (all(na.omit(x) %in% c("1", "2"))) {
    return(as.integer(x == "2"))
  }
  
  ifelse(x %in% c("yes", "true", "dead", "event", "e"), 1L, 0L)
}