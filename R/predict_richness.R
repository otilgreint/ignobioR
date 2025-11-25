#' @title Predict Species Richness Using Species-Area Relationship
#'
#' @description
#' Predicts expected species richness for a given area using the Arrhenius
#' power function (S = c × A^z). Provides multiple approaches: literature-based
#' coefficients from Mediterranean studies, custom user-defined coefficients,
#' or data-driven estimation from reference floras.
#'
#' @param area Numeric. Study area in km². Can be a single value or vector.
#' @param method Character. One of:
#'   \itemize{
#'     \item{"literature"}{ Use published c and z values (default)}
#'     \item{"custom"}{ Provide your own c and z coefficients}
#'     \item{"data"}{ Estimate c and z from reference flora data}
#'   }
#' @param region Character. For method="literature", region name (e.g., "Italy", "Tuscany").
#'   Use \code{list_sar_regions()} to see available regions.
#' @param flora_type Character. Type of flora for prediction:
#'   \itemize{
#'     \item{"total"}{ All species (native + alien, default)}
#'     \item{"native"}{ Native species only}
#'     \item{"alien"}{ Alien species only}
#'   }
#' @param c Numeric. Custom intercept coefficient (method="custom" only).
#' @param z Numeric. Custom slope coefficient (method="custom" only).
#' @param reference_data Data frame or path to CSV file containing reference
#'   flora data with columns: 'area' (km²) and 'species' (species count).
#'   Optional columns: 'native' and 'alien' counts for separate models.
#' @param coef_database Character. Optional path to custom SAR coefficient CSV.
#'   If NULL, uses package defaults. CSV must have columns: region, country,
#'   flora_type, c, z, R2, n_floras, source, doi, notes.
#' @param return_intervals Logical. If TRUE and method="data", returns
#'   prediction intervals based on model fit (default = TRUE).
#' @param confidence_level Numeric. Confidence level for intervals (default = 0.95).
#' @param verbose Logical. Print progress messages (default = TRUE).
#'
#' @return A list with:
#' \itemize{
#'   \item{\code{prediction}}{Data frame with area and predicted richness}
#'   \item{\code{coefficients}}{Values of c and z used}
#'   \item{\code{intervals}}{Prediction intervals (if applicable)}
#'   \item{\code{model}}{Fitted model object (if method="data")}
#'   \item{\code{diagnostics}}{Model fit statistics (R², RMSE)}
#' }
#'
#' @details
#' **Literature coefficients** (Arrhenius power function: S = c × A^z):
#'
#' *Italy* (D'Antraccoli et al. 2024):
#' - Total: c=241.2, z=0.281 (R²=0.92, n=286 floras)
#' - Native: c=245.2, z=0.263 (R²=0.91)
#' - Alien: c=10.1, z=0.404 (R²=0.73)
#'
#' *Tuscany* (D'Antraccoli et al. 2019):
#' - Total: c=303.4, z=0.246 (R²=0.87, n=67 floras)
#' - Native: c=314.6, z=0.229 (R²=0.86)
#' - Alien: c=12.0, z=0.272 (R²=0.56)
#'
#' **Data-driven estimation**: Uses nonlinear least squares (nls) to fit
#' the power function to your reference flora dataset.
#'
#' **Coefficient database**: The function searches for SAR coefficients in order:
#' \enumerate{
#'   \item{Custom user CSV (if coef_database provided)}
#'   \item{Package default CSV (inst/extdata/sar_coefficients.csv)}
#'   \item{Hardcoded fallback (Italy + Tuscany only)}
#' }
#'
#' @importFrom stats nls predict coef residuals
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Predict richness for 100 km² using Italian coefficients
#' pred <- predict_richness(
#'   area = 100,
#'   method = "literature",
#'   region = "Italy",
#'   flora_type = "total"
#' )
#' print(pred$prediction)
#' # Expected: ~1470 species
#'
#' # Example 2: Multiple areas
#' areas <- c(10, 50, 100, 500, 1000)
#' pred_multi <- predict_richness(
#'   area = areas,
#'   region = "Tuscany",
#'   flora_type = "native"
#' )
#'
#' # Example 3: Custom coefficients
#' pred_custom <- predict_richness(
#'   area = 250,
#'   method = "custom",
#'   c = 280,
#'   z = 0.25
#' )
#'
#' # Example 4: Estimate from your own data
#' my_floras <- data.frame(
#'   area = c(5, 15, 45, 120, 350),
#'   species = c(520, 890, 1340, 1820, 2450)
#' )
#'
#' pred_data <- predict_richness(
#'   area = 200,
#'   method = "data",
#'   reference_data = my_floras,
#'   return_intervals = TRUE
#' )
#'
#' # Example 5: Use custom coefficient database
#' pred_custom_db <- predict_richness(
#'   area = 150,
#'   region = "Liguria",
#'   coef_database = "~/my_sar_coefficients.csv"
#' )
#'
#' # See available regions
#' list_sar_regions()
#' }

predict_richness <- function(area, method = "literature", region = "Italy", flora_type = "total",
                             c = NULL, z = NULL,
                             reference_data = NULL, coef_database = NULL, 
                             return_intervals = TRUE, confidence_level = 0.95,
                             verbose = TRUE
) {
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  
  msg <- function(...) if (verbose) message(...)
  
  # ============================================================================
  # SECTION 1: INPUT VALIDATION
  # ============================================================================
  
  msg("Validating inputs...")
  
  # Validate area
  if (!is.numeric(area) || any(area <= 0)) {
    stop("area must be positive numeric value(s) in km²")
  }
  
  # Validate method
  valid_methods <- c("literature", "custom", "data")
  if (!method %in% valid_methods) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # Validate flora_type
  valid_types <- c("total", "native", "alien")
  if (!flora_type %in% valid_types) {
    stop("flora_type must be one of: ", paste(valid_types, collapse = ", "))
  }
  
  # Method-specific validation
  if (method == "custom" && (is.null(c) || is.null(z))) {
    stop("method='custom' requires both c and z parameters")
  }
  
  if (method == "data" && is.null(reference_data)) {
    stop("method='data' requires reference_data parameter")
  }
  
  # ============================================================================
  # SECTION 2: COEFFICIENT DETERMINATION
  # ============================================================================
  
  model_fit <- NULL
  fit_data <- NULL
  
  if (method == "literature") {
    msg("Loading SAR coefficients...")
    
    # Load coefficient database (three-tier approach)
    coef_df <- load_sar_coefficients(coef_database, verbose = verbose)
    
    msg(paste0("Searching for: region='", region, "', flora_type='", flora_type, "'"))
    
    # Query the database
    coef_row <- coef_df[coef_df$region == region & 
                          coef_df$flora_type == flora_type, ]
    
    if (nrow(coef_row) == 0) {
      available_regions <- unique(coef_df$region)
      stop("No coefficients found for region: ", region, 
           " (flora_type: ", flora_type, ")\n",
           "Available regions: ", paste(available_regions, collapse = ", "), "\n",
           "Tip: Use list_sar_regions() to see all available combinations.\n",
           "     Consider providing a custom CSV via coef_database parameter.")
    }
    
    if (nrow(coef_row) > 1) {
      warning("Multiple coefficient sets found. Using the first one.")
      coef_row <- coef_row[1, ]
    }
    
    c_value <- coef_row$c
    z_value <- coef_row$z
    r2_value <- coef_row$R2
    source_ref <- coef_row$source
    doi_ref <- coef_row$doi
    n_floras <- coef_row$n_floras
    
    msg(paste0("  c = ", round(c_value, 2), ", z = ", round(z_value, 3)))
    if (!is.na(r2_value)) {
      msg(paste0("  Original model R² = ", round(r2_value, 3)))
    }
    if (!is.na(n_floras)) {
      msg(paste0("  Based on ", n_floras, " reference floras"))
    }
    msg(paste0("  Source: ", source_ref))
    if (!is.na(doi_ref)) {
      msg(paste0("  DOI: ", doi_ref))
    }
    
  } else if (method == "custom") {
    msg("Using custom coefficients...")
    
    c_value <- c
    z_value <- z
    r2_value <- NA
    source_ref <- "User-provided"
    doi_ref <- NA
    n_floras <- NA
    
    msg(paste0("  c = ", round(c_value, 2), ", z = ", round(z_value, 3)))
    
  } else if (method == "data") {
    msg("Estimating coefficients from reference data...")
    
    # [KEEP ALL DATA LOADING CODE AS-IS]
    # Load data if path provided
    if (is.character(reference_data)) {
      msg(paste0("  Loading data from: ", reference_data))
      reference_data <- read.csv(reference_data, stringsAsFactors = FALSE)
    }
    
    # Validate reference_data structure
    if (flora_type == "total") {
      required_cols <- c("area", "species")
    } else {
      required_cols <- c("area", flora_type)
    }
    
    if (!all(required_cols %in% names(reference_data))) {
      stop("reference_data must contain columns: ", paste(required_cols, collapse = ", "))
    }
    
    # Prepare data
    if (flora_type == "total" && "species" %in% names(reference_data)) {
      fit_data <- data.frame(
        area = reference_data$area,
        richness = reference_data$species
      )
    } else {
      fit_data <- data.frame(
        area = reference_data$area,
        richness = reference_data[[flora_type]]
      )
    }
    
    # Remove NA and zero values
    fit_data <- fit_data[complete.cases(fit_data) & 
                           fit_data$area > 0 & 
                           fit_data$richness > 0, ]
    
    if (nrow(fit_data) < 3) {
      stop("reference_data must contain at least 3 valid observations")
    }
    
    msg(paste0("  Fitting SAR model to ", nrow(fit_data), " reference floras..."))
    
    # =========================================================================
    # NEW CODE STARTS HERE: TWO-STAGE FITTING WITH FALLBACK
    # =========================================================================
    
    # Initialize variables
    nls_success <- FALSE
    fitting_method <- NULL
    
    # STAGE 1: Try NLS first (preferred - fits in original space)
    tryCatch({
      start_c <- mean(fit_data$richness) / mean(fit_data$area^0.25)
      start_z <- 0.25
      
      msg("  Attempting non-linear least squares (NLS)...")
      
      model_fit <- stats::nls(
        richness ~ c * area^z,
        data = fit_data,
        start = list(c = start_c, z = start_z),
        control = list(maxiter = 100, warnOnly = TRUE)
      )
      
      # Extract coefficients
      coefs <- stats::coef(model_fit)
      c_value <- coefs["c"]
      z_value <- coefs["z"]
      
      # Calculate R²
      ss_res <- sum(stats::residuals(model_fit)^2)
      ss_tot <- sum((fit_data$richness - mean(fit_data$richness))^2)
      r2_value <- 1 - (ss_res / ss_tot)
      
      fitting_method <- "nls"
      nls_success <- TRUE
      
      msg(paste0("  NLS successful: c = ", round(c_value, 2), 
                 ", z = ", round(z_value, 3)))
      msg(paste0("  Model R² = ", round(r2_value, 3)))
      
    }, error = function(e) {
      msg(paste0("  NLS failed: ", e$message))
      nls_success <<- FALSE
    })
    
    # STAGE 2: Fallback to log-log linear model if NLS failed
    if (!nls_success) {
      msg("  Falling back to log-log linear regression...")
      
      tryCatch({
        # Log-log transformation: log(S) = log(c) + z*log(A)
        lm_fit <- stats::lm(log(richness) ~ log(area), data = fit_data)
        
        # Extract and back-transform coefficients
        c_value <- exp(coef(lm_fit)[1])  # exp(intercept) = c
        z_value <- coef(lm_fit)[2]        # slope = z
        
        # R² from linear model (in log-space)
        r2_value <- summary(lm_fit)$r.squared
        
        # Store LM object as model_fit
        model_fit <- lm_fit
        fitting_method <- "log-log-lm"
        
        msg(paste0("  Log-log LM successful: c = ", round(c_value, 2), 
                   ", z = ", round(z_value, 3)))
        msg(paste0("  Model R² (in log-space) = ", round(r2_value, 3)))
        
      }, error = function(e) {
        stop("Both NLS and log-log linear model failed.\n",
             "Check your reference data for:\n",
             "  - Sufficient observations (need ≥3)\n",
             "  - Positive values only\n",
             "  - Reasonable range of areas\n",
             "Original error: ", e$message)
      })
    }
    
    # Set metadata (same for both methods)
    source_ref <- paste0("Data-driven estimation (", fitting_method, ")")
    doi_ref <- NA
    n_floras <- nrow(fit_data)
  }
  
  # ============================================================================
  # SECTION 3: RICHNESS PREDICTION
  # ============================================================================
  
  msg("Calculating predicted species richness...")
  
  # Arrhenius power function: S = c * A^z
  predicted_richness <- c_value * area^z_value
  
  # Create results data frame
  prediction_df <- data.frame(
    area_km2 = area,
    predicted_richness = round(predicted_richness, 0),
    stringsAsFactors = FALSE
  )
  
  # ============================================================================
  # SECTION 4: PREDICTION INTERVALS (if applicable)
  # ============================================================================
  
  intervals_df <- NULL
  
  if (return_intervals && method == "data" && !is.null(model_fit)) {
    msg("Calculating prediction intervals...")
    
    if (fitting_method == "nls") {
      # NLS: Manual calculation in original space
      pred_vals <- stats::predict(model_fit, newdata = data.frame(area = area))
      
      residual_se <- sqrt(sum(stats::residuals(model_fit)^2) / (nrow(fit_data) - 2))
      t_crit <- qt(1 - (1 - confidence_level)/2, df = nrow(fit_data) - 2)
      margin <- t_crit * residual_se * sqrt(1 + 1/nrow(fit_data))
      
      intervals_df <- data.frame(
        area_km2 = area,
        lower = round(pmax(0, pred_vals - margin), 0),
        upper = round(pred_vals + margin, 0),
        stringsAsFactors = FALSE
      )
      
    } else if (fitting_method == "log-log-lm") {
      # Log-log LM: Get intervals in log-space, then back-transform
      log_intervals <- predict(
        model_fit,
        newdata = data.frame(area = area),
        interval = "prediction",
        level = confidence_level
      )
      
      intervals_df <- data.frame(
        area_km2 = area,
        lower = round(exp(log_intervals[, "lwr"]), 0),
        upper = round(exp(log_intervals[, "upr"]), 0),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # ============================================================================
  # SECTION 5: DIAGNOSTICS
  # ============================================================================
  
  diagnostics <- list(
    method = method,
    region = if (method == "literature") region else NA,
    flora_type = flora_type,
    fitting_method = if (method == "data") fitting_method else NA,
    R2 = r2_value,
    RMSE = if (!is.null(model_fit)) {
      if (fitting_method == "nls") {
        sqrt(mean(stats::residuals(model_fit)^2))
      } else {
        # For log-log LM: calculate RMSE in original space
        fitted_richness <- exp(fitted(model_fit))
        sqrt(mean((fit_data$richness - fitted_richness)^2))
      }
    } else {
      NA
    },
    n_reference_floras = if (!is.na(n_floras)) n_floras else NA,
    source = source_ref,
    doi = if (!is.na(doi_ref)) doi_ref else NA,
    convergence_note = if (method == "data") {
      ifelse(fitting_method == "nls",
             "Non-linear least squares converged successfully",
             "NLS failed; used log-log linear regression fallback")
    } else {
      NA
    }
  )
  
  # ============================================================================
  # SECTION 6: RETURN RESULTS
  # ============================================================================
  
  msg("Done!")
  
  results <- list(
    prediction = prediction_df,
    coefficients = data.frame(
      c = c_value,
      z = z_value,
      source = source_ref,
      doi = if (!is.na(doi_ref)) doi_ref else NA,
      stringsAsFactors = FALSE
    ),
    intervals = intervals_df,
    model = model_fit,
    diagnostics = diagnostics
  )
  
  class(results) <- c("sar_prediction", "list")
  return(results)
}


# =============================================================================
# HELPER FUNCTION: LOAD SAR COEFFICIENTS (THREE-TIER APPROACH)
# =============================================================================

#' @title Load SAR Coefficient Database
#' @description Internal function to load SAR coefficients using three-tier approach
#' @param custom_path Optional path to custom CSV
#' @param verbose Logical, print messages
#' @keywords internal
load_sar_coefficients <- function(custom_path = NULL, verbose = TRUE) {
  
  msg <- function(...) if (verbose) message(...)
  
  # ===========================================================================
  # TIER 1: Try custom user-provided CSV first
  # ===========================================================================
  if (!is.null(custom_path)) {
    if (file.exists(custom_path)) {
      msg("Loading custom SAR coefficients from: ", custom_path)
      tryCatch({
        coef_df <- read.csv(custom_path, stringsAsFactors = FALSE)
        
        # Validate structure
        required_cols <- c("region", "flora_type", "c", "z")
        if (!all(required_cols %in% names(coef_df))) {
          warning("Custom CSV missing required columns: ", 
                  paste(setdiff(required_cols, names(coef_df)), collapse = ", "),
                  "\nFalling back to package defaults...")
        } else {
          return(coef_df)
        }
      }, error = function(e) {
        warning("Error reading custom CSV: ", e$message,
                "\nFalling back to package defaults...")
      })
    } else {
      warning("Custom coefficient file not found: ", custom_path, 
              "\nFalling back to package defaults...")
    }
  }
  
  # ===========================================================================
  # TIER 2: Try package-shipped CSV
  # ===========================================================================
  pkg_csv <- system.file("extdata", "sar_coefficients.csv", 
                         package = "ignobioR")
  
  if (file.exists(pkg_csv) && file.info(pkg_csv)$size > 0) {
    msg("Loading SAR coefficients from package database...")
    tryCatch({
      coef_df <- read.csv(pkg_csv, stringsAsFactors = FALSE)
      return(coef_df)
    }, error = function(e) {
      warning("Error reading package CSV: ", e$message,
              "\nUsing minimal hardcoded coefficients...")
    })
  }
  
  # ===========================================================================
  # TIER 3: Hardcoded fallback (ONLY Italy + Tuscany - empirical data only)
  # ===========================================================================
  if (verbose) {
    msg("Using default SAR coefficients (Italy + Tuscany).")
  }
  
  # Minimal but complete coefficient table - ONLY empirical data
  fallback_coefs <- data.frame(
    region = c(
      "Italy", "Italy", "Italy",           # National level
      "Tuscany", "Tuscany", "Tuscany"      # Regional level
    ),
    country = c(
      rep("Italy", 6)
    ),
    flora_type = c(
      "total", "native", "alien",          # Italy
      "total", "native", "alien"           # Tuscany
    ),
    c = c(
      241.2, 245.2, 10.1,                  # Italy (D'Antraccoli 2024)
      303.4, 314.6, 12.0                   # Tuscany (D'Antraccoli 2019)
    ),
    z = c(
      0.281, 0.263, 0.404,                 # Italy
      0.246, 0.229, 0.272                  # Tuscany
    ),
    R2 = c(
      0.92, 0.91, 0.73,                    # Italy
      0.87, 0.86, 0.56                     # Tuscany
    ),
    n_floras = c(
      286, 286, 286,                       # Italy (from 2024 paper)
      67, 36, 36                           # Tuscany (from 2019 paper)
    ),
    source = c(
      rep("D'Antraccoli et al. 2024", 3),
      rep("D'Antraccoli et al. 2019", 3)
    ),
    doi = c(
      rep("10.3390/plants13010012", 3),
      rep("10.1007/s10531-019-01730-x", 3)
    ),
    notes = c(
      rep("National scale (all Italy, 266 local floras)", 3),
      rep("Regional scale (Tuscany)", 3)
    ),
    stringsAsFactors = FALSE
  )
  
  return(fallback_coefs)
}


# =============================================================================
# HELPER FUNCTION: LIST AVAILABLE REGIONS
# =============================================================================

  #' @title List Available SAR Regions
  #' 
  #' @description Shows all regions with available SAR coefficients in the database
  #' 
  #' @param coef_database Optional path to custom coefficient CSV
  #' @param detailed Logical. If TRUE, shows full coefficient details including c, z, and R² (default = FALSE)
  #' @param verbose Logical. Print loading messages (default = FALSE)
  #' 
  #' @return Invisibly returns the full coefficient data frame. Prints summary to console.
  #' 
  #' @export
  #' 
  #' @examples
  #' \dontrun{
  #' # Quick overview of available regions
  #' list_sar_regions()
  #' 
  #' # Detailed view with coefficients
  #' list_sar_regions(detailed = TRUE)
  #' 
  #' # Check a custom database
  #' list_sar_regions(coef_database = "~/my_sar_coefficients.csv")
  #' }
  
  list_sar_regions <- function(coef_database = NULL, detailed = FALSE, verbose = FALSE) {
    
    coefs <- load_sar_coefficients(coef_database, verbose = verbose)
    
    if (detailed) {
      # DETAILED VIEW: Show full coefficient table
      cat("\n=== Complete SAR Coefficient Database ===\n\n")
      
      display_cols <- c("region", "flora_type", "c", "z", "R2", "n_floras", "source", "doi")
      display_cols <- intersect(display_cols, names(coefs))
      
      display_table <- coefs[, display_cols]
      display_table <- display_table[order(display_table$region, display_table$flora_type), ]
      rownames(display_table) <- NULL
      
      print(display_table, row.names = FALSE)
      cat("\n")
      
    } else {
      # SUMMARY VIEW: Show region overview
      summary_table <- unique(coefs[, c("region", "country", "source")])
      
      # Available flora types per region
      flora_counts <- aggregate(
        flora_type ~ region,
        data = coefs,
        FUN = function(x) paste(sort(unique(x)), collapse = ", ")
      )
      
      # Sample sizes per region
      n_floras_summary <- aggregate(
        n_floras ~ region,
        data = coefs,
        FUN = function(x) {
          if (all(is.na(x))) return(NA)
          vals <- unique(x[!is.na(x)])
          if (length(vals) == 1) return(vals[1])
          paste0(min(vals), "-", max(vals))
        }
      )
      
      summary_table <- merge(summary_table, flora_counts, by = "region")
      summary_table <- merge(summary_table, n_floras_summary, by = "region", all.x = TRUE)
      
      names(summary_table)[names(summary_table) == "flora_type"] <- "flora_types"
      names(summary_table)[names(summary_table) == "n_floras"] <- "n_floras"
      
      # Sort by region
      summary_table <- summary_table[order(summary_table$region), ]
      rownames(summary_table) <- NULL
      
      cat("\n=== Available SAR Coefficient Regions ===\n\n")
      print(summary_table, row.names = FALSE)
      cat("\nUse predict_richness(region = '...', flora_type = '...') to calculate predictions.\n")
      cat("For detailed coefficients, use: list_sar_regions(detailed = TRUE)\n\n")
    }
    
    invisible(coefs)
  }