#' Linear Regression using QR Decomposition
#'
#' Fits a linear regression model using QR decomposition.
#'
#' @param formula A formula specifying the model (e.g., y ~ x).
#' @param data A data frame containing the variables in the formula.
#' @return An object of class \code{linreg} containing model results.
#' @importFrom stats model.frame model.matrix pt
#' @export
linreg <- function(formula, data) {
  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # Create model frame and matrices
  model_frame <- stats::model.frame(formula, data)
  y <- model_frame[[1]]  # Response variable
  X <- stats::model.matrix(formula, model_frame)  # Design matrix

  # QR decomposition approach
  qr_decomp <- qr(X)
  coefficients <- qr.coef(qr_decomp, y)
  fitted <- as.vector(X %*% coefficients)
  residuals <- y - fitted

  # Degrees of freedom and variance estimates
  n <- nrow(X)
  p <- ncol(X)
  df <- n - p

  # Residual variance using QR (more stable)
  residual_var <- sum(residuals^2) / df

  # Variance-covariance matrix using QR decomposition
  R_inv <- solve(qr.R(qr_decomp))
  coef_var <- residual_var * (R_inv %*% t(R_inv))

  # Standard errors and hypothesis tests
  std_errors <- sqrt(diag(coef_var))
  t_values <- coefficients / std_errors
  p_values <- 2 * pt(-abs(t_values), df)

  result <- list(
    coefficients = coefficients,
    fitted = fitted,
    residuals = residuals,
    formula = formula,
    data = data,
    df = df,
    residual_var = residual_var,
    coef_var = coef_var,
    t_values = t_values,
    p_values = p_values,
    std_errors = std_errors,
    call = match.call()
  )
  class(result) <- "linreg"
  return(result)
}


#' Print Method for linreg Objects
#'
#' Displays the formula and estimated coefficients.
#'
#' @param x An object of class \code{linreg}.
#' @param ... Additional arguments (currently unused).
#' @export
print.linreg <- function(x, ...) {
  cat("linreg(formula = ")
  cat(deparse(x$formula))
  cat(", data = iris)\n\n")
  cat("Coefficients:\n")
  print(x$coefficients)
}



#' Plot Method for linreg Objects
#'
#' Creates diagnostic plots: residuals vs fitted and scale-location.
#'
#' @param x An object of class \code{linreg}.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns a list of ggplot objects.
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth labs
#' @export
plot.linreg <- function(x, ...) {
  plot_data <- data.frame(
    fitted = x$fitted,
    residuals = x$residuals,
    std_resid = sqrt(abs(x$residuals))
  )

  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fitted, y = residuals)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = "Residuals vs Fitted",
                  x = "Fitted values",
                  y = "Residuals")

  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = fitted, y = std_resid)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(se = FALSE, method = "loess", span = 0.5) +
    ggplot2::labs(title = "Scale-Location",
                  x = "Fitted values",
                  y = "sqrt(|Standardized residuals|)")

  print(p1)
  print(p2)

  invisible(list(residuals_vs_fitted = p1, scale_location = p2))
}

# Declare global variables outside any function (at the top-level)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("fitted", "residuals", "std_resid"))
}




#' Residuals Method for linreg Objects
#'
#' Returns the residuals from the fitted model.
#'
#' @param object An object of class \code{linreg}.
#' @param ... Additional arguments (currently unused).
#' @return A numeric vector of residuals.
#' @export
residuals.linreg <- function(object, ...) {
  return(object$residuals)
}


#' Prediction Generic Function
#'
#' A generic function for obtaining predictions from model objects.
#'
#' @param object A model object (can be any object for which a method has been defined).
#' @param ... Additional arguments passed to methods.
#' @export
pred <- function(object, ...) {
  UseMethod("pred")
}

#' Prediction Method for linreg Objects
#'
#' Returns the fitted values from a linear regression model.
#'
#' @param object An object of class \code{linreg}.
#' @param ... Additional arguments (currently unused).
#' @return A numeric vector of fitted values.
#' @export
pred.linreg <- function(object, ...) {
  return(object$fitted)
}




#' Coefficients Method for linreg Objects
#'
#' Returns the estimated regression coefficients.
#'
#' @param object An object of class \code{linreg}.
#' @param ... Additional arguments (currently unused).
#' @return A named numeric vector of coefficients.
#' @export
coef.linreg <- function(object, ...) {
  return(object$coefficients)
}


#' Summary Method for linreg Objects
#'
#' Displays a table of coefficients, standard errors, t-values, and p-values.
#'
#' @param object An object of class \code{linreg}.
#' @param ... Additional arguments (currently unused).
#' @importFrom stats symnum
#' @export
summary.linreg <- function(object, ...) {
  cat("Call:\n")
  cat(deparse(object$formula), "\n\n")  # Print formula without environment info

  # Compute significance stars based on p-values
  stars <- symnum(object$p_values, corr = FALSE, na = FALSE,
                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))

  # Prepare summary table
  summary_table <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = sqrt(diag(object$coef_var)),
    `t value` = object$t_values,
    `Pr(>|t|)` = object$p_values,
    Signif. = stars
  )

  print(summary_table, digits = 5)

  cat("\nResidual standard error:",
      format(sqrt(object$residual_var), digits=5),
      "on", object$df, "degrees of freedom\n")
}





