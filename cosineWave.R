#' Fitting cosine-wave equation
#' @description Fitting a cosine-wave equation: y = baseline+(amplitude∙cos(2∙π∙((x−[phase shift]/24))) or 
#' a double harmonic cosine-wave equation: 
#' y = baseline+([amplitude A]∙cos(2∙π∙((x−[phase shift A])/24))) + ([amplitude B]∙cos(4∙π∙((x−[phase shift B])/24)))
#' with a fixed 24-h period
#' @param x x
#' @param y y
#' @param model either "single" (consine wave) or "double" (double harmonic cosine-wave)
#' @return a list object with two elements:
#'  - fit: the fitted model
#'  - summary: a summary of fitted parameters
#' @importFrom minpack.lm nlsLM
#' @examples 
#' # consine-wave equation example
#' x <- seq(0, 25, length.out = 100)
#' A_phase_shift <- 5
#' baseline <- 0
#' A <- 10
#' y <- baseline + A * cos(2 * pi * ((x - A_phase_shift)/24)) + rnorm(length(x), sd = 0.1)
#' plot(x, y)
#' cosineWave(x, y, model = "single")
#' 
#' ## double harmonic consine-wave equation example
#' x <- seq(0, 50, length.out = 100)
#' baseline <- 0
#' A <- 10
#' A_phase_shift <- 5
#' B <- 3
#' B_phase_shift <- 3
#' y <- baseline + A * cos(2 * pi * ((x - A_phase_shift)/24)) + B * cos(4 * pi * ((x - B_phase_shift)/24)) +rnorm(length(x), sd = 0.1)
#' y <- sample(y)
#' plot(x, y)
#' s <- cosineWave(x, y, model = "double")
#' s$summary

cosineWave <- function(x, y, model = c("single", "double")[1]) {
  
  i <- !is.na(x) & !is.na(y)
  x <- x[i]
  y <- y[i]
  
  par.init <-  list( A_phase_shift = 0, baseline = mean(x), A = max((max(y)-min(y))/2, 1e-2) )
  par.lower <- c( A_phase_shift = 0, baseline = min(x)/1.5, A = 1e-3)
  par.upper <- c( A_phase_shift = 24, baseline = max(x)*1.5, A = max(max(y)-min(y), 10) )
  if (model == "single") {
    form <- y ~ baseline + A * cos(2 * pi * ((x - A_phase_shift)/24))
  } else if (model == "double") {
    form <- y ~ baseline + A * cos(2 * pi * ((x - A_phase_shift)/24)) + B * cos(4 * pi * ((x - B_phase_shift)/24))
    par.init$B <- par.init$A
    par.init$B_phase_shift <- par.init$A_phase_shift
    par.lower["B"] <- par.lower["A"]
    par.lower["B_phase_shift"] <- par.lower["A_phase_shift"]
    par.upper["B"] <- par.upper["A"]
    par.upper["B_phase_shift"] <- par.upper["A_phase_shift"]
  }
  
  goodNLS <- function(fit) {
    if (!inherits(fit, c("nls", 'nlrob')))
      return(FALSE)
    if (inherits(fit, "nls"))
      r <- fit$convInfo$isConv else
        r <- fit$status == "converged"
    r
  }
  
  cin <-  paste(rep(names(par.init), each = 2), rep(c("2.5%", "97.5%"), time = 3), sep = "_")
  ci <- structure(rep(NA, length(cin)), names = cin)
  cf <- structure(rep(NA, length(par.init)), names = names(par.init))
  ret_invalid <- c(cf, pseudoRsq = NA, MSE = NA, n = NA, CI = ci)
  
  if (length(x) <= length(par.init))
    return(  
      list(fit = NA, summary = ret_invalid)
          )
  
  suppressWarnings({
    fit <- try(
      nls(form, start = par.init, algorithm="port", lower = par.lower, upper = par.upper, 
          control = list(warnOnly = TRUE, maxiter = 100)),
      silent = TRUE)
    alg <- "port"
    if (!goodNLS(fit))  {
      fit <- try(nls(form, start = par.init, algorithm="default"), silent = TRUE)
      alg <- "default"
    }
    if (!goodNLS(fit))  {
      fit <- try(nls(form, start = par.init, algorithm="plinear"), silent = TRUE)
      alg <- "plinear"
    }
    if (!goodNLS(fit))  { 
      fit <- minpack.lm::nlsLM(
        form, data = data.frame(x = x, y = y, pi = pi), start = par.init, lower = par.lower, upper = par.upper
      )
      alg <- "nlsLM"
    }
  })
  
  if (goodNLS(fit)) {
    ci <- suppressWarnings(suppressMessages(try(confint(fit), silent = TRUE)))
    if (inherits(ci, "try-error")) {
      cin <-  paste(rep(names(coef(fit)), each = 2), rep(c("2.5%", "97.5%"), time = 3), sep = "_")
      ci <- structure(rep(NA, length(cin)), names = cin)
    } else 
      ci <- structure(c(t(ci)), names = paste(rep(rownames(ci), each = 2), rep(colnames(ci), time = 3), sep = "_"))
    mse <- mean(residuals(fit)^2)
    rsq <- max(1 - var(residuals(fit))/var(y), 0)
    res <- c(coef(fit), pseudoRsq = rsq, MSE = mse, n = length(x), CI = ci)
    attr(res, "algrithm") <- alg
  } else {    
    res <- ret_invalid
  }
  list(
    fit = fit, summary = res
  )
}






