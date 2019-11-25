#' plot fittedChromatogram
#' @param x a fittedChromatogram object
#' @export 
plot.fittedChromatogram <- function(x) {
  
  pdcnorm <- function(x, sigma=1, mu=0) {
    exp(-((x-mu)^2)/(2*sigma^2)) / sqrt(2*pi*sigma^2)
  }
  
  m1 <- sprintf(
    "mz=%s, rt=%s, sn=%s, rsq=%s",
    round(x$peak["mz"], digits = 4),
    round(x$peak["rt"], digits = 4),
    x$peak["sn"], 
    signif(x$peak["rsq"], digits = 3))
  
  plot(MSnbase::rtime(x$chromatogram), MSnbase::intensity(x$chromatogram), 
       xlab = "Retention time", ylab = "Intensity", main = m1)
  xx <- seq(min(x$rtint$rt), max(x$rtint$rt), length.out = 100)
  if (!is.na(x$model)[1]) {
    est <- summary(x$model)
    est <- est$parameters[, "Estimate"]
    lines(xx, est["k"]*pdcnorm(xx, sigma = est["sigma"], mu = est["mu"]) + est["b"], type = "l")
  }
}

#' fit chromatogram with a gaussian curve using nonlinear least square
#' @param x XCMSnSet after peak detection
#' @param mc.cores passded to mclapply
#' @importFrom parallel mclapply
#' @importFrom minpack.lm nlsLM
#' @export
fitChromatogram <- function(x, mc.cores = 1) {
  
  ## define internal functions 
  nlsrsq <- function (mdl, y, param) {
    adj <- (sum(!is.na(y)) - 1)/(sum(!is.na(y)) - param)
    sum.sq <- (sum(!is.na(y)) - 1) * var(y, na.rm = TRUE)
    rsq <- 1 - (adj * (deviance(mdl)/sum.sq))
    return(rsq)
  }
  
  ## peaks
  peaks <- chromPeaks(x)
  # mclapply(1:3, function(i) {
  mclapply(1:nrow(peaks), function(i) {
    ii <<- i
    
    ch <- xcms::chromatogram(
      x, 
      rt=c(peaks[i, "rtmin"], peaks[i, "rtmax"]), 
      mz=c(peaks[i, "mzmin"], peaks[i, "mzmax"]), 
      aggregationFun = "sum")
    
    ch <- ch[[1]]
    tab <- data.frame(
      rt = MSnbase::rtime(ch), 
      int = MSnbase::intensity(ch))
    tab <- na.omit(tab)
    
    if (nrow(tab) < 5) 
      res <- NA else {
        res <- try(
          nlsLM(
            int ~ k*exp(-((rt-mu)^2)/(2*sigma^2)) / sqrt(2*pi*sigma^2)+b, 
            start=c(mu=mean(tab$rt),sigma=sd(tab$rt), k = max(tab$int), b = 0) , 
            data = tab, control = nls.lm.control(maxiter = 500)
          ), silent = TRUE)
        if (inherits(res, "try-error")) 
          res <- NA
      }
    
    rsq <- rtgap <- intgap <- rtintgap  <- truncated <- b <- NA
    
    if (!is.na(res)[1]) {
      rsq <- nlsrsq(res, tab$int, param  = 4)
      
      d <- dist(tab$rt)
      rtgap <- max(d)/median(d)
      
      d <- dist(tab$int)
      intgap <- max(d)/median(d)
      
      tabn <- apply(tab, 2, function(x) {
        x <- x-min(x)
        x/max(x)
      })
      d <- dist(tabn)
      rtintgap <- max(d)/median(d)
      
      estmu <- summary(res)$parameters["mu", "Estimate"]
      if (estmu < min(tab$rt) || estmu > max(tab$rt)) 
        truncated <- Inf else
          truncated <- log10((max(tab$rt) - estmu)/(estmu - min(tab$rt)))
      
      b <- summary(res)$parameters["b", "Estimate"]
    }
    
    addinfo <- peaks[i, ]
    addinfo["rsq"] <- rsq
    addinfo["rtgap"] <- rtgap
    addinfo["intgap"] <- intgap
    addinfo["rtintgap"] <- rtintgap
    addinfo["truncated"] <- truncated
    addinfo["b"] <- b
    
    l <- list(chromatogram = ch,
              rtint = tab, 
              model = res, 
              peak = addinfo)
    class(l) <- "fittedChromatogram"
    l
  }, mc.cores = mc.cores)
  
}


#' Comparison of MS2 spectra with reference mass table
#' @param x an object of class MSnExp
#' @param refmass the reference mass table, it has at least two columns named "mass" and "cpd"
#' @param ppm.tol the mass tolerence in ppm comparing spectra
#' @param mc.cores passed to mclapply
#' @import parallel mclapply
#' @export
#' 
compareMSnExp <- function(x, refmass, ppm.tol=10, mc.cores=1) {
  
  msl <- MSnbase::msLevel(x)
  msl <- which(msl == 2)
  x <- mclapply(msl, function(i, raw) {
    x <- raw[[i]]
    if (MSnbase::msLevel(x) != 2)
      stop("Only MS2 spectrum allowed!")
    
    if (length(x@mz) == 0)
      return()
    data.frame(
      precScanNum = x@precScanNum,
      precMz = x@precursorMz,
      precInt = x@precursorIntensity,
      mz = x@mz,
      int = x@intensity, 
      scanNum = x@scanIndex,
      raw = x@fromFile,
      stringsAsFactors = FALSE)
  }, raw = x, mc.cores = mc.cores)
  sx <- do.call(rbind, x)
  
  mappedSpec <- mclapply(1:nrow(refmass), function(j) {
    ppm <- 1e6*(abs(refmass$mass[j] - sx$mz))/refmass$mass[j]
    i <-  ppm < ppm.tol
    dfx <- sx[i, ]
    dfx$queryMass <- refmass$mass[j]
    dfx$queryCpd <- refmass$cpd[j]
    dfx$dppm <- ppm[i] 
    dfx
  }, mc.cores = mc.cores)
  
  do.call(rbind, mappedSpec)
}

#' Comparison of a MS2 spectrum with reference mass table
#' @param x an object of class Spectrum2
#' @param refmass the reference mass table, it has at least two columns named "mass" and "cpd"
#' @param ppm.tol the mass tolerence in ppm comparing spectra
#' 
compareMS2 <- function(x, refmass, ppm.tol=10) {
  
  if (MSnbase::msLevel(x) != 2)
    stop("Only MS2 spectrum allowed!")
  
  if (length(x@mz) == 0)
    return()
  
  df <- data.frame(
    precScanNum = x@precScanNum,
    mz = x@mz,
    int = x@intensity, 
    scanNum = x@scanIndex,
    raw = x@fromFile,
    stringsAsFactors = FALSE)
  
  q <- dist(c(refmass$mass, df$mz))
  q <- as.matrix(q)
  q <- q[1:nrow(refmass), -(1:nrow(refmass))]
  q_ppm <- q/refmass$mass * 1e6
  q_idx <- which(q_ppm < ppm.tol, arr.ind = TRUE)
  
  if (length(q_idx) == 0)
    return()
  
  cbind(refmass[q_idx[, "row"], ], df[q_idx[, "col"], ])
}

