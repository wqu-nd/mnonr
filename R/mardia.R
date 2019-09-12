#' Univariate and Multivariate skewness and kurtosis checker
#'
#' @param x A data matrix
#' @param na.rm An indication of the missing data, the default value is True
#'
#' @return Data information: sample size and number of variables. The marginal and multivariate test (Mardia's Test) of skewness and kurtosis.
#' @export mardia
#' @importFrom stats na.omit cov pchisq pnorm
#'
#'
mardia <- function(x, na.rm = TRUE){
  if (na.rm)
    x = na.omit(x)
  n = dim(x)[1]
  p = dim(x)[2]

  uni = function(x){
    n = length(x)
    xbar = mean(x)
    m2 = sum((x-xbar)^2)/n
    m3 = sum((x-xbar)^3)/n
    m4 = sum((x-xbar)^4)/n

    skewness = sqrt(n*(n-1))/(n-2)*m3/m2^1.5
    kurtosis = (n-1)/((n-2)*(n-3))*((n+1)*(m4/m2^2-3)+6)

    skew.se = sqrt(6*n*(n-1) / ((n-2)*(n+1)*(n+3)))
    kurt.se = sqrt(4*(n^2-1)*skew.se^2 / ((n-3)*(n+5)))
    c(skewness, skew.se=skew.se, kurtosis, kurt.se)
  }

  univariate = apply(x, 2, uni)
  rownames(univariate) = c('Skewness', 'SE_skew', 'Kurtosis', 'SE_kurt')

  x = scale(x, scale = FALSE)
  S = cov(x)*(n-1)/n
  S.inv = MASS::ginv(S)
  D = x %*% S.inv %*% t(x)
  b1p = sum(D^3)/n^2
  b2p = sum(diag(D^2))/n
  chi.df = p * (p + 1) * (p + 2)/6
  k = (p + 1) * (n + 1) * (n + 3)/(n * ((n + 1) * (p + 1) - 6))
  small.skew = n * k * b1p/6
  M.skew = n * b1p/6
  M.kurt = (b2p - p * (p + 2)) * sqrt(n/(8 * p * (p + 2)))
  p.skew = 1 - pchisq(M.skew, chi.df)
  p.small = 1 - pchisq(small.skew, chi.df)
  p.kurt = 2 * (1 - pnorm(abs(M.kurt)))

  multivariate = rbind(c(b1p, M.skew, p.skew), c(b2p, M.kurt, p.kurt))
  rownames(multivariate) = c('Skewness', 'Kurtosis')
  colnames(multivariate) = c('b', 'z', 'p-value')

  results = list(n.obs = n, n.var = p, univariate = univariate, multivariate=multivariate)

  cat('Sample size: ', n, "\n")
  cat('Number of variables: ', p, "\n\n")

  cat("Marginal skewness and kurtosis\n")
  print(t(univariate))

  cat("\nMardia's multivariate skewness and kurtosis\n")
  print(multivariate)
  cat("\n")

  invisible(results)
}


