
#' Multivariate Non-normal Random Number Generator based on Multivariate Measures
#'
#' @param n Sample size
#' @param p Number of variables
#' @param ms A value of multivariate skewness
#' @param mk A value of multivariate kurtosis
#' @param Sigma A covariance matrix (In this function, the generated data are standarized. A correlation matrix is equal to its corresponding covariance matrix.)
#' @param initial A vector with 3 numbers for initial polynominal coefficients' (b,c,d). The default setting is (0.9,0.4,0).
#' @return A data matrix (multivariate data)
#' @importFrom stats optim rnorm
#' @export mnonr
#'
#' @examples mnonr(n=10000,p=2,ms=3,mk=61,Sigma=matrix(c(1,0.5,0.5,1),2,2),initial=NULL)
mnonr <- function(n, p, ms, mk, Sigma, initial=NULL){
  ##
  #Generates n random variables with specified multivariate skewness and kurtosis
  # using revised Fleishman's power method. Note that not all combinations
  # of skew and kurtosis are possible
  # p independent variable
  ##
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  if (!isSymmetric(Sigma))
    stop(cat("Covariance matrix is must be symmetric."))

  eS = eigen(Sigma, symmetric = TRUE)
  ev = eS$values

  if (!all(ev >= -1e-06 * abs(ev[1L])))
    stop(cat("Covariance matrix is not positive definite."))

  if (is.null(initial))
    initial = c(0.9, 0.4, 0)

  #Multivariate skewness and kurtosis range
  sug_mk = 1.641 * ms + p * (p + 0.774)
  sug_ms = (mk - p * (p + 0.774)) / 1.641
  sug_min_mk = p * (p + 0.774)
  if (ms < 0)
    stop(cat('The multivariate skewness must be non-negtive.'))
  if(mk < p * (p + 0.774))
    stop(stop(cat('The minimun multivariate kurtosis in your setting is',sub_min_mk)))
  else if(!(mk >= 1.641 * ms + p * (p + 0.774)))
    stop(cat("The multivariate skewness and kurtosis must follow the range of:  MK>= 1.641*MS + p*(p + 0.774) and MS cannot be negative, where p is the number of variables.\n For your reference:\n For the given p and multivariate skewness, the kurtosis must be no less than",sug_mk,".\n For the given p and multivariate kurtosis, the skewness must be no more than",sug_ms,'.'))

  multouni = function(p, ms, mk){
    beta1 = ms
    beta2 = mk

    #Fourth moment of ksi
    f = beta2 / p - (p - 1)

    #third moment of ksi
    t = sqrt(beta1 / p)

    out = c(t, f)
    return(out)
  }

  fleishtarget = function(x, a){
    ##Revised Fleishman method
    b = x[1]
    cc = x[2]
    d = x[3]
    g1 = a[1]
    g2 = a[2]
    (1 - (b^2+2*cc^2+6*b*d+15*d^2))^2 +
      (g1 - (6*b^2*cc+8*cc^3+72*b*cc*d+270*cc*d^2))^2 +
      (g2 - (3*b^4+60*b^2*cc^2+60*cc^4+60*b^3*d+936*b*cc^2*d+630*b^2*d^2+4500*cc^2*d^2+3780*b*d^3+10395*d^4))^2
  }

  findcoe = function(tf,initial){
    ##
    #Uses the built in minimization function to solve for b, c, and d#
    # if the third and fourth moment of ksi are given.
    ##
    output = optim(initial,fleishtarget,a = tf,method = "BFGS",
                 control = list(ndeps = rep(1e-10, 3),reltol = 1e-10,maxit = 1e8))$par

    return(output)
  }

  Z = matrix(rnorm(n*p, mean=0, sd=1), nrow=n, ncol=p)

  bcd = findcoe(multouni(p, ms, mk),initial)
  b = bcd[1]
  cc = bcd[2]
  d = bcd[3]
  a = -1*cc
  xi = a+b*Z+cc*Z^2+d*Z^3

  x = matrix(0, nrow = n, ncol = p)
  r = chol(Sigma)
  for (j in 1:n){
    for(m in 1:p){
      for(i in 1:p){
        x[j,m] = x[j,m] + r[i,m] * xi[j,i]}}}

  return(x)
}


