#' Regularized regression incorporating external information for summary statistics analysis
#'
#' @param betas.Gx A list of numeric vectors that contains the summary statistics from each ethnic population
#' @param N.Gx A list of numeric that contains sample sizes in all original GWAS studies
#' @param Geno A list of matrices with individual-level SNP dosage data in each study/population as reference data. Each column corresponds to a SNP. Note that the order of the study/population should be matched in all lists.
#' @param eaf.Gy The effect allele frequency of the SNPs in betas.Gx
#' @param Z Prior information matrix about the predictors/SNPs (\eqn{p} rows, each corresponding to a SNP in X; \eqn{q} columns of external information about the predictors, such as prior biological importance). If Z is the grouping of predictors, it is best if user codes it as a dummy variable (i.e. each column indicating whether predictors belong to a specific group).
#' @param c The elastic-net mixing parameter ranging from 0 to 1. When  \eqn{c} = 1, the model corresponds to Lasso. When \eqn{c} is set to 0, it corresponds to Ridge. For values between 0 and 1 (with a default of 0.5), the model corresponds to the elastic net.
#' @param sigma.square A user-supplied noise variance estimate. Typically, this is left unspecified, and the function automatically computes an estimated sigma square values using R package \code{selectiveinference}.
#' @param alpha.est.init
#' @param maxstep
#' @param tolerance
#' @param maxstep_inner
#' @param tolerance_inner
#' @param compute.likelihood
#' @param verbosity
#' @param standardize
#' @param intercept
#' @param epsilon
#' @import glmnet crayon selectiveInference lbfgs
#'
#' @return
#' @export
#'
#' @examples
xtune_mJAM <- function(betas.Gx, N.Gx, Geno, eaf.Gy = NULL, Z = as.matrix(rep(1, ncol(Geno[[1]]))),
                     c = 0.5,
                     sigma.square = NULL,
                     alpha.est.init = rep(0,ncol(Z)+1),
                     maxstep = 100,
                     tolerance = 0.001,
                     maxstep_inner = 50,
                     tolerance_inner = 0.1,
                     compute.likelihood = FALSE,
                     verbosity = FALSE,
                     standardize = TRUE,
                     intercept = TRUE,
                     epsilon = 3){

  X = 0; Y = 0; sigma2 = 0
  for (i in 1:length(betas.Gx)) {
  X_i = get_zL_X(betas.Gx[[i]], N.Gx[[i]], Geno[[i]], eaf.Gx = NULL)$L
  X_i = scale(X_i, center = TRUE, scale = F)
  X = X + X_i
  Y_i = get_zL_X(betas.Gx[[i]], N.Gx[[i]], Geno[[i]], eaf.Gx = NULL)$zL
  Y = Y + Y_i
  sigma2_i = ifelse(is.null(sigma.square), estimateVariance(X_i,Y_i), sigma.square)
  sigma2 = sigma2 + sigma2_i
  }

  n = nrow(X)
  p = ncol(X)

  ##------------ Elastic-net regression
  ## Initialize
  alpha.old = alpha.est.init
  likelihood.score = c()
  k = 1

  ## calculate alpha.max
  alpha.max = max(abs(colSums(X%*%Y)))/ (c * n)

  ## reparameterize Z
  Z = sweep(Z,2,colMeans(Z))
  Z = cbind(rep(1,p), Z)

  while(k < maxstep){
    # Given alpha, update theta
    lambda = exp(Z%*%alpha.old)
    gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
    Sigma_y = sigma2 * diag(n) + (t(t(X) * c(gamma))) %*% t(X)
    theta = colSums(X * solve(Sigma_y, X))

    # Compute likelihood
    if (compute.likelihood == TRUE) {
      likelihood.score = c(likelihood.score, approx_likelihood.EN(alpha.old,c,
                                                                  X, Y, Z, sigma2))
    }

    # Given theta, update alpha
    update.result <- update_alpha.EN(X, Y, Z,c=c, alpha.old = alpha.old, alpha.max = alpha.max, epsilon = epsilon,
                                     sigma.square = sigma2, theta = theta, maxstep_inner = maxstep_inner,
                                     tolerance_inner = tolerance_inner)
    alpha.new <- update.result$alpha.est

    # Check convergence
    if (sum(abs(alpha.new - alpha.old)) < tolerance) {
      break
    }
    alpha.old <- alpha.new

    # Track iteration progress
    if (verbosity == TRUE) {
      cat("#-----------------Iteration ", k, " Done -----------------#\n",
          sep = "")
    }
    k <- k + 1
  }

  tauEst = exp(Z%*%alpha.old)
  pen_vec = tauEst * sigma2/n

  pen_vec[pen_vec>1e6] <- 1e6
  C = sum(pen_vec)/p

  obj <- glmnet(X, Y, alpha = c,family="gaussian",standardize = standardize, intercept = intercept)
  cus.coef <- coef(obj,x=X,y=Y,alpha = c, exact=TRUE,s= C, penalty.factor = pen_vec,standardize=standardize,intercept = intercept)

  return(list(beta.est = cus.coef, penalty.vector = pen_vec, lambda = C, alpha.est = alpha.old,
              n_iter = k - 1, sigma.square = sigma2, likelihood.score = likelihood.score))

}

estimateVariance <- function(X, Y, n_rep = 5) {
  Y <- as.double(drop(Y))
  dimY = dim(Y)
  nrowY = ifelse(is.null(dimY), length(Y), dimY[1])
  if (nrowY < 10) {
    stop("Need at least 10 observations to estimate variance")
  }

  temp = array(NA, n_rep)
  for (i in 1:n_rep) {
    c = suppressWarnings(estimateSigma(X, Y)$sigmahat^2)
    temp[i] = ifelse(is.infinite(c), NA, c)
  }
  return(mean(temp, na.rm = T))
}

get_zL_X <- function(betas.Gx, N.Gx, Geno, eaf.Gx = NULL, ridgeTerm = T) {
  # Check the dimension of betas.Gy, Geno and A
  dim_betas = length(betas.Gx)
  dim_Geno = ncol(Geno)

  if(dim_betas == dim_Geno){

    # Remove rows with zero in Genotype file
    if(sum(is.na(Geno))>0){
      Geno = Geno[complete.cases(Geno), ]
    }


    if(!is.null(eaf.Gx)){
      p_D = eaf.Gx
    }else{
      p_D = apply(Geno, 2, mean)/2
    }

    # Obtain the JAM variables: zL and L
    n0 = N.Gx*(1-p_D)^2
    n1 = N.Gx*2*p_D*(1-p_D)
    n2 = N.Gx*p_D^2

    y0 = -(n1*betas.Gx+2*n2*betas.Gx)/(n0+n1+n2)
    y1 = y0+betas.Gx
    y2 = y0+2*betas.Gx
    z = n1*y1 + 2*n2*y2

    ## Compute G0'G0
    G0 = scale(Geno, center=T, scale=F)
    G0_t_G0 = t(G0)%*%G0

    ## Modify G0'G0 if the sample sizes of Geno and Gy are different
    Dj = 2*p_D*(1-p_D)*N.Gx
    D_sqrt = diag(sqrt(Dj))
    Dw_sqrt_inv = diag(1/sqrt(diag(G0_t_G0)))
    G0_t_G0.scaled = D_sqrt %*% Dw_sqrt_inv  %*% G0_t_G0 %*% Dw_sqrt_inv %*% D_sqrt

    ## Add a ridge term in case G0'G0 is singular
    ridgeValue = ifelse(ridgeTerm, min(1, min(diag(G0_t_G0.scaled)*.001)), 0)
    G0_t_G0.ridge = G0_t_G0.scaled + ridgeValue*diag(length(betas.Gx))

    # Perfrom Cholesky decompostion and construct zL
    L = chol(G0_t_G0.ridge)
    zL = solve(t(L))%*%z

    return(list(zL = zL, L=L))
  }
}

approx_likelihood.EN <- function(to_estimate,c, X, Y, Z, sigma.square.est) {
  n = nrow(X)
  lambda = exp(Z %*% to_estimate)
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
  K = sigma.square.est * diag(n) + X %*% diag(c(gamma)) %*% t(X)
  logdetK = determinant(K)$modulus[1]
  part1 = t(Y) %*% solve(K, Y)
  normapprox = 1/2 * (part1 + logdetK)
  return(-as.numeric(normapprox))
}


update_alpha.EN <- function(X, Y, Z,c,alpha.old, alpha.max, epsilon, sigma.square, theta, maxstep_inner,
                            tolerance_inner) {
  ## initial
  alpha.inner.old = alpha.old
  k_inner = 1
  n = nrow(X)
  p = ncol(X)
  while (k_inner < maxstep_inner) {
    # given alpha update delta
    lambda = exp(Z %*% alpha.inner.old)
    gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))

    sd_y <- sqrt(var(Y) * (n - 1)/n)
    C = sum(1/gamma)/p * sd_y * sigma.square/n
    delta.est = coef(glmnet(X, Y, alpha = 0, penalty.factor = 1/gamma, lambda = C,
                            standardize = F, intercept = FALSE))[-1]

    ## given delta update alpha
    alpha.inner.new <- optim(alpha.old, likelihood.alpha.theta.EN,likelihood.alpha.theta.gradient.EN,
                             c =c,Z = Z, theta = theta, delta = delta.est,method = "BFGS", upper = c(alpha.max*epsilon, rep(Inf, length(alpha.old)-1)))$par
    if (sum(abs(alpha.inner.new - alpha.inner.old)) < tolerance_inner) {
      break
    }
    k_inner = k_inner + 1
    alpha.inner.old <- alpha.inner.new
  }
  return(list(alpha.est = alpha.inner.old, inner_iter = k_inner))
}

likelihood.alpha.theta.EN <- function(Z,c, alpha, theta, delta) {
  lambda = exp(Z %*% alpha)
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))
  L = as.numeric(t(theta) %*% gamma + delta^2 %*% (1/gamma))

  return(ifelse(is.infinite(L), 1e6, L))
}


likelihood.alpha.theta.gradient.EN <- function(Z,c, alpha, theta, delta) {
  lambda = exp(Z %*% alpha)
  gamma = 2/(lambda^2*c^2 + 2*lambda*(1-c))

  dev_gamma = theta - delta^2/(gamma^2)
  dev_gamma_alpha = as.vector((-2*(2*(1-c) + 2*c^2*lambda))/(2*lambda*(1-c) + (c*lambda)^2)^2 *lambda) *Z
  return(crossprod(dev_gamma, dev_gamma_alpha))
}



relative_change<-function(a,b){
  return((a-b)/b)
}
