
# MetaXcan function
metaXcanRef = function(Gx, Gl, Z){
    G_mean = matrix(rep(apply(Gx, 2, mean),each=nrow(Gx)), 
                  nrow = nrow(Gx), ncol=ncol(Gx))
    Gamma_G = t(Gx-G_mean)%*%(Gx-G_mean)/nrow(Gx)
    sigma_g = sqrt(t(Z)%*%Gamma_G%*%Z)
    varG.gwas = apply(Gl, 2, var)
    sigma_l = sqrt(varG.gwas)
    return(list(Gamma_G=Gamma_G, sigma_g=sigma_g, sigma_l=sigma_l))
}
metaXcan = function(betas, se.betaG, Gamma_G, Z, sigma_g, sigma_l){
    if(dim(sigma_g)[1] != 1){
      cat('dim(sigma_g)[1] != 1 --> More than one X in the model \n')
      break
    }else{
      beta_metaXcan = sum(Z*betas*(sigma_l^2))/(sigma_g^2) #beta to return
      vec.sigma_G = rep(sigma_g, times=length(sigma_l))
      z.score = sum(Z*(sigma_l/vec.sigma_G)*(betas/se.betaG))
      p = 2*pnorm(-abs(z.score))
      se_metaXcan = beta_metaXcan/z.score
      return(list(betas=beta_metaXcan, ses=se_metaXcan, p=p))
    }
}

# Simulation function build up -----------------------------------------------------
createDistAR <- function(n, order=NULL, diag.v =1) {
    order <- ifelse(is.null(order), n, order)
    order <- ifelse(order>=n,n-1,order)
    A <- diag(n)*diag.v
    for(d in 1:order) {
        v <- 1/(d*d)
        for(t in 1:(n-d)) {
            A[t+d, t] <- v
            A[t,t+d] <- v
        }
    }    
    A
}
createSIGMA <- function(n, order=NULL, diag.v =1) {
    order <- ifelse(is.null(order), n, order)
    order <- ifelse(order>=n,n-1,order)
    A <- diag(n)*diag.v
    for(d in 1:order) {
        v <- runif(1,0.3,0.5)
        for(t in 1:(n-d)) {
            A[t+d, t] <- v
            A[t,t+d] <- v
        }
    }
}

# G: correlation structure --> set the cut-off point for G = 0,1,2
createGeno = function(N, numSNP, MAF, SIGMA){
    if(numSNP != length(MAF)){
        cat('numSNP != length(MAF)\n')
        break
    }else{
        c0 = qnorm((1-MAF)^2) # (1-p)^2 = q^2 --> G=0
        c2 = qnorm(1-MAF^2) # p^2 --> G=2
        
        G = rmvnorm(N, sigma=SIGMA) # default --> mean = 0
        G = lapply(1:numSNP, FUN=function(m) {ifelse(G[,m] > c2[m], 2, ifelse(G[,m]<c0[m], 0,1))})
        G = do.call(cbind, G)
        return(G)
    }
}

# simulate beta based on r^2 in y~x
simBetabyR2 = function(r2, Nrow, Ncol, X, varY, sd){
    beta = matrix(NA, ncol=Ncol, nrow = Nrow)
    if(Ncol>1){
        for(i in 1:Ncol){
            beta[,i] = rnorm(Nrow, mean=sqrt(r2*varY/(sum(var(X)))), sd)
        }
    }else{
        beta = rnorm(Nrow, mean=sqrt(r2*varY/(sum(var(X)))), sd)
    }
    return(abs(beta))
}
simCorrEbyGZ = function(G, alpha, numGenes, sigmaE, N){
    E = mat.or.vec(nrow(G), numGenes)
    meanE = G%*%alpha
    for(i in 1:nrow(G)){
        E[i,] = rmvnorm(1, mean=t(meanE[i,]), sigmaE)
    }
    return(E)
}

# legend function for grid arrange
g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# Simulation scenarios build up -----------------------------------------------------
# Get variable list for simulation
getSIGMA = function(numSNP, numBlock, order, presetcorrG){
    SIGMA = matrix(0, nrow = numSNP*numBlock, ncol = numSNP*numBlock)
    for(blocknumber in 1:numBlock){
      blockstart = (blocknumber-1)*numSNP+1
      blockstop = blocknumber*numSNP
      SIGMA[blockstart:blockstop, blockstart:blockstop] = createDistAR(numSNP, order)*presetcorrG
    }
    diag(SIGMA) = 1 
    return(SIGMA)
}

# elnetsum
elnetsum <- function(betas.Gx, N.Gx, Geno, eaf.Gx = NULL, c = 0.5){
  
  X = 0; Y = 0
  for (i in 1:length(betas.Gx)) {
    X_i = get_zL_X(betas.Gx[[i]], N.Gx[[i]], Geno, eaf.Gx = NULL)$L
    X = X + X_i
    Y_i = get_zL_X(betas.Gx[[i]], N.Gx[[i]], Geno, eaf.Gx = NULL)$zL
    Y = Y + Y_i
  }
 
  fit = cv.glmnet(X, Y, alpha = c, family = "gaussian") 
  
  return(fit)
}

