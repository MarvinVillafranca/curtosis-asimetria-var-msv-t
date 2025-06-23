
# LIBRERIAS  =====================================================================
# ================================================================================

library(MASS)

# ================================================================================
# ================================================================================





# GENERACION DE PARAMETROS =======================================================
# ================================================================================

generar_parametros_simulacion <- function(phi_diag,var_diag,var_out,p, k)
{
  # ===== Phi ===== 
  phi <- diag(phi_diag, p)
  
  # ===== v | A1 | ... | Ak =====
  VAR <- matrix(0, nrow = p, ncol = 1 + k*p)
  
  for (l in 1:k)
  {
    Al <- matrix(0, p, p)
    for (i in 1:p)
    {
      for (j in 1:p)
      {
        if (i == j)
        {
          Al[i, j] <- var_diag - 0.02 * l
        } 
        else if (l <= 2)
        {
          Al[i, j] <- var_out - 0.02 * l
        }
      }
    }
    VAR[, (l*p - (p - 1)):(l*p) + 1] <- Al
  }
  
  # ===== Sigma =====
  
  # See
  See <- matrix(0.864, p, p)
  diag(See) <- 1.44
  
  # Snn
  Snn <- matrix(0.028, p, p)
  diag(Snn) <- 0.04
  
  # Sen
  Sen <- matrix(-0.072, p, p)
  diag(Sen) <- -0.096
  
  # Ensamblar Sigma
  Sigma <- matrix(0, 2*p, 2*p)
  Sigma[1:p, 1:p] <- See
  Sigma[(p+1):(2*p), (p+1):(2*p)] <- Snn
  Sigma[1:p, (p+1):(2*p)] <- Sen
  Sigma[(p+1):(2*p), 1:p] <- t(Sen)
  
  return(list(phi = phi, VAR = VAR, Sigma = Sigma))
}

# ================================================================================
# ================================================================================





# GENERACION DE DATOS ============================================================
# ================================================================================

sigma0_VAR_MSV_t <- function(phi,Snn,p)
{
  sigma0 <- Snn
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      sigma0[i,j] <- (1/(1 - phi[i,i]*phi[j,j]))*Snn[i,j]
    } 
  }
  return(sigma0)
}

simulacion_VAR_MSV_t <- function(phi, sigma, VAR, nu, k, p, n)
{
  
  Snn <- sigma[(p+1):(2*p), (p+1):(2*p)]
  sigma0 <- sigma0_VAR_MSV_t(phi,Snn,p)
  
  # ===== Errores et y nt ===== 
  media <- rep(0, 2*p)
  #errores <- t(mvrnorm(n = n, mu = media, Sigma = sigma))
  errores <- matrix(NA,nrow = 2*p,ncol = n,byrow = T)
  media <- matrix(0L,nrow = 2*p,ncol = 1,byrow = T)
  for(t  in 1:n)
  {
    errores[1:(2*p),t] <- mvrnorm(1,media,sigma)
  }
  
  # ===== variables no observables alpha =====
  alpha <- matrix(NA, nrow = p, ncol = n)
  alpha[, 1] <- mvrnorm(1, rep(0, p), sigma0)
  for(i in 2:n)
  {
    alpha[, i] <- phi %*% alpha[, (i - 1)] + errores[(p+1):(2*p),(i - 1)]
  }
  
  # ===== ut = Vt^{1/2}et =====
  Un <- matrix(NA, nrow = p, ncol = n)
  for(i in 1:n)
  {
    Vt <- diag(exp(alpha[, i] / 2))  
    Un[, i] <- Vt %*% errores[1:p, i]
  }
  
  # ===== wt = lambda_t^{-1/2}ut =====
  #lambdas_N <- rgamma(N, shape = nu/2, rate = nu/2)  
  #lambdas_I <- 1 / sqrt(lambdas_N)                   
  #Wn <- Un * matrix(rep(lambdas_I, each = p), nrow = p)
  Wn <- Un
  lambdas_N <- matrix(NA, nrow = 1,ncol=N)
  lambdas_I <- matrix(NA, nrow = 1,ncol=N)
  for (i in 1:N)
  {
    lambdas_N[1,i] <- rgamma(1,nu/2,nu/2)
    lambdas_I[1,i] <- (1/sqrt(lambdas_N[1,i]))
    Wn[1:p,i] <- lambdas_I[i]*Un[1:p,i]  
  }
  
  # ===== yt = v + A1*y(t-1) + ... + + Aky*(t-k) + wt =====
  v <- VAR[, 1]
  A1k <- lapply(1:k, function(j) VAR[, ((j-1)*p + 2):((j*p) + 1)])
  Yn <- matrix(0L, nrow = p, ncol = N)
  # Inicialización de los primeros k valores
  Yn[, 1] <- v + Wn[, 1]
  if(1<k)
  {
    for (i in 2:k)
    {
      Yn[, i] <- v
      for (j in 1:(i-1))
      {
        Yn[, i] <- Yn[, i] + A1k[[j]] %*% Yn[, (i-j)]
      }
      Yn[, i] <- Yn[, i] + Wn[, i]
    }
  }
  
  # Simulación recursiva a partir de t = k+1
  for (i in (k + 1):n)
  {
    Yn[, i] <- v
    for (j in 1:k)
    {
      Yn[, i] <- Yn[, i] + A1k[[j]] %*% Yn[, i - j]
    }
    Yn[, i] <- Yn[, i] + Wn[, i]
  }
  
  return(list(sigma0 = sigma0, errores = errores, alpha = alpha, Un = Un, Wn = Wn, Yn = Yn))
}

# ================================================================================
# ================================================================================






# PRUEBA DE DATOS SIMULADOS ======================================================
# ================================================================================

# Generar parametros Sigma, Phi, v , A1,...,Ak.
phi_diag = 0.96
var_diag <- 0.13
var_out <- 0.09
p_nu <- 12
p <- 4
k <- 5
gps_temp <- generar_parametros_simulacion(phi_diag,var_diag,var_out,p,k)
p_phi <- gps_temp$phi
p_sigma <- gps_temp$Sigma
p_var <- gps_temp$VAR

# Datos Simulados
N <- 850
p_nu <- 12
sim_temp <- simulacion_VAR_MSV_t(p_phi, p_sigma, p_var, p_nu, k, p, N)
sim_temp$sigma0
sim_temp$errores  
sim_temp$Un       # Vt^{-1/2}et
sim_temp$Wn 
sim_temp$Yn       # observaciones

# ================================================================================
# ================================================================================


       