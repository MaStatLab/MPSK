# Mode = function(x) {
#   ux = unique(x)
#   ux[which.max(tabulate(match(x, ux)))]
# }

calibrate <- function(x)
{
  C = x$data$C - 1
  Z = x$chain$t - 1
  ns  = dim(x$chain$xi0)[3]
  K = dim(x$chain$xi0)[2]

  # S = matrix(0, nrow=ns, ncol=K)
  # for (i in 1:ns) {
  #   S[i,] = apply(x$chain$S[,,i], 2, Mode)
  # }
  output = calib(x$data$Y,
                 matrix(C,ncol=1),
                 Z,
                 x$chain$xi, dim(x$chain$xi),
                 x$chain$xi0, dim(x$chain$xi0),
                 x$chain$S )
  colnames(output$Y_cal) = colnames(x$data$Y)
  return(output)

}

relabelChain = function(res) {
  relabeled_chain = relabel(res)
  res$chain = relabeled_chain 
  res
}

# Recover stochastic representation parameters from distribution parameters:
transform_params = function(Omega, alpha) {
  n = nrow(Omega)
  m = ncol(Omega)
  if (n!=m) stop("Omega is not sqaure")
  if (length(alpha)!=n) stop("alpha is of wrong length")
  omega = sqrt(diag(diag(Omega)))
  omega_inv = solve(omega)
  OmegaBar = omega_inv %*% Omega %*% omega_inv
  alpha = matrix(alpha, ncol=1)
  numer = OmegaBar %*% alpha
  denom = as.numeric(sqrt(1 + t(alpha) %*% OmegaBar %*% alpha))
  delta = numer/denom
  psi = sqrt(diag(Omega)) * delta
  alpha.tilde = delta / sqrt(1-delta^2)
  inv.Delta = diag(1/c(sqrt(1-delta^2)))
  Omega.epsilon = inv.Delta %*% OmegaBar %*% inv.Delta - alpha.tilde%*%t(alpha.tilde)
  return(list(delta=delta, omega=omega, Omega.epsilon=Omega.epsilon,
              psi=psi, Sigma=(Omega - psi%*%t(psi))))
}