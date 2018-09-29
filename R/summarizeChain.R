Mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

summarizeChain = function( res ) {

  chainSummary = list()
  K = res$prior$K
  chain = res$chain
  p = ncol(res$data$Y)
  J = length(unique(res$data$C))
  ns = res$pmc$nsave

  xi_raw = matrix(0, nrow=J, ncol=p*K)
  chainSummary$xi0 = matrix(0, nrow=p,ncol=K)
  chainSummary$psi = matrix(0, nrow=p,ncol=K)
  Omega_raw = matrix(0, nrow=p, ncol=p*K)
  Sigma_raw = matrix(0, nrow=p, ncol=p*K)
  E_raw = matrix(0, nrow=p, ncol=p*K)
  chainSummary$alpha = matrix(0, nrow=p,ncol=K)
  chainSummary$W = matrix(0, nrow=J, ncol=K)
  for (i in 1:ns) {
    xi_raw = xi_raw + chain$xi[,,i]
    chainSummary$xi0 = chainSummary$xi0 + chain$xi0[,,i]
    chainSummary$psi = chainSummary$psi + chain$psi[,,i]
    Omega_raw = Omega_raw + chain$Omega[,,i]
    Sigma_raw = Sigma_raw + chain$G[,,i]
    E_raw = E_raw + chain$E[,,i]
    chainSummary$alpha = chainSummary$alpha + chain$alpha[,,i]
    chainSummary$W = chainSummary$W + chain$W[,,i]
  }

  chainSummary$W = chainSummary$W/ns

  xi_raw = xi_raw/ns
  chainSummary$xi = array(0, dim=c(J,p,K))
  for (k in 1:K) {
    chainSummary$xi[,,k] = xi_raw[,(1+p*(k-1)):(p*k)]
  }
  
  chainSummary$xi0 = chainSummary$xi0/ns
  chainSummary$psi = chainSummary$psi/ns

  Omega_raw = Omega_raw/ns
  chainSummary$Omega = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$Omega[,,k] = Omega_raw[,(1+p*(k-1)):(p*k)]
  }

  Sigma_raw = Sigma_raw/ns
  chainSummary$Sigma = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$Sigma[,,k] = Sigma_raw[,(1+p*(k-1)):(p*k)]
  }

  E_raw = E_raw/ns
  chainSummary$E = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$E[,,k] = E_raw[,(1+p*(k-1)):(p*k)]
  }
  
  chainSummary$alpha = chainSummary$alpha/ns

  chainSummary$meanvec = array(0, c(J, p, K))
  chainSummary$meanvec0 = matrix(0, p, K)
  for (k in 1:K) {
    del.om = transform_params(chainSummary$Omega[,,k], chainSummary$alpha[,k])
    # chainSummary$psi[,k] = del.om$psi
    chainSummary$Sigma[,,k] = del.om$Sigma
    for (j in 1:J) {
      chainSummary$meanvec[j,,k] = chainSummary$xi[j,,k] + del.om$omega %*% del.om$delta*sqrt(2/pi)
    }
    chainSummary$meanvec0[,k] = chainSummary$xi0[,k] + del.om$omega %*% del.om$delta*sqrt(2/pi)
  }

  chainSummary$t = apply(chain$t,2,Mode)
  chainSummary$S = apply(chain$S, 2, Mode)
  chainSummary$varphi = mean(chain$varphi)
  chainSummary$a0 = mean(chain$a0)
  
  chainSummary
}