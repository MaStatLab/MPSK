Mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Calibrate <- function(x)
{
  C = x$data$C - 1
  Z = x$chain$t - 1
  ns  = dim(x$chain$S)[3]
  K = dim(x$chain$S)[2]

  S = matrix(0, nrow=ns, ncol=K)
  for (i in 1:ns) {
    S[i,] = apply(x$chain$S[,,i], 2, Mode)
  }
  output = calib(x$data$Y,
                 matrix(C,ncol=1),
                 Z,
                 x$chain$xi, dim(x$chain$xi),
                 x$chain$xi0, dim(x$chain$xi0),
                 S )
  colnames(output$Y_cal) = colnames(x$data$Y)
  return(output)

}
