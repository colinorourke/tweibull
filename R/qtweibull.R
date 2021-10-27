qtweibull = function(p, a=0, b=Inf, shape, scale, log.p=FALSE, ...){
  log_p = if(isFALSE(log.p)) log(p) else p
  res = a^shape - scale^shape * log1mpexp(log_p, a, b, shape, scale)

  res^(1/shape)
}
