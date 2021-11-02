# Point to point test of permutation test

p2pPerm <- function (x, p1 = NA, p2 = NA, p3 = NA, focalm=1, type='spd', plot=FALSE){
  if (!is.na(p1) & !is.na(p2)) {
    if (p2 > p1) {
      stop("the end point should be more recent than the start point")
    }
  }
  if (!type%in%c('spd','roc','res','ros'))
  {
    stop("The argument 'type' should be either 'spd', 'roc', 'res', or 'ros'.")
  }
  if (type=='spd'){
    index1 = match(p1, x$observed[[focalm]]$calBP)
    index2 = match(p2, x$observed[[focalm]]$calBP)
    p1.y = x$observed[[focalm]][index1, 2]
    p2.y = x$observed[[focalm]][index2, 2]
    if (plot) {
      plot(x, type=type, focalm=focalm)
      points(p1, p1.y, pch = 20, col="red", cex=1.2)
      points(p2, p2.y, pch = 20, col="red", cex=1.2)
      lines(x$observed[[focalm]]$calBP[index1:index2], 
            x$observed[[focalm]]$PrDens[index1:index2], 
            lwd = 2, col="red")
    }
    obs.diff = p1.y - p2.y
    sim.diff = x$raw[[focalm]][index1, ] - x$raw[[focalm]][index2, ]
    nsim = x$nsim
    lo = sum(obs.diff < sim.diff)
    hi = sum(obs.diff > sim.diff)  
    eq = sum(obs.diff == sim.diff)
    pvalHi = (lo + eq + 1)/c(nsim + 1)
    pvalLo = (hi + eq + 1)/c(nsim + 1)
    pval = ifelse(pvalHi < pvalLo, pvalHi, pvalLo) * 2 
    return(list(p1 = p1, p2 = p2, pval = pval))
  }
  if (type=='roc'){
    index1 = match(p1, x$observed.roc[[focalm]]$calBP)
    index2 = match(p2, x$observed.roc[[focalm]]$calBP) 
    p1.y = x$observed.roc[[focalm]][index1, 2] 
    p2.y = x$observed.roc[[focalm]][index2, 2] 
    if (plot) {
      plot(x, type=type, focalm=focalm)
      points(p1, p1.y, pch = 20, col="red")
      points(p2, p2.y, pch = 20, col="red")
      lines(x$observed.roc[[focalm]]$calBP[index1:index2],
            x$observed.roc[[focalm]]$roc[index1:index2], 
            lwd = 2, col="red")
    }
    obs.diff = abs(p1.y - p2.y)
    sim.diff = na.omit(abs(x$raw.row[[focalm]][index1, ] - x$raw.row[[focalm]][index2, ]))
    nsim = x$nsim
    lo = sum(obs.diff < sim.diff)
    hi = sum(obs.diff > sim.diff)
    eq = sum(obs.diff == sim.diff)
    pvalHi = (lo + eq + 1)/c(nsim + 1)
    pvalLo = (hi + eq + 1)/c(nsim + 1)
    pval = ifelse(pvalHi < pvalLo, pvalHi, pvalLo) * 2
    return(list(p1 = p1, p2 = p2, pval = pval))
  }
  if (type=='res'){
    
    index1 = match(p1, x$observed[[focalm]]$calBP)# index of beginning
    index2 = match(p2, x$observed[[focalm]]$calBP)  # index of end
    p1.y = x$observed[[focalm]][index1, 2] # PrDens beginning
    p2.y = x$observed[[focalm]][index2, 2] # PrDens end
    
    if (is.na(p3)){
      print("p3 unspecified, defaulting to minimum")
      p3 = x$observed[[focalm]][index1:index2, 1][which.min(x$observed[[focalm]][index1:index2, 2])] # Find year for min value of PrDens in range index1:index2
    }
    
    index3 = match(p3, x$observed[[focalm]]$calBP) # index of middle
    p3.y = x$observed[[focalm]][index3, 2] # PrDens middle
    
    if (plot) {
      plot(x, type="spd", focalm=focalm)
      points(p1, p1.y, pch = 20, col="#e66101", cex=1.5)
      points(p2, p2.y, pch = 20, col="#e66101", cex=1.5)
      lines(x$observed[[focalm]]$calBP[index1:index2], 
            x$observed[[focalm]]$PrDens[index1:index2],
            lwd = 2, col="#e66101")
      points(p3, p3.y, pch = 20, col="#5e3c99", cex=1.5)
      
    }
    d <- abs(p1.y - p3.y) #resistance
    dd <- 2 * d
    p <- abs(p1.y) + d # abs() because negative growth rates
    dp <- dd/p
    obs.rt <- 1-dp
    
    d.sim <- abs(x$raw[[focalm]][index1, ] - x$raw[[focalm]][index3, ])
    dd.sim <- 2 * d.sim
    p.sim <- abs(x$raw[[focalm]][index1, ]) + d.sim
    dp.sim <- d.sim/p.sim
    sim.rt <- 1-dp.sim
    
    e <- abs(p1.y - p2.y) #resilience
    e <- d + e
    e <- dd/e
    obs.rs <- e - 1
    
    e.sim <- abs(x$raw[[focalm]][index1, ] - x$raw[[focalm]][index2, ])
    e.sim <- dd.sim - e.sim
    e.sim <- dd.sim/e.sim
    sim.rs <- e.sim - 1
    
    nsim = x$nsim 
    
    lo.rt = sum(obs.rt < sim.rt) # resistance
    hi.rt = sum(obs.rt > sim.rt)
    eq.rt = sum(obs.rt == sim.rt)
    pvalHi.rt = (lo.rt + eq.rt + 1)/c(nsim + 1)
    pvalLo.rt = (hi.rt + eq.rt + 1)/c(nsim + 1)
    
    lo.rs = sum(obs.rs < sim.rs) # resilience
    hi.rs = sum(obs.rs > sim.rs)
    eq.rs = sum(obs.rs == sim.rs)
    pvalHi.rs = (lo.rs + eq.rs + 1)/c(nsim + 1)
    pvalLo.rs = (hi.rs + eq.rs + 1)/c(nsim + 1)
    
    return(list(pval.resistance=(ifelse(pvalHi.rt < pvalLo.rt, pvalHi.rt, pvalLo.rt) * 2),
                resistance = obs.rt,
                lag.resistance = p1-p3,
                pval.resilience=(ifelse(pvalHi.rs < pvalLo.rs, pvalHi.rs, pvalLo.rs) * 2),
                resilience = obs.rs) )
  }
  if (type=='ros'){
    
    index1 = match(p1, x$observed.roc[[focalm]]$calBP)# index of beginning
    index2 = match(p2, x$observed.roc[[focalm]]$calBP)  # index of end
    p1.y = x$observed.roc[[focalm]][index1, 2] # PrDens beginning
    p2.y = x$observed.roc[[focalm]][index2, 2] # PrDens end
    
    if (is.na(p3)){
      print("p3 unspecified, defaulting to minimum")
      p3 = x$observed.roc[[focalm]][index1:index2, 1][which.min(x$observed.roc[[focalm]][index1:index2, 2])] # Find year for min value of PrDens in range index1:index2
    }
    
    index3 = match(p3, x$observed.roc[[focalm]]$calBP) # index of middle
    p3.y = x$observed.roc[[focalm]][index3, 2] # PrDens middle
    
    if (plot) {
      plot(x, type="roc", focalm=focalm)
      points(p1, p1.y, pch = 20, col="#e66101", cex=1.5)
      points(p2, p2.y, pch = 20, col="#e66101", cex=1.5)
      lines(x$observed.roc[[focalm]]$calBP[index1:index2], 
            x$observed.roc[[focalm]]$roc[index1:index2],
            lwd = 2, col="#e66101")
      points(p3, p3.y, pch = 20, col="#5e3c99", cex=1.5)
      
    }
    
    d <- abs(p1.y - p3.y) #resistance
    dd <- 2 * d
    p <- abs(p1.y) + d # abs() because negative growth rates
    dp <- dd/p
    obs.rt <- 1-dp
    
    d.sim <- abs(x$raw.row[[focalm]][index1, ] - x$raw.row[[focalm]][index3, ])
    dd.sim <- 2 * d.sim
    p.sim <- abs(x$raw.row[[focalm]][index1, ]) + d.sim
    dp.sim <- d.sim/p.sim
    sim.rt <- 1-dp.sim
    
    e <- abs(p1.y - p2.y) #resilience
    e <- d + e
    e <- dd/e
    obs.rs <- e - 1
    
    e.sim <- abs(x$raw.row[[focalm]][index1, ] - x$raw.row[[focalm]][index2, ])
    e.sim <- dd.sim - e.sim
    e.sim <- dd.sim/e.sim
    sim.rs <- e.sim - 1
    
    nsim = x$nsim 
    
    lo.rt = sum(obs.rt < sim.rt) # resistance
    hi.rt = sum(obs.rt > sim.rt)
    eq.rt = sum(obs.rt == sim.rt)
    pvalHi.rt = (lo.rt + eq.rt + 1)/c(nsim + 1)
    pvalLo.rt = (hi.rt + eq.rt + 1)/c(nsim + 1)
    
    lo.rs = sum(obs.rs < sim.rs) # resilience
    hi.rs = sum(obs.rs > sim.rs)
    eq.rs = sum(obs.rs == sim.rs)
    pvalHi.rs = (lo.rs + eq.rs + 1)/c(nsim + 1)
    pvalLo.rs = (hi.rs + eq.rs + 1)/c(nsim + 1)
    
    return(list(pval.resistance=(ifelse(pvalHi.rt < pvalLo.rt, pvalHi.rt, pvalLo.rt) * 2),
                resistance = obs.rt,
                lag.resistance = p1-p3,
                pval.resilience=(ifelse(pvalHi.rs < pvalLo.rs, pvalHi.rs, pvalLo.rs) * 2),
                resilience = obs.rs) )
  }
}



