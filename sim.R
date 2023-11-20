gendata = function(max.size, 
                   n.study, 
                   prob, 
                   ind = NULL) {
  if (is.null(ind)) ind = rep(TRUE, n.study)
  size = rep(0, n.study)
  names(size) = 1:n.study
  a = 1
  while (length(unique(a)) == 1) {
    a = sort(sample((1:n.study)[ind], max.size, replace = TRUE))
  }
  size[names(size) %in% names(table(a))] = table(a)
  resp = rbinom(length(size), size, prob)
  
  data.frame(size,resp)
}

sampler = function(response, 
                   stage.size, 
                   niter, 
                   nburn, 
                   nchain = 2, 
                   nthin = 1) {
  idst = stage.size != 0
  jags_data = list( "y" = response[idst], 
                    "N" = length(response[idst]), 
                    "n" = stage.size[idst] ) 
  jags_inits = function() list("mu" = 1, "tau2" = rep(1, jags_data$N))
  params = c("prob","mu","tau2")
  jags_fit = jagsUI::jags(data = jags_data, inits = jags_inits, 
                          parameters.to.save = params, 
                          model.file = model.file, 
                          n.chains = nchain, n.iter = niter, verbose = FALSE,
                          n.burnin = nburn, n.thin = nthin)
  ww = jags_fit$sims.list$tau2
  pooled = rowSums((ww*jags_fit$sims.list$prob))/rowSums(ww)
  
  list(post.p = pooled, 
       max.rhat = max(abs(unlist(jags_fit$Rhat) - 1)))
}

library(rjags)
t1err = 0.078
t2err = 0.15 

nsim = 2000
niter = 40000
nburn = 20000

model.file = "rehb1.jags"
mu = qlogis(0.1)
maxnum.stg1 = 8
maxnum.stg2 = 24 - maxnum.stg1

pooled.se = max.rhat = efficacy = pooled.ci = pooled.prob = avg.samp = exd.samp = futi.ind = NA
for (i in 1:nsim) {
  set.seed(i)
  n1 = 28; n2 = 24; n3 = 24
  n = c(n1,n2,n3)
  p1 = plogis(mu + rnorm(1, sd = 0.1))
  p2 = plogis(mu + rnorm(1, sd = 0.1))
  p3 = plogis(mu + rnorm(1, sd = 0.1))
  p0 = c(p1,p2,p3)
  cutoff.U = 1 - t2err
  cutoff.L = t1err
  s1dat = gendata(max.size = maxnum.stg1, n.study = 3, prob = p0)
  stg1.size = s1dat$size
  y1 = s1dat$resp
  avg.samp[i] = sum(stg1.size)
  st1_jags = sampler(response = y1,
                     stage.size = stg1.size,
                     niter = niter, nburn = nburn)
  st1_post_p = st1_jags$post.p
  futi.ind[i] = mean(st1_post_p > 0.3) < cutoff.L
  exd.samp[i] = maxnum.stg1 + mean(st1_post_p > 0.1)*maxnum.stg2
  
  s2data = gendata(maxnum.stg2, n.study = 3, prob=p0, ind = (stg1.size != 0))
  stg2.size = s2data$size
  y2 = s2data$resp
  y = y1 + y2
  stg.size = stg1.size + stg2.size
  if (futi.ind[i] == FALSE) {
    avg.samp[i] = avg.samp[i] + sum(stg2.size)
    st2_jags = sampler(response = y,
                       stage.size = stg.size,
                       niter = niter, nburn = nburn)
    st2_post_p = st2_jags$post.p
    max.rhat[i] = st2_jags$max.rhat
    pooled.prob[i] = mean(st2_post_p)
    pooled.se[i] = sd(st2_post_p)
    pooled.ci[i] = quantile(st2_post_p, 1 - t1err)
    efficacy[i] = quantile(st2_post_p, t1err) > 0.1
  }
}
mean(pooled.prob, na.rm = TRUE)
sd(pooled.prob, na.rm = TRUE)
mean(pooled.se, na.rm = TRUE)
mean(pooled.ci, na.rm = TRUE)
mean(futi.ind, na.rm = TRUE)
mean(efficacy, na.rm = TRUE)
mean(avg.samp, na.rm = TRUE)
mean(exd.samp, na.rm = TRUE)
max(max.rhat, na.rm = TRUE)