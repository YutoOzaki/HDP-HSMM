### Embedded HMM message passing ###
## load data ##
datalist = readRDS("simdata_rseed-88.rds")

## setup ##
y = datalist$y
T = length(y)

u = 1
v = 1
rmax = 10
nu = rep(1, rmax)/rmax

## initialization ##
N = 20
Abar = datalist$Abar[1:N, 1:N]
mu = datalist$mu[1:N]
tau = datalist$tau[1:N]
theta = cbind(mu, 1/sqrt(tau))

p_pos = rbeta(N, u, v)
r_pos = sample(1:rmax, size=N, replace=TRUE, nu)

## embedded HMM message passing ##
print(Sys.time())
lnB = matrix(0, nrow=T, ncol=N)
lnBlik = matrix(0, nrow=T, ncol=N)

lnBbar = sapply(1:N, function(i){matrix(0, nrow=T, ncol=r_pos[i])})
for(i in 1:N) lnBbar[[i]][T, ] = 0

lnc = sapply(1:N, function(i){dbinom((r_pos[i] - 1):0, r_pos[i] - 1, p_pos[i], log=TRUE)})
lnAbar = log(Abar)
lnp_pos = log(p_pos)
ln1mp_pos = log1p(-p_pos)
lnf = function(y, theta) {dnorm(y, mean=theta[1], sd=theta[2], log=TRUE)}

for(t in (T-1):1) {
  for(i in 1:N) {
    lnQ = lnc[[i]] + lnBbar[[i]][t+1, ]
    lnC = max(lnQ)
    lnBlik[t, i] = lnf(y[t+1], theta[i, ]) + (lnC + log(sum(exp(lnQ - lnC))))
  }
  
  for(i in 1:N) {
    lnQ = lnBlik[t, ] + lnAbar[i, ]
    lnC = max(lnQ)
    lnB[t, i] = lnC + log(sum(exp(lnQ - lnC)))
  }
  
  for(i in 1:N) {
    lnf_i = lnf(y[t+1], theta[i, ])
    
    lnQ = c(
      lnp_pos[i] + lnf_i + lnBbar[[i]][t+1, r_pos[i]],
      ln1mp_pos[i] + lnAbar[i, ] + lnBlik[t, ]
    )
    lnC = max(lnQ)
    lnBbar[[i]][t, r_pos[i]] = lnC + log(sum(exp(lnQ - lnC)))
    
    if(r_pos[i] > 1) {
      for(j in (r_pos[i]-1):1) {
        lnQ = c(lnp_pos[i] + lnBbar[[i]][t+1, j], ln1mp_pos[i] + lnBbar[[i]][t+1, j+1])
        lnC = max(lnQ)
        lnBbar[[i]][t, j] = lnf_i + lnC + log(sum(exp(lnQ - lnC)))
      }
    }
  }
}

t = 0
lnBlik_0 = vector(mode="numeric", length=N)
for(i in 1:N) {
  lnQ = lnc[[i]] + lnBbar[[i]][t+1, ]
  lnC = max(lnQ)
  lnBlik_0[i] = lnf(y[t+1], theta[i, ]) + (lnC + log(sum(exp(lnQ - lnC))))
}
print(Sys.time())

## HSMM message passing ##
print(Sys.time())
lnB_hsmm = matrix(0, nrow=T, ncol=N)
lnBlik_hsmm = matrix(0, nrow=T, ncol=N)
lnp_D = function(i, d){dnbinom(d-1, size=r_pos[i], prob=1-p_pos[i], log=TRUE)} # the range of d is 1, 2... so shift it by 1

for(t in (T-1):1) {
  for(i in 1:N) {
    lnQ = vector(mode="numeric", length=T-t)
    for(d in 1:(T-t)) {
      lnQ[d] = lnB_hsmm[t+d, i] + lnp_D(i, d) + sum(lnf(y[(t+1):(t+d)], theta[i, ]))
    }
    
    # censoring term P(D>T-t | x_{t+1}=i) is survival function
    lnp_D_censored = pnbinom(T-t-1, size=r_pos[i], prob=1-p_pos[i], lower.tail=FALSE, log.p=TRUE)
    lnR = lnp_D_censored + sum(lnf(y[(t+1):T], theta[i, ]))
    lnC = max(c(lnQ, lnR))
    lnBlik_hsmm[t, i] = lnC + log(sum(exp(c(lnQ, lnR) - lnC)))
  }
  
  for(i in 1:N) {
    lnQ = lnBlik_hsmm[t, ] + lnAbar[i, ]
    lnC = max(lnQ)
    lnB_hsmm[t, i] = lnC + log(sum(exp(lnQ - lnC)))
  }
}

t = 0
lnBlik_0_hsmm = vector(mode="numeric", length=N)
for(i in 1:N) {
  lnQ = vector(mode="numeric", length=T-t)
  for(d in 1:(T-t)) {
    lnQ[d] = lnB_hsmm[t+d, i] + lnp_D(i, d) + sum(lnf(y[(t+1):(t+d)], theta[i, ]))
  }
  
  lnp_D_censored = pnbinom(T-t-1, size=r_pos[i], prob=1-p_pos[i], lower.tail=FALSE, log.p=TRUE)
  lnR = lnp_D_censored + sum(lnf(y[(t+1):T], theta[i, ]))
  lnC = max(c(lnQ, lnR))
  lnBlik_0_hsmm[i] = lnC + log(sum(exp(c(lnQ, lnR) - lnC)))
}
print(Sys.time())