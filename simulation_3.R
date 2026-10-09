### Posterior inference of Recurrent HDP-HSMM with HMM embedding ###
## load library ##
library(MCMCpack)
library(MASS)
library(posterior)

## load data ##
datalist = readRDS("simdata_rseed-623.rds")

## setup ##
y = datalist$y
T = length(y)

N = 20

gam = 2
alp = 3

a = 0.8
b = 2
m = c(0, 0)
Lmd = diag(c(0.5, 0.5))

u = 1
v = 1
rmax = 15
nu = rep(1, rmax)/rmax

## Logarithm of the Stirling numbers of the first kind
lnstrl = matrix(-Inf, nrow=1000, ncol=1000)
lnstrl[1, 1] = 0

for (n in 2:dim(lnstrl)[1]) {
  log_n_minus_1 = log(n - 1)
  max_j = min(n, dim(lnstrl)[2])
  
  for (k in 1:max_j) {
    lnQ = c(
      log_n_minus_1 + lnstrl[n-1, k],
      if(k > 1) lnstrl[n-1, k-1] else -Inf
    )
    
    lnC <- max(lnQ)
    if (lnC != -Inf) {
      lnstrl[n, k] <- lnC + log(sum(exp(lnQ - lnC)))
    }
  }
}

## initialization ##
bet_ = rbeta(N, 1, gam)
bet_pos = c(bet_[1], sapply(2:N, function(l){bet_[l]*prod(1 - bet_[1:(l-1)])}))

A = matrix(0, nrow=N, ncol=N)
Abar = matrix(0, nrow=N, ncol=N)
for(j in 1:N) {
  A_ = sapply(1:N, function(k){rbeta(1, alp*bet_pos[k], alp*(1 - sum(bet_pos[1:k])))})
  A[j, ] = c(A_[1], sapply(2:N, function(l){A_[l]*prod(1 - A_[1:(l-1)])}))
  Abar[j, ] = A[j, ]/(sum(A[j, ]) - A[j, j])
  Abar[j, j] = 0
}

tau_pos = rgamma(N, shape=a, rate=b)
mu_pos = sapply(1:N, function(i){mvrnorm(1, mu=m, Sigma=solve(Lmd*tau_pos[i]))})
theta_pos = t(rbind(mu_pos, 1/sqrt(tau_pos)))

p_pos = rbeta(N, u, v)
r_pos = sample(1:rmax, size=N, replace=TRUE, nu)

## Posterior inference ##
numsampling = 200
burnin = 100

lnf = function(y, y_, theta) {dnorm(y, mean=y_*theta[1] + theta[2], sd=theta[3], log=TRUE)} # AR(1)
lng = function(d, r, p){dnbinom(d-1, size=r, prob=1-p, log=TRUE)} # the range of d is 1, 2... so shift it by 1 and match Johnson's parameterization
y_0 = 1 # define a value for the 0th sequence

z_pos_seq = vector(mode="list", length=numsampling)
d_pos_seq = vector(mode="list", length=numsampling)
x_pos_seq = matrix(0, nrow=T, ncol=numsampling)
mu_pos_seq = matrix(data=0, nrow=N*2, ncol=numsampling)
tau_pos_seq = matrix(data=0, nrow=N, ncol=numsampling)
r_pos_seq = matrix(data=0, nrow=N, ncol=numsampling)
p_pos_seq = matrix(data=0, nrow=N, ncol=numsampling)
A_pos_seq = matrix(data=0, nrow=N*N, ncol=numsampling)
bet_pos_seq = matrix(data=0, nrow=N, ncol=numsampling)
y_pos_seq = matrix(data=0, nrow=T, ncol=numsampling)

for(cntsampling in 1:numsampling) {
  if(cntsampling%%10 == 1) {
    cat(sprintf("%s - Number of sampling: %d/%d\n", Sys.time(), cntsampling, numsampling))
  }
  
  ## embedded HMM message passing ##
  lnB = matrix(0, nrow=T, ncol=N)
  lnBlik = matrix(0, nrow=T, ncol=N)
  lnBbar = sapply(1:N, function(i){matrix(0, nrow=T, ncol=r_pos[i])})
  lnc = sapply(1:N, function(i){dbinom((r_pos[i] - 1):0, r_pos[i] - 1, p_pos[i], log=TRUE)})
  lnAbar = log(Abar)
  lnp_pos = log(p_pos)
  ln1mp_pos = log1p(-p_pos)
  
  for(t in (T-1):1) {
    for(i in 1:N) {
      lnQ = lnc[[i]] + lnBbar[[i]][t+1, ]
      lnC = max(lnQ)
      lnBlik[t, i] = lnf(y[t+1], y[t], theta_pos[i, ]) + (lnC + log(sum(exp(lnQ - lnC))))
    }
    
    for(i in 1:N) {
      lnQ = lnBlik[t, ] + lnAbar[i, ]
      lnC = max(lnQ)
      lnB[t, i] = lnC + log(sum(exp(lnQ - lnC)))
    }
    
    for(i in 1:N) {
      lnf_i = lnf(y[t+1], y[t], theta_pos[i, ])
      
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
    lnBlik_0[i] = lnf(y[t+1], y_0, theta_pos[i, ]) + (lnC + log(sum(exp(lnQ - lnC))))
  }
  
  ## sampling z and d
  z_pos = c()
  d_pos = c()
  d_cens = c()
  s = 1
  t = 1
  lnAini = log(bet_pos)
  
  while(t <= T) {
    ## sampling z
    if(s == 1) {
      lnQ = lnAini + lnBlik_0
    } else {
      lnQ = lnAbar[z_pos[s-1], ] + lnBlik[t-1, ]
    }
    lnC = max(lnQ)
    p_z = exp(lnQ - lnC)/sum(exp(lnQ - lnC))
    z_pos[s] = sample(1:N, size=1, replace=TRUE, prob=p_z)
    
    ## sampling d
    lnp_d_pos = vector(mode="numeric", length=T-t+1)
    for(d in 1:(T-t+1)) {
      if(t == 1) {
        lnf_d = sum(lnf(y[t:(t+d-1)], c(y_0, y[(t-1):(t+d-2)]), theta_pos[z_pos[s], ]))
      } else {
        lnf_d = sum(lnf(y[t:(t+d-1)], y[(t-1):(t+d-2)], theta_pos[z_pos[s], ]))
      }
      
      lnp_d_pos[d] = 
        lng(d, r_pos[z_pos[s]], p_pos[z_pos[s]]) + 
        lnf_d + 
        lnB[t+d-1, z_pos[s]]
    }
    
    if(t == 1) {
      lnf_d = sum(lnf(y[t:T], c(y_0, y[(t-1):(T-1)]), theta_pos[z_pos[s], ]))
    } else {
      lnf_d = sum(lnf(y[t:T], y[(t-1):(T-1)], theta_pos[z_pos[s], ]))
    }
    
    lnp_d_cens_pos = 
      pnbinom(T-t-1, size=r_pos[z_pos[s]], prob=1-p_pos[z_pos[s]], lower.tail=FALSE, log.p=TRUE) +
        lnf_d # recall B(T > t) = 1 
    
    lnQ = c(lnp_d_pos, lnp_d_cens_pos)
    lnC = max(lnQ)
    p_d = exp(lnQ - lnC)/sum(exp(lnQ - lnC))
    d_pos_draw = sample(1:(T-t+2), size=1, replace=TRUE, prob=p_d)
    
    if(d_pos_draw <= (T-t+1)) {
      d_pos[s] = d_pos_draw
      d_cens[s] = FALSE
    } else {
      d_pos[s] = T-t+1
      d_cens[s] = TRUE
    }
    
    ## updating indices
    t = t + d_pos[s]
    s = s + 1
  }
  
  z_pos_seq[[cntsampling]] = z_pos
  d_pos_seq[[cntsampling]] = d_pos
  x_pos_seq[, cntsampling] = unlist(sapply(1:(s-1), function(i){rep(z_pos[i], d_pos[i])}))
  
  ## sampling mu, tau, r, p, and A
  z_pos_set = unique(z_pos)
  n_z = matrix(0, nrow=N, ncol=N)
  for(i in 1:N) {
    if(i %in% z_pos_set) {
      ## posterior draws of mu and sgm
      idx_i = which(x_pos_seq[, cntsampling] == i)
      n = length(idx_i)
      y_it = y[idx_i]
        
      if(1 %in% idx_i) {
        y_it_ = c(y_0, y[idx_i - 1])
      } else {
        y_it_ = y[idx_i - 1]
      }
      
      S_x = sum(y_it_)
      S_y = sum(y_it)
      S_xx = sum(y_it_^2)
      S_yy = sum(y_it^2)
      S_xy = sum(y_it_*y_it)
      XX = matrix(data=c(S_xx, S_x, S_x, n), nrow=2, ncol=2)
      Xy = c(S_xy, S_y)
      
      Lmd_n = Lmd + XX
      m_n = solve(Lmd_n)%*%(Lmd%*%m + Xy)
      a_n = a + n/2
      b_n = b + 0.5*(S_yy + t(m)%*%Lmd%*%m - t(m_n)%*%Lmd_n%*%m_n)
      
      tau_pos[i] = rgamma(1, shape=a_n, rate=b_n)
      mu_pos[, i] = mvrnorm(1, mu=m_n, Sigma=solve(Lmd_n*tau_pos[i]))
      theta_pos[i, ] = c(mu_pos[, i], 1/sqrt(tau_pos[i]))
      
      ## posterior draws of r and p
      idx_i = z_pos ==i & (!d_cens) # duration of censoring term is not identifiable so remove it
      D = sum(idx_i)
      if(D == 0) { # censored state case
        r_pos[i] =  sample(1:rmax, size=1, replace=TRUE, nu)
        p_pos[i] = rbeta(1, u, v)
      } else {
        d_pos_i = d_pos[idx_i]
        u_n = u + sum(d_pos_i - 1)
        
        lnp_r = vector(mode="numeric", length=rmax)
        for(r in 1:rmax) {
          lnp_r[r] = 
            log(nu[r]) + 
            sum(sapply(1:D, function(k){lchoose(d_pos_i[k] + r - 2, d_pos_i[k] - 1)})) +
            lbeta(u_n, v + r*D)
        }
        lnC = max(lnp_r)
        p_r = exp(lnp_r - lnC)/sum(exp(lnp_r - lnC))
        r_pos[i] = sample(1:rmax, size=1, replace=TRUE, prob=p_r)
        
        v_n = v + r_pos[i]*D   # n_s is reduced because duration is interpreted as the number of trials + 1
        p_pos[i] = rbeta(1, u_n, v_n)
      } 
      
      ## posterior draws of A
      j = z_pos[which(z_pos[-length(z_pos)] == i) + 1]
      if(length(j) == 0) {
        A[i, ] = rdirichlet(1, alp*bet_pos) # appeared only in the final state and did no enter into other states
      } else {
        for(l in 1:length(j)) n_z[i, j[l]] = n_z[i, j[l]] + 1
        rho = rgeom(length(j), min(1, 1 - A[i, i] + 1e-15)) # to avoid 0 due to underflow
        n_z[i, i] = n_z[i, i] + sum(rho)
        A[i, ] = rdirichlet(1, alp*bet_pos + n_z[i, ])
      }
    } else {
      # redraw from prior
      tau_pos[i] = rgamma(1, shape=a, rate=b)
      mu_pos[, i] = mvrnorm(1, mu=m, Sigma=solve(Lmd*tau_pos[i]))
      theta_pos[i, ] = c(mu_pos[, i], 1/sqrt(tau_pos[i]))
      r_pos[i] = sample(1:rmax, size=1, replace=TRUE, nu)
      p_pos[i] = rbeta(1, u, v)
      A[i, ] = rdirichlet(1, alp*bet_pos)
    }
    
    Abar[i, ] = A[i, ]/(1 - A[i, i] + 1e-15) # to avoid 0 division due to underflow
    Abar[i, i] = 0
    Abar[i, ] = Abar[i, ]/sum(Abar[i, ])
  }
  
  ## sampling beta
  m_pos = matrix(0, nrow=N, ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      if(n_z[i, j] > 0) {
        m_ij = 1:n_z[i, j]
        lnp_m = lnstrl[n_z[i, j], m_ij] + m_ij*log(alp*bet_pos[j])
        lnC = max(lnp_m)
        m_pos[i, j] = sample(m_ij, 1, replace=TRUE, exp(lnp_m - lnC)/sum(exp(lnp_m - lnC)))
      }
    }
  }
  
  bet_pos = rdirichlet(1, gam/N + colSums(m_pos))
  
  ## record samples
  mu_pos_seq[, cntsampling] = c(mu_pos)
  tau_pos_seq[, cntsampling] = tau_pos
  r_pos_seq[, cntsampling] = r_pos
  p_pos_seq[, cntsampling] = p_pos
  A_pos_seq[, cntsampling] = c(A)
  bet_pos_seq[, cntsampling] = bet_pos
  
  ## plot
  if(cntsampling%%10 == 1) {
    mu_pos_w_1 = mu_pos[1, x_pos_seq[, cntsampling]]
    mu_pos_w_2 = mu_pos[2, x_pos_seq[, cntsampling]]
    sgm_pos_w = 1/sqrt(tau_pos[x_pos_seq[, cntsampling]])
    y_pos_seq[1, cntsampling] = rnorm(1, y_0*mu_pos_w_1[1] + mu_pos_w_2[1], sgm_pos_w[1])
    for(t in 2:T)  y_pos_seq[t, cntsampling] = rnorm(1, y_pos_seq[t-1, cntsampling]*mu_pos_w_1[t] + mu_pos_w_2[t], sgm_pos_w[t])
    
    plot(y_pos_seq[, cntsampling], type="l")
    lines(y, col="blue", lty="dotted")
  }
}

## rhat ##
# Needs to implement relabeling for label-switching and label birth-death

## Plot (y)
for(w in (burnin+1):numsampling) {
  idx_i = x_pos_seq[, w]
  mu_pos_w = mu_pos_seq[, w]
  mu_pos_w_1 = mu_pos_w[seq(from=1, to=N*2, by=2)]
  mu_pos_w_2 = mu_pos_w[seq(from=2, to=N*2, by=2)]
  mu_pos_w_1 = mu_pos_w_1[idx_i]
  mu_pos_w_2 = mu_pos_w_2[idx_i]
  sgm_pos_w = 1/sqrt(tau_pos_seq[idx_i, w])
  
  y_pos_seq[1, w] = rnorm(1, y_0*mu_pos_w_1[1] + mu_pos_w_2[1], sgm_pos_w[1])
  
  for(t in 2:T)  y_pos_seq[t, w] = rnorm(1, y_pos_seq[t-1, w]*mu_pos_w_1[t] + mu_pos_w_2[t], sgm_pos_w[t])
}
y_seq_hat = sapply(1:T, function(t){mean(y_pos_seq[t, (burnin+1):numsampling])})

plot(y_seq_hat, type="l")
lines(y, col="blue", lty="dotted")

plot(y_seq_hat - y, type="l")

hist(y_seq_hat - y)

print(var(y_seq_hat - y))