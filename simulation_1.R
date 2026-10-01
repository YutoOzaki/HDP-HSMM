### Posterior inference of z, D, mu, tau, p, pi, and beta ###
## load libraries ##
library(MCMCpack)

## read simulation data
datalist = readRDS("simdata_rseed-25.rds")
y = datalist$y
T = length(y)

## parameters
numsampling = 100

# Gaussian emission distribution (likelihood)
lnf = function(y, mu, sgm) {dnorm(y, mu, sgm, log=TRUE)}

# parameters of normal-gamma prior for the likelihood
m = 0
a = 0.1
b = 0.1
lmd = 0.1

# parameters of beta prior for the duration distribution (negative binomial)
c = 1
d = 1

# parameters of HDP
K = 25
gam = 2
al = 2

## initialization
tau_pos = rgamma(K, a, b)
mu_pos = rnorm(K, m, 1/sqrt(lmd*tau_pos))
sgm_pos = 1/sqrt(tau_pos)
p_pos = rbeta(K, c, d)

bet = rbeta(K, 1, gam)
bet[2:K] = sapply(2:K, function(i){bet[i]*prod(1 - bet[1:(i-1)])})
bet = bet/sum(bet)

pi = matrix(0, nrow=K, ncol=K)
for(i in 1:K) {
  for(j in 1:K) {
    pi[i, j] = rbeta(1, al*bet[j], al*(1 - sum(bet[1:j])))
  }
  pi[i, 2:K] = sapply(2:K, function(j){pi[i, j]*prod(1 - pi[i, 1:(j-1)])})
}
pi = pi/rowSums(pi)

pi_bar = pi*0
for(i in 1:K) {
  pi_bar[i, ] = pi[i, ]/(1 - pi[i, i])
  pi_bar[i, i] = 0
}

lnbet = log(bet)
lnpi = log(pi_bar)

## Stirling numbers of the first kind
s = matrix(0, nrow=50, ncol=50)
s[1, 1] = 1
for(i in 1:49) {
  for(j in 1:50) {
    if(j > 1) {
      s[i+1, j] = s[i, j-1] + i*s[i, j]
    } else {
      s[i+1, j] = i*s[i, j]
    }
  }
}

## posterior sampling
z_pos_list = vector(mode="list", length=numsampling)
z_pos_seq = matrix(data=0, nrow=numsampling, ncol=T)
mu_pos_seq = matrix(data=0, nrow=numsampling, ncol=K)
sgm_pos_seq = matrix(data=0, nrow=numsampling, ncol=K)
p_pos_seq = matrix(data=0, nrow=numsampling, ncol=K)
pi_pos_seq = matrix(data=0, nrow=numsampling, ncol=K*K)
bet_pos_seq = matrix(data=0, nrow=numsampling, ncol=K)

for(cntsampling in 1:numsampling) {
  if(cntsampling%%10 == 1) {
    cat(sprintf("%s - Number of sampling: %d/%d\n", Sys.time(), cntsampling, numsampling))
  }
  
  ## message passing
  # t = 0, 1, 2, ..., T for B and B* (B_ is used as a variable name for the latter)
  # t = 1, 2, ..., T for y
  lnB = matrix(0, nrow=K, ncol=T+1) 
  lnB_ = matrix(0, nrow=K, ncol=T+1)
  lnB[, T+1] = log(1)
  
  # the support of the distribution starts with 0, 1, 2... but it is interpreted as duration of 1, 2, 3...
  lnp_D = sapply(0:(T-1), function(d){dnbinom(d, 1, prob=p_pos[1:K], log=TRUE)})
  
  for(t in (T-1):0) {
    for(i in 1:K) {
      lnQ = vector(mode="numeric", length=T-t)
      for(d in 1:(T-t)) {
        lnQ[d] = lnB[i, t+1+d] + lnp_D[i, d] + sum(lnf(y[(t+1):(t+d)], mu_pos[i], sgm_pos[i]))
      }
      
      # survival function of the censoring term P(D>T-t | x_{t+1}=i) based on NB(1, p): 1 - (1-p)^(T-t)
      lnp_D_censored = (T-t)*log1p(-p_pos[i])
      lnR = lnp_D_censored + sum(lnf(y[(t+1):T], mu_pos[i], sgm_pos[i]))
      
      lnC = max(c(lnQ, lnR))
      
      lnB_[i, t+1] = log(sum(exp(c(lnQ, lnR) - lnC))) + lnC
    }
    
    for(i in 1:K) {
      lnQ = lnB_[, t+1] + lnpi[i, ]
      lnC = max(lnQ[!is.infinite(lnQ)])
      lnB[i, t+1] = log(sum(exp(lnQ - lnC))) + lnC
    }
  }
  
  ## sampling z and D
  t_pos = c()
  z_pos = c()
  D_pos = c()
  D_cens = c()
  i = 1
  t = 1
  
  while(t <= T) {
    t_pos[i] = t
    
    ## sampling z
    if(t == 1) {
      lnQ = lnbet + lnB_[, t]
    } else {
      lnQ = lnpi[z_pos[i-1], ] + lnB_[, t]
    }
    lnC = max(lnQ)
    p_x = exp(lnQ - lnC)/sum(exp(lnQ - lnC))
    z_pos[i] = sample(1:K, size=1, replace=TRUE, prob=p_x)
    
    ## sampling D
    lnp_D_pos = vector(mode="numeric", length=T-t+1)
    for(d in 1:(T-t+1)) {
      lnp_D_pos[d] = lnp_D[z_pos[i], d] + sum(lnf(y[t:(t+d-1)], mu_pos[z_pos[i]], sgm_pos[z_pos[i]])) + lnB[z_pos[i], t+d]
    }
    lnp_D_cens_pos = (T-t+1)*log1p(-p_pos[z_pos[i]]) + sum(lnf(y[t:T], mu_pos[z_pos[i]], sgm_pos[z_pos[i]]))  #   recall B(T > t) = 1 
    lnp_D_all_pos = c(lnp_D_pos, lnp_D_cens_pos)
    lnC = max(lnp_D_all_pos)
    p_D_pos = exp(lnp_D_all_pos - lnC)/sum(exp(lnp_D_all_pos - lnC))
    D_pos_draw = sample(1:(T-t+2), size=1, replace=TRUE, prob=p_D_pos)
    
    if(D_pos_draw <= (T-t+1)) {
      D_pos[i] = D_pos_draw
      D_cens[i] = FALSE
    } else {
      D_pos[i] = T-t+1
      D_cens[i] = TRUE
    }
    
    ## updating indices
    t = t + D_pos[i]
    i = i + 1
  }
  
  z_pos_list[[cntsampling]] = z_pos
  z_pos_seq[cntsampling, ] = unlist(sapply(1:length(z_pos), function(i){rep(z_pos[i], D_pos[i])}))[1:T]
  cat(sprintf(" z: "))
  cat(z_pos)
  cat("\n")
  
  ## sampling mu, tau, p, and pi
  z_pos_set = unique(z_pos_seq[cntsampling, ])
  n_z = matrix(0, nrow=K, ncol=K)
  for(i in 1:K) {
    if(i %in% z_pos_set) {
      ## posterior draws of mu and sgm
      n = sum(z_pos_seq[cntsampling, ] == i)
      y_i = y[z_pos_seq[cntsampling, ] == i]
      y_bar = mean(y_i)
      
      m_n = ((lmd*m) + n*y_bar)/(lmd + n)
      lmd_n = lmd + n
      a_n = a + n/2
      b_n = b + 0.5*sum((y_i - y_bar)^2) + (lmd*n*(y_bar - m)^2)/(2*(lmd + n))
      
      tau_pos[i] = rgamma(1, a_n, b_n)
      mu_pos[i] = rnorm(1, m_n, 1/sqrt(lmd_n*tau_pos[i]))
      sgm_pos[i] = 1/sqrt(tau_pos[i])
      
      ## posterior draws of p
      n_s = sum(z_pos == i & !D_cens)   # duration of censoring term is not identifiable so remove it
      c_n = c + n_s
      d_n = d + n - n_s   # n_s is reduced because duration is interpreted as the number of trials + 1
      p_pos[i] = rbeta(1, c_n, d_n)
      
      ## posterior draws of pi
      j = z_pos[which(z_pos == i) + 1]
      if(length(j) == 1 && is.na(j)) {
        pi[i, ] = rdirichlet(1, al*bet)
      } else {
        for(l in 1:length(j)) n_z[i, j[l]] = n_z[i, j[l]] + 1
        n_zi = n_z[i, ]
        rho = rgeom(length(j), min(1, pi[i, i] + 1e-15)) # to avoid 0 due to underflow
        n_zi[i] = n_zi[i] + sum(rho)
        pi[i, ] = rdirichlet(1, al*bet + n_zi)
      }
    } else {
      # redraw from prior
      tau_pos[i] = rgamma(1, a, b)
      mu_pos[i]  = rnorm(1, m, 1/sqrt(lmd*tau_pos[i]))
      sgm_pos[i] = 1/sqrt(tau_pos[i])
      p_pos[i]   = rbeta(1, c, d)
      pi[i, ] = rdirichlet(1, al*bet)
    }
    
    pi_bar[i, ] = pi[i, ]/(1 - pi[i, i] + 1e-15) # to avoid 0 division due to underflow
    pi_bar[i, i] = 0
    pi_bar[i, ] = pi_bar[i, ]/sum(pi_bar[i, ])
  }
  
  lnpi = log(pi_bar)
  
  ## sampling beta
  m_pos = matrix(0, nrow=K, ncol=K)
  for(i in 1:K) {
    for(j in 1:K) {
      if(n_z[i, j] > 0) {
        m_ij = 1:n_z[i, j]
        p_m = s[n_z[i, j], m_ij]*(al*bet[j])^m_ij
        m_pos[i, j] = sample(m_ij, 1, replace=TRUE, p_m/sum(p_m))
      }
    }
  }
  
  bet = rdirichlet(1, gam/K + colSums(m_pos))
  
  lnbet = log(bet)
  
  ## record samples
  mu_pos_seq[cntsampling, ] = mu_pos
  sgm_pos_seq[cntsampling, ] = sgm_pos
  p_pos_seq[cntsampling, ] = p_pos
  pi_pos_seq[cntsampling, ] = c(pi)
  bet_pos_seq[cntsampling, ] = bet
  
  ## plot
  mu_pos_seq_hat = mu_pos[z_pos_seq[cntsampling, ]]
  plot(mu_pos_seq_hat, type="l")
  lines(y, col="blue", lty="dotted")
}

## plot
mu = datalist$mu
z = datalist$z
D = datalist$D

# mu
mu_pos_seq_hat = mu_pos[sapply(1:T, function(t){median(z_pos_seq[(numsampling-4):numsampling, t])})]
mu_seq = unlist(sapply(1:length(z), function(i){rep(mu[z[i]], D[i])}))

plot(mu_pos_seq_hat, type="l")
lines(mu_seq, col="blue", lty="dotted")
lines(mu_pos_seq_hat - mu_seq, col="red", lty="dashed")