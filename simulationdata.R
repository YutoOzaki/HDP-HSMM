### Simulation of the generative process of HDP-HSMM ###
## configuration
save_result = TRUE
rseed = 25
set.seed(rseed)

## parameters
gam = 2.5
alp = 2.3
a = 3.0
b = 2.0
c = 3.0
d = 8.0
m = 0.0
lmd = 0.06

K = 50
L = 400

## β ~ GEM(γ)
b = rbeta(K, 1, gam)
bet = vector(mode="numeric", length=K)
C_sbp = 1
for(k in 1:K) {
  bet[k] = b[k]*C_sbp
  C_sbp = C_sbp*(1 - b[k])
}

## π_i ~ DP(α, β)
pi_dp = matrix(0, nrow=K, ncol=K)
b_dp = matrix(0, nrow=K, ncol=K)
C_dp = vector(mode="numeric", length=K) + 1
for(i in 1:K) {
  for(k in 1:K) {
    b_dp[i, k] = rbeta(1, alp*bet[k], alp*(1 - sum(bet[1:k])))
    pi_dp[i, k] = b_dp[i, k]*C_dp[i]
    C_dp[i] = C_dp[i]*(1 - b_dp[i, k])
  }
}

pi_bar = pi_dp*0
for(i in 1:K) {
  pi_bar[i, ] = pi_dp[i, ]/(1 - pi_dp[i, i])
  pi_bar[i, i] = 0
}

## (θ_i, ω_i) ~ H x G
tau <- rgamma(K, shape=a, rate=b)
mu <- rnorm(K, mean=m, sd=1/sqrt(lmd*tau))
p <- rbeta(K, c, d)

## sampling y
y = vector(mode="numeric", length=L)
t = c()
z = c()
D = c()
s = 1
t[s] = 1
z[s] <- sample(1:K, size=1, replace=TRUE, prob=bet)
D[s] <- rnbinom(n=1, size=1, prob=p[z[s]]) + 1

while((t[s] + D[s] - 1) < L) {
  y[t[s]:(t[s] + D[s] - 1)] = rnorm(D[s], mean=mu[z[s]], sd=1/sqrt(tau[z[s]]))
  s = s + 1
  t[s] = t[s-1] + D[s-1]
  z[s] = sample(1:K, size=1, replace=TRUE, prob=pi_bar[z[s-1], ])
  D[s] = rnbinom(n=1, size=1, prob=p[z[s]]) + 1
}

y[t[s]:L] = rnorm(L - t[s] + 1, mean=mu[z[s]], sd=1/sqrt(tau[z[s]]))
D[s] = L - t[s] + 1

# plot
y_mu = unlist(sapply(1:length(z), function(i){rep(mu[z[i]], D[i])}))

plot(y, type="l", lty="dotted")
lines(y_mu, col="blue")

# save
if(save_result) {
  saveRDS(
    list(y=y, z=z, D=D, t=t, tau=tau, mu=mu, p=p, bet=bet, pi_bar=pi_bar, pi_dp=pi_dp),
    file=paste("./simdata_rseed-", rseed, ".rds", sep=""))
}