### Simulation data for embedded HMM message passing ###
## configuration
save_result = TRUE
rseed = 88
set.seed(rseed)

## setup ##
T = 300
N = 50
gam = 2
alp = 3
a = 3
b = 2
kap = 0.5
m = 0
c = 1
d = 1
rmax = 6
nu = rep(1, rmax)/rmax

## simulation data sampling ##
bet_ = rbeta(N, 1, gam)
bet = c(bet_[1], sapply(2:N, function(l){bet_[l]*prod(1 - bet_[1:(l-1)])}))

A = matrix(0, nrow=N, ncol=N)
Abar = matrix(0, nrow=N, ncol=N)
for(j in 1:N) {
  A_ = sapply(1:N, function(k){rbeta(1, alp*bet[k], alp*(1 - sum(bet[1:k])))})
  A[j, ] = c(A_[1], sapply(2:N, function(l){A_[l]*prod(1 - A_[1:(l-1)])}))
  Abar[j, ] = A[j, ]/(sum(A[j, ]) - A[j, j])
  Abar[j, j] = 0
}

P_ = t(A) - diag(N)
P_[N, ] = rep(1, N)
Aini = solve(P_, c(rep(0, N-1), 1))
Aini[Aini < 0] = 0

tau = rgamma(N, shape=a, rate=b)
mu = rnorm(N, mean=m, sd=1/sqrt(kap*tau))

p = rbeta(N, c, d)
r = sample(1:rmax, size=N, replace=TRUE, nu)

s = 1
t = 1
z = sample(1:N, size=1, replace=TRUE, Aini)
d = rnbinom(n=1, size=r[z[s]], prob=p[z[s]]) + 1
x = rep(z[s], d[s])
y = rnorm(d[s], mu[z[s]], 1/sqrt(tau[z[s]]))
while((t[s] + d[s]) <= T) {
  s = s + 1
  t[s] = t[s-1] + d[s-1]
  z[s] = sample(1:N, size=1, replace=TRUE, Abar[z[s-1], ])
  d[s] = rnbinom(n=1, size=r[z[s]], prob=p[z[s]]) + 1
  x = c(x, rep(z[s], d[s]))
  y = c(y, rnorm(d[s], mu[z[s]], 1/sqrt(tau[z[s]])))
}
x = x[1:T]
y = y[1:T]

# plot
plot(y, type="l", lty="dotted")
lines(mu[x], col="blue")

# save
if(save_result) {
  saveRDS(
    list(y=y, s=s, x=x, z=z, d=d, t=t, tau=tau, mu=mu, p=p, r=r, bet=bet, A=A, Abar=Abar, Aini=Aini),
    file=paste("./simdata_rseed-", rseed, ".rds", sep=""))
}