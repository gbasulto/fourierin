## -- Example 1: Recovering std. normal from its characteristic function ------

out <- fourierin(f = function(t) exp(-t^2/2),
                 a = -5, b = 5, c = -3, d = 3,
                 r = -1, s = -1, resol = 64)
grid <- out$w
values <- Re(out$values)

plot(grid, values, type = "l", col = 3)
lines(grid, dnorm(grid), col = 4)

## -- Example 2: Computing characteristic function of a gamma r. v. -----------

out <- fourierin(f = function(t) dgamma(t, shape = 5, scale = 3),
                 a = -5, b = 5, c = -3, d = 3,
                 r = -1, s = -1, resol = 64)
grid <- out$w
values <- Re(out$values)

plot(grid, values, type = "l", col = 3)
lines(grid, dnorm(grid), col = 4)
