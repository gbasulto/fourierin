## ------ Example 1 -----------------------------------------------------------
## -- Example 1: Recovering std. normal from its characteristic function ------

out <- fourierin(f = function(t) exp(-t^2/2),
                 a = -5, b = 5, c = -3, d = 3,
                 r = -1, s = -1, resol = 128)
grid <- out$w
values <- Re(out$values)

plot(grid, values, type = "l", col = 3)
lines(grid, dnorm(grid), col = 4)
