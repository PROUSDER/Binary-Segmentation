## diagram script for the partition of the interval [k_1, k_2] by changepoints.

library(latex2exp)
library(calibrate)

line.x <- c(10, 50)*10
line.y <- c(1,1)

plot(line.x, line.y, type = 'l', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", main=TeX('Changepoints in between $\\lbrack \\k_1, \\k_2 \\rbrack$'))

interval.points <- c(15, 45)*10

points(interval.points[1], 1, pch = '[')
points(interval.points[2], 1, pch = ']')

changepoints.x <- c(13.5, 17, 21, 27, 33, 37, 43, 47)*10

n <- length(changepoints.x)
changepoints.y <- rep(1, n)

points(changepoints.x, changepoints.y)

i <- 1:(n-1)
changepoints.labels<- TeX(sprintf("$\\nu_{i_0 + %d}$", i))

textxy(c(interval.points,changepoints.x), c(c(1,1),changepoints.y), c(TeX("$k_1$"), TeX("$k_2$"), TeX("$\\nu_{i_0}$"),changepoints.labels), cex = 0.75)
# use the above line for the setup.

#textxy(c(interval.points), c(1,1), c(TeX("$k_1$"), TeX("$k_2$")), cex = 0.75)
# use the above line for the distance conditions for changepoints

enclosure.points.x_red <- c(changepoints.x[2:(n-1)] -1.5*10, changepoints.x[2:(n-1)] + 1.5*10)
enclosure.points.x_blue <- c(changepoints.x[2:n] -1.5*2*10, changepoints.x[1:(n-1)] + 1.5*2*10)


points(enclosure.points.x_red[c(1,(n-2))], changepoints.y[1:2], pch = "[", col = "red")
points(enclosure.points.x_red[c(n-1, (2*(n-2)))], changepoints.y[1:2], pch = "]", col = "red")

points(enclosure.points.x_blue[1:n-1], rep(1, length(enclosure.points.x_blue)/2), pch = "[", col = "blue")
points(enclosure.points.x_blue[n:(2*(n-1))], rep(1, length(enclosure.points.x_blue)/2), pch = "]", col = "blue")