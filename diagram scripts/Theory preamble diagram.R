## preamble diagrams

par(mfrow = c(1,3))

data <- rep(c(4,-4), each = 100) + rnorm(200, 0.5)

mean1 <- c(mean(data[1:50]), mean(data[51:200]))
mean2 <- c(mean(data[1:100]), mean(data[101:200]))
mean3 <- c(mean(data[1:150]), mean(data[151:200]))

y1 <- rep(mean1, c(50, 150))
y2 <- rep(mean2, each = 100)
y3 <- rep(mean3, c(150, 51))

plot(data, type = 'l' , ylab= TeX("$X_i$"), xlab = TeX("$i$"))
lines(y1, col = "red")
plot(data, type = 'l' , ylab= TeX("$X_i$"), xlab = TeX("$i$"))
lines(y2, col = "red")
plot(data, type = 'l' , ylab= TeX("$X_i$"), xlab = TeX("$i$"))
lines(y3, col = "red")