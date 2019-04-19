## Dice roll example

D1.distribution <- rep(1/6, 6)
D1 <- sample(1:6, 200, replace = TRUE ,prob = D1.distribution)

D2.distribution <- c(0.1, 0.1, 0.1, 0.1, 0.3, 0.3)
D2 <- sample(1:6, 200, replace = TRUE ,prob = D2.distribution)

dicerolls <- c(D1, D2)

signal <- rep(c(100,130), c(1300,3068)) + rnorm(4368, sd = 50)
plot(signal, type = 'l', xlab = "Hours after the start of the year", ylab = "units produced", xaxp = c(0,0,1))
axis(side = 1, at = c(0, 1000, 2000, 3000, 4000), labels = c(0, 2000, 4000, 6000, 8000))

s <- Binary.Segmentation(signal)
lines(s$est.signal, col = "red")

print(s[c(2,3)])
