## Louis Forbes Wright

#using _-_- as a basis for the algorithm created here


# we start by creating some data which contains change points.

means <- c(0,3,4,5,1,4,0)/1.5

changePoints <- c(0,100,200,300,400,500,600,700)

X_data <- NULL
U_data <- NULL
Theta_data <- NULL


for (i in 1:length(changePoints)-1){
  U_data <- c(U_data, rnorm(n = (changePoints[i+1] - changePoints[i]), mean = 0, sd = 1))
} 

for (i in 1:length(changePoints)-1){
  Theta_data <- c(Theta_data, rnorm(n = (changePoints[i+1] - changePoints[i]), mean = means[i], sd = 0))
}
X_data <- U_data + Theta_data
par(mfrow = c(1,3))

plot(X_data, type = 'l', ylab = 'X_i', xlab = 'i')

plot(U_data, type = 'l', ylab = 'U_i', xlab = 'i')

plot(Theta_data, type = 'l', ylab = 'Theta_i', xlab = 'i')
