library(tidyverse)

data <- read.csv(file = "C:/Users/User/Documents/programming/Stats R work/datasets/ukhpi-comparison-all-pmc-london-from-2006-01-01-to-2019-04-01.csv")
names(data)[7] <- "Percentage.change"

#From Jan 2010
data <- data[193:nrow(data),]

london_data <- data[data$Name == "London",]
london_data <- london_data[!is.na(london_data$Percentage.change),]

NEast_data <- data[data$Name == "North East",]
Neast_data <- Neast_data[!is.na(NEast_data$Percentage.change),]

Scotland_data <- data[data$Name == "Scotland",]
Scotland_data <- Scotland_data[!is.na(Scotland_data$Percentage.change),]

SEast_data <- data[data$Name == "South East",]
SEast_data <- SEast_data[!is.na(SEast_data$Percentage.change),]

avg <- (london_data$Percentage.change + NEast_data$Percentage.change
        + SEast_data$Percentage.change + Scotland_data$Percentage.change)/4


vec1 <- c(london_data$Percentage.change, SEast_data$Percentage.change, NEast_data$Percentage.change, Scotland_data$Percentage.change)

plot(avg, type = 'l', yaxp = c(min(vec1), max(vec1)), 10)



points(london_data$Percentage.change, pch = ".", cex = 2)
points(SEast_data$Percentage.change, pch = ".", col = "red", cex = 2)
points(NEast_data$Percentage.change, pch = ".", col = "green", cex = 2)
points(Scotland_data$Percentage.change, pch = ".", col = "blue", cex = 2)

vec1 <- c(london_data$Percentage.change, SEast_data$Percentage.change, NEast_data$Percentage.change, Scotland_data$Percentage.change)

ts.plot(Binary.Segmentation(london_data$Percentage.change)$est.signal, ylab = "")
#lines(Binary.Segmentation(SEast_data$Percentage.change)$est.signal, lty = 2)
lines(Binary.Segmentation(NEast_data$Percentage.change)$est.signal, lty = 2)
lines(Binary.Segmentation(Scotland_data$Percentage.change)$est.signal, lty =3)

abline(v = c(which(london_data$Period == "2016-06")), col = "red")


london_AR <- ar(london_data$Percentage.change, order.max = 12)
NEast_AR <- ar(NEast_data$Percentage.change, order.max = 12)
Scotland_AR <- ar(Scotland_data$Percentage.change, order.max = 12)

london_GARCH <- garch(london_AR$resid[!is.na(london_AR$resid)])
NEast_GARCH <- garch(NEast_AR$resid[!is.na(NEast_AR$resid)])
Scotland_GARCH <- garch(Scotland_AR$resid[!is.na(Scotland_AR$resid)])

x <- Binary.Segmentation(london_GARCH$residuals[!is.na(london_GARCH$residuals)])$est.signal
y <- Binary.Segmentation(NEast_GARCH$residuals[!is.na(NEast_GARCH$residuals)])$est.signal
z <- Binary.Segmentation(Scotland_GARCH$residuals[!is.na(Scotland_GARCH$residuals)])$est.signal

comb <- c(x, y, z)

ts.plot(NA, xlim = c(0, 105), ylim = c(1.2*min(comb),1.2*max(comb)), ylab = "")
lines(x)
lines(y, lty = 2)
lines(z, lty = 3)
abline(v = c(which(london_data$Period == "2016-06") - 13), col = "red")




