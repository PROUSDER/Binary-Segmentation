# Real data example
# Data from https://catalog.data.gov/dataset/air-quality-measures-on-the-national-environmental-health-tracking-network/resource/8c3d19a7-d13e-446b-abb4-cc5a96c6d7a7

# aim: exploratory data analysis

myData <- read.csv(file="file:///C:/Users/User/Documents/programming/Stats R work/datasets/AMZN.csv",
                   header = TRUE, sep = ',')

head(data.frame(myData))
        