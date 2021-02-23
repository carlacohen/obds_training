###Distributions###

#Generate a vector of 1000 normally distributed values with mean 10 and standard deviation 5
norm_dist <- rnorm(1000, mean = 10, sd = 5)

?rnorm #access the help

#Inspect the output of the summary() function for that vector
summary(norm_dist)

#Compute the mean and standard deviation for those values
mean(norm_dist)
sd(norm_dist)

#Compute the deciles (i.e. 10 evenly spaced quantiles) for those values
quantile(norm_dist, probs = seq(0, 1, 0.1))

#Visualise the distribution of those values as a histogram
hist(norm_dist, breaks = 50)
#breaks tells us the number of bins
#later we will learn ggplot do to this

#Visualise as vertical lines on the histogram: 
#Don't need to re-generate the histogram, it will add on top of the current plot
#the mean (red solid), 
abline(v = mean(norm_dist), col = "red")
#median (red dashed), 
abline(v = median(norm_dist), col = "red", lty = 2)
#one standard deviation from the mean (blue solid), 
abline(v = mean(norm_dist) + sd(norm_dist), col = "blue")
abline(v = mean(norm_dist) - sd(norm_dist), col = "blue")
#one median absolute deviation from the median (blue dashed).
abline(v = median(norm_dist) + mad(norm_dist), col = "blue", lty = 2)
abline(v = median(norm_dist) - mad(norm_dist), col = "blue", lty = 2)
#This is the median of the difference between each value and the median
#good for finding outliers

#Generate a new vector with a lot more values (e.g., one million). 
#Draw again a histogram. 
#How does the distribution compare with more data points?
norm_dist_2 <- rnorm(1000000, mean = 10, sd = 5)
hist(norm_dist_2,  breaks =50)

#Visualise both histograms together in a plotting grid
par(mfrow = c(2, 1))
hist(norm_dist, breaks = 50)
hist(norm_dist_2,  breaks =50)
par(mfrow = c(1, 1))#set it back to original




