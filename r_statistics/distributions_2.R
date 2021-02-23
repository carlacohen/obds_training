###Query distributions and probabilities### slide 12
#For the standard normal distribution mean = 0, sd = 1
#Plot the cumulative distribution function in the range -5, 5
q = seq(-5, 5, 0.1) # make a vector from -5 to 5 in increments of 1
vector_probs <- pnorm(q, mean = 0, sd = 1) #creates a vector of probabilities
plot(x = q, y = vector_probs)
#this shows you for any value of q, what is the probability of any value up to q occurring
#e.g. a t-test uses this info because you are asking what is the probability of q being between these two values
#e.g. > 2 SDs (the packages will do this for you)

#Plot the inverse cumulative distribution function for quantiles in 0.01 increment.
#p is a vector of probabilities (so must be between 0 and 1)
p = seq(0,1,0.01)
vector_values <- qnorm(p, mean = 0, sd = 1)
plot (x = p, y = vector_values)

#Plot the density function in the range -5, 5
vector_density <- dnorm(q, mean = 0, sd = 1)
plot(x = q, y = vector_density)
#less intuitive to understand!

#For more explanations see https://en.wikipedia.org/wiki/Normal_distribution

#What is the probability of observing a value greater than 2

prob <- 1 - pnorm(2, mean = 0, sd = 1)
prob 
#What is the probability of observing a value between -2 and 2?
prob_2 <- pnorm(2, mean = 0, sd = 1) - pnorm(-2, mean = 0, sd = 1)
prob_2

#What is the probability of observing a value more extreme than -2 or 2?
prob_3 <- 1 - (pnorm(2, mean = 0, sd = 1) - pnorm(-2, mean = 0, sd = 1))
prob_3

