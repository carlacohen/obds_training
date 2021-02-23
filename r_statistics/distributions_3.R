###slide 13
###Compute an Empirical Cumulative Distribution Function###

#This is useful if you don't know if your data are normally distributed
#The mean and SD may not make sense e.g. if they are binomially distributed

#Use the ecdf() function to compute the empirical cumulative distribution function 
#for the variable Sepal.Length in the iris data set.
#iris data is builr into R so we can just call it immediately

head(iris)
iris_ecdf <-  ecdf(iris$Sepal.Length)
str(iris_ecdf)
#Tells us that ecdf is a function

#Use the plot() function to visualise the empirical cumulative distribution function.
plot(iris_ecdf)
#The dot is the sepal length, the line is the gap until the next sepal length


#Use the knots() function on the ecdf output and compare this with the list of 
#unique values for the variable Sepal.Length
knots(iris_ecdf)
#this tells us the unique values of sepal length
#you can also do
sort(unique(iris$Sepal.Length))
#and here I have sorted it as well

#What is the probability of any value up to 6?
iris_ecdf(6)
