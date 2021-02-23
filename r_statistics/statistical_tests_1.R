###Slide 31###
###Statistical tests###
#Iris data set

#Use the summary() function to view some information about each column.
#Since "Species" is a character vector we just get a count of each species
summary(iris)

#Visualise the distribution of Sepal.Length , stratified by species
hist(iris$Sepal.Length)

#This way makes 3 histograms next to each other
par(mfrow=c(2,2))
hist(iris$Sepal.Length[iris$Species == "setosa"], breaks = 20, col = "red")
hist(iris$Sepal.Length[iris$Species == "versicolor"], breaks = 20, col = "blue")
hist(iris$Sepal.Length[iris$Species == "virginica"], breaks = 20, col = "green")
par(mfrow=c(1,1))

#This way makes a plot and then uses lines to add the additional data
plot.new() #Create a new empty plotting window
#Calculate the range of Sepal.Length
range(iris$Sepal.Length)
#Set appropriate size of the axis
plot.window(xlim = c(4, 8), ylim = c(0, 1.4))
lines(density(iris$Sepal.Length[iris$Species == "setosa"]), col = "red")
lines(density(iris$Sepal.Length[iris$Species == "versicolor"]), col = "blue")
lines(density(iris$Sepal.Length[iris$Species == "virginica"]), col = "green")
#draw the axes and specify the tick marks
axis(side = 1, at = seq(4, 8))
axis(side = 2, at = seq(0, 1.4, 0.2))

#Is Sepal.Length length normally distributed? Overall? Within each species
#First get an idea of the data by plotting
plot(density(iris$Sepal.Length))
#Perform the shapiro test
shapiro.test(iris$Sepal.Length)
#p<0.05
#reject the null hypothesis that the data are normally distributed

shapiro.test(iris$Sepal.Length[iris$Species == "setosa"])
shapiro.test(iris$Sepal.Length[iris$Species == "versicolor"])
shapiro.test(iris$Sepal.Length[iris$Species == "virginica"])
#for all these, p>0.05 so we accept the hypothesis that the data are normally distributed

#Is there a significant variation of Sepal.Length between the various species?
#We can use ANOVA since the three datasets are normally distributed
#Kruskal_wallace test first (to show why not to do it??)

#iris data is already a df
iris$Sepal.Length
#explain the sepal length for each species
#SO in the first part you specify the columns, then in the second part you put what is the df (i.e. iris)
anova_iris <- aov(Sepal.Length ~ Species, data = iris)
summary(anova_iris)
#This tells us that p<2e-16
#It doesn't tell us where the difference is, so we may also want to peform a post-hoc analysis
#Run t-tests between each pair of species
t.test((iris$Sepal.Length[iris$Species == "setosa"]), (iris$Sepal.Length[iris$Species == "versicolor"]))
t.test((iris$Sepal.Length[iris$Species == "setosa"]), (iris$Sepal.Length[iris$Species == "virginica"]))
t.test((iris$Sepal.Length[iris$Species == "versicolor"]), (iris$Sepal.Length[iris$Species == "virginica"]))
#We could do this using a Tukey test as well
tukey_iris <- TukeyHSD(anova_iris)
View(tukey_iris$Species)
#Tukey output is a list, but we can view just the Species column to make it a better view
#These show that all comparisons are significant


