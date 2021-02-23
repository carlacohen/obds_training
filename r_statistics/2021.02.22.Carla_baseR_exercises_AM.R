###22 Feb 2021 morning talk on Base R###
x = 4 
print (x)
y <- 4
print (y)
typeof(x)
typeof(joined_matrix)
#Create vector
num_1 <- c(1, 4, 7, 3, 5)
num_1 * 5
num_2 <- c(4, 6, 6, 12, 13)
num_2
num_1+ num_2

#Create a vector
animals <- c("dog", "cat", "rabbit")
#Assign attributes
names(animals) <- c("Harvey", "Max", "Jessica")
names(animals)
#Alternatively
attr(animals, "names") <- c("Harvey", "Max", "Jessica")
animals <- c(Harvey = "dog", Max = "cat", Jessica = "rabbit")

#Subsetting vectors (NB position is 1-based)
animals[2]
animals [1:2]
animals["Jessica"]

#Create a factor (represent categorical data)
my_colours <- factor(c("blue", "yellow", "red", "blue", "green"))
levels(my_colours)
min(my_colours)
#Create an ordered factor
quality <- factor(c("low", "high", "medium", "low", "high"), 
                  levels = c("low", "medium", "high"),
                  ordered = TRUE)
levels(quality)
min(quality)
#Add another level
levels(quality) <- c("low", "average", "medium", "high")

#Lists
animals <- c("dog", "cat", "rabbit")
bool <- c(TRUE, FALSE, TRUE, TRUE)
num <-c(6, 3, 2, 18)
my_list <- list(animals, bool, num)
View(my_list)
my_list[[2]] #access 2nd item in list
my_list[[2]][[1]] #access 1st element of 2nd item in list
my_list[2:3] #access 2nd and 3rd item in list

###Exercises###
#1. Generate a character vector of length 5. 
#Assign names to the elements
foods <- c("apple", "banana", "carrot", "date", "egg")
typeof(foods)
names(foods) <- c("Andy", "bob", "Charlie", "Dave", "Eric")
foods


#2. Access elements 1 and 3 by name. 
#Access the last element of the vector. 
#Replace element 4 with a new element.

foods["Andy"]
foods ["Charlie"]
foods [5]
foods [4] <- "dog"
foods

#3. num <- c(6, 3, 2, 18), bool <- c(TRUE, FALSE, TRUE, TRUE)
#Predict the result of sum(num[bool]) and num * bool
#Check whether you are correct
sum(num[bool])
num * bool

#4. Predict the result of adding num to num_2 <- c(2, 4, 12)
#Check whether you are correct
num_2 <- c(2, 4, 12)
num + num_2
#It doesn't work because they are not the same length

#5. missing <-  c(9, 4, NA, 2, 43, NA, 4, 2)
#How would you count the number of NA values in missing?
missing <-  c(9, 4, NA, 2, 43, NA, 4, 2)
bool_missing <- is.na(missing)
sum(bool_missing)
#or
sum(is.na(missing))

#6. Generate an ordered factor of length 6
my_factor <- factor(c("one", "two", "three", "one", "one", "two"), 
                    levels = c("one", "two", "three"),
                    ordered = TRUE)
View(my_factor)

#Check the levels of the factor
levels(my_factor)
#Replace element 2 with a value not listed in the levels
#Note the warning
my_factor[2] <- "five"
my_factor
#Now there is an NA in position 2

#7. Make a list containing a character vector, a boolean vector, 
#a factor and a numerical vector

#character vector
foods
#boolean vector
my_boolean <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#factor
my_factor 
#numerical vector
num
my_list <- list(foods, my_boolean, my_factor, num)
View(my_list)

#Use [[]] or [] to access:
#1. the 2nd and 3rd items in the list
my_list[2:3]
#2. the second element of the third item in the list
my_list[[3]][[2]]

#8. Add names to the list elements
my_list_names <- c("foods", "boolean", "factor", "numerical")
names(my_list) <- my_list_names
View(my_list)

#Access list elements using the $notation
my_list$boolean
my_list$boolean[2]

#9. Change the order of the factor levels inside the list 
#add a new level to the factor
levels(my_factor) <- c("four", "three", "two", "one")
my_factor

###lapply/sapply###

num
sapply(num, function(x) x + (x-1))

#1. Create a numeric vector of length 10.
num_10 <- c(1:10)
num_10

#Write an lapply and sapply statement to square each element. 
#Compare the two outputs

sapply(num_10, function(x) x*x)
lapply(num_10, function(x) x*x)

#2. Generate a list of length 4 containing both numeric and logical vectors
num
num_2
num_1
num_3
bool
bool_2 <- c(FALSE, FALSE, TRUE, FALSE)
num_3 <- c(5, 6, 7, 8)

my_list_2 <- list(num, num_3, bool, bool_2)
# Write an lapply or sapply statement to calculate the sum of
# the elements in each vector
sapply(my_list_2, sum)
lapply(my_list_2, sum)

#Write an sapply statement to repeat each element of each vector in
#your list three times
#Assign the output to a new list.
my_list_3 <- sapply(my_list_2, function(x) rep(x, each =3))
my_list_3

#4. Convert your new list into a single vector.
my_vector <- as.vector(my_list_3)



