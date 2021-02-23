#Slide 35
###Testing & Multiple testing correction###
#For each gene (i.e. row) in logcounts.csv , use cell_metadata.csv and a statistical test of your
#choice to identify gene differentially expressed in cells infected with Salmonella relative to the
#control uninfected cells

#Suggestion: write the code to test one gene, refactor the code into a function that returns
#the p-value, and use vapply to apply that function to all the genes

logcounts <- read.csv("data/logcounts.csv", row.names = 1)
logcounts = as.matrix(logcounts)
#using row.names = 1 assignes the first column (in this case the gene name) as the row name
#This is single cell data where each row is a gene and each column is a cell
#NB this is the other way round from Python!
#Also turn this into a matrix as it will be a df automatically

cell_metadata <- read.csv("data/cell_metadata.csv", row.names = 1)
#This tells us for each cell, whether it is mock or infected with salmonella (STM-LT2)

#Make a function in order to be able to run the wilcox test on every gene
#name the function "test_row" and tell it that you need to input an index (i.e. row) and matrix (i.e. data)
test_row <- function(index, matrix) {
  #create a temporary df consisting of a value and a group
  test_data <- data.frame(
    #this gives you the value of the gene expression for the specified row (using index to say which row)
    value = as.numeric(matrix[index, ]),
    #the group is mock or infected
    group = cell_metadata$Infection
  )
  #perform a wilcox test on the df you have made
  #we don't know if the data are normally distributed hence using wilcox.test not t.test
  out <- wilcox.test(value ~ group, test_data)
  #output of the wilcox test is a list, we want to extract the p value only
  out$p.value
  }

#run the function on the first row (i.e. one gene only)
test_row(1, logcounts)
#so for the first row, we get a p value of 6.72e-10

#now use vapply to run this function over the whole matrix
#vapply(X, FUN, FUN.VALUE, ..., USE.NAMES = TRUE)
#X is the row index
#FUN is the function
#seq(1, nrow(logcounts)) means make a sequence from 1 to the number of rows that exist in logcounts
# numeric(1) creates a placeholder for a vector of the specified length
#...means add additional parameters that are required by your function, whihc in this case means specify the matrix
#the first parameter is automatically passed to the first input of your function so we don' need to specify index in this case
rowselect <- seq(1, nrow(logcounts))
p_values <- vapply(rowselect, test_row, FUN.VALUE = numeric(1), matrix = logcounts)

#Visualise a bar plot of the p-values.
hist(p_values)

#Correct p-values for multiple testing. 
#We could use Bonferoni, BH or FDR (depends on the expt)
corrected_p_values <- p.adjust(p_values, method = "bonferroni")
hist(corrected_p_values)

plot(p_values, corrected_p_values)
abline(a = 0, b = 1)
abline(a = 0.05, b =0)
abline(v=0.05)
corrected_p_values < 0.05


#How many genes remain before and after multiple testing?

#Need to save code from this morning to shared repo 
