#From slide 54 onwards
###Matrix exercises###
#in R you can use a ? to get help e.g. ?matrix ?seq

#Generate a matrix with 5 rows containing the numbers 2:100 in increments of 2, fill by row
seq(from = 2, to = 100, by = 2)
exercise_1_matrix <- matrix(seq(from = 2, to = 100, by = 2), nrow = 5, byrow = TRUE)
View(exercise_1_matrix)

#Calculate the sum of squares of each row
apply(exercise_1_matrix, 1, function(x) sum (x*x))
#means apply this function to every row in this matrix

#Calculate the min and max of each column

apply(exercise_1_matrix, 2, min)
apply(exercise_1_matrix, 2, max)

#In the same line of code
#Returns a matrix where min is first row and max is bottom row
apply(exercise_1_matrix, 2, function(x) c(min(x), max(x)))

#Generate a second matrix containing 10 rows and 6 columns
#Calculate the transpose of this matrix
#Join the transposed matrix to the matrix from #1 (join by row)
exercise_3_matrix  <- matrix (seq(from=1, to = 60), nrow = 10, ncol = 6, byrow = TRUE)
View(exercise_3_matrix)
dim (exercise_3_matrix)    #tells you the dimensions of your matrix

#Checking to see what happens if you don't specifcy rows and columns
exercise_3_matrix_1 <- matrix(seq(from = 1, to = 70), byrow = TRUE)
dim (exercise_3_matrix_1)

transposed_exercise_3_matrix <- t(exercise_3_matrix)
View (transposed_exercise_3_matrix)
dim (transposed_exercise_3_matrix) #6 rows and 10 columns
dim (exercise_1_matrix) #5 rows and 10 colunns

#Join the transposed matrix to the matrix from #1 (join by row)
joined_matrix <- rbind (exercise_1_matrix, transposed_exercise_3_matrix)
dim(joined_matrix)
class(joined_matrix)
typeof(joined_matrix)

#Convert your joined matrix into a data frame
joined_df <- as.data.frame((joined_matrix))
View(joined_df)
class(joined_df) #tells us what the type of data is
typeof (joined_df) #tells you the internal structure of the object

#A df can have multiple data types whereas a matrix can only have one

#Convert the data frame to a list
joined_list <- as.list(joined_df)
View(joined_list)
class(joined_list)
typeof(joined_list)

###Data frame exercises###
#1. Load the coding_gene_region.bed file into R 
#(in /ifs/obds-training/jan21/shared/r/baseR).
#Copy file to local drive using Filezilla to Documents/R/OBDS
setwd("~/R/OBDS")
bed_file <- read.table("coding_gene_region.bed", 
                       header = FALSE, 
                       sep = "\t")
#don't actually need to specify the tab
#Use stringsAsFactors = FALSE,if you don't want character vectors to be converted to factors

#Check the dimensions of the data frame and the class of each variable
dim(bed_file)
class(bed_file)

#Add column names
colnames(bed_file) <- c("chr", "start", "end", "name", "score", "strand")

#Add a new column containing the length of each genomic interval 
#sort this column from largest to smallest using a base R function
bed_file$length = bed_file$end - bed_file$start
#NB Could also make a new object and use the order function e.g.
sorted_table <- bed_file[order(bed_file$length, decreasing = TRUE), ]
head(sorted_table)

#Extract the element at row 30, column 3
bed_file[30,3]
sorted_table[30,3]

#Extract the second column by index and by name
bed_file[,2]
bed_file$start

#On which chromosome is the largest interval?
#Output just the chromosome value and store in the variable max_chrom
max_length <- max(bed_file$length) #This finds the max length
max_length
max_columns <- apply(bed_file, 2, max) #so does this
max_columns

#Since you have already ordered the table, you can in fact just do:
max_chrom <- bed_file[1,1]
max_chrom

#Alternatively you can do
max_chr <- bed_file$chr[bed_file$length == max(bed_file$length)]
max_chr

#Subset the data frame to contain only regions with a length from
#100,001-200,000 bp - assign to a new variable. 

data_subset <- subset(bed_file, bed_file$length >= 100001 & bed_file$length <= 200000)
data_subset

#Could also do using square brackets, but this is more computationally intensive
data_subset_2 <- bed_file[bed_file$length %in% c(100001:200000),]
#%in% you give it a vector and it will give you the length within that range true or false
head(data_subset_2)
bed_file$length %in% c(100001:200000) #gives you a list of TRUE FALSE etc

#Or else
data_subset_3 <- bed_file[bed_file$length >= 100001 & bed_file$length <= 20000, ]
data_subset_3


#Write your subset data frame to a tab separated file 
#include column names but not row names
write.table(data_subset, "data_subset.txt", 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

#In the original data frame, replace the score value with 100 for
#genomic intervals on chr4 or chr17 that are on the + strand and longer than 200,000 bp. 

bed_file$score[bed_file$chr %in% c("chr4", "chr17") 
               & bed_file$strand == "+" 
               & bed_file$length > 200000] <- 100
View(bed_file)

#Count the number of regions that have a score of 100.
sum(bed_file$score == 100)

#Add a new row to the original data frame - you can make up the values. 
bed_file_2 <- rbind(bed_file, list("chrM", 10321321, 65165132, "test_gene", 0, "+", 12354354))
View(bed_file_2)

#Make sure the class of each variable in the data frame is correct.
# This is done by including "" or not, according to whether it is a number or a character

#Remove the score variable from the data frame


