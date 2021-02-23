###This was used by Alina in testing.R but we ended up doing it a different way

#transform the logcounts
logcounts_transform <- t(logcounts)

#Check if the row names are identical in both tables (sanity check)
all(row.names(logcounts_transform) == row.names(cell_metadata))
#If not true then you can tell it to match the row names up.

cell_names <- row.names(logcounts_transform) #make a reference of row names
cell_metadata <- cell_metadata[cell_names, ] #tell the second table to use that order of names

#now use cbind to add a new column from cell_metadata to the transformed matrix
#Also rename the new column "Infection" otherwise it would be called "cell_metatdata"
combined_data <- cbind(Infection = cell_metadata, as.data.frame(logcounts_transform))

#Find a gene name
colnames(logcounts_transform)
#Run a t test on one gene
t.test(ENSG00000131203 ~ Infection, data = combined_data)
