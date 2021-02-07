# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 13:40:55 2021
Find the complementary DNA strand and reverse it.
@author: ccohen
"""

#Write a function to return the complement of a DNA sequence

def complementarynucleotide(nucleotide):
    #input is nucleotide A, T, C or G
    #output is the complement of A, T, C or G
    output = None
    if nucleotide == "A":
        output = "T"
    elif nucleotide == "T":
        output = "A"
    elif nucleotide == "C":
        output = "G"
    elif nucleotide == "G":
        output = "C"
    else: 
        output = "*"
    #print(output)
    return(output)
    
#Run on single nucleotide
#complementarynucleotide("A")

#Run on a sequence
sequence = "ATGNCCT"

# My first go that doesn't completely work
#for nucleotide in sequence:
 #   complementarynucleotide(nucleotide)

#Niamh's function
#Make a function to get the complementary strand
def complementaryDNAstrand(strand):
   #print(strand)
   #Tell it that it will need to use a string called complementstrand
   complementstrand = ""
   #for every letter in the string, do the complementarynucleotide function defined above
   for n in strand:
      x = complementarynucleotide(n)
      #print(x)
      #then move on to the next letter in the string
      complementstrand = complementstrand + x #or could use complementstrand += x I think
   print(complementstrand) #finally print the string once all the loop is finished
   return(complementstrand)#the end of the function and can be used to define a variable

#Run the complementaryDNAstrand function on the sequence
complement_seq = complementaryDNAstrand(sequence)

#Makea function to reverse the string
def reversecomplement(string):
     return string[::-1]

#Run the function on the output of the previous function, and print it
print(reversecomplement(complement_seq))
