# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:30:38 2021

@author: ccohen
"""
#Accept a number from the user and calculate the sum of all numbers between 1 and the given number


guess = int(input("Enter an integer : "))
y = 1
total = 0
while y < guess:
    #print (y)
    y = y + 1
    total = total + y
    #print(total)
print (total)

#Aneesha's solution
num = int(input("Enter a number : "))
sum = 0
for num in range(0, num, 1):
    sum = sum+num
    print ("Sum of first", num, "numbers is", sum)
    
#Niamh's solution
running = True


num = int(input("Enter a number : "))
y=1
for x in range (0, num):
    print (x)
    y = y + x
print (y)
print ("The sum of all the numbers up to your number")