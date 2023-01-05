# Lesson 7
[Class Slides](Slides7.pdf)

## Assignment
Your goal is to implement the Needleman-Wunsch algorithm. You can read about the Needleman-Wunsch algorithm on [Wikipedia](https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm). The Wikipedia page holds psuedo-code which you might find helpful.

### A
Write a function that takes two sequences as input and returns a matrix of scores as we saw in class. No need for backtracing here.</br>
You can use this code skeleton written in python to get you started
```python
# Use these values to calculate scores
gap_penalty = -1
match_award = 1
mismatch_penalty = -1

# Make a score matrix with these two sequences
seq1 = "ATTACA"
seq2 = "ATGCT"

# A function for making a matrix of zeroes
def zeros(rows, cols):
    # Define an empty list
    retval = []
    # Set up the rows of the matrix
    for x in range(rows):
        # For each row, add an empty list
        retval.append([])
        # Set up the columns in each row
        for y in range(cols):
            # Add a zero to each column in each row
            retval[-1].append(0)
    # Return the matrix of zeros
    return retval

# A function for determining the score between any two bases in alignment
def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty

# The function that actually fills out a matrix of scores
def needleman_wunsch(seq1, seq2):
    
    # length of two sequences
    n = len(seq1) 
    m = len(seq2)
    
    # Generate matrix of zeros to store scores
    score = zeros(m+1, n+1)
    
    ########################
    # Your code starts here
    ########################
    
    # Use the following steps as a guide to calculate the score matrix
    
    # 1. Fill out first column
    
    # 2. Fill out first row
    
    # 3. Fill out all other values in the score matrix

    # Return the final score matrix
    return score

# Test out the needleman_wunsch() function
print_matrix(needleman_wunsch(seq1, seq2))
```
### B
Do the backtracing through the score matrix and print out the final alignment</br></br></br>


## Assignment - BWT
Read about the BWT and inverse BWT on [Wikipedia](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform). 

**A** Write functions that perform the BWT and inverse BWT on a given string. The BWT function should take a string as input and return the transformed string as output, while the inverse BWT function should take a transformed string as input and return the original string as output. Include the output of *bwt('hello world')* </br>
For each of your functions, explain the time and space complexity of the algorithm?

**B** How do suffix arrays are used in the BWT? How they are constructed and used in the BWT algorithm? What are the benefits and drawbacks of using suffix arrays in the BWT, and how they can be used to improve the performance of the algorithm. </br></br></br>


You are free to use whatever programing language you feel comfortable in. The only requirement is to provide a detailed explanation on how to run your code on a 64 bit Linux.

