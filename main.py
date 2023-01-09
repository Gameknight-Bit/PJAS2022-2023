# Main Python file for PJAS Project #
# Dev start date: 12/20/2022 #
# Dev end date:              #

#Does image compression through multiple different methods and compares#

import numpy as np
import math
import matplotlib.pyplot as plt
from skimage.color import rgb2gray
import myMath #Custom math library >:)

#SDV Compression#

a = np.array([[3, 2, 1], #A == matrix of pixle values (b&w)
              [2, 1, 4]])
aTransposed = np.array(myMath.invertMatrix(a))

ata = np.dot(aTransposed, a)

eigenvalues,vt = np.linalg.eig(ata) #Gets Eigenvalues and Vt

print(vt)

# Get valid non-zero eigenvalues for sigma matrix
sig = []
eigenvalues = np.sort(eigenvalues)[::-1] #sorts in descending order
for val in eigenvalues:
    for i in range(len(sig)):
        sig[i].append(0)

    if val > 0.0000001: #keep eigenval if not 0
        #add sqrt of value to sig_arr#
        arrToAppend = []
        for _ in range(len(sig)):
            arrToAppend.append(0)
        arrToAppend.append(math.sqrt(val))
        sig.append(arrToAppend)
    if val < -0.00001:
        print("Found negative eigenval????: "+str(val))    
    
sigma = np.array(sig)

print(sigma)

colVectors = []

for i in range(len(sigma)):
    vtThing = []
    for values in vt:
        vtThing.append([values[i]])
    #vtThing[len(vtThing)-1][0] *= -1
    print(vtThing)

    matrix = 1/math.sqrt(sigma[i][i])*np.dot(a, np.array(vtThing))
    print(matrix)