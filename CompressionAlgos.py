#Custom Module that holds all the img compression algorithms >:)#

import numpy as np
import math

#Custom Libraries
import myMath #Custom math library >:)
import Picture


def SDV(a, kVal: int):
    maxKVal = max(kVal, 1) #has to be greater than 0!!! (Compression amount)

    aTransposed = np.array(myMath.invertMatrix(a), dtype=np.float64)

    ata = np.dot(aTransposed, a)

    #print(ata)

    #TODO: Replace with awesome eigenvector finding algorithm in myMath.py >:)
    eigenvalues, vt = np.linalg.eig(ata) #Gets Eigenvalues and Vt

    #EigenSort (Stupid sorting stuff >:() 
    #Set EigenRelations
    eigenRelations = []
    for i in range(len(eigenvalues)):
        eigenRelations.append([i, eigenvalues[i]])

    eigenvalues = np.sort(eigenvalues)[::-1] #sorts in descending order
    newvt = []
    for i in range(len(eigenvalues)):
        for o in range(len(eigenRelations)):
            if eigenvalues[i] == eigenRelations[o][1]:
                
                vectorTemp = [] #Contains values of a col in vt
                for vt_length in range(len(vt)):
                    vectorTemp.append(vt[vt_length][eigenRelations[o][0]])
                newvt.append(vectorTemp)

    newvt = np.array(newvt, dtype=np.float64)
    vt = newvt

    #print(vt)

    # Get valid non-zero eigenvalues for sigma matrix
    sig = []
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

    sigma = np.array(sig, dtype=np.float64)

    #print(len(sigma))

    #Calculating U
    colVectors = []

    for i in range(len(sigma)):
        vtThing = vt[i]

        matrix = 1/sigma[i][i]*np.dot(a, np.array(vtThing, dtype=np.float64))
        for i in range(len(matrix)):
            if len(colVectors) <= i:
                colVectors.append([])
            colVectors[i].append(matrix[i])

    U = np.array(colVectors, dtype=np.float64)
    #print(U)

    #testA = np.rint(np.dot(np.dot(U, sigma), vt)).astype("uint8") #Successfully separated >:)

    ####### Next part is the actual compression!!! #######
    #Divide U(columns), Sigma(Singular Values), and VT(rows) into groups
    #Only take certain number of groups (min compression num 'k')
    #  Sum all multiples for USVt less than min compress
    #  There is the new compressed form >:)

    minSig = min(sigma.shape)

    #grab Cols from U
    uFinal = None
    if minSig <= maxKVal:
        uFinal = U
    else:
        uFinal = []
        for i in range(len(U)):
            uFinal.append(U[i][:maxKVal])
        uFinal = np.array(uFinal)

    #Grab Rows from Vt
    vtFinal = None
    if minSig <= maxKVal:
        vtFinal = vt
    else:
        vtFinal = []
        for i in range(maxKVal):
            vtFinal.append(vt[i])
        vtFinal = np.array(vtFinal)

    #Grab sigma stuff >:)
    sigFinal = []
    if maxKVal >= minSig:
        sigFinal = sigma
    else:
        for i in range(minSig):
            if i >= maxKVal:
                break
            sigFinal.append(sigma[i][:maxKVal])
        sigFinal = np.array(sigFinal)
        
    #print(vtFinal.shape)
    #print(uFinal.shape)
    #print(sigFinal.shape)
    #print("----")

    Afinal = np.rint(np.dot(np.dot(uFinal, sigFinal), vtFinal)).astype("uint8")

    return Afinal