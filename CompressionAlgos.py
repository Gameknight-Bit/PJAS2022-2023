#Custom Module that holds all the img compression algorithms >:)#

import numpy as np
import math

#Custom Libraries
import myMath #Custom math library >:)
import Picture

DATATYPE = np.float64

def SDV(a, kVal: int):
    maxKVal = max(kVal, 1) #has to be greater than 0!!! (Compression amount)
    #print(a.shape)

    a = a.astype("float64")

    #a = np.array([ ##### Testing Stuff #####
    #    [3, 2],
    #    [2, 1],
    #    [3, 3]
    #], dtype=DATATYPE)#####################################################

    n = a.shape[1]
    m = a.shape[0]

    aTransposed = np.array(myMath.invertMatrix(a), dtype=DATATYPE)

    ata = np.dot(aTransposed, a)

    #print(ata)

    #TODO: Replace with awesome eigenvector finding algorithm in myMath.py >:)
    eigenvalues, vt = np.linalg.eigh(ata) #Gets Eigen values and Vt

    #print(vt.shape)

    #EigenSort (Stupid sorting stuff >:() 
    #Set EigenRelations
    eigenRelations = []
    for i in range(len(eigenvalues)):
        eigenRelations.append([i, eigenvalues[i]])

    eigenvalues = np.sort(eigenvalues)[::-1] #sorts in descending order
    newvt = [ [0]*len(vt) for _ in range(len(vt)) ]
    #print("eigenRelations: "+str(len(eigenRelations)))
    #print("eigenValueNum: "+str(len(eigenvalues)))

    for i in range(len(eigenvalues)):
        for o in range(len(eigenvalues)):
            if eigenvalues[i] == eigenRelations[o][1]:
                
                #vectorTemp = [] #Contains values of a col in vt
                for vt_length in range(len(vt)):
                #    vectorTemp.append(vt[vt_length][eigenRelations[o][0]])
                    newvt[vt_length][i] = vt[vt_length][eigenRelations[o][0]]
                #newvt.append(vectorTemp)
    
    newvt = np.array(myMath.invertMatrix(newvt), dtype=DATATYPE)
    #for row in range(len(newvt)):
    #    if np.all(newvt[row] < 0):
    #        for i in range(len(newvt[row])):
    #            newvt[row][i] *= -1

    vt = np.array(newvt, dtype=DATATYPE)
    #print(vt.shape)
    #print(vt)

    # Get valid non-zero eigenvalues for sigma matrix
    sig = []

    #print(len(eigenvalues))
    for val in eigenvalues:
        if val > 0.0000000000001: #keep eigenval if not 0
        #    #add sqrt of value to sig_arr#
            if len(sig) < min(a.shape):
                arrToAppend = []                   #
                for _ in range(len(sig)):          #
                    arrToAppend.append(0)          #
                arrToAppend.append(math.sqrt(val)) #
                sig.append(arrToAppend)            #
        #print(val)
        if val < -0.0001:
            print("Found negative eigenval????: "+str(val))    

        for i in range(len(sig)): #Add zeros to make matrix orthonagonal
            sig[i].append(0)

    #Now sigma is orthonormal, however, it needs to have sides that correspond with each U and VT

    sigma = np.array(sig, dtype=DATATYPE)
    sigma = np.delete(sigma, len(sigma[0])-1, 1) #delete last column (Was 0)

    #print(sigma)

    #Calculating U
    colVectors = []

    for i in range(len(sigma)):
        vtThing = vt[i]

        matrix = (1/sigma[i][i])*np.dot(a, np.array(vtThing, dtype=DATATYPE))
        for i in range(len(matrix)):
            if len(colVectors) <= i:
                colVectors.append([])
            colVectors[i].append(matrix[i])

    U = np.array(colVectors, dtype=DATATYPE)
    #print(U)

    #print(vt.shape)
    #print(U.shape)
    #print(sigma.shape)
    #print("----")

    #print(np.rint(np.dot(np.dot(U, sigma), vt)).astype("uint8")) #Successfully separated >:)

    ####### Next part is the actual compression!!! #######
    #Divide U(columns), Sigma(Singular Values), and VT(rows) into groups
    #Only take certain number of groups (min compression num 'k')
    #  Sum all multiples for USVt less than min compress
    #  There is the new compressed form >:)

    #minSig = min(sigma.shape)

    #grab Cols from vt
    #UFinal = [] #np.matrix(U[:, :kVal])#
    #UVal = maxKVal
    #if len(U) <= maxKVal:
    #    UVal = len(U)
    #for i in range(len(U)):
    #    UFinal.append(U[i][:UVal])
    #UFinal = np.array(UFinal)
    #print(uFinal)

    #Grab Rows from Vt
    #vtFinal = [] #np.matrix(vt[:kVal, :])#
    #vtVal = kVal
    #if len(vt) <= maxKVal:
    #    vtVal = len(vt)
    #for i in range(len(vt)):
    #    if i >= vtVal:
    #        break
    #    vtFinal.append(vt[i])
    #vtFinal = np.array(vtFinal)
    #print(vtFinal)

    #Grab sigma stuff >:)
    #sigFinal = [] #np.matrix(sigma[:kVal,:kVal])#
    #for i in range(minSig):
    #    if i >= maxKVal:
    #        break
    #    sigFinal.append(sigma[i][:maxKVal])
    #sigFinal = np.array(sigFinal)
    #print(sigFinal)

    ##################Testing out new method ####################

    Afinal = np.array([ [0]*len(list(vt)) for _ in range(len(U)) ], dtype=DATATYPE)

    for kindex in range(min(kVal, min(a.shape))):

        #Get sigma val
        sigVal = sigma[kindex][kindex]

        #Get U col
        uCol = U[:, kindex]
        newCol = []
        for val in uCol:
            newCol.append([val])
        uCol = np.array(newCol, dtype=DATATYPE)
        #print(uCol)

        #Get vt row
        vtRow = np.array([vt[kindex, :]], dtype=DATATYPE) # ':' is reduntant :)

        #if Afinal == None:
        #    Afinal = sigVal*np.dot(uCol, vtRow)
        #else:
        prod = np.dot(uCol, vtRow)
        #print(prod.shape)
        #print(Afinal.shape)
        Afinal = Afinal + sigVal*prod

    #print(vtFinal.shape)
    #print(UFinal.shape)
    #print(sigFinal.shape)
    #print("----")

    Afinal = np.rint(Afinal)#.astype("uint8") #np.rint(np.dot(np.dot(U, sigma), vt)).astype("uint8")

    for x in range(len(Afinal)): #Reformatting image
        for y in range(len(Afinal[x])):
            if Afinal[x][y] < 0:
                Afinal[x][y] = 0
            if Afinal[x][y] > 255:
                Afinal[x][y] = 255

    return Afinal.astype('uint8')