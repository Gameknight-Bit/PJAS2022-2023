#Custom Module that holds all the img compression algorithms >:)#

import numpy as np
import math
import matplotlib.pyplot as plt

#Custom Libraries
import myMath #Custom math library >:)
import Picture

np.set_printoptions(precision=10, suppress=True)

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

    ###EigenSort (Stupid sorting stuff >:( )###
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

    ### Get valid non-zero eigenvalues for sigma matrix ###
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

    ###Calculating U###
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

# - - - - - - - - - - Fourier Calculations - - - - - - - - - - - - - #

def DFTransform(row): #Calculates a DFT on the given array
    
    N = len(row)
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp((-2.0j * math.pi * k * n )/ N)

    retArr = np.dot(e, row)

    return retArr

def InvDFT(row): #Calculates an Inverse DFT on the given array

    N = len(row)
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp((2.0j * math.pi * k * n )/ N) #Removed negative for j

    retArr = np.dot(e, row)

    return retArr*(1/N) #Just added 1/N

def FFTransform(row): #Calculates a Fast Fourier transform on the Array
    #Using Cooley-Turkey Variation#
    N = len(row)
    
    if N % 2 > 0:
        return DFTransform(row)
    else:
        even = FFTransform(row[::2])
        odd = FFTransform(row[1::2])
        factor = np.exp(-2.0j * np.pi * np.arange(N) / N)
        return np.concatenate([even + factor[:N // 2] * odd,
                               even + factor[N // 2:] * odd])

def IFFTransform(row): #Calculates an Inverse Fast Fourier transform on the Array
    #Using Cooley-Turkey Variation#

    N = len(row)
    return (1/N)*np.conj(FFTransform(np.conj(row)))

def FT(a, TVal: int): 
    '''a = [[1, 20, 8],
        [2, 11, 13],
        [1, 100, 4],
        [0, 42, 55]]'''

    threshold = TVal

    #print(np.array(DFTransform(a)))

    x, y = a.shape

    ### 2d DFT ###
    #Run fourier transform against rows of a
    rowTrans = []
    for row in a:
        rowTrans.append(DFTransform(row))
    
    #Run Fourier transform against cols of a
    fourierSum = [] #End product
    for col in myMath.invertMatrix(rowTrans):
        fourierSum.append(DFTransform(col))
    fourierSum = np.array(myMath.invertMatrix(fourierSum))

    #print(fourierSum)

    magnitudeSpec = 20*np.log10(abs(np.fft.fftshift(fourierSum))) #For fourier analysis
    #return magnitudeSpec.astype('uint8')

    print(abs(fourierSum))
    threshold = 0.1*threshold * abs(fourierSum).max()
    print(threshold)

    indecies = (abs(fourierSum)>threshold)*1
    print(indecies)
    fourierSumFiltered = fourierSum*indecies

    count = x*y - sum(sum(indecies))
    percentLoss = 100-count/(x*y)*100
    print("Keeping "+str(percentLoss)+"% of the original image...")
    
    ### Inverse 2d DFT ###
    rowTrans = []
    for row in fourierSumFiltered: #Done on rows
        rowTrans.append(InvDFT(row))

    fourierEndSum = []
    for col in myMath.invertMatrix(rowTrans):
        fourierEndSum.append(InvDFT(col))
    compressedImgData = np.rint(np.array(myMath.invertMatrix(fourierEndSum)))

    for x in range(len(compressedImgData)): #Reformatting image
        for y in range(len(compressedImgData[x])):
            if compressedImgData[x][y] < 0:
                compressedImgData[x][y] = 0
            if compressedImgData[x][y] > 255:
                compressedImgData[x][y] = 255

    return compressedImgData.astype('uint8')

def FFT(a, TVal: int): 
    threshold = TVal

    x, y = a.shape

    ### 2d DFT ###
    #Run fourier transform against rows of a
    rowTrans = []
    for row in a:
        #print(len(row))
        rowTrans.append(FFTransform(row))
    
    #Run Fourier transform against cols of a
    fourierSum = [] #End product
    for col in myMath.invertMatrix(rowTrans):
        #print(len(col))
        fourierSum.append(FFTransform(col))
    fourierSum = np.array(myMath.invertMatrix(fourierSum))

    #print(fourierSum)

    magnitudeSpec = 20*np.log10(abs(np.fft.fftshift(fourierSum))) #For fourier analysis
    #return magnitudeSpec.astype('uint8')

    print(abs(fourierSum))
    threshold = 0.1*threshold * abs(fourierSum).max()
    print(threshold)

    indecies = (abs(fourierSum)>threshold)*1
    print(indecies)
    fourierSumFiltered = fourierSum*indecies

    count = x*y - sum(sum(indecies))
    percentLoss = 100-count/(x*y)*100
    print("Keeping "+str(percentLoss)+"% of the original image...")
    
    ### Inverse 2d DFT ###
    rowTrans = []
    for row in fourierSumFiltered: #Done on rows
        rowTrans.append(IFFTransform(row))

    fourierEndSum = []
    for col in myMath.invertMatrix(rowTrans):
        fourierEndSum.append(IFFTransform(col))
    compressedImgData = np.rint(np.array(myMath.invertMatrix(fourierEndSum)))

    for x in range(len(compressedImgData)): #Reformatting image
        for y in range(len(compressedImgData[x])):
            if compressedImgData[x][y] < 0:
                compressedImgData[x][y] = 0
            if compressedImgData[x][y] > 255:
                compressedImgData[x][y] = 255

    return compressedImgData.astype('uint8')