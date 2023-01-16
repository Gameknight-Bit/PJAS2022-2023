# Custom Module that Converts images into readable formats #
# Does things like output BW tables and exports pixel matricies to imgs #

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

PATH = "testimgs/"
SAVE_PATH = "testimgs/newimgs/"


### Image Processing Things ###
def convGreyscale(imgData, encode: bool = True): #only needed for rgba images already converted
    #Converts to matrix of rows and cols of just pixel strength
    ar = []
    if encode == True: #Encodes image Data
        #print(imgData)
        for x in imgData: #Rows
            rowData = []
            for y in x: #Cols
                rowData.append(y[0])
            ar.append(rowData)

        return np.array(ar)
    else: #Decodes greyscale image data (Encodes into )
        for x in imgData: #Rows
            rowData = []
            for y in x: #Cols
                pixelData = [y, y, y, 255] #R, G, B, Alpha Channel
                rowData.append(pixelData)
            ar.append(rowData)

        return np.array(ar)
        
def convRGBA(imgData, encode:bool = True): #returns r, g, b, a as own img
    r, g, b, a = [], [], [], []

    if encode == True:
        for x in imgData: #Rows
            rRow, gRow, bRow, aRow = [], [], [], []
            for y in x: #Cols
                rRow.append(y[0])
                gRow.append(y[1])
                bRow.append(y[2])
                aRow.append(y[3])
            r.append(rRow)
            g.append(gRow)
            b.append(bRow)
            a.append(aRow)
        
        return (np.array(r, dtype="uint8"), np.array(g, dtype="uint8"), np.array(b, dtype="uint8"), np.array(a, dtype="uint8"))
    else:
        r = imgData[0]
        g = imgData[1]
        b = imgData[2]
        a = imgData[3]
        arr = []
        for x in range(len(r)): #Rows of r
            rowData = []
            for y in range(len(r[x])): #Cols
                rowData.append([r[x][y], g[x][y], b[x][y], a[x][y]])
            arr.append(rowData)

        return np.array(arr, dtype="uint8")
            

#################################

#   Global Functions   #
def getImg(imgName: str, greyscale: bool = True):
    """
    Gives Numpy Matrix of pixel values
    """
    img = Image.open(PATH+imgName, 'r')
    img.show()
    if greyscale == True:
        img = img.convert("L")
    else:
        img = img.convert("RGBA")
    pix_vals = np.array(img) #Contains for each pixel 

    #print(pix_vals)

    return pix_vals

def savImg(newImgName: str, data, greyscale: bool = True): #MUST BE SAVED AS PNG UNLESS ALPHA IS GONE!!!!

    img = None
    if greyscale == True:
        img = Image.fromarray(data, "L")
    else:
        img = Image.fromarray(data, "RGBA")
        #img = img.convert("P", palette=Image.ADAPTIVE, colors=8)
    #img.show()
    img.save(SAVE_PATH+newImgName)
    
#data = getImg("flowerbw.jpg")
#savImg("NewImg.jpg", data, True)
