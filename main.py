# Main Python file for PJAS Project #
# Dev start date: 12/20/2022 #
# Dev end date:              #

#Does image compression through multiple different methods and compares#

import numpy as np
import math

#Custom Libraries
import Picture
import CompressionAlgos



##### SDV Compression Example #####
KValue = 100

#- - - - - - - Greyscale - - - - - - #
a = Picture.getImg("pjaslogo.jpg", greyscale=True)

compressiondata = CompressionAlgos.SDV(a, KValue) #Compress the image data

#Save the image!!!
#print(np.array_equal(a, compressiondata))

Picture.savImg("Tester.jpg", compressiondata)

#- - - - - - - - RGBA - - - - - - - -#
'''rgbaValues = Picture.getImg("steinsgate.png", greyscale=False)

#Separate each channel (Red, Green, Blue, and Alpha)
r, g, b, a = Picture.convRGBA(rgbaValues, encode=True)

#Compress Image Data (for each channel)
rcomp = CompressionAlgos.SDV(r, KValue)
gcomp = CompressionAlgos.SDV(g, KValue)
bcomp = CompressionAlgos.SDV(b, KValue)
#acomp = CompressionAlgos.SDV(a, KValue)

#Recombine each channel (Red, Green, Blue, and Alpha)
pictureData = Picture.convRGBA([rcomp, gcomp, bcomp, a], encode=False)

Picture.savImg("SteinsGATECOLOR.png", pictureData, greyscale=False)'''