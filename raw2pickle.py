import MSFileReader
import pickle 
import numpy
import os
import sys

r = MSFileReader.ThermoRawfile(sys.argv[1])

rval = []
for i in range(r.GetFirstSpectrumNumber(),r.GetLastSpectrumNumber()):
    if r.GetMSOrderForScanNum(i) == 1:
        rl = r.GetLabelData(i)
        rval.append(numpy.array([rl[0][0],rl[0][1]])) 
        
print(len(rval))
pickle.dump( rval, open(os.path.splitext(sys.argv[1])[0]+".spectra.p", "wb" ) )
