import sys
import numpy as np
import array
import matplotlib.pyplot as plt

excludelist = ["<positionref","<rotationref"]
def corrfile(infile,outfile):
    file1 = open(infile)
    outlines = ""
    while(True):
        line1 = file1.readline()
        if (not line1):
            break
        rmv_flag=False
        for i in excludelist:
            if (line1[:10]==i[:10]):
                rmv_flag=True
        if (not rmv_flag):
            outlines+=line1
    print(outlines)
        

corrfile(sys.argv[1],"test")