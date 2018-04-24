import numpy as np
import sys as system

def openFile():
    try:
        inputFile = open(system.argv[1],'r')
        return inputFile
    except Exception as error:
        print('Failed to open file => %s' %error)
        system.exit(2)
meshfile = openFile()
meshfile.close()

""" nodes = np.fromfile(fname, dtype=int, count=3, sep=" ")
point = nodes[0]
element = nodes[1]
boundary = nodes[2]
print(type(point))
 """
fname = system.argv[1]
meshfile = open(fname, mode='r')
header = meshfile.readline()
data = meshfile.readlines()
meshfile.close
header = header.split(' ')
point = int(header[2])
print(point)
data[0] = str(data[0])
data[0] = data[0].split(' ')
print(data[0])

point = int(data[0])
element = int(data[1])
border = int(data[2])
print(point, element, border)
print(type(element))
print(data.shape)
i = 3
while (i <= point*3):
    x1 = float(data[i])
    x2 = float(data[i+1])
    pointlabel = float(data[i+2])
    print(x1,x2,pointlabel)
    i = i+3
""" for i in range(3,point*3,3):
    x1 = float(data[i])
    x2 = float(data[i+1])
    pointlabel = float(data[i+2])
    print(x1,x2,pointlabel) """