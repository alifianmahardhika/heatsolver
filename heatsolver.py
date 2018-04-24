import numpy as np
import sys as system

fname = system.argv[1]
data = np.fromfile(fname, sep="\n")
first_line = data[0:3:1]
point = int(first_line[0])
element = int(first_line[1])
border = int(first_line[2])
npxylim = point*3+3
elnplim = element*4+npxylim
bdnplim = border*3+elnplim
print(npxylim,elnplim,bdnplim)
npxy = data[3:npxylim:1]
elnp = data[npxylim:elnplim:1]
bdnp = data[elnplim:bdnplim:1]

""" print(point, element, border)
print(npxy)
print(npxy.shape)
print(type(npxy))
npxy1 = np.reshape(npxy,(point,3))
print(npxy1)
print(npxy1.shape) """

""" print(elnp)
print(elnp.shape)
print(type(elnp))
elnp1 = np.reshape(elnp,(element,4))
print(elnp1)
print(elnp1.shape) """

print(bdnp)
print(bdnp.shape)
print(type(bdnp))
bdnp1 = np.reshape(bdnp,(border,3))
print(bdnp1)
print(bdnp1.shape)
