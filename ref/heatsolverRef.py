#3D heat equation
 
 
import numpy as np
from numpy import pi,exp,sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import matplotlib.animation as animation
 
fig = plt.figure()
fig.set_dpi(100)
ax1 = fig.gca(projection='3d')
 
x = np.linspace(-5,5,30)
y = np.linspace(-5,5,30)
 
X,Y = np.meshgrid(x,y)
 
#Initial time
t0 = 0
 
#Time increment
dt = 0.05
 
#Initial temperature at (0,0) at t0=0
T = 5
 
#Sigma squared 
s = 2
 
#Temperature function
def u(x,y,t):
    return (T/sqrt(1+4*t/s))*exp(-(x**2+y**2)/(s+4*t))
 
 
a = []
 
for i in range(500):
    z = u(X,Y,t0)
    t0 = t0 + dt
    a.append(z)
     
m = plt.cm.ScalarMappable(cmap=plt.cm.jet)
m.set_array(a[0])
cbar = plt.colorbar(m)
 
k = 0
def animate(i):
    global k
    temp = a[k]
    k += 1
    ax1.clear()
    ax1.plot_surface(X,Y,temp,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
 
    #ax1.contour(x,y,temp)
    ax1.set_zlim(0,T)
    ax1.set_xlim(-5,5)
    ax1.set_ylim(-5,5)
     
     
anim = animation.FuncAnimation(fig,animate,frames=220,interval=20)
plt.show()