import numpy as np
import matplotlib.pyplot as plt
import time, sys


nx = 41
domain = np.linspace(0,2,nx)
dx = 2 / (nx-1)
nt = 25
dt = 0.025
c = 1

u = np.ones(nx)
u[int(0.5 / dx): int(1 / dx +1)] = 2
plt.plot(domain,u)


un = np.ones(nx)
for n in range(nt):
    un = u.copy()
    for i in range(1,nx):
        #for i in range(nx):
            u[i] = un[i] - c*dt/dx*(un[i] - un[i-1])
    
    #if(n%2 ==0):
        #plt.plot(domain,u)
        
    

        
        
        
#plt.plot(domain,u)
plt.grid()
