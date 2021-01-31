from electron_rk4_2 import *
import numpy as np
import matplotlib.pyplot as plt

Ex = -0.001   # 
Ey = 0     # 
B  = 1e-6  #

t0  = 0.0     # initial time
dt  = 1e-7    # integration time step
xmax= 0.02    # end of integration

x0= -0.01999              # initial x position
y0=0.0              # initial y position
vx0=100             # initial speed along the x axis
vy0=0.0             # initial speed along the y axis

ode = ElectronRK4([x0, y0, vx0, vy0], dt, t0, Ex, Ey, B) 

ode.iterate_xmax(xmax)
print(ode.t,ode.V[0],ode.V[1])

plt.xlabel("x", fontsize=22) # Set horizontal (x) figure label to "x"     
plt.ylabel("y", fontsize=22) # Set vertical (y) figure label to "y"   
plt.axis('equal')  # The x and y scales are identical. To see a circle 
ode.plot(2, 1, "g*") # plot y(x) using green stars
#plt.savefig("5_4_c_1.pdf")
plt.show()

plt.xlabel("t", fontsize=22) # Set horizontal (x) figure label to "t"     
plt.ylabel("x", fontsize=22) # Set vertical (y) figure label to "x"   
ode.plot(1, 0, "b-") # plot x(t) in blue
#plt.savefig("5_4_c_2.pdf")
plt.show()

plt.xlabel("t", fontsize=22) # Set horizontal (x) figure label to "t"     
plt.ylabel("y", fontsize=22) # Set vertical (y) figure label to "y"   
ode.plot(2, 0, "r-") # plot y(t) in red
#plt.savefig("5_4_c_3.pdf")
plt.show()

