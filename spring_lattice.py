import numpy as np
import matplotlib.pyplot as plt
from ode_rk4 import ODE_RK4
import timeit
import time

class springs(ODE_RK4):
   def __init__(self, M, k, V0, dt, t0=0):
     """ : param M : mass of each item
       : param k : spring constant
       : param V0: Initial value [x0,x1,...xN, v0,v1,...vN] 
       : param dt: integration time step
       : param t0: intial time
     """
     self.M = M
     self.k = k
     self.N = len(V0)//2 # number of masses. must be an integer
     self.reset(V0, dt, t0)
     
   def F(self, t, V):
      f = np.zeros(self.N*2, dtype='float64')
      for i in range(self.N):
          f[i] = V[self.N+i]   # dx[i]/dt = v[i]
          if (i==0):
            # dv[0]/dt = k/M(x[1]-2x[0])
            f[self.N+i] = self.k/self.M*(V[i+1]-2*V[i])
          elif (i== self.N-1):
            # dv[N-1]/dt = k/M(x[N-2]-2x[N-1])
            f[self.N+i] = self.k/self.M*(V[i-1]-2*V[i])
          else:
            # dv[i]/dt = k/M(x[i+1]+x[i-1]-2x[i])
            f[self.N+i] = self.k/self.M*(V[i+1]+V[i-1]-2*V[i])
      return(f)
  
   def plot_last(self):
      #plt.axis('off')
      amp = 0.2
      plt.axis((-1,self.N,0,2*amp))
      plt.plot(range(0,self.N),amp+self.V_list[-1][0:self.N],"b.")
      plt.show()

if __name__ == "__main__":
      N = 1000
      V0 = np.zeros(N*2, dtype='float64')
      V0[N//2] = 0.5 # we move 1 mass away from its rest position
      sp = springs(M=1, k=1, V0=V0, dt=0.1)

      t1 = time.clock() # we time how long iterate takes to run
      sp.iterate(100,1)
      t2 = time.clock()
      print("Running time: {}s".format(t2-t1))
      sp.plot_last()

