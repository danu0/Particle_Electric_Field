import ode_rk4
import numpy as np
import matplotlib.pyplot as plt

class ElectronRK4(ode_rk4.ODE_RK4):
  """A class derived from the class Ode_rk4.
  """
  def __init__(self, V0=[], dt=0.1, t0=0, Ex=0, Ey=0, B=0, r_circle=0.01):
      """Set the ODE parameters
      : param V0 : initial value of V 
      : param t0 : initial time (usualy 0)
      : param dt : integration time step
      : param a  : equation parameter
      """  
      super().__init__(V0, dt=dt, t0=t0) 
      self.Ex=Ex
      self.Ey=Ey
      self.B=B
      self.q=-1.602e-19
      self.m=9.109e-31
      self.qom=self.q/self.m
      self.r_circle = r_circle
      self.t = t0

  def F(self, t, V): #input is (x,y,vx,vy)
      """ equation to solve: 
          dx/dt = 
          dy/dt =
          dg/dt = self.q(self.Ex + h(x)self.B)
          dh/dt = 
          v[0] is f  and v[1] is g
      : param t : current time 
      : param v : current function as a vector 
      """
      if(np.hypot(V[0],V[1])>self.r_circle):
        B = 0
        Ex = 0
        Ey = 0
      else:
        B = self.B
        Ex = self.Ex
        Ey = self.Ey
      eq1 = V[2]
      eq2 = V[3]
      eq3 = self.qom*(Ex + V[3]*B)
      eq4 = self.qom*(Ey - V[2]*B)
      return(np.array([eq1,eq2,eq3,eq4]))

  def iterate_xmax(self, xmax, fig_dt=-1):
      """ Solve the system of equations DN/dt = F(N) until t=tmax.
          Save N and t in lists N_list and t_list every fig_dt.

      :param xmax   : integration upper bound
      :param fig_dt : interval between data point for figures (use dt if < 0)
      """
       
      if(fig_dt < 0) : fig_dt = self.dt*0.99 # save all data
      
      next_fig_t = self.t*(1-1e-15) # ensure we save the initial values

      xmax -= self.dt*0.1   # stop as close to tmax as possible
      while(abs(self.V[0]) < xmax): # integrate until tmax
        self.one_step()
        if(self.t >= next_fig_t): # save fig when next_fig_t is reached       
          self.V_list.append(np.array(self.V)) # force a copy of V!
          self.t_list.append(self.t)
          next_fig_t += fig_dt # set the next figure time

