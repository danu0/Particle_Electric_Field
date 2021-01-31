import ode_rk4
import numpy as np
import matplotlib.pyplot as plt

class ElectronRK4(ode_rk4.ODE_RK4):
  """A class derived from the class Ode_rk4.
  """
  def __init__(self, V0=[], dt=0.1, t0=0, Ex=0, Ey=0, B=0):
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
      eq1 = V[2]
      eq2 = V[3]
      eq3 = self.qom*(self.Ex + V[3]*self.B)
      eq4 = self.qom*(self.Ey - V[2]*self.B)
      return(np.array([eq1,eq2,eq3,eq4]))


