from __future__ import division
from pylab import *
from scipy.optimize import fmin, fmin_cg, fmin_powell
from scipy.integrate import odeint
from numpy.random import *
ion()

def ActivationCurve(V, Vm, k):
  return 1/(1+exp(-(V-Vm)/k))
  
def Kinetic(V, Vmax, sigma, Camp, Cbase):
  return Cbase + Camp*exp(-(Vmax-V)**2/sigma**2)

def Gate(x,t):
  global pars
  xinf = ActivationCurve(V[:,1], pars[0], pars[1])
  tau = Kinetic(V[:,1], pars[2], pars[3], pars[4], pars[5])
  return (xinf-x)/tau

#figure(facecolor=[1,1,1])  

def FullTrace(p,time,y):
  global pars
  pars = p[0:6]
  rinit1 = ActivationCurve(V[:,0], pars[0], pars[1]) # initial values
  r1 = odeint(Gate,rinit1,time)
  pars = p[7:13]
  rinit2 = ActivationCurve(V[:,0], pars[0], pars[1]) # initial values
  r2 = odeint(Gate,rinit2,time)
  I = p[6]*r1*(V[:,1]-E)+p[13]*r2*(V[:,1]-E)
  plot(time,y,'b'); hold(True)
  plot(time,I,'--r');hold(False)
  draw()
  A = 0
  for a in arange(len(p)):
    if p[a] > HB[a] :
      A = A + (p[a] - HB[a])*1e8
    if p[a] < LB[a] :
      A = A + (LB[a]-p[a])*1e8
  return sum(array(y - I)*array(y - I))+A

global pars
  
V = array([[-50, -60],[-50, -70],[-50, -80],[-50, -90],[-50, -100],\
[-50, -110],[-50, -120],[-50, -130],[-50, -140],[-50, -150]])

par1 = [-120, -7, -110, 60, 1000, 10, 4]
par2 = [-100, -5, -90, 40, 300, 5, 4]

E = -35
time = arange(0.0,1000,10)

pars = par1
xinit1 = ActivationCurve(V[:,0], pars[0], pars[1]) # initial values
x1 = odeint(Gate,xinit1,time)
pars = par2
xinit2 = ActivationCurve(V[:,0], pars[0], pars[1]) # initial values
x2 = odeint(Gate,xinit2,time)
y = par1[6]*x1*(V[:,1]-E)+par2[6]*x2*(V[:,1]-E)+10*randn(len(time),len(V))

for NbRandTries in arange(0,10):
  P = par1+par2
  tol = 0.01
  LB = [(1-tol)*x for x in P]
  HB = [(1+tol)*x for x in P]
  p0 = LB+rand(len(P))*(array(HB)-array(LB))

  
  figure(0)
  parId = fmin_powell(FullTrace, p0, args=(time,y),maxfun = 5000)
  print parId

  figure(1,facecolor=[1,1,1])
  Vplot = arange(-180,-50)

  subplot(2,1,1)
  hold(True)
  plot(Vplot,ActivationCurve(Vplot, parId[0], parId[1]),'k',linewidth = 0.5)
  plot(Vplot,ActivationCurve(Vplot, parId[7], parId[8]),'k',linewidth = 0.5)

  subplot(2,1,2)
  hold(True)
  plot(Vplot,Kinetic(Vplot, parId[2], parId[3], parId[4], parId[5]),'k',linewidth = 0.5)
  plot(Vplot,Kinetic(Vplot, parId[9], parId[10], parId[11], parId[12]),'k',linewidth = 0.5)
  
  draw()


subplot(2,1,1)
hold(True)
plot(Vplot,ActivationCurve(Vplot, par1[0], par1[1]),'r')
plot(Vplot,ActivationCurve(Vplot, par2[0], par2[1]),'r')

subplot(2,1,2)
hold(True)
plot(Vplot,Kinetic(Vplot, par1[2], par1[3], par1[4], par1[5]),'r')
plot(Vplot,Kinetic(Vplot, par2[2], par2[3], par2[4], par2[5]),'r')

show()


