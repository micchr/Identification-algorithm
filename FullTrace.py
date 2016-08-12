from __future__ import division
from pylab import *
from scipy.optimize import fmin, fmin_cg, fmin_powell
from scipy.integrate import odeint
from numpy.random import *

# This algorithm has been tested on hyperpolarisation activated cation currents for to evaluate independently 
# the effects of both Ih components in the vestibular nerve
# see Michel CB et al., Eur J Neurosci. 2015 Nov;42(10):2867-77. doi: 10.1111/ejn.13021. Epub 2015 Aug 6.
# Identification and modelling of fast and slow Ih current components in vestibular ganglion neurons.
# Python 2.7+ is required

ion()

def ActivationCurve(V, Vh, k):
  # these two curves are called characteristics curves and model the activation...
  # Vh is the half activation, k is the slope
  return 1/(1+exp(-(V-Vh)/k))

def Kinetic(V, Vmax, sigma, Camp, Cbase): 
  # ... and the kinetics of the ionic channel against voltage
  # bell-shaped curve centred in Vmax, with amplitude Camp, width sigma, and offset Cbase
  return Cbase + Camp*exp(-(Vmax-V)**2/sigma**2)

def Gate(x,t):
  # this equation modelize the opening of a channel, given the membrane potential and the two precedent equations
  global pars
  xinf = ActivationCurve(V[:,1], pars[0], pars[1])
  tau = Kinetic(V[:,1], pars[2], pars[3], pars[4], pars[5])
  return (xinf-x)/tau

  global pars
  # the following function is the optimisation part of the algorithm. The two chanels are modeled with initial conditions 
  # and the result is compared to the data, taking account the parameter physiological ranges for identificationdef FullTrace(p,time,y):
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
  # this plot is to see in real time the reconstructed curves trying to fit the data (very fun)
  
  A = 0
  # this part is added to constrain the algorithm with physiological search ranges
  for a in arange(len(p)):
    if p[a] > HB[a] :
      A = A + (p[a] - HB[a])*1e8
    if p[a] < LB[a] :
      A = A + (LB[a]-p[a])*1e8
  return sum(array(y - I)*array(y - I))+A

global pars
  
V = array([[-50, -60],[-50, -70],[-50, -80],[-50, -90],[-50, -100],\
[-50, -110],[-50, -120],[-50, -130],[-50, -140],[-50, -150]])
# The first part of the code is for the data modelling to fit (the algorithm is firt tested on simulated data) 
# voltage clamp protocol given rise to the currents

par1 = [-120, -7, -110, 60, 1000, 10, 4]
par2 = [-100, -5, -90, 40, 300, 5, 4]
# Initial parmeters the algorithm is supposed to find
# first two correspond to the activation curve
# following to the kinetic curve

E = -35
# reversal potential of the Ih currents

time = arange(0.0,1000,10)
# simulation duration in ms

# construction of the simulated data
pars = par1
xinit1 = ActivationCurve(V[:,0], pars[0], pars[1]) # initial values
x1 = odeint(Gate,xinit1,time)
# current 1

pars = par2
xinit2 = ActivationCurve(V[:,0], pars[0], pars[1]) # initial values
x2 = odeint(Gate,xinit2,time)
y = par1[6]*x1*(V[:,1]-E)+par2[6]*x2*(V[:,1]-E)+10*randn(len(time),len(V))
# current 2


for NbRandTries in arange(0,10):
  # the algorithm comprise NbRandomTries tries with random initial conditions, in order to evaluate and overcome
  # eventual local minima
  P = par1+par2
  
  tol = 0.50
  LB = [(1-tol)*x for x in P]
  HB = [(1+tol)*x for x in P]
  p0 = LB+rand(len(P))*(array(HB)-array(LB))
  # here the search ranges are around the actual parameter vetor with a tolerance of +/-50%

  figure(0)
  parId = fmin_powell(FullTrace, p0, args=(time,y),maxfun = 5000)
  print parId
  # optimisation function with real time display of the minimization on the...

  figure(1,facecolor=[1,1,1])
  Vplot = arange(-180,-50)
  # results are plotted for each individual try

  subplot(2,1,1)
  # for the activation curve
  hold(True)
  # for activation 1...
  plot(Vplot,ActivationCurve(Vplot, parId[0], parId[1]),'k',linewidth = 0.5)
  # ... and activation 2
  plot(Vplot,ActivationCurve(Vplot, parId[7], parId[8]),'k',linewidth = 0.5)
  
  # and for the kinetic curve
  subplot(2,1,2)
  hold(True)
  # for kinetic 1...
  plot(Vplot,Kinetic(Vplot, parId[2], parId[3], parId[4], parId[5]),'k',linewidth = 0.5)
  # ... and kinetic 2
  plot(Vplot,Kinetic(Vplot, parId[9], parId[10], parId[11], parId[12]),'k',linewidth = 0.5)
  
  draw()

# the good results (determined from first mode of the goodness of fit histogram) are averaged and plotted
subplot(2,1,1)
hold(True)
plot(Vplot,ActivationCurve(Vplot, par1[0], par1[1]),'r')
plot(Vplot,ActivationCurve(Vplot, par2[0], par2[1]),'r')

subplot(2,1,2)
hold(True)
plot(Vplot,Kinetic(Vplot, par1[2], par1[3], par1[4], par1[5]),'r')
plot(Vplot,Kinetic(Vplot, par2[2], par2[3], par2[4], par2[5]),'r')

show()


