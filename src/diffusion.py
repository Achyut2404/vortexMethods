# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombay

# Functions for applying vorticity diffusion

import numpy
import definations as dfn
import math
import random
import matplotlib.pyplot as plt

def applyRVM(dt,nu=0.1,toMod=[dfn.vortexList(),dfn.traceList()],BCList=[]):
	"""Apply RVM diffusion based on kinematic viscocity nu on "toMod" list of lists
	Also, reflect any vortex which is getting diffused into the boundary"""
	
	#Calculate mu and sigma for random number generation gaussian distribution
	mu=0.0
	sigma=(2*nu*dt)**0.5
	
	#Modify position for each point
	for eachList in toMod:
		#Diffusion should not be applied on tracers. This is only for vortex lists. check that condition
		if eachList.__class__.__name__=='vortexList':
			for eachPoint in eachList.allV:
				rad=random.gauss(mu,sigma)
				theta=random.uniform(-math.pi,math.pi)
				dPos=numpy.array([rad*math.cos(theta),rad*math.sin(theta)])
				newPos=eachPoint.position+dPos
				#Reflect the new position if it is in boundary
				for eachBC in BCList:
					if (eachBC.inBoundary(newPos)):
						newPos=eachBC.reflect(eachPoint.position,dPos)
				eachPoint.modifyPos(newPos)          #This will also change fieldGen point positions
	return()

def testRVM(cpos=numpy.array([1.,1.]),Np=100):
	""" Test RVM diffusion function"""
	
	#Initiate vortex list and tracer list and toMod
	V1=dfn.vortexList()
	T1=dfn.traceList()
	toMod=[V1,T1]
		
	#Add points at the same location
	for i in range(Np):
		V1.addVortex(cpos,0.1)
	
	#Add arbitary tracers, just for check
	T1.addTracer(numpy.array([.5,.5])+cpos)
	T1.addTracer(numpy.array([.5,-.5])+cpos)
	T1.addTracer(numpy.array([-.5,.5])+cpos)
	T1.addTracer(numpy.array([-.5,-.5])+cpos)
	
	#Apply RVM
	dt=0.1
	nu=0.1
	applyRVM(dt,nu,toMod)

	#Positions of all points in toMod
	pos=[]
	for eachList in toMod:
		pos=pos+eachList.posit()
	pos=numpy.array(pos)

	#Plot all points. Red is the original vortex, Blue is new distribution
	plt.figure()
	plt.title("Random vortex method diffusion test")
	for i in range(len(pos)):
		plt.plot(pos[i][0],pos[i][1],'bo',markersize=2.0)
	plt.plot(cpos[0],cpos[1],'ro',markersize=4.0)
	plt.show()
	
	return()
