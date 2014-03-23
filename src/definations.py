# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombay

# Class and function definations for basic structure

import numpy
import math
import cmath
from copy import deepcopy
import matplotlib.pyplot as plt
import threading
NODETOL=1e-08

#### Define Elements ####

class vortex:
	"""A single vortex. It can be either a point vortex or a blob. 
	It has following properties. Location of vortex([x,y,z]), 
	Strength of vortex, Blob type (0 for no blob, 1 for Krasny blob, 2 for Chorin blob), 
	delta for blobs, Trace flag"""
	
	def __init__(self,position, strength,blobType=0,delta=0.0,traceFlag=0):
		self.position=position
		self.strength=strength
		self.pointid=0
		self.blobType=blobType
		self.delta=delta
		self.traceFlag=traceFlag
		if self.traceFlag==1:
			self.traceHist=[self.position]
	
	def fieldEffect(self,pos):
		"""Calculates field at "pos" due to the vortex"""
		z=pos[0]+((1j)*pos[1])
		z0=self.position[0]+((1j)*self.position[1])
		gamma=self.strength
		if abs(z-z0)<NODETOL:
			return([0.,0.])
		if self.blobType==0:
			complexVelocity=gamma*(1j)/2/math.pi/(z-z0)
			return([complexVelocity.real,-complexVelocity.imag])
		elif self.blobType==1:				#Krasny Blob
			delta=self.delta
			complexVelocity=gamma*(1j)/2/math.pi/(z-z0)*(abs(z-z0)**2/(abs(z-z0)**2+(delta**2)))
			return([complexVelocity.real,-complexVelocity.imag])
		elif self.blobType==2:				#Chorin Blob
			delta=self.delta
			if abs(z-z0)>delta:
				complexVelocity=gamma*(1j)/2/math.pi/(z-z0)
				return([complexVelocity.real,-complexVelocity.imag])
			else:
				complexVelocity=gamma*(1j)/2/math.pi/(z-z0)*abs(z-z0)/delta
				return([complexVelocity.real,-complexVelocity.imag])
				
	def modifyPos(self,newPos):
		"""Modify position of the vortex"""
		self.position=newPos
		if self.traceFlag==1:
			self.traceHist.append(self.position)
			
def testVortex(strength=1.0,blobType=0,delta=0.0):
	"""Test function for a point vortex"""
	position=numpy.array([0.,0.])
	V=vortex(position, strength,blobType,delta,traceFlag=0)
	X,Y = numpy.meshgrid(numpy.arange(-1.,1.,0.1),numpy.arange(-1,1,0.1))
	u=numpy.copy(X)
	v=numpy.copy(X)
	
	for i in range(len(X.flatten())):
		vel=V.fieldEffect([X.flatten()[i],Y.flatten()[i]])
		u.ravel()[i]=vel[0]
		v.ravel()[i]=vel[1]
	
	plt.figure()
	plt.quiver(X,Y,u,v)
	plt.title('Vector field due to a single vortex')
	plt.scatter(position[0],position[1])
	plt.show()
	return()
	
class tracer:
	"""A tracer point : Points with zero strength"""
	
	def __init__(self,position,traceFlag=0):
		self.position=position
		self.traceid=0
		self.traceFlag=traceFlag
		if self.traceFlag==1:
			self.traceHist=[self.position]
	
	def modifyPos(self,newPos):
		"""Modify the position of tracer particle"""
		self.position=newPos
		if self.traceFlag==1:
			self.traceHist.append(self.position)
	
	def fieldEffect(self,pos):
		"""field effect due to a tracer point will be zero"""
		return(vector.zeroValue)

class linVortex:
	"""Linear vortex sheet"""
	
	def __init__(self,gamma1,gamma2,x1,x2):
		self.gamma1=gamma1
		self.gamma2=gamma2
		self.x1=x1
		self.x2=x2
		self.diff=x2-x1
		self.lemda=(self.diff.dot(self.diff))**0.5
		self.linid=0
		zdiff=complex(self.diff[0],self.diff[1])
		self.theta=(cmath.log(zdiff/abs(zdiff))/1.j).real
	
	def fieldEffect(self,pos):
		"""Velocity field due to this sheet on a pos"""
		z1=complex(self.x1[0],self.x1[1])
		z2=complex(self.x2[0],self.x2[1])
		zd=complex(pos[0],pos[1])
		if (abs(z1-zd)<NODETOL) or (abs(z2-zd)<NODETOL):
			return(numpy.array([0.,0.]))
		z=(zd-z1)*(math.e**(-1j*self.theta))
		pa1=self.gamma1*((((z/self.lemda)-1)*cmath.log((z-self.lemda)/z))+1)
		pa2=self.gamma2*((z/self.lemda*cmath.log((z-self.lemda)/z))+1)
		vel=1.j/2./math.pi*(pa1-pa2)
		vel=vel*(math.e**(-1j*self.theta))
		return(numpy.array([vel.real,-vel.imag]))

def testLinVort(p1=numpy.array([0.,0.]),p2=numpy.array([1.,1.]),g1=0.1,g2=0.1):
	L=linVortex(g1,g2,p1,p2)
	X,Y = numpy.meshgrid(numpy.arange(-4.,4.,0.4),numpy.arange(-3.,3.,0.4) )
	u=numpy.copy(X)
	v=numpy.copy(X)

	for i in range(len(X.flatten())):
		vel=L.fieldEffect([X.flatten()[i],Y.flatten()[i]])
		u.ravel()[i]=vel[0]
		v.ravel()[i]=vel[1]

	plt.figure()
	Q = plt.quiver(X,Y,u,v)
	plt.title('Vector field due to Linear Vortex')
	plt.plot([p1[0],p2[0]],[p1[1],p2[1]])
	plt.show()
	return()

class vector:
	"""Vector: Properties: 1. Location, 2. Value . This is a 2 dimensional vector"""
	zeroValue=numpy.zeros(2)
	def __init__(self,position,value=[]):
		self.value=value
		self.position=position

#### Define element lists ####

class vortexList:
	"""List of vortices. allV is list of all the vortices, nP is total number of vortices"""
	
	def __init__(self):
		self.allV=[]
		self.nPoints=0
	
	def addVortex(self,position, strength,blobType=0,delta=0.0,traceFlag=0):
		"""add a vortex to this vortex list"""
		self.allV.append(vortex(position, strength,blobType,delta,traceFlag))
		self.allV[-1].pointid=self.nPoints
		self.nPoints=self.nPoints+1
	
	def addLists(self,list1):
		"""add two vortex lists"""
		for i in range(list1.nPoints):
			position=list1.allV[i].position
			strength=list1.allV[i].strength
			blobType=list1.allV[i].blobType
			delta=list1.allV[i].delta
			traceFlag=list1.allV[i].traceFlag
			self.addVortex(position, strength,blobType,delta,traceFlag)
	
	def fieldEffect(self,pos):
		"""velocity field due to the whole vortex list at a position"""
		newValue=vector.zeroValue
		activeVortices=self.allV
		for eachPointActive in activeVortices:
			newValue=newValue+eachPointActive.fieldEffect(pos)
		return(newValue)
	
	def posit(self):
		"""returns a list containing positions of all the points"""
		pos=[eachV.position for eachV in self.allV]
		return(pos)

def testVortexList(strength=1.0,blobType=0,delta=0.0):
	"""Test function for a point vortex list"""
	V=vortexList()
	Npoints=50
	lowLim=numpy.array([0.,0.])
	upLim=numpy.array([1.,1.])
	for i in range(Npoints):
		[x,y]=(lowLim+(float(i)/Npoints*(upLim-lowLim)))
		V.addVortex([x,y], strength,blobType,delta)
	
	X,Y = numpy.meshgrid(numpy.arange(-2.,2.,0.2),numpy.arange(-2,2,0.2))
	u=numpy.copy(X)
	v=numpy.copy(X)
	
	for i in range(len(X.flatten())):
		vel=V.fieldEffect([X.flatten()[i],Y.flatten()[i]])
		u.ravel()[i]=vel[0]
		v.ravel()[i]=vel[1]
	
	plt.figure()
	plt.quiver(X,Y,u,v)
	plt.title('Vector field due to multiple vortices')
	plt.scatter([lowLim[0],upLim[0]],[lowLim[1],upLim[1]])
	plt.show()
	return()

class traceList:
	"""List of tracers. allV is list of all the vortices, nP is total number of vortices"""
	
	def __init__(self):
		self.allV=[]
		self.nPoints=0
	
	def addTracer(self,position,traceFlag=0):
		""" Add a tracer point to this tracer list"""
		self.allV.append(tracer(position,traceFlag))
		self.allV[-1].pointid=self.nPoints
		self.nPoints=self.nPoints+1
	
	def addLists(self,list1):
		"""Add two tracer lists"""
		for i in range(list1.nPoints):
			position=list1.allV[i].position
			traceFlag=list1.allV[i].traceFlag
			self.addTracer(position,traceFlag)
	
	def fieldEffect(self,pos):
		"""Field effect at any position due to a tracer field will be zero"""
		return(vector.zeroValue)
	
	def posit(self):
		"""returns a list containing positions of all the points"""
		pos=[eachV.position for eachV in self.allV]
		return(pos)

class linVortList:
	"""List of linear vortex sheets"""
	
	def __init__(self):
		self.allLin=[]
		self.nPoints=0
	
	def addLinVortex(self,gamma1,gamma2,x1,x2):
		"""Add a linear vortex sheet to the list"""
		self.allLin.append(linVortex(gamma1,gamma2,x1,x2))
		self.allLin[-1].linid=self.nPoints
		self.nPoints=self.nPoints+1
	
	def addLists(self,list1):
		""""Add two linear vortex sheet lists"""
		for i in range(list1.nPoints):
			gamma1=list1.allLin[i].gamma1
			gamma2=list1.allLin[i].gamma2
			x1=list1.allLin[i].x1
			x2=list1.allLin[i].x2
			self.addLinVortex(gamma1,gamma2,x1,x2)
	
	def fieldEffect(self,pos):
		"""Velocity field at pos due to all linear Vortex sheets in this list"""
		newValue=vector.zeroValue
		linSheets=self.allLin
		for eachLin in linSheets:
			newValue=newValue+eachLin.fieldEffect(pos)
		return(newValue)

def testLinVortList():
	"""Test function for a linear vortex sheet list"""
	V=linVortList()
	Npoints=50
	lowLim=numpy.array([0.,0.])
	upLim=numpy.array([1.,1.])
	for i in range(Npoints):
		p1=(lowLim+(float(i)/Npoints*(upLim-lowLim)))
		p2=(lowLim+(float(i+1)/Npoints*(upLim-lowLim)))
		V.addLinVortex(0.1,0.1,p1,p2)
	
	X,Y = numpy.meshgrid(numpy.arange(-2.,2.,0.2),numpy.arange(-2,2,0.2))
	u=numpy.copy(X)
	v=numpy.copy(X)
	
	for i in range(len(X.flatten())):
		vel=V.fieldEffect([X.flatten()[i],Y.flatten()[i]])
		u.ravel()[i]=vel[0]
		v.ravel()[i]=vel[1]
	
	plt.figure()
	plt.quiver(X,Y,u,v)
	plt.title('Vector field due to multiple linear vortex sheets')
	plt.plot([lowLim[0],upLim[0]],[lowLim[1],upLim[1]])
	plt.show()
	return()

#### Define additional Functions ####
		
def linearTimeMesh(startTime, timeStep, endTime):
	"""create a linear time mesh"""
	timeMesh=[]
	currentTime=startTime
	while (currentTime<endTime):
		timeMesh.append(currentTime)
		currentTime=currentTime+timeStep
	timeMesh.append(endTime)
	return(timeMesh)

def eulerInt(oldPos,fVal,dt):
	"""Returns new Position from old Position using fVal as the field value and using Euler integration"""
	oldPos=numpy.array(oldPos)
	newPos=numpy.copy(oldPos)
	for k in range(len(oldPos)):
		rhs=fVal[k]
		old=oldPos[k]
		new=old+(dt*(rhs))
		newPos[k]=new
	return(newPos)

def findParticle(pos,vList):
	"""find particle in vortex list or tracer list"""
	for eachVortex in vList.allV:
		distance=numpy.array(eachVortex.position)-numpy.array(pos)
		if distance.dot(distance)<=NODETOL and gotParticle==0:
			tpid=eachVortex.pointid
			return(pid)
	return()

def velField(pos,fieldGen,vinf=0.0):
	"""Define velocity field at set of points due to vinf and fieldGen( a list containing various velocity generator lists)"""
	field=numpy.zeros(pos.shape)
	print "number of Positions = %i" %len(pos)
	for i in range(len(pos)):
		vel=vector.zeroValue
		for eachGen in fieldGen:
			vel=vel+eachGen.fieldEffect(pos[i])
		field[i]=vel+vinf
	return(field)
