# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombayimport definations

# Boundary condition implementations
# Can apply no-penetration and no-slip boundary conditions
# Need to modify for moving wall boundary conditions

import math
import cmath
import numpy
import definations as dfn
import matplotlib.pyplot as plt
NODETOL=1e-08

#### Define boundary condition classes ####

class wallNPNSBC:
	""" wall boundary with no penetration and no-slip condition"""
	
	def __init__(self,points,cp,normals,inFunction=0,reflectFunction=0):
		self.points=points
		self.cp=cp
		self.normals=normals
		self.A=linNPA(points,cp,normals)
		self.inBoundary=inFunction
		self.reflect=reflectFunction
		# Find cps points, that are control points for no Slip conditions.
		# located just above the normal control points for no penetration
		self.cps=self.cp+(100*NODETOL)*(self.normals)
	
	def findVcp(self,fieldGens=[],vinf=0.0):
		"""find V at control points.
		This step is to be executed at starting of each step"""
		self.vcp=definations.velField(self.cp,fieldGens,vinf)
	
	def findVcps(self,fieldGens,vinf=0.0):
		"""For slip boundary conditions Vcp is to be find slightly above the control point
		that is done by using this function. """
		self.vcps=definations.velField(self.cps,fieldGens,vinf)
	
	
	def closeNPBC(self,fieldGens):
		"""This step is to remove all the no penetration field generators introduced.
		This is to be applied after each time step is over"""
		fieldGens.pop(self.fieldGenIndex)	
	
	def applyNPBC(self,fieldGens=[]):
		"""modifies fieldGens in order to apply boundary conditions"""
		self.fieldGenIndex=len(fieldGens)
		linList=NPLinList(self.points,self.vcp,self.A,self.cp,self.normals)		
		fieldGens.append(linList)
		self.initFlag=1
	
	def applyNSBC(self,fieldGens,toMod,gmin=0.1,delta=0.03):
		""" adds vortices near the boundary to both toMod and fieldGens.
		In order to satisfy no-slip condition"""
		
		#Find vSlips
		N=len(self.points)
		vSlip=numpy.zeros(N)
		tange=numpy.array([numpy.array([-self.normals[i][1],self.normals[i][0]]) for i in range(N)])
		vSlip =numpy.array([-self.vcps[i].dot(tange[i]) for i in range(N)])
		#Calculate parameters
		lemda=delta*math.pi
		nBlobs=vSlip/gmin
		direc=numpy.sign(nBlobs)
		nBlobs=numpy.vectorize(int)(abs(nBlobs))
		blobLoc=self.cp+delta*(self.normals)
		blobStren=gmin*lemda
		
		#introduce Chorin blobs
		BlobList=dfn.vortexList()
		for i in range(N):
			for j in range(nBlobs[i]):
				BlobList.addVortex(position=blobLoc[i], strength=direc[i]*blobStren, blobType=2, delta=delta, traceFlag=0)
		toMod.append(BlobList)
		fieldGens.append(BlobList)
		return()
		
#### Non penetration wall boundary functions ####

def NPLinList(points,vcp,A,cp,normals):
	"""Gives back a linear vortex sheet list, which will satisfy non penetration boundary conditions on wall defined by points. Points should be in the order of wall boundary line"""
	B=NPB(points,vcp,normals)
	gammas=findGamma(A,B)
	linList=definations.linVortList()
	for i in range(len(points)):
		if i==len(points)-1:
			linList.addLinVortex(gammas[i],gammas[0],points[i],points[0])
		else:
			linList.addLinVortex(gammas[i],gammas[i+1],points[i],points[i+1])
	return(linList)

def linNPA(points,cp,normals):
	"""returns the A matrix required for calculating wall Non penetration boundary conditions based on linear vortex panels"""
	N=len(points)
	A=numpy.zeros([N+1,N])
	gamma=numpy.zeros(N)
	for i in range(N):
		for j in range(N):
			l1=j
			if j==N-1:
				l2=0
				l0=j-1
			elif j==0:
				l2=j+1
				l0=-1
			else:
				l2=j+1
				l0=j-1
			linVo1=definations.linVortex(1.,0.,points[l1],points[l2])
			linVo2=definations.linVortex(0.,1.,points[l0],points[l1])
			vel=linVo1.fieldEffect(cp[i])+linVo2.fieldEffect(cp[i])
			vdotN=vel.dot(normals[i])
			A[i][j]=vdotN
	
	#Keep circulation zero
	A[N][:]=1.0
	return(A)

def findGamma(A,B):
	"""Solve Ax=B to give gammas for boundary conditions"""
	gamma=numpy.linalg.lstsq(A,B)
	return(gamma[0])


def NPB(points,vcp,normals):
	"""Give the B matrx(right hand side) for boundary condition matrix for wall non penetration boundary conditions"""
	N=len(points)
	B=numpy.zeros(N+1)
	for i in range(N):
		B[i]=	-vcp[i].dot(normals[i])
	return(B)

#### Functions for reading or defining boundary points ####	
def cylBCPoints(rad=1.,NPoints=20):
	"""Create cylinder boundary points with radius rad, and number of points NPoints
	Also it is the boundary creaters respolsibility to provide normals and control points"""
	
	points=[]
	
	#Find points
	for i in range(NPoints):
		theta=2*math.pi/NPoints*i
		point=numpy.array([rad*math.cos(theta),rad*math.sin(theta)])
		points.append(point)
	
	# Find control points and normals
	normals=numpy.zeros([NPoints,2])
	cpl=numpy.zeros([NPoints,2])
	for i in range(NPoints):
		k1=i
		if i==NPoints-1:
			k2=0
		else:
			k2=i+1
		cp=(points[k1]+points[k2])/2.0
		normals[i]=cp/(cp.dot(cp))**0.5		#For circle with origin center		
		cpl[i]=cp
	
	def inFunction(pos):
		"""A function which says whether a point is inside the boundary or nor"""
		return(pos[0]**2+pos[1]**2<rad**2)
	
	def reflectFunction(pos,dPos):
		"""A function that reflects a point which otherwise could have got inside the boundary"""
	
		#Find intersection point
		A=dPos.dot(dPos)
		B=2*pos.dot(dPos)
		C=pos.dot(pos)-rad**2
		L1=(-B+(B**2-4*A*C)**0.5)/2/A
		L2=(-B-(B**2-4*A*C)**0.5)/2/A
	
		# use a point with minimum L
		if abs(L1)<abs(L2):
			L=L1
		else:
			L=L2  
		iPoint=pos+L*dPos
	
		#Find normal at intersection point
		iNormal=iPoint/(iPoint.dot(iPoint))**0.5
	
		#Find reflected position
		newPos=pos+dPos-((dPos.dot(iNormal))*iNormal)
		return(newPos)
	return(points,cpl,normals,inFunction,reflectFunction)

#### Tests ####
def testCylBC1(Np=50):
	"""Test cylinder-constant velocity-no penetration boundary condition"""
	vinf=numpy.array([10.0,0.0])
	fieldGens=[]
	toMod=[]
	
	#Create cylinder points
	rad=1.0
	[points,cpl,normals]=cylBCPoints(rad,Np)
		
	#Add boundary condnitions
	BC=wallNPNSBC(points,cpl,normals)
	BC.findVcp(fieldGens,vinf)
	BC.applyNPBC(fieldGens)
	
	#Define mesh
	X,Y = numpy.meshgrid(numpy.arange(-3.,3.,10./Np),numpy.arange(-3.,3.,10./Np))
	u=numpy.copy(X)
	v=numpy.copy(X)
	
	#Create position array
	pos=[]
	for i in range(len(X.flatten())):
		pos.append([X.flatten()[i],Y.flatten()[i]])
	pos=numpy.array(pos)
	
	#Calculate field at these positions
	F=definations.velField(pos,fieldGens,vinf)
	
	# Arrange field values in u and v in order to plot them
	for i in range(len(X.flatten())):
		u.ravel()[i]=F[i][0]
		v.ravel()[i]=F[i][1]
		if (X.flatten()[i]**2+Y.flatten()[i]**2-rad**2)<=0.2:
			u.ravel()[i]=0.0
			v.ravel()[i]=0.0
	BC.closeNPBC(fieldGens)
	plt.figure()
	plt.quiver(X,Y,u,v)
	plt.title('Vector field due to boundary conditions')
	plt.plot(numpy.array(points).transpose()[0],numpy.array(points).transpose()[1],'o')
	plt.show()
