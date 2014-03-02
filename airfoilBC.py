# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombayimport numpy

# Functions needed for creation of airfoil boundary conditions
# reflect function is to be modified

import math
NODETOL = 1e-08
def readPoints(fileName):
	"""Read boundary points from a file. Example: naca0012.txt"""
	f=open(fileName)
	a=f.readlines()
	a=a[2:]
	points=numpy.zeros([len(a),2])
	for i in range(len(a)):
		b=a[i]
		b=b.split(' ')
		k=0
		for j in range(len(b)):
			if b[j]!='':
				points[i][k]=float(b[j])
				k=k+1
	return(points)

def AOA(alpha,points):
	"""Gives an positive angle of attack alpha to the airfoil or any body specified by points"""
	newPoints=points.copy()
	for i in range(len(points)):
		z=complex(points[i][0],points[i][1])
		z=z*(math.e**(-1j*alpha/180*math.pi))
		newPoints[i][0]=z.real
		newPoints[i][1]=z.imag
	return(newPoints)

def findNormal(x1,x2):
	"""Find a normal which connects points x1 and x2"""
	if x2[1]==x1[1]:
		n=numpy.array([0.,1.])
		return(n)
	if x2[0]==x1[0]:
		nSlope=0.0
	else:
		slope=(x2[1]-x1[1])/(x2[0]-x1[0])
		nSlope=-1./slope
	ny=nSlope
	nx=1.
	n=numpy.array([nx,ny])
	n=n/(nx*nx+ny*ny)**0.5
	return(n)

def airfoilBC(filename,alpha=0.0):
	"""Returns required properties for a boundary condition.
	Reads filename in order to get airfoil points
	Applies angle of attach alpha"""
	
	# Read points
	points=readPoints(filename)
	points=AOA(alpha,points)
	NPoints=len(points)
	
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
		normals[i]=findNormal(points[k1],points[k2])		
		cpl[i]=cp
	
	# Check whether a point is inside or outside the airfoil
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
	# Find normals at control points
	# Find whether a point is inside or outside the airfoil
	# Find a reflect function
	return(BCpoints,BCcp,BCNormals,BCinFunction,BCreflectFunction)

def intersect(BCpoints,pos1,pos2):
	""" Find intersection point of a line connecting pos1, pos2 and 
	a polygon with BCpoints. The function returns, the sheet id, 
	and point location(pos) where they intersect"""
	NPoints=len(BCpoints)
	for i in range(NPoints):
		k1=i
		if i==NPoints-1:
			k2=0
		else:
			k2=i+1
		[inters,pos]=intersecrLine(pos1,pos2,BCpoints[k1],BCpoints[k2])
		if inters==True:
			return(i,pos)
	return(0)

def intersectLine(pos1,pos2,pos3,pos4):
	""" Checks whether lines [pos1,pos2] and [pos3,pos4] intersect.
	Gives the point of interaction"""
	s1=(pos2[1]-pos1[1])/(pos2[0]-pos1[0])
	s2=(pos4[1]-pos3[1])/(pos4[0]-pos3[0])
	
	# Find is lines are parallel
	if abs(s2)-abs(s1)<NODETOL:
		return(False,0)
	
	# Find intersection point
	A=numpy.array([pos2-pos1,pos4-pos3]).transpose()
	B=numpy.transpose(pos3-pos1)
	paras=numpy.linalg.solve(A,B)
	
	# Check whether the point actually lies between the given points
	if paras.all>=0 and paras.all<=1:
		return(True,pos1+paras[0]*(pos2-pos1))
	return(False,0)
