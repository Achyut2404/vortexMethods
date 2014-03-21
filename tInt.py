# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombay

# RK-2 implementation
# Provided test cases can solve viscous incompressible N-S around a circular cylinder
# python -c "import tInt; tInt.test4RK2()"      will give a flow around circular cylinder

import numpy
import math
import cmath
import wallBC
import definations as dfn
import diffusion
import matplotlib.pyplot as plt
import plots

def advectRK2(dt,toMod=[dfn.vortexList(),dfn.traceList()],fieldGens=[dfn.vortexList(),dfn.linVortList()],vinf=0.0,BCList=[]):
	"""Modify "toMod(List format)" objects according to RK 2 advection based on fieldGens(List format) for advection"""
	
	#Calculate total number of toMod Points and append their positions
	N=0
	pos=[]
	for eachList in toMod:
		N=N+eachList.nPoints
		pos=pos+eachList.posit()
	pos=numpy.array(pos)
	oldPos=pos.copy()
	
	# Satisfy no penetration
	[eachBC.findVcp(fieldGens,vinf) for eachBC in BCList]
	[eachBC.applyNPBC(fieldGens) for eachBC in BCList]
	
	#Calculate no slip velocities
	[eachBC.findVcps(fieldGens,vinf) for eachBC in BCList]

	# Calculate field at all positions 
	field=dfn.velField(pos,fieldGens,vinf)

	# Close NP BC
	[eachBC.closeNPBC(fieldGens) for eachBC in BCList]

	#Take first step of RK 2 and modify positions
	n=0
	for eachList in toMod:
		for eachPoint in eachList.allV:	
			newPos=dfn.eulerInt(pos[n],field[n],dt/2.0)
			eachPoint.modifyPos(newPos)     # Note that this step will also modify appropriate positions in fieldGens, because of pointers
			n=n+1

	#Calculate total number of toMod Points and append their positions
	newN=0
	pos=[]
	for eachList in toMod:
		newN=newN+eachList.nPoints
		pos=pos+eachList.posit()
	pos=numpy.array(pos)
	
	if newN!=N:
		print "number of points cannot change in an advection step"
	
	# Satisfy no penetration
	[eachBC.findVcp(fieldGens,vinf) for eachBC in BCList]
	[eachBC.applyNPBC(fieldGens) for eachBC in BCList]
	
	#Calculate field with new positions
	field=dfn.velField(pos,fieldGens,vinf)

	# Close NP BC
	[eachBC.closeNPBC(fieldGens) for eachBC in BCList]

    # Take second step of RK 2 and modify final positions
	n=0
	for eachList in toMod:
		for eachPoint in eachList.allV:	
			newPos=dfn.eulerInt(oldPos[n],field[n],dt)
			dPos=newPos-oldPos[n]
			#Reflect the new position if it is in boundary
			for eachBC in BCList:
				if (eachBC.inBoundary(newPos)):
					newPos=eachBC.reflect(oldPos[n],dPos)
			eachPoint.modifyPos(newPos)
			n=n+1

#### Test functions ####
def test1RK2(endTime=100.):
	"""A test case for advection RK-2 integrator without any boundaries . Take two vortices of same strength. They should rotate in a circle"""
	
	#Initiate vortex list
	VList=dfn.vortexList()
	
	# Define vortices
	p1=numpy.array([-0.5,0.0])
	p2=numpy.array([0.5,0.0])
	rad=0.5
	VList.addVortex(position=p1, strength=1.0,blobType=0,delta=0.0,traceFlag=1)
	VList.addVortex(position=p2, strength=1.0,blobType=0,delta=0.0,traceFlag=1)
	
	# Initiate time mesh
	startTime=0.0
	timeStep=0.01
	timeMesh=dfn.linearTimeMesh(startTime,timeStep,endTime)
	
	#Start time loop
	for i in range(1,len(timeMesh)):
		print "time=%f"%timeMesh[i]
		dt=timeMesh[i]-timeMesh[i-1]
		advectRK2(dt,toMod=[VList],fieldGens=[VList],vinf=0.0)
	
	#Plot traces of particles
	plt.figure(1)
	plt.title('test case for RK-2 two particle motion in a circle')
	for i in range(VList.nPoints):
		plt.plot(numpy.array(VList.allV[i].traceHist).transpose()[0],numpy.array(VList.allV[i].traceHist).transpose()[1])
	
	#Calculate Error history
	newPos1=VList.allV[0].traceHist
	newPos2=VList.allV[0].traceHist
	err1=[abs(newPos1[i].dot(newPos1[i])-rad**2) for i in range(len(newPos1))]
	err2=[abs(newPos2[i].dot(newPos2[i])-rad**2) for i in range(len(newPos2))]
	
	#Plot errors
	plt.figure(2)
	plt.title('errors with time')
	plt.plot(err1)
	plt.plot(err2)
	return()

def test2RK2(Npoints=30,Npanels=50):
	"""Test RK-2 time integrator with boundary conditions. Simulate a non-viscous flow around a cylinder.
	NPoins being number of tracer points and NPanels being number of panels on the cylinder surface"""
	
	#Define V_Infinite, field generators and modifiable points
	vinf=numpy.array([1.0,0.0])
	fieldGens=[]
	toMod=[]
	
	#Define tracer points
	TList=dfn.traceList()
	lowLim=numpy.array([-5.0,-2.0])
	upLim=numpy.array([-5.0,2.0])
	for i in range(Npoints):
		pos=(lowLim+(float(i)/Npoints*(upLim-lowLim)))
		TList.addTracer(pos,traceFlag=1)
	toMod.append(TList)
	
	#Create cylinder points for boundary conditions
	rad=1.0
	[BCpoints,BCcp,BCNormals,iF,rF]=wallBC.cylBCPoints(rad,Npanels)
		
	#Initiate boundary condnitions
	BCList=[]
	BC=wallBC.wallNPNSBC(BCpoints,BCcp,BCNormals,iF,rF)
	BCList.append(BC)

	# Initiate time mesh
	startTime=0.0
	timeStep=0.08
	endTime=10.0
	timeMesh=dfn.linearTimeMesh(startTime,timeStep,endTime)
	
	#Start time loop
	for i in range(1,len(timeMesh)):
		print "time=%f"%timeMesh[i]
		dt=timeMesh[i]-timeMesh[i-1]
		advectRK2(dt,toMod,fieldGens,vinf,BCList)
	#Plot traces of particles
	plt.figure(1)
	plt.title('test case for RK-2 non-viscous flow around a cylinder')
	plt.plot(numpy.array(BCpoints).transpose()[0],numpy.array(BCpoints).transpose()[1],linewidth=3.0)
	for i in range(TList.nPoints):
		plt.plot(numpy.array(TList.allV[i].traceHist).transpose()[0],numpy.array(TList.allV[i].traceHist).transpose()[1])
	
	return()

def test3RK2(Npanels=50):
	"""Compute the motion of a single point vortex at a distance of 1.5 units from the center of the cylinder and plot its trajectory"""
	
	#Define V_Infinite, field generators and modifiable points
	vinf=numpy.array([0.0,0.0])
	fieldGens=[]
	toMod=[]
	
	#Define vortex points
	VList=dfn.vortexList()
	pos=numpy.array([0.0,1.1])
	VList.addVortex(position=pos,strength=1.0,blobType=1,delta=0.1,traceFlag=1)
	toMod.append(VList)
	fieldGens.append(VList)
	
	#Create cylinder points for boundary conditions
	rad=1.0
	[BCpoints,BCcp,BCNormals,iF,rF]=wallBC.cylBCPoints(rad,Npanels)
		
	#Initiate boundary condnitions
	BCList=[]
	BC=wallBC.wallNPNSBC(BCpoints,BCcp,BCNormals,iF,rF)
	BCList.append(BC)

	# Initiate time mesh
	startTime=0.0
	timeStep=0.08
	endTime=100.0
	timeMesh=dfn.linearTimeMesh(startTime,timeStep,endTime)
	
	#Start time loop
	for i in range(1,len(timeMesh)):
		print "time=%f"%timeMesh[i]
		dt=timeMesh[i]-timeMesh[i-1]
		advectRK2(dt,toMod,fieldGens,vinf,BCList)
	
	#Plot traces of particles
	plt.figure(1)
	plt.title('Assignment 3 Problem 3')
	plt.plot(numpy.array(BCpoints).transpose()[0],numpy.array(BCpoints).transpose()[1],linewidth=3.0)
	for i in range(VList.nPoints):
		plt.plot(numpy.array(VList.allV[i].traceHist).transpose()[0],numpy.array(VList.allV[i].traceHist).transpose()[1])
	plt.show()
	return()

def test4RK2(Npanels=50):
	"""Test RK-2 time integrator with boundary conditions. Simulate a Viscous flow around a cylinder.
	NPoins being number of tracer points and NPanels being number of panels on the cylinder surface"""
	
	#Define parameters
	Re=1000
	delta=(1./Re)**0.5
	lemda=delta*math.pi
	gmin=0.2
	#Define V_Infinite, field generators and modifiable points
	vinf=numpy.array([1.0,0.0])
	fieldGens=[]
	toMod=[]
	
	#Create cylinder points for boundary conditions
	rad=1.0
	[BCpoints,BCcp,BCNormals,BCinFunction,BCreflectFunction]=wallBC.cylBCPoints(rad,Npanels)
	
	# Viscocity
	nu=((vinf.dot(vinf))**0.5)*2*rad/Re
	
	#Initiate boundary condnitions
	BCList=[]
	BC=wallBC.wallNPNSBC(BCpoints,BCcp,BCNormals,BCinFunction,BCreflectFunction)
	BCList.append(BC)

	# Initiate time mesh
	startTime=0.0
	CFL=1.0
	timeStep=CFL*lemda/((vinf.dot(vinf))**0.5)
	endTime=6.0
	timeMesh=dfn.linearTimeMesh(startTime,timeStep,endTime)

	#Start time loop
	for i in range(1,len(timeMesh)):
		print "time=%f"%timeMesh[i]
		dt=timeMesh[i]-timeMesh[i-1]
		
		# Advect particles
		advectRK2(dt,toMod,fieldGens,vinf,BCList)

		# Introduce no slip boundary conditions 
		# Vslips are calculated in the advection step. That is to be changed for better architecture
		[eachBC.applyNSBC(fieldGens,toMod,gmin,delta) for eachBC in BCList]

		#Diffuse all particles using RVM
		diffusion.applyRVM(dt,nu,toMod,BCList)

		# Post Processing: to be done after 5 time steps
		if i%5==0:	
			# Plot Vorticity Particles			
			plots.plotVort(toMod,BCList,timeMesh[i])
			# Plot Velocity Field
			plots.plotVel(fieldGens,vinf,BCList,timeMesh[i])
	return()
