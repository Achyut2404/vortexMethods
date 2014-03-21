# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombay

# Plots and saves figures based on testcase-4 RK-2 requirement
# Change as per requirement

import numpy
import matplotlib.pyplot as plt
import definations

## Plot Velocity Field
def plotVel(fieldGens,vinf,BCList,tim):
	"""Plots velocity vector plot
	Designed to be used with tInt test case 4.
	Change for different cases as required"""
	
	plt.figure()
	[plt.plot(numpy.array(eachBC.points).transpose()[0],numpy.array(eachBC.points).transpose()[1],linewidth=3.0) for eachBC in BCList]
	plt.title('Test Case for RK-2 Viscous Flow Around a Cylinder \n Velocity Plot: Time=%f'%tim)
	X,Y = numpy.meshgrid(numpy.arange(0.0,5.0,0.2),numpy.arange(0.0,2.5,0.1))
	plt.xlim([-0.1,5.0])
	plt.ylim([-0.1,2.5])
	pos=numpy.array([X.flatten(),Y.flatten()]).transpose()
	
	# Apply no penetration boundary condition in order to find velocity field
	[eachBC.findVcp(fieldGens,vinf) for eachBC in BCList]
	[eachBC.applyNPBC(fieldGens) for eachBC in BCList]
	Vel = definations.velField(pos,fieldGens,vinf)
	[eachBC.closeNPBC(fieldGens) for eachBC in BCList]
	
	# If any point is inside the boundary give it a zero velocity
	for j in range(len(pos)):
		for eachBC in BCList:
			if eachBC.inBoundary(pos[j]):
				Vel[j]=numpy.zeros(2)
	u=Vel.transpose()[0].reshape(X.shape)
	v=Vel.transpose()[1].reshape(X.shape)
	plt.quiver(X,Y,u,v)
	plt.savefig('Velocity_%f.png'%tim)
	return()

#Plot Vorticities of points
def plotVort(toMod,BCList,tim):
	"""Plots velocity vector plot
	Designed to be used with tInt test case 4.
	Change for different cases as required"""
	
	plt.figure()
	plt.axis('equal')
	plt.title('Test Case for RK-2 Viscous Flow Around a Cylinder \n Vorticity Plot: Time=%f'%tim)
	[plt.plot(numpy.array(eachBC.points).transpose()[0],numpy.array(eachBC.points).transpose()[1],linewidth=3.0) for eachBC in BCList]
	for eachList in toMod:
		for eachV in eachList.allV:
			if numpy.sign(eachV.strength)==1:
				plt.plot(eachV.position[0],eachV.position[1],'ro',markersize=3.0)
			else:
				plt.plot(eachV.position[0],eachV.position[1],'bo',markersize=3.0)
	plt.savefig('Vorticity%f.png'%tim)
