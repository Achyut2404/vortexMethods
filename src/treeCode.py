# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombay

# Tree code algorithm to achieve N.logN speed up
# Incomplete code

import definations as dfn
import numpy
import random
import matplotlib.pyplot as plt
import time
class tree():
	"""A tree structure that contains info aboout various clusters"""
	def __init__(self):
		"""Initiate zeroth level box and decide lowest number of particles in a cluster"""
		self.allBoxes
		self.levelMap 
	def mkTree(self,fieldGens):
		"""Generate tree in a downward direction, based on fieldGens"""
	def findCoefs(self):
		"""Calculates coffecients in downward direction"""
	def findVel(self,pos):
		"""Find velocity induced due to the tree structure at 'pos'; velocity induced by vortex sheets and any other field Generators that are not included in the tree structure, needs to be calculated seperately""" 

def velFieldFMM(pos,fieldGen,vinf = 0.0):
	"""Final function that calculates velocity field at 'pos' positions"""
	
	field=numpy.zeros(pos.shape)
	print "number of Positions = %i" %len(pos)
	
	# Initiate tree
	myTree = tree()
	
	# Create the tree
	myTree.mkTree(fieldGens)

	for i in range(len(pos)):
		vel=vector.zeroValue
		# Seperately calculate velocity is it is linear vortex sheet
		for eachGen in fieldGen:
			if type(eachGen) == linVortList or type(eachGen) == linVortex:
				vel=vel+eachGen.fieldEffect(pos[i])
		# Add vinf
		field[i]=vel+vinf
		
		# Calculate velocity from tree struct
		vel = vel + myTree.findVel(pos)
	return(field)

class clustr():
	"""A rectengular cluster which contains specified number of particles, and stores coefficients"""
	def __init__(self,xmin,xmax,ymin,ymax):
		"""Initiate a cluster with  appropriate xmin, xmax, ymin, ymax values"""
		self.glbId
		self.parentId = []
		self.childrenId = []
		self.levelId
		self.Coeffs
	def findClustrFieldGens(fieldGens):
		"""Find fieldGens that reside inside this cluster"""
		self.fieldGensInsideId
	def findVel(self,pos):
		"""Find velocity due to this cluster on 'pos'"""
		# Check whether 'pos' is sufficiently far away
	def transCoefs(self,allBoxes):
		"""Transfers coefficients from childs to parents"""
	def findCoef(self,fieldGens):
		"""Find coefficients for childless clusters"""
	def isFar(self,pos):
		"""Find whether 'pos' is far from the cluster"""

#### Test functions ####
def randomFieldGens(N = 200):
	"""Generate N random fieldGens in (0,0), (1,1) square box"""
	R1 = dfn.vortexList()
	for i in range(N/4):
		p1 = numpy.array([random.random(),random.random()])
		R1.addVortex(position = p1, strength = 1.0, blobType = 0.0, delta = 0.0)
	R2 = dfn.vortexList()
	for i in range(N/4):
		p1 = numpy.array([random.random(),random.random()])
		R2.addVortex(position = p1, strength = 1.0, blobType = 0.0, delta = 0.0)
	R3 = dfn.vortexList()
	for i in range(N/4):
		p1 = 0.5*numpy.array([random.random(),random.random()])
		R3.addVortex(position = p1, strength = 1.0, blobType = 0.0, delta = 0.0)
	R4 = dfn.vortexList()
	for i in range(N/4):
		p1 = 0.01*numpy.array([random.random(),random.random()])
		R4.addVortex(position = p1, strength = 1.0, blobType = 0.0, delta = 0.0)
	fieldGens = [R1,R2,R3,R4]
	return(fieldGens)

def checkTreeGen():
	"""Checks whether tree generation is proper or not"""
	# Generate random fieldGens and put them in random lists
	fieldGens = randomFieldGens()
	
	# Initiate tree
	myTree = tree()

	# Create the tree
	myTree.mkTree(fieldGens)

	# Plot the tree
	plt.figure()
	# Plot all the field Gen points
	for eachList in fieldGens:
		if type(eachList) != linVortList:
			for eachV in eachList:
				plt.plot(eachV.position[0],eachV.position[1],'ro',markersize = 3.0)
	
	#Plot tree structure
	for eachClustr in myTree.allBoxes:
		xmin = eachClustr.xmin
		xmax = eachClustr.xmax
		ymin = eachClustr.ymin
		ymax = eachClustr.ymax
		plt.plot([xmin,xmin,xmax,xmax,xmin],[ymin,ymax,ymax,ymin,ymin])
	
	plt.show()

def checkTreeCoeffs():
	"""Check whether tree coefficients are properly generated or not. Generates random fieldGens. Find velocities using FMM and regular method, and compares accuracies at different positions"""
	# Generate random fieldGens and put them in random lists
	fieldGens = randomFieldGens()
		
	# Create random 20 positions; where velocity field will be evaluated
	pos = numpy.zeros([20,2])
	for i in range(20):
		pos[i][0] = random.random()	
		pos[i][1] = random.random()
	
	# Find Velocity with regular method
	fieldReg = dfn.velField(pos,fieldGens,vinf = 0.0)
	
	# Find Velocity with FMM
	fieldFMM = velFieldFMM(pos,fieldGens,vinf = 0.0)
	
	# Find and print errors between the two
	for i in range(20):
		err = numpy.norm(fieldReg - fieldFMM)
		print err
		if err > 1e-06:
			print "Error is greater than 1e-06"
	
def checkFMMTime():
	"""Compare the time taken to calculate velocity fields on each other with varying number of fieldGens"""
 	
	Np = [100,300,600,1000,2000]
	tReg = []
	tFMM = []
	for i in Np:
		# Generate random fieldGens and put them in random lists
		fieldGens = randomFieldGens(i)
	
		# Calculate pos matrix
		N = 0
		pos = []
		for eachList in fieldGens:
			N=N+eachList.nPoints
			pos=pos+eachList.posit()
		pos=numpy.array(pos)
	
		# Calculate time for regular velocity calculation
		initT = time.time()
		fieldReg = dfn.velField(pos,fieldGens,vinf = 0.)
		dt = time.time() - initT
		tReg.append(dt)
	
		# Calculate time fot FMM velocity calculation
		initT = time.time()
		fieldFMM = velFieldFMM(pos,fieldGens,vinf = 0.)
		dt = time.time() - initT
		tFMM.append(dt)
	
	# Plot trends
	plt.figure()
	plt.plot(Np,tReg,label = 'Regular')
	plt.plot(Np,tFMM,label = 'FMM')
	plt.legend(loc = 'upper center')
	plt.show()
