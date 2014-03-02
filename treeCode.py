# Code for 2-dimensional vortex methods
# Written by: Achyut Panchal
# Aerospace Engineering, Indian Institute of Technology Bombay
# Inspired by lectures from Prof. Prabhu Ramachandran, IIT Bombay

# Tree code algorithm to achieve N.logN speed up
# Incomplete code

import definations as dfn
import numpy
def velField(pos,fieldGen,vinf=0.0):
	"""Define velocity field at set of points due to vinf and 
	fieldGen( a list containing various velocity generator lists)
	Based on tree code clustering approach"""
	field=numpy.zeros(pos.shape)
	print "number of Positions = %i" %len(pos)
	for i in range(len(pos)):
		vel=vector.zeroValue
		for eachGen in fieldGen:
			vel=vel+eachGen.fieldEffect(pos[i])
		field[i]=vel+vinf
	return(field)
def mkCluster(fieldGen):
	"""Creates a box,tree based cluster of particles.
	Currently, this only works for vortices/Blobs.
	Need to be coded for Sheets"""
	clusters=[]
	# Initiate level 0 Box
	xlim=fieldGens.limx
	ylim=fieldGens.limy
	clusters.append([clustBox(0,xlim,ylim,points,fieldGen)])
	# If number of particles in the box is large enough divide it into subBoxes
	for eachClust in clusters[-1]:
		if !eachClust.isFull:
			clusters.append(eachClust.subClusts)
	return(clusters)

class clustBox():
	""" Rectengular cluster"""
	def __init__(self,parent,xlim,ylim,stren):
		""" Initiate a rectengular cluster.
		parent is the id of parent box, xlim and ylim
		specify the boundaries of rectangle.
		stren is the total strenght of a particular box"""
		[self.x1,self.x2]=xlim
		[self.y1,self.y2]=ylim
		self.centr=numpy.array([numpy.average(xlim),numpy.average(ylim)])
		self.stren=stren
		self.parent=numpy.array([parent])
		self.child = numpy.array([]).astype(int)
	def hasChild(self):
		if len(self.child)==0:
			return(False)
		else:
			return(True)
	def subClusts
	def vel(pos)
