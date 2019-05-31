"""
/***********************************************************************************************************\

This Python script is parameters parser for EnergyChaser.pyx library associated with paper:                                                        
 Ruben Tikidji-Hamburyan , Tarek El-Ghazawi , Jason Triplett
    Novel Models of Visual Topographic Map Alignment into Superior Colliculus

 Copyright: Ruben Tikidji-Hamburyan <rath@gwu.edu> Apr.2016 - Sep.2016

\************************************************************************************************************/    
"""
import os,sys,csv
from numpy import *
from numpy import random as rnd
#from multiprocs import multiprocs as mps
import subprocess as sbp
try:
	import cPickle as pkl
except:
	import pickle as pkl

import time
import pyximport; pyximport.install()
from EnergyChaser import chaser
			
########## DEFAULTs HERE ##########
#Size
xsize,  ysize	= 100, 100

#Chemistri
Aaf, Baf		= 60., 90.
#Activity
Bca, Gca		= 11., 20.
Rca, Vca, Sca	= 5., 3.,1.
A2Pca			= 0.625
#Competition
Apr, Bpr, Dpr	= 5., 1., 1.

E12				= 1e5
Nstep			= 150 #number steps per neuron
ParentDir		= ""

#Flags
Knocked			= False

Report			= False
TotalEnergy		= False
Init			= True
StartRec		= False
StopRec			= True
Indicator		= False
Graphs			= False
Log				= True
RunDB			= True
Model			= 3 #"ScaledCor" #"Correlation" or "V1Int"
Norm			= True
ModelID			= time.strftime("%Y%m%d%H%M%S")+"%03d"%rnd.randint(1000)
##################################
	
if __name__ == "__main__":
	for arg in sys.argv[1:]:
		if arg[0] != "/": continue
		if not "="in arg[1:] : continue
		print "Applay: ",arg[1:],
		try:
			exec arg[1:]
		except:
			sys.stderr.write("ERROR in parameter {}\n\n".format(arg))
		print " Done"

if type(Model) is str:
	if   Model == "ScaledCor"   : Model = 2
	elif Model == "Correlation" : Model = 1
	elif Model == "V1Int"       : Model = 3
	elif Model == "V1Int-D"     : Model = 4
	else                        : Model = 3
if Model != 1 and Model != 2 and Model != 3 and Model != 4 : Model = 3

p = chaser(xsize, ysize, Nstep,
			E12 = E12,
			Aaf = Aaf, Baf = Baf, 
			Bca = Bca, Gca = Gca, Rca=Rca, Vca = Vca, Sca=Sca, A2Pca = A2Pca,
			Apr = Apr, Bpr = Bpr, Dpr = Dpr,
			TotalEnergy	= TotalEnergy, Init      = Init,
			StartRec    = StartRec   , StopRec   = StopRec,
			Knocked     = Knocked    , Indicator = Indicator,
			Log         = Log        , Reporting = Report,
			Graphs      = Graphs     , RunDB     = RunDB,
			Model       = Model      , Norm      = Norm,
			ParentDir = ParentDir    , timestemp = ModelID)
