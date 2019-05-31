import os,sys,csv
from numpy import *
from numpy import random as rnd
import scipy as sp
import scipy.fftpack as spfft
import scipy.signal as spsignal
from scipy import stats
import matplotlib
matplotlib.rcParams["savefig.directory"] = ""
from matplotlib.pyplot import *
try:
	import cPickle as pkl
except:
	import pickle as pkl

import gzip
#cmap = matplotlib.cm.get_cmap('jet')
#cmap = matplotlib.cm.get_cmap('plasma')
#cmap = matplotlib.cm.get_cmap('autumn')
cmap = matplotlib.cm.get_cmap('rainbow')



rec=[]
recname=[]
##==##
D1      = None
recs    = None
Xax     = None
Nraw    = None
FigName = None
Sml     = None
for arg in sys.argv[2:]:
	sub = arg.split("=")
	if len(sub) != 2: continue
	if   sub[0] == "/D1"     : exec "D1     = "+sub[1]
	elif sub[0] == "/recs"   : exec "recs   = "+sub[1]
	elif sub[0] == "/Nraw"   : exec "Nraw   = "+sub[1]
	elif sub[0] == "/Xax"    : exec "Xax    = "+sub[1]
	elif sub[0] == "/FigName": exec "FigName= "+sub[1]
	elif sub[0] == "/Sml"    : exec "Sml    = "+sub[1]
	else:
		print "Cannot recognized key %s"%sub[0]
		print "USAGE:"
		exit(1)	

if D1   == None: D1="NT"
elif type(D1) is int or type(D1) is float:
	if float(D1) > 0.5: D1="DV"
	else: D1="NT"

if Xax  == None: Xax=(0,100)
if Nraw == None: Nraw=3


#DB>>
print "File  = ",sys.argv[1]
print "D1    = ",D1
print "Recs  = ",recs
print "Xax   = ",Xax
print "Nraw  = ",Nraw
print "Fig   = ",FigName
#<<DB

#ax = linspace(0,100,401)
ax = linspace(0,100,101)

#with open(sys.argv[1]) as fd:
with gzip.open(sys.argv[1], 'rb') as fd:
	for ri,r in enumerate(fd.readlines()):
		print "Read Record #",ri,"...",
		if recs != None:
			if type(recs) is int:
				if ri != recs:  print "Skip"; continue
			elif type(recs) is list or type(recs) is tuple:
				if not ri in recs: print "Skip"; continue
		fld = r[:-1].split(":")
		n = int(fld[3])
		k=[[] for x in xrange(Xax[0],Xax[1])]
		for x in fld[4:]:
			exec "m=["+x+"]"
			if D1 == 'NT':
				if m[2] >= Xax[1] or m[2] < Xax[0]:continue
				else:
					for q in xrange(m[4]):
						k[m[2]-Xax[0]].append(m[0])
			if D1 == 'DV' :
				if m[3] >= Xax[1] or m[3] < Xax[0]:continue
				else:
					for q in xrange(m[4]):
						k[m[3]-Xax[0]].append(m[1])
		for kid in xrange(Xax[1]-Xax[0]):
			k[kid] = array(k[kid]).astype(float)
			kde = stats.gaussian_kde(k[kid], bw_method='silverman')
			k[kid] = kde(ax)
			#print k[kid]
			#exit(0)
			#k,ax = histogram(k,bins=100,range=(0,100))
			k[kid] = k[kid].astype(float)/float(max(k[kid]))
		rec.append(array(k))
		recname.append(fld[0])
		print "Done"

def plotit(rec):
	hotmap=pcolormesh(Xax,ax,rec.T, cmap=cmap)	
def plotsome(rec):
	for idx,dst in enumerate(rec[Sml/2::Sml,:]):
		plot(ax,dst,c=cmap((float(idx)/float(smlXax.shape[0]-1))) )#,label="{}".format(smlXax[idx]))

if FigName != None:
	figure("Parameters ID:"+sys.argv[1][:-7]+D1+" View",figsize=(17,17))
				
suptitle("Parameters ID:"+sys.argv[1][:-7])

n=len(rec)
Xax = array(range(Xax[0],Xax[1]))
if  n == 0 :
	print "Nothing to plot"
	exit(0)
elif n == 1:
	plotit(rec[0])
	ylabel(recname[0])
elif n <= Nraw:
	Sax = None
	for i,r in enumerate(rec):
		if Sax == None:
			Sax = subplot(1,n,i+1)
		else:
			subplot(1,n,i+1,sharex=Sax,sharey=Sax)
		plotit(r)
		ylabel(recname[i])
elif n > Nraw:
	Sax = None
	nr = n/Nraw +1 if n%Nraw else n/Nraw
	for i,r in enumerate(rec):
		print n,Nraw,"-->",nr,Nraw,i+1
		if Sax == None:
			Sax = subplot(nr,Nraw,i+1)
		else:
			subplot(nr,Nraw,i+1,sharex=Sax,sharey=Sax)
		plotit(r)
		title(recname[i])

#hotmap=pcolormesh(array(Xax),(ax[1:]+ax[:-1])/2,rec.T)
xlim(0,100)
ylim(0,100)

if Sml != None:
	figure(2)
	smlXax = Xax[Sml/2::Sml]
	if n == 1:
		plotsome(rec[0])
		ylabel(recname[0])
	elif n <= Nraw:
		Sax = None
		for i,r in enumerate(rec):
			if Sax == None:
				Sax = subplot(1,n,i+1)
			else:
				subplot(1,n,i+1,sharex=Sax,sharey=Sax)
			plotsome(r)
			ylabel(recname[i])
	elif n > Nraw:
		Sax = None
		nr = n/Nraw +1 if n%Nraw else n/Nraw
		for i,r in enumerate(rec):
			print n,Nraw,"-->",nr,Nraw,i+1
			if Sax == None:
				Sax = subplot(nr,Nraw,i+1)
			else:
				subplot(nr,Nraw,i+1,sharex=Sax,sharey=Sax)
			plotsome(r)
			title(recname[i])
	

if FigName != None:
	#savefig(FigName+".svg")
	savefig(FigName+".jpg")
else:
	show()
