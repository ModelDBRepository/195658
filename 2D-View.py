import os,sys,csv
from numpy import *
from numpy import random as rnd
import scipy as sp
import scipy.fftpack as spfft
import scipy.signal as spsignal
import matplotlib
matplotlib.rcParams["savefig.directory"] = ""
from matplotlib.pyplot import *
try:
	import cPickle as pkl
except:
	import pickle as pkl

import gzip
cmap = matplotlib.cm.get_cmap('jet')
#cmap = matplotlib.cm.get_cmap('plasma')
#cmap = matplotlib.cm.get_cmap('autumn')
#cmap = matplotlib.cm.get_cmap('gist_rainbow')



rec=[]
recname=[]

pname=["recs","FigName","Nraaw","CS","R"]
for p in pname:
	exec p+"=None"
for arg in sys.argv[2:]:
	sop = arg.split("=")
	if len(sop) != 2:
		print "Cannot find = in parameter",arg
		continue
	if sop[0][0] != "/" : 
		print "Cannot find / in parameter name",arg
		continue
	if sop[0][1:] in pname:
		exec arg[1:]
	else:
		continue
		#print "Cannot read parameter",arg
		#print "STOP!!!!"
		#exit(0)
		

if CS == None:	CS = [25,50,75]
elif CS == "ALL" or CS == "all":CS = range(100)

if Nraaw == None :Nraaw = 3

#DB>>
print "File  = ",sys.argv[1]
print "R     = ",R
print "Recs  = ",recs
print "CS    = ",CS
print "Nraaw = ",Nraaw
print "Fig   = ",FigName
#<<DB


ml = ['*', '.', ',', '+', 'o', 'D', "-"]
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
		k=[]
		for x in fld[4:]:
			#if ri == 9: print x,
			exec "m=["+x+"]"
			#if ri == 9: print "=>",m
			k.append(m)
		if len(k) != n:
			print "\nError sise of rec #%d isn't fit db marker: %d vs %d"%(ri,len(k),n)
			#exit(0)
		else:
			rec.append(k)
			recname.append(fld[0])
			print "Done"

def plotit(REC):
	REC = array(REC)
	samples = [ (sx,sy) for sx in CS for sy in CS ]
	
	nscale = 1./float(len(samples))
	cpl=[]
	for nidx,ncor in enumerate(samples):
		sx,sy = ncor
		c = cmap(float(nidx)*nscale)
		cpl.append(c)
		idx = where( (REC[:,2] == sx)*(REC[:,3] == sy) )[0]
		plot(REC[idx,0],REC[idx,1],'.',c=c,mfc=c,mec=c,ms=5)#REC[idx,4])
		if R != None:
			for xp in range(sx-R,sx+R):
				for yp in range(sy-R,sy+R):
					if sqrt(float(xp-sx)**2+float(yp-sy)**2) > float(R): continue
					idx = where( (REC[:,2] == xp)*(REC[:,3] == yp) )[0]
					plot(REC[idx,0],REC[idx,1],'.',c=c,mfc=c,mec=c,ms=5)#REC[idx,4])
					
	for c,nor in zip(cpl,array(samples)):
		plot(nor[0],nor[1],"x",ms=8,mew=4, mfc=c,mec=c)
	#plot(samples[:,0],samples[:,1],"kx",ms=5,mew=2, mfc='k',mec='k')
			
			

if FigName != None:
	figure(1,figsize=(17,17))
suptitle("Parameters ID:"+sys.argv[1][:-7])
n=len(rec)
if  n == 0 :
	print "Nothing to plot"
	exit(0)
elif n == 1:
	plotit(rec[0])
	ylabel(recname[0])
elif n <= Nraaw:
	ax = None
	for i,r in enumerate(rec):
		if ax == None:
			ax = subplot(1,n,i+1)
		else:
			subplot(1,n,i+1,sharex=ax,sharey=ax)
		plotit(r)
		ylabel(recname[i])
elif n > Nraaw:
	ax = None
	nr = n/Nraaw +1 if n%Nraaw else n/Nraaw
	for i,r in enumerate(rec):
		print n,Nraaw,"-->",nr,Nraaw,i+1
		if ax == None:
			ax = subplot(nr,Nraaw,i+1)
		else:
			subplot(nr,Nraaw,i+1,sharex=ax,sharey=ax)
		plotit(r)
		title(recname[i])
xlim(0,100)
ylim(0,100)
if FigName != None:
#	savefig(FigName+".svg")
	savefig(FigName+".jpg")
else:
	show()
	
