"""
/***********************************************************************************************************\

This Cython script associated with paper:                                                        
 Ruben Tikidji-Hamburyan , Tarek El-Ghazawi , Jason Triplett
    Novel Models of Visual Topographic Map Alignment into Superior Colliculus

 Copyright: Ruben Tikidji-Hamburyan <rath@gwu.edu> Apr.2016 - Sep.2016

\************************************************************************************************************/    
"""
from __future__ import division
import  numpy as np
cimport numpy as np
import sys,os
from libc.math cimport sqrt,exp,floor
from libc.stdlib cimport rand,srand,RAND_MAX,malloc
from libc.stdio cimport *
import time


PTYPE = np.uint
DTYPE = np.float64
ctypedef np.uint_t    PTYPE_T
ctypedef np.float64_t DTYPE_T

cdef extern from "stdio.h":
	#FILE * fopen ( const char * filename, const char * mode )
	FILE *fopen(const char *, const char *)
	#int fclose ( FILE * stream )
	int fclose(FILE *)
	#int fprintf(FILE *stream, const char *format, ...)
	int fprintf(FILE *, const char *, ...)
	#void * malloc(
	
#Afinity
cdef double cAaf = 20.,	cBaf = 30.
#Competition
cdef double cApr = 5.,	cBpr = 1.,	cDpr = 1.
#Activity
cdef double cBca = 11.,	cGca = 200., cRca = 3, cVca=11, cSca=1., cA2Pca=1.5 #???

cdef unsigned int xsize, ysize, Nstep, ic_sizeX, ic_sizeY
cdef double xscl = 0, yscl = 0

cdef unsigned int cNorm, Don=1, Den3sigma=0

cdef double *pDf=NULL, *pLf=NULL


cdef inline double FdEpr(PTYPE_T na, unsigned long nd):
	return cBpr*na**2-cApr*sqrt(na) + cDpr*nd**2


cdef inline double FdEaf(long xs, long ys, long xl, long yl, PTYPE_T p ):
	return (   cAaf*exp(float(xs)*xscl-1.)*exp(-float(xl)*xscl   ) \
		     - cBaf*exp(float(ys)*yscl-1.)*exp( float(yl)*yscl-1.) )*float(p)

cdef inline double Lfactor(long xs, long ys, long xl, long yl, double sigma):
	return         exp( - sqrt( float( (xs-xl)**2 + (ys-yl)**2 ) ) / sigma )

cdef inline double Dfactor(long xs, long ys, long xl, long yl, double sigma):
	if Don: return exp( -       float( (xs-xl)**2 + (ys-yl)**2 ) / (2. * sigma**2) )
	else:   return 1.0

cdef inline double FdEcaW(long xs, long ys, long xl, long yl, float add, np.ndarray[PTYPE_T, ndim=4] p):
	cdef double norm = Dfactor(xs,ys,xl,yl,cVca) * add,  res = norm, df
	cdef unsigned int xls, yls
	for xls in xrange(xsize):
		for yls in xrange(ysize):
			if p[xs,ys,xls,yls] == 0: continue
			df    = Dfactor(xs,ys,xls,yls,cVca) * float(p[xs,ys,xls,yls])
			res  += Lfactor(xl,yl,xls,yls,cBca) * df
			norm += df
	if norm > 0. and cNorm:
		res /= norm
	res += cRca*exp( - sqrt( (float(xs-xl)/cBca/cSca)**2 + (float(ys-yl)/cBca)**2 ) )
	return - cGca * res

cdef inline double FdEcaM(long xs, long ys, long xl, long yl, float add, np.ndarray[PTYPE_T, ndim=4] p):
	cdef double norm = Dfactor(xs,ys,xl,yl,cVca) * add,  res = norm, df
	cdef unsigned int xls, yls
	for xls in xrange(xsize):
		for yls in xrange(ysize):
			if p[xs,ys,xls,yls] == 0: continue
			df    = Dfactor(xs,ys,xls,yls,cVca) * float(p[xs,ys,xls,yls])
			res  += Lfactor(xl,yl,xls,yls,cBca) * df 
			norm += df
	if norm > 0. and cNorm:
		res /= norm
	res += cRca * (\
		     cA2Pca  * exp( - sqrt( (float(xs-xl*0.5        )/cBca/cSca)**2 + (float(ys-yl)/cBca)**2 ) )\
		             + exp( - sqrt( (float(xs-xl*0.5-xsize/2)/cBca/cSca)**2 + (float(ys-yl)/cBca)**2 ) )\
				  ) 
	return - cGca * res


cdef inline double FdEcaW_JC(long xs, long ys, long xl, long yl, float add, np.ndarray[PTYPE_T, ndim=4] p):
	return - cGca * exp( - sqrt(float((xs-xl)**2+(ys-yl)**2))/cBca /cSca )

cdef inline double FdEcaM_JC(long xs, long ys, long xl, long yl, float add, np.ndarray[PTYPE_T, ndim=4] p):
	return - cGca * (\
			 cA2Pca  * exp( - sqrt( (float(xs-xl*0.5        )/cBca/cSca)**2 + (float(ys-yl)/cBca)**2 ) )\
		             + exp( - sqrt( (float(xs-xl*0.5-xsize/2)/cBca/cSca)**2 + (float(ys-yl)/cBca)**2 ) )\
					 )

cdef inline double FdEcaW_SC(long xs, long ys, long xl, long yl, float add, np.ndarray[PTYPE_T, ndim=4] p):
	cdef double res = 0.
	cdef unsigned int xls, yls
	for xls in xrange(xsize):
		for yls in xrange(ysize):
			res += exp( - sqrt( (float(xs-xls)/cBca/cSca)**2 + (float(ys-yls)/cBca)**2 ) ) * Dfactor(xs,ys,xls,yls,cRca)
	return - cGca * res

cdef inline double FdEcaM_SC(long xs, long ys, long xl, long yl, float add, np.ndarray[PTYPE_T, ndim=4] p):
	cdef double res = 0.
	cdef unsigned int xls, yls
	for xls in xrange(xsize):
		for yls in xrange(ysize):
			res +=     cA2Pca  * exp( - sqrt( (float(xs-xl*0.5        )/cBca/cSca)**2 + (float(ys-yl)/cBca)**2 ) ) * Dfactor(xs,ys,int(xl*0.5        ),yl,cRca)
			res +=               exp( - sqrt( (float(xs-xl*0.5-xsize/2)/cBca/cSca)**2 + (float(ys-yl)/cBca)**2 ) ) * Dfactor(xs,ys,int(xl*0.5-xsize/2),yl,cRca)
	return - cGca * res 

cdef inline double random():
	return float(rand())/float(RAND_MAX)

cdef inline unsigned int randint( unsigned int maxi):
	cdef unsigned int x = int( floor( float(rand()) * float(maxi) / float(RAND_MAX) ) )
	if x >= maxi: return maxi-1
	if x <   0  : return 0
	return x
	
def chaser(ixsize, iysize, iNstep,
			E12 = 1e5,
			Aaf = 20., Baf = 30., 
			Bca = 11., Gca = 200., Rca=3., Vca=11., Sca=1., A2Pca=1.5, 
			Apr =  5., Bpr =   1., Dpr = 1.,
			TotalEnergy	= False, Init      = 'RANDOM',
			StartRec    = False, StopRec   = True,
			Knocked     = False, Indicator = False,
			Log         = True,  Reporting = False,
			Graphs      = False, RunDB     = True,
			Model       = 3    , Norm      = True,
			ParentDir	= ""   , timestemp = "" ):
	global xscl, yscl
	global cAaf,cBaf,cApr,cBpr,cDpr,cBca,cGca,cRca,cVca,cSca,cA2Pca
	global xsize, ysize, Nstep, ic_sizeX, ic_sizeY
	global cNorm, Don, Den3sigma
	global pDf,pLf
	
	cAaf,cBaf,cApr,cBpr,cDpr,cBca,cGca,cRca,cVca,cSca, cA2Pca = Aaf,Baf,Apr,Bpr,Dpr,Bca,Gca,Rca,Vca,Sca,A2Pca
	xsize, ysize, Nstep = ixsize, iysize, iNstep
	cNorm = Norm
	cdef unsigned int fTotalEnergy = int(TotalEnergy),\
					  fIndicator   = int(Indicator),\
					  fLog         = int(Log),\
					  fReporting   = int(Reporting),\
					  nGraphs      = 0
	
	cdef double cE12 = E12

	timestemp_t = timestemp.encode("UTF-8")
	cdef char *timestemp_c = timestemp_t
	filenameprefix = ParentDir + timestemp
	
	xscl = 1./float(xsize)
	yscl = 1./float(ysize)

	cdef FILE *log = NULL
	logfilename_t = filenameprefix+"-log.csv"
	logfilename_t = logfilename_t.encode("UTF-8")
	cdef char *logfilename = logfilename_t

	cdef FILE *recdb = NULL
	redbfilename_t = filenameprefix+"-rec.db"
	redbfilename_t = redbfilename_t.encode("UTF-8")
	cdef char *redbfilename = redbfilename_t

	cdef long xs, ys, xl, yl, Niter, itr_c = 0, npos=0
	cdef np.ndarray[long, ndim=1] pXs, pYs, pXl, pYl
	cdef double Epr=0., Eaf=0., Eca=0., E=0.
	
	if Model == 5:
		Den3sigma = int( np.ceil(cVca*3) )
		pDf = <double *>malloc(Den3sigma**2*sizeof(double))
		for xl in xrange(Den3sigma):
			for yl in xrange(Den3sigma):
				pDf[xl+yl*Den3sigma] = exp( -       float( xl**2 + yl**2 )   / (2. * cVca**2) )
		pLf = <double *>malloc(xsize*ysize*sizeof(double))
		for xl in xrange(xsize):
			for yl in xrange(ysize):
				pLf[xl+yl*xsize    ] = exp( - sqrt( float( xl**2 + yl**2 ) ) /       cBca     )

	#Backword compability
	if type(Init) is bool:
		Init = "RANDOM" if Init else None
	
	if Graphs:
		if type(Graphs) is bool:
			nGraphs = 20
		elif type(Graphs) is float:
			nGraphs = int(xsize*ysize*Graphs)
		elif type(Graphs) is int:
			nGraphs = int(Graphs)
	
	if RunDB:
		if not type(RunDB) is str:
			RunDB = timestemp + "-rundb.csv"		
	print "#############################"
	print "#       :MODEL ID:          #"
	print "  > Time Stamp    :",timestemp
	print "#      :PARAMETERS:         #"
	print "  > Aaf           :",cAaf
	print "  > Baf           :",cBaf
	print "  > Apr           :",cApr
	print "  > Bpr           :",cBpr
	print "  > Dpr           :",cDpr
	print "  > Bca           :",cBca
	print "  > Gca           :",cGca
	print "  > Rca           :",cRca
	print "  > Vca           :",cVca
	print "  > Sca           :",cSca
	print "  > A2Pca         :",cA2Pca
	print "  > XxY           :",xsize,"x",ysize
	print "  > Nstep         :",Nstep
	print "  > E1/2          :",cE12
	print "#         :FLAGS:           #"
	print "  > Model         :",Model,
	if   Model == 1: print "(Just Correlation)"
	elif Model == 2: print "(Sclaed Correlation)"
	elif Model == 3: print "(V1 and Retinal Inputs Interation)"
	elif Model == 4: print "(V1 and Retinal Inputs Interation withput convolution)"
	print "  > Normalization :",bool(cNorm)
	print "  > Knocked       :",Knocked
	print "  > TotalEnergy   :",bool(fTotalEnergy)
	print "  > Init          :",Init
	print "  > StartRec      :",StartRec
	print "  > StopRec       :",StopRec
	print "  > Log           :",bool(fLog)
	print "  > Reporting     :",bool(fReporting)
	if fReporting:
		print "    > every iter  :",fReporting
	print "  > Graphs        :",bool(nGraphs)
	if nGraphs:
		print "    > rendom cells:",nGraphs
	print "  > Indicator     :",bool(fIndicator)
	print "#         :PATHS:           #"
	print "  > ParentDir     :",ParentDir
	print "  > Prefix        :",filenameprefix
	print "#         :FILES:           #"
	if fLog:
		print "  > Log           :",logfilename
	if StartRec or StopRec or bool(fReporting):
		print "  > Record        :",redbfilename
	print "  > RunDB         :",RunDB
	print "#############################\n"
	
	
	if type(RunDB) is str:
		with open(RunDB,"a") as fd:
			#TIMESTEMP,xsize,ysize,Aaf,Baf,Bca,Gca,Rca,Vca,Sca,A2Pca,Apr,Bpr,Dpr,E12,Nstep,Knocked,TotalEnergy,Init,Record,Log,StartRec,StopRec,Report,Indicator,Model,Normalization
			fd.write("%s,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"%(\
				timestemp,xsize,ysize,\
				Aaf,Baf,Bca,Gca,Rca,Vca,Sca,A2Pca,Apr, Bpr, Dpr,\
				E12,Nstep,\
				"Y" if Knocked else "N",
				"Y" if bool(fTotalEnergy) else "N",
				Init,
				filenameprefix+"-rec.db" if StartRec or StopRec or bool(fReporting) else "N",
				filenameprefix+"-log.csv" if bool(fLog) else "N",
				
				"Y" if StartRec else "N",
				"Y" if StopRec else "N",
				str(fReporting) if bool(fReporting) else "N",
				"Y" if bool(fIndicator) else "N",
				"JC" if Model == 1 else "SC" if Model == 2 else "INT" if Model == 3 else "INT-D" if Model == 4 else "DI",
				"Y" if bool(cNorm) else 'N' ) )

	### Squared distances
	#cBca,cSca  = cBca**2,cSca**2
	
	print "#############################"
	print "# Create synaptic densities #"
	cdef np.ndarray[PTYPE_T, ndim=2] ax = np.zeros((xsize,  ysize),dtype=PTYPE)
	cdef np.ndarray[PTYPE_T, ndim=2] dn = np.zeros((xsize,  ysize),dtype=PTYPE)
	cdef np.ndarray[PTYPE_T, ndim=4] p  = np.zeros([xsize,ysize,xsize,ysize], dtype=PTYPE)
	print "#           DONE            #"
	print "#############################\n"

	if Knocked:
		if Model == 1:
			FdEca = FdEcaM_JC
		elif Model == 2:
			FdEca = FdEcaM_SC
		elif Model == 3:
			FdEca = FdEcaM
		elif Model == 4:
			FdEca = FdEcaM
			Don = 0
		else:
			FdEca = FdEcaM
	else:
		if Model == 1:
			FdEca = FdEcaW_JC
		elif Model == 2:
			FdEca = FdEcaW_SC
		elif Model == 3:
			FdEca = FdEcaW
		elif Model == 4:
			FdEca = FdEcaW
			Don = 0
		else:
			FdEca = FdEcaW
	
	srand( np.random.randint( RAND_MAX ) )
		
	if Init == "RANDOM":
		print "#############################"
		print "#   Random Initialization   #"
		Niter = int( xsize * ysize * 50) # half of population
		for i in xrange(Niter):
			xs = randint(xsize)
			ys = randint(ysize)
			xl = randint(xsize)
			yl = randint(ysize)
			ax[xs,ys] += 1
			dn[xl,yl] += 1
			p[xs,ys,xl,yl] += 1
		print "#           DONE            #"
		print "#############################\n"
	elif Init == None:
		pass
	else:
		if type(Init) is str:
			Init = (Init,-1)
		if len(Init) < 2: Init = (Init, -1)
		print "#############################"
		print "# Init conditions from File #"
		print "  > File          :", Init[0]
		print "  > Fentch record :", Init[1],type(Init[1])
		lastrec = None
		with open(Init[0],"r") as fd:
			for nl,l in enumerate(fd.readlines()):
				if len(l) <= 5: continue
				if type(Init[1]) is str:
					if l.split(":")[0] == Init[1] : lastrec = l
				elif type(Init[1]) is int:					
					if Init[1] >= 0 and nl == Init[1] : lastrec = l
					elif Init[1] < 0: lastrec = l
		if lastrec == None:
			print " ERROR: Cannon find required record in ",Init, "... Skip"
			exit(1)
		else:
			lastrec = lastrec.split(":")
			if len(lastrec) < 5:
				print " ERROR: Last line in ",Init, "Has not records ....Skip"
				exit(1)
			else:
				for xs,ys,xl,yl,pn in [ ( int(xyxyp) for xyxyp in r.split(",") ) for r in lastrec[4:] ]:
					if xs >= xsize or ys >= ysize or xl >= xsize or yl >= ysize:
						print " ERROR: dementions {} doesn't match ....Skip".format(xs,ys,xl,yl,p)
						continue
					p[xs,ys,xl,yl] = pn				
		print "#           DONE            #"
		print "#############################\n"
	
	if Init != None and bool(fTotalEnergy):
		print "#############################"
		print "# Calculate initial Energy  #"
		for a,d in zip(ax,dn):
			for na,nd in zip(a,d):
				Epr += FdEpr(na, nd)
		#NP
		#Epr = np.sum(Bpr*ax**2-Apr*np.sqrt(ax)) + np.sum(Dpr*dn**2)
		pXs, pYs, pXl, pYl = np.where( p  > 0. )
		for xs, ys, xl, yl in zip(pXs, pYs, pXl, pYl):
			Eaf += FdEaf(xs, ys, xl, yl, p[xs, ys, xl, yl])
			Eca += FdEca(xs, ys, xl, yl, 0., p) * p[xs, ys, xl, yl]
			#DB>>
#			import time
#			tb= time.time()
#			print "TB:",tb
#			Eca += FdEca(xs, ys, xl, yl, 0., p) * p[xs, ys, xl, yl]
#			te= time.time()
#			print te-tb
#			exit(0)
			#<<DB
		#NP
#		Eaf = sum( ( cAaf*np.exp(1.-pXs.astype(float)*xscl)*np.exp(1.-pXl.astype(float)*xscl)\
#			       - cBaf*np.exp(1.-pYs.astype(float)*yscl)*np.exp(1.-pYl.astype(float)*yscl))*p[pXs,pYs,pXl,pYl] )
#		Eca = sum( -cGca * np.exp( - np.sqrt((pXs-pXl)**2+(pYs-pYl)**2)/cBca )*p[pXs, pYs, pXl, pYl] )
		print "   > Chemoaffinity          :", Eaf
		print "   > Competition            :", Epr
		print "   > Activity               :", Eca
		E   = Eaf + Eca + Epr
		print "   > Total                  :", E
		print "#############################\n"
		
		#DB>>
#		exit(0)
		#<<DB
	if Graphs and Init != None:
		from matplotlib import pyplot as plt
		from matplotlib import cm
		cmap = cm.get_cmap('jet')
		Gax = plt.subplot(121)
		Gcells=[ [randint(xsize),randint(ysize),cmap(float(i)/float(nGraphs))] for i in xrange(nGraphs)]
		for xs,ys,c in Gcells:
			pXl, pYl = np.where( p[xs,ys,:,:].astype(int)  > 0 )
			for xl,yl in zip(pXl, pYl):
				plt.plot([xs,xl],[ys,yl],"-*",c=c,mfc=c,mec=c,lw=int(p[xs,ys,xl,yl]) )

	if fLog:
		log = fopen(logfilename,"w")
		if log == NULL:
			sys.stderr.write("Cannot open Log file: '%s' for writing" %(filenameprefix+"-log.csv"))
			raise
		fprintf(log,"iter,+/-,xs,ys,xl,yl,dEpr,dEaf,dEca,dE,E,P0\n")
	
	if StartRec or StopRec or bool(fReporting):
		recdb = fopen(redbfilename,"w")
		if recdb == NULL:
			sys.stderr.write("Cannot open Record file: '%s' for writing" %(filenameprefix+"-rec.db"))
			raise
	
	if StartRec:
		fprintf(recdb,"%s:%010d",timestemp_c,itr_c)
		fprintf(recdb,":%g",E)
		pXs, pYs, pXl, pYl = np.where( p.astype(int)  > 0 )
		npos = int(pXs.shape[0])
		fprintf(recdb,":%d",npos)
		for xs,ys,xl,yl in zip(pXs, pYs, pXl, pYl):
			fprintf(recdb,":%d,%d,%d,%d,%d",xs,ys,xl,yl,p[xs,ys,xl,yl])
		fprintf(recdb,"\n")
	
	
	cdef double P0, dEpr, dEaf, dEca, dE
	Niter = xsize*ysize*Nstep
	#DB>>
	#Niter = Nstep
	#<<DB
	for itr in xrange(Niter):
		itr_c = itr
		# Add new synapse
		xs, ys, xl, yl  = randint(xsize),randint(ysize),randint(xsize),randint(ysize)
		dEpr = FdEpr(ax[xs,ys]+1,dn[xl,yl]+1) - FdEpr(ax[xs,ys],dn[xl,yl])
		dEaf = FdEaf(xs,ys,xl,yl,1)
		dEca = FdEca(xs,ys,xl,yl,1.,p)
		dE   = dEpr + dEaf + dEca
		if dE > 100:
			P0 = 0.
		elif dE < -100:
			P0 = 1.
		else:
			P0 = 1./(1.+exp(4.*dE) )
		if TotalEnergy:
			P0 *= E/(cE12+E)
		Pplus = False
		if np.random.random() < P0:
			Pplus = True
			p[xs,ys,xl,yl] += 1
			ax[xs,ys] += 1
			dn[xl,yl] += 1
			Epr += dEpr
			Eaf += dEaf
			Eca += dEca
			E   += dE
			if fLog:
				fprintf(log,"%d,+,%d,%d,%d,%d,%g,%g,%g,%g,%g,%g\n",itr_c,xs,ys,xl,yl,dEpr,dEaf,dEca,dEpr + dEaf + dEca,E,P0)
			
		if E <= 0. and fTotalEnergy: break
		if not (itr%1000) and fIndicator: 
			sys.stderr.write("%d (+) %g  %g  [%g,%g,%g]  %g ||"%(itr,E,dEpr + dEaf + dEca,dEpr,dEaf,dEca,P0))

		#Remove connection
		cnt=0
		xs, ys, xl, yl  = randint(xsize),randint(ysize),randint(xsize),randint(ysize)
		while p[xs,ys,xl,yl] == 0 and cnt < 1000 :
			xs, ys, xl, yl  = randint(xsize),randint(ysize),randint(xsize),randint(ysize)
			cnt += 1
		if cnt < 1000:
			dEpr = FdEpr(ax[xs,ys]-1,dn[xl,yl]-1) - FdEpr(ax[xs,ys],dn[xl,yl])
			dEaf = -FdEaf(xs,ys,xl,yl,1)
			dEca = -FdEca(xs,ys,xl,yl,1.,p)
			dE   = dEpr + dEaf + dEca
			if dE > 100:
				P0 = 0.
			elif dE < -100:
				P0 = 1.
			else:
				P0 = 1./(1.+exp(4.*dE) )
			if TotalEnergy:
				P0 *= E/(cE12+E)
			Pminus = False
			if np.random.random() < P0:
				Pminus = True
				p[xs,ys,xl,yl] -= 1
				ax[xs,ys] -= 1
				dn[xl,yl] -= 1
				Epr += dEpr
				Eaf += dEaf
				Eca += dEca
				E   += dE
				if fLog:
					fprintf(log,"%d,-,%d,%d,%d,%d,%g,%g,%g,%g,%g,%g\n",itr_c,xs,ys,xl,yl,dEpr,dEaf,dEca,dEpr + dEaf + dEca,E,P0)
				
		
			if not (itr%1000) and fIndicator: 
				sys.stderr.write(" (-)  %g   %g   [%g,%g,%g]   %g              \r[%s%s] "%(E,dEpr + dEaf + dEca,dEpr,dEaf,dEca,P0,"+" if Pplus else " ","-" if Pminus else " "))
		else:
			if not (itr%1000) and fIndicator: 
				sys.stderr.write(" (-) Not found              \r[%s%s] "%("+" if Pplus else " "," " ))
		if E <= 0. and fTotalEnergy: break
		
		if fReporting and (not (itr%fReporting)):
			fprintf(recdb,"%s:%010d",timestemp_c,itr_c)
			fprintf(recdb,":%g",E)
			pXs, pYs, pXl, pYl = np.where( p.astype(int)  > 0 )
			npos = int(pXs.shape[0])
			fprintf(recdb,":%d",npos)
			for xs,ys,xl,yl in zip(pXs, pYs, pXl, pYl):
				fprintf(recdb,":%d,%d,%d,%d,%d",xs,ys,xl,yl,p[xs,ys,xl,yl])
			fprintf(recdb,"\n")
	print
	print "#############################"
	print "#        STATISTIC          #"
	print "  > Mean Number of Axons     :",np.mean(ax)
	print "  > Mean Number of Dendrites :",np.mean(dn)
	print "#############################\n"

	if StopRec:
		fprintf(recdb,"%s:%010d",timestemp_c,itr_c)
		fprintf(recdb,":%g",E)
		pXs, pYs, pXl, pYl = np.where( p.astype(int)  > 0 )
		npos = int(pXs.shape[0])
		fprintf(recdb,":%d",npos)
		for xs,ys,xl,yl in zip(pXs, pYs, pXl, pYl):
			fprintf(recdb,":%d,%d,%d,%d,%d",xs,ys,xl,yl,p[xs,ys,xl,yl])
		fprintf(recdb,"\n")

	if StartRec or StopRec or bool(fReporting):
		fclose(recdb)
	if fLog:
		fclose(log)

	if Graphs:
		if Init == None:
			from matplotlib import pyplot as plt
			from matplotlib import cm
			cmap = cm.get_cmap('jet')
			Gcells=[ [randint(xsize),randint(ysize),cmap(float(i)/float(nGraphs))] for i in xrange(nGraphs)]
		else:
			plt.subplot(122,sharex=Gax,sharey=Gax)
		for xs,ys,c in Gcells:
			pXl, pYl = np.where( p[xs,ys,:,:]  > 0 )
			for xl,yl in zip(pXl, pYl):
				plt.plot([xs,xl],[ys,yl],"-*",c=c,mfc=c,mec=c,lw=int(p[xs,ys,xl,yl]) )

		plt.title("Model: "+timestemp)
		plt.xlim(0,xsize)
		plt.ylim(0,ysize)
		plt.show()

	return p
