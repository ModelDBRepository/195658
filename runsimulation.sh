#! /bin/bash

#setup bigger number of step Nstep for better convergence.
#This parameter is 15000 for any simulation in the paper.
#For NSTEP=15000 each simulation takes ~ 1 day
[ -z "$NSTEP" ] && NSTEP=2000

#Setup INICATOR=True if you want see progress report.
[ -z "$INDICATOR" ] && INDICATOR=False




if [[ "$NSTEP" != '0' ]]; then
	[ -d 'local' ] || mkdir local
	python steup.py install --install-lib=local ||{ echo "Cannot compile EnergyChaser module"; exit 1; }
	
	pushd local

	#Random model ID for all runs
	[ -z "$MODELID" ] && MODELID=`python -c 'import time; from numpy import random as rnd; print time.strftime("%Y%m%d%H%M%S") + "-%03d"%rnd.randint(1000)'`
	
	#Correlational Model WT
	python EnergyModel.py /Model=1 /Gca=20.                          /Knocked=False /Apr=450 /ModelID=\"$MODELID"-CM-WT"\"   /Graphs=False /Log=False /RunDB=\"$MODELID"-CM-WT.csv"\"   /Report=False /StopRec=True /Nstep=$NSTEP /Indicator=$INDICATOR &

	#Correlational Model KI/KI
	python EnergyModel.py /Model=1 /Gca=20.     /Sca=1.              /Knocked=True  /Apr=450 /ModelID=\"$MODELID"-CM-KI"\"   /Graphs=False /Log=False /RunDB=\"$MODELID"-CM-KI.csv"\"   /Report=False /StopRec=True /Nstep=$NSTEP /Indicator=$INDICATOR &

	#Correlational Model KI/KI / B2-/-
	python EnergyModel.py /Model=1 /Gca=20./2.3 /Sca=4. /Bca=Bca*2.1 /Knocked=True  /Apr=450 /ModelID=\"$MODELID"-CM-KIB2"\" /Graphs=False /Log=False /RunDB=\"$MODELID"-CM-KIB2.csv"\" /Report=False /StopRec=True /Nstep=$NSTEP /Indicator=$INDICATOR &

	wait

	#Integrational Model WT
	python EnergyModel.py /Model=3 /Vca=15 /Rca=3. /Gca=20.                         /Knocked=False /Apr=450 /ModelID=\"$MODELID"-IM-WT"\"   /Graphs=False /Log=False /RunDB=\"$MODELID"-IM-WT.csv"\"   /Report=False /StopRec=True /Nstep=$NSTEP /Indicator=$INDICATOR &

	#Integrational Model KI/KI
	python EnergyModel.py /Model=3 /Vca=15 /Rca=3. /Gca=20.     /Sca=1.             /Knocked=True  /Apr=450 /ModelID=\"$MODELID"-IM-KI"\"   /Graphs=False /Log=False /RunDB=\"$MODELID"-IM-KI.csv"\"   /Report=False /StopRec=True /Nstep=$NSTEP /Indicator=$INDICATOR &

	#Integrational Model KI/KI / B2-/-
	python EnergyModel.py /Model=3 /Vca=15 /Rca=3. /Gca=20./2.3 /Sca=4. /Bca=Bca*2.1 /Knocked=True  /Apr=450 /ModelID=\"$MODELID"-IM-KIB2"\" /Graphs=False /Log=False /RunDB=\"$MODELID"-IM-KIB2.csv"\" /Report=False /StopRec=True /Nstep=$NSTEP /Indicator=$INDICATOR &

	wait
	cat $MODELID"-CM-WT-rec.db" $MODELID"-CM-KI-rec.db" $MODELID"-CM-KIB2-rec.db" $MODELID"-IM-WT-rec.db" $MODELID"-IM-KI-rec.db" $MODELID"-IM-KIB2-rec.db" >$MODELID"-tot.db" &&\
	gzip -9 $MODELID"-tot.db" &&\
	rm $MODELID"-CM-WT-rec.db" $MODELID"-CM-KI-rec.db" $MODELID"-CM-KIB2-rec.db" $MODELID"-IM-WT-rec.db" $MODELID"-IM-KI-rec.db" $MODELID"-IM-KIB2-rec.db"  &&\
	cat $MODELID"-CM-WT.csv" $MODELID"-CM-KI.csv" $MODELID"-CM-KIB2.csv" $MODELID"-IM-WT.csv" $MODELID"-IM-KI.csv" $MODELID"-IM-KIB2.csv" >$MODELID"-tot.csv" &&\
	rm $MODELID"-CM-WT.csv" $MODELID"-CM-KI.csv" $MODELID"-CM-KIB2.csv" $MODELID"-IM-WT.csv" $MODELID"-IM-KI.csv" $MODELID"-IM-KIB2.csv" &&\
	python ../2D-View.py        $MODELID"-tot.db.gz" /Nraaw=3 /R=4 & python ../1D-DensityView.py $MODELID"-tot.db.gz" /Nraw=3  /Sml=20 &&\
	popd
else
	[ -z "$MODELID" ] && MODELID="20160910104211-402"
	python 2D-View.py        $MODELID"-tot.db.gz"  /Nraaw=3 /R=4 & python 1D-DensityView.py $MODELID"-tot.db.gz"  /Nraw=3  /Sml=20
fi

echo 'done'

