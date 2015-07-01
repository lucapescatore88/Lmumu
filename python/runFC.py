#! /usr/bin/env python

import os
import sys
from string import find

maxAfb = 0.75
minAfb = -0.75
maxfL = 1
minfL = 0
stepAfb = 0.035
stepfL = 0.035

fLcut = 0.2
afbminCut = -1
afbmaxCut = 1
obs = "afb"

args = 1
restore = False
if sys.argv[1] == "-r" :
	restore = True
	args+=1;

namexec = 'AnafitAfb_3D_new.out'
dim = '3D'
if sys.argv[2] == "-1D" :
	print "Performing 1D fit"
	namexec = 'AnafitAfb_nosW_new.out'
	dim = '1D'
	args+=1;

if len(sys.argv) < args+2 :
	print "Not enough arguments"
	sys.exit()

bin = sys.argv[args]
obs = sys.argv[args+1]

print bin, obs, dim, restore

if obs == "afbfL" :
	if bin == '2' :
		afbmaxCut = 0.75
		fLcut = 0.
		stepfL = 0.02
		stepAfb = 0.02
	if bin == '4' :
		fLcut = 0.35
		stepAfb = 0.02
		stepfL = 0.02
	if bin == '3' :
		stepAfb = 0.025
		stepfL = 0.025
		fLcut = 0.05
		afbminCut = -0.5
		afbmaxCut = 0.35
	if bin == '5' :
		fLcut = 0.25
		stepAfb = 0.025
		stepfL = 0.025
		afbminCut = -0.3
		afbmaxCut = 0.3
	if bin == '6' :
		stepAfb = 0.025
		stepfL = 0.025
		fLcut = 0.
elif obs == "afb":
	minAfb = -0.75
	maxAfb = 0.75
	stepAfb = 0.01
elif obs == "afbB":
	minAfb = -0.5
	maxAfb = 0.2
	stepAfb = 0.01
elif obs == "fL":
	minAfb = 0.
	maxAfb = 1.
	stepAfb = 0.01

#print stepAfb, stepfL, fLcut

np_per_job = 1

def drange(start, stop, step):
	r = start
	while r <= stop:
		yield r
		r += step



fs = []
As = []

if obs == "afbfL" :
	for a in drange(minAfb,maxAfb,stepAfb) :
		for f in drange(minfL,maxfL,stepfL) :
			if not ( (f-1)*3./4. > a or a > -(f-1)*3./4. ) and f >= fLcut and a > afbminCut and a < afbmaxCut :
				fs.append(f)
				As.append(a)
else :
	for a in drange(minAfb,maxAfb,stepAfb) :
		fs.append(0)
		As.append(a)



ntot_jobs = len(fs)//np_per_job + int(len(fs)%np_per_job!=0)
print ntot_jobs, "jobs,", len(fs), "points"
#if(len(sys.argv) > 4 ) :
#	sys.exit()

namebin = "Fit_"+obs+"_"+dim+"_B" + str(bin)
fold = "/afs/cern.ch/work/p/pluca/jobs/"+namebin

os.system("mkdir -p " + fold)
stepfile = open(fold + "/steps.txt","w")
stepfile.write("stepAfb = " + str(stepAfb)+"\n")
stepfile.write("stepfL = " + str(stepfL))
stepfile.close()

sfs = "["
sAs = "["
i=0
while i < len(fs) :

	sfs += str(fs[i])+","
	sAs += str(As[i])+","
	
	if((i+1)%np_per_job==0 or i==(len(fs)-1) ) :
		sfs = sfs[:-1]
		sAs = sAs[:-1]
		sfs += "]"
		sAs += "]"
		jobID = str((i+1)//np_per_job + int((i+1)%np_per_job!=0))
	
		#print "/afs/cern.ch/work/p/pluca/jobs/"+namebin+"/"+jobID+"/run.sh", not os.path.isfile("/afs/cern.ch/work/p/pluca/jobs/"+namebin+"/"+jobID+"/run.sh")
		if not restore or not os.path.isfile("/afs/cern.ch/work/p/pluca/jobs/"+namebin+"/"+jobID+"/run.sh"):
			print "---- Job " + jobID + " / " + str(ntot_jobs) + " : submit -uexe -d " + namebin + " -n " + jobID + " '" + namexec + " " + sAs + " " + sfs + " -E" + obs + " -B" + str(bin) + "'"
			os.system( "~/python/submit.py -uexe -d " + namebin + " -n " + jobID + " '" + namexec + " " + sAs + " " + sfs + " -E" + obs + " -B" + str(sys.argv[1]) + "'" ) 
		
		sfs = "["
		sAs = "["
	
	i+=1

