import sys
import os
import subprocess

rootdir = "/afs/cern.ch/work/p/pluca/jobs/"+sys.argv[1]

count = 0
total = 0

def run(command) :
	return subprocess.Popen(command,stdout=subprocess.PIPE).communicate()

for subdir, dirs, files in os.walk(rootdir):
	for dir in dirs :
		total+=1
		dirname = rootdir+"/"+dir
		#if not os.path.isfile(dirname+"/out") :
		if run(["/bin/grep", "Best", dirname+"/out"])[0] == "" :
			count += 1
			print "Resubmitting job ID ", dir 
			os.system("bsub -R 'pool>30000' -o " + dirname + "/out -e " + dirname + "/err -q 8nh -J " + sys.argv[1]+"_"+dir + "  < " + dirname + "/run.sh")

print count, " jobs resubmitted (of ", total, " total jobs)"
