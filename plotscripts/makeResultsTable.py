import math

normal = open("sum_yields_normal.txt")
lines = normal.readlines()
nosig = open("sum_yields_nosig.txt")
lines_nosig = nosig.readlines()

LogL = []
LogL_nosig = []
q2 = []
y = []

for l in lines :
	el = l.split()
	q2.append( el[0] )
	y.append( el[2]+el[3]+el[4] ) 
	LogL.append( float(el[6].replace("\\","")) )
for l in lines_nosig :
	el = l.split()
	LogL_nosig.append( float(el[6].replace("\\","")) )

for q,yy,s,ns in zip(q2,y,LogL,LogL_nosig) :
	print q, " & ", yy, " & ", '{:4.3}'.format(math.sqrt(2*abs(ns-s))), " \\\\"

