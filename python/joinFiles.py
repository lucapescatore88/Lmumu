import sys

lines1 = [line.strip() for line in open(sys.argv[1])]
lines2 = [line.strip() for line in open(sys.argv[2])]

for l in range  (0,len(lines1)):
	print lines1[l]+"\t"+lines2[l]
