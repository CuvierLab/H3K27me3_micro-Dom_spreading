import sys
f=open(sys.argv[1])
s=sys.stdin.readline()
line =f.readline()
count=0
while s != '':
  s=s.split()
  rank = line.split()[0]
  print rank, s[1]
  count=count+1
  s=sys.stdin.readline()
  line =f.readline()	 
