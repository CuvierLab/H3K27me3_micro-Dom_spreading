import sys
import string

miny = float(sys.argv[1])
maxy = float(sys.argv[2])
Nbin = int(sys.argv[3])
dy = (maxy-miny)/Nbin

seq = []
s = sys.stdin.readline()
while s != '':
	s = s.split()
	y = float(s[1])
	if y <= miny:
		y = miny+0.001
	if y >= maxy:
		y = maxy - 0.001
	lett = `int((y-miny)/dy)+1`
	seq.append(lett)
	s = sys.stdin.readline()

print 'T=', len(seq)
print string.join(seq)
