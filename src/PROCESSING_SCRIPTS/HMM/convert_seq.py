import sys
import string

win = int(sys.argv[1])
shift = int(sys.argv[2])
state = sys.argv[3]
chrom = sys.argv[4]
label = sys.argv[5]

s = sys.stdin.readline() # token line
s = sys.stdin.readline()
s = s.split()

xi = win/2
xf = xi + shift
flag = 0
ii = 0
for x in s:

    if x == state:
        if flag == 0:  # start new domain
            xL = xi
            xR = xi
            flag = 1
        else:
            xR = xf
    else:
        if flag == 1:
            ii = ii+1
            print string.join([chrom,'domain',label,`xL`, `xR`, `xR - xL + 1`,'+','.','Feature '+chrom+'_'+label+'_'+`ii`],'\t')
            flag = 0

    xi = xi + shift
    xf = xi + shift
        
if flag == 1:
    ii = ii + 1
    print string.join([chrom,'domain',label,`xL`, `xR`, `xR - xL + 1`,'+','.','Feature '+chrom+'_'+label+'_'+`ii`],'\t')
    
