# -*- coding: utf-8 -*-
import sys
import string

type = sys.argv[1]
label = sys.argv[2]
#prfx = sys.argv[3]

ii = 0
s = sys.stdin.readline()
while s != '':
  if (s[0] != '#') and (s.count('track')==0) and (s.count('browser')==0):
    ii = ii + 1
    s = s.split()

    c = s[0]
#    if len(c) > 3:
#	c = c[3:]
    b = int(s[1])
    e = int(s[2])
    n =s[3]
#    if len(s) > 3:
#	    n = s[3]
#    else:
#  	    n = c+'_'+prfx+'_'+`ii`
    print string.join([c,type,label,`b`,`e`,`e-b+1`,'+','.','Feature '+n],'\t')

  s = sys.stdin.readline()
