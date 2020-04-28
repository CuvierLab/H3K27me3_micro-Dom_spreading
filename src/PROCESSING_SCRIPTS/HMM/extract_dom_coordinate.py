# -*- coding: utf-8 -*-
import sys

gap_size=int(sys.argv[1])
label=sys.argv[2]
#step=int(sys.argv[2])

chrom = ['2L','2R','3L','3R','2LHet','3LHet','2RHet','3RHet','4','X','XHet','YHet']



for c in chrom:
  count_dom=0
  begin=0
  f=open(c+'_'+label,'r')
  shift=0
  line=f.readline()
  while (line !=''):
   line=line.split()
   x=int(line[0])
   y=float(line[1])
   if (y>0.5 or begin==1):
    if begin==0:
      begin_pos=x
    line=f.readline()
    line=line.split()
    y=float(line[1])
    while (line!='' and y>0.5):
      line=f.readline() 
      if line!='':
	line=line.split()
	x=int(line[0])
	y=float(line[1])         
    pos_inter=x
    if line !='':
     line=f.readline() 
     line=line.split()
     x=int(line[0])
     y=float(line[1])
     while (line!='' and y<=0.5):
      line=f.readline()
      if line !='':
	line=line.split()
	x=int(line[0])
	y=float(line[1])
    if (x - pos_inter) > gap_size and count_dom==0 :
      print c, str(begin_pos-100), str(pos_inter-100), 'big_domain'
      begin_pos=x
      count_dom=0
    elif (x - pos_inter) > gap_size and count_dom!=0 :
      print c, str(begin_pos-100), str(pos_inter-100), 'small_domain_end'
      begin_pos=x
      count_dom=0
    elif (x - pos_inter) < gap_size and count_dom==0:
      print c, str(begin_pos-100), str(pos_inter-100), 'small_domain'+str(count_dom)
      begin_pos=x
      count_dom=1
    elif (x - pos_inter) < gap_size and count_dom!=0:
      print c, str(begin_pos-100), str(pos_inter-100), 'small_domain_'+str(count_dom)
      begin_pos=x
      count_dom+=1     
    begin=1
   line=f.readline()
      
      
      
    
  

