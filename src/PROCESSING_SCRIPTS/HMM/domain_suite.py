# -*- coding: utf-8 -*-

import sys

dist=int(sys.argv[1]) #dist la distance avec le next boundary pour le fusionner si cette distance est inférieure à dist

long_min=int(sys.argv[2]) #la longueur minimale d'un boundary

f = open(sys.argv[3],'r') #typiquement ~/Documents/gff/WTK27-5000bpdom.gff

line = f.readline()

t = {}

#après avoir extrait les domaines à partir du fichier python extract_domain_coordinate.py 
#faire en sorte de fusionner les domaines trop petits en fonction de la distance avec le prochain

while line!='':
	
	line = line.split()
	chrom = line[0]
	pos1 = int(line[1])
	pos2 = int(line[2])
	
	if not t.has_key(chrom):
		t[chrom]={}
		
	if not t[chrom].has_key(`pos1`):
		t[chrom][`pos1`]=[pos2]
		
	line = f.readline()
	
vect = {}

for chrom in t.keys():
	
	vect[chrom] = map(int,t[chrom].keys())
	vect[chrom].sort()
	
t2 = {}
	
for chrom in vect.keys():
	
	if not t2.has_key(chrom):
		t2[chrom]={}
	
	count=0
	
	i=0
	while (i <( len(vect[chrom])-1)):
		
		pos1_1 = int(vect[chrom][i])
		
		pos2_1 = int(vect[chrom][i+1])
		
		pos1_2 = int(t[chrom][`pos1_1`][0])
		
		if count==0 and not t2[chrom].has_key(pos1_1):
			
			t2[chrom][`pos1_1`]=0
			count=1
			pos_tmp=`pos1_1`
			
		if (pos2_1-pos1_2)<dist :
			i+=1
			
		else :
			
			t2[chrom][pos_tmp]=pos1_2
			count=0
			i+=1
			
for chrom in t2.keys():
	
	sorted=map(int,t2[chrom].keys())
	sorted.sort()
	
	for pos1 in sorted:
		
	#	print t2[chrom][pos1]
		if int(t2[chrom][`pos1`])-int(`pos1`)>long_min :
		
			print chrom, pos1, t2[chrom][`pos1`]		
		
		
		
		




	
	
