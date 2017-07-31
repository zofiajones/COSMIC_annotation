import sys
import os
from pybedtools import BedTool
import tabix
import re
import pysam
tb=tabix.open('CosmicCoding_sort.bed.gz')
file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

file=open('exceptions.txt','w')

for line in sys.stdin:
	print(line)
	m=line.strip().split()
	aa=m[4].split(';')[1].replace('p.','')
	genomic=m[4].split(';')[2]
	records=tb.query(m[0],int(m[1])-1 ,int(m[1]))
	aminos=[]
	lines=[]
	for record in records:
		print('record',record)
		index=('@').join(record[-1].split('@')[0:4])
		amino=record[-1].split('@')[-1].split(';')[3].replace('p.','').replace('AA=','')
		strand=record[-1].split('@')[-1].split(';')[1].replace('STRAND=','')
		#print(m[5],index)
		if m[5] == index:
#			print('here')
			aminos.append(amino)
			found=1;
			all_ref=list(set(m[2]))
			all_alt=list(set(m[3]))
			if ( all_ref[0]=='O' )  & ( len(all_ref)==1  ):
#				loc=m[0]+':'+str(int(m[1])-1)+'-'+str(int(m[1])-1)
#				base=file.fetch(region=loc)
#				m[1]=str(int(m[1])-1)
#				m[2]=base.strip().upper()
#				m[3]=base.strip().upper()+m[3]
				m[2]='.'
			elif ( all_alt[0]=='O' )  & ( len(all_alt)==1  ):
#				loc=m[0]+':'+str(int(m[1])-1)+'-'+str(int(m[1])-1)
#				base=file.fetch(region=loc)
#				m[1]=str(int(m[1])-1)
#				m[3]=base.strip().upper()
#				m[2]=base.strip().upper()+m[2]
				m[3]='.'
			elif ('O' in all_ref ) & ('O' in all_alt )  :
				#print('else')
				m[1]='NA'
				m[2]='NA'
				m[3]='NA'
			if '-' in strand :
#				print('here')
				new_ref=''
				ref=m[2]
				if 'O' not in ref:
					for i in ref:
						if 'C' in i:
							new=i.replace('C','G')
						if 'A' in i:
							new=i.replace('A','T')
						if 'G' in i:
							new=i.replace('G','C')
						if 'T' in i:
							new=i.replace('T','A')
						new_ref = new + new_ref
					ref=new_ref
					new_alt=''
				alt=m[3]
				if 'O' not in alt:
					for i in alt:
						if 'C' in i:
							new=i.replace('C','G')
						if 'A' in i:
							new=i.replace('A','T')
						if 'G' in i:
							new=i.replace('G','C')
						if 'T' in i:
							new=i.replace('T','A')
						new_alt = new + new_alt
					alt=new_alt
			fields=record[-1].split('@')[-1].split(';')
			#print(fields)
			for field in fields:
				if 'CDS' in field:
					cds=field.replace('CDS=','')
				if 'AA' in field:
					aa=field.replace('AA=','').replace('p.','')
				if 'STRAND' in field:
					strand=field.replace('STRAND=','')
				if 'COSM' in field:
					cosm=field
				if 'GENE' in field:
					gene=field.replace('GENE=','').split('_')[0]
			info=(',').join([gene,cosm,aa,strand,'NA',genomic])
			lines.append(m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+m[3]+'\t'+info)
#	print(lines)
	indices=[]
	if len(aminos)>1:
		for i in range(len(aminos)):
			if aa in aminos[i]:
				indices.append(i)
#		print(indices,len(indices))
		if len(indices)==0:
			for i in range(len(aminos)):
				if 'COSM' in lines[i]:
					print(lines[i])
	elif len(aminos)==1:
		indices=[0]
	else:
		file.write(line)
		#print('here1',line)
		#		print(aa,aminos[i])
	#print(indices)
	#print(lines)
	for i in indices:
		print(lines[i])

file.close()
