import sys
import os
from pybedtools import BedTool
import tabix
import re
import pysam
tb=tabix.open('CosmicCoding_sort.bed.gz')
file=pysam.FastaFile('/home/zjones/fb_genome/genome.fa')

#file=open('exceptions.txt','w')

for line in sys.stdin:
	#print(line)
	m=line.strip().split()
	genomic=m[4].split(';')[2]
	aa=m[4].split(';')[1].replace('p.','')
	records=tb.query(m[0],int(m[1])-1 ,int(m[1]))
	aminos=[]
	lines=[]
	for record in records:
		#print(record)
		index=('@').join(record[-1].split('@')[0:2])
		m_index=('@').join(m[5].split('@')[0:2])
		m_base=('@').join(m[5].split('@')[2:4])
		cos_base=('@').join(record[-1].split('@')[2:4])
		#print(m_base,m[5])
		amino=record[-1].split('@')[-1].split(';')[3].replace('p.','').replace('AA=','')
		#print(m[5],index)
		if m_index == index:
			#print('here')
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
			info=(',').join([gene,cosm,aa,strand,'NA',genomic,cos_base])
			lines.append(m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+m[3]+'\t'+info)
	for i in range(len(lines)):
		print(lines[i])
	#if len(lines)==0:
	#	print(line)
