import subprocess
import os
import re
import MySQLdb
import pandas as pd

#Initialise mysql
db=MySQLdb.connect(db="mysql",read_default_file="~/config.cnf")
c=db.cursor()

cmd0='select gene,protein_change,genomic_change,id from mutations;'
c.execute(cmd0)
results=c.fetchall()
print(type(results))
gene=[]
protein=[]
genomic=[]
id=[]

file2=open('full.vcf','w')

for entry in results:
	gene0=entry[0]
	protein0=entry[1].replace('p.','')
	if not entry[2]:
		loc='NA'
	else:
		loc=entry[2]
	if not loc:
		loc='NA'
	if not protein0:
		prot='NA'
	loc0=loc;
	genomic.append(loc0)
	protein.append(protein0)
	gene.append(gene0)
	id.append(entry[3])
	if 'NA' not in loc:
		m=loc.replace('g.','').split(':')
		print(m)
		chr=m[0]
		m0=re.search('(\d+)',m[1])
		loc=m0.group(1)
		if '>' in m[1]:
			n=m[1].split('>')
			m0=re.search('([A-Z]+)',n[0])
			#print(n)
			ref=m0.group(1)
			#print(ref)
			alt=n[1]
		elif 'del' in m[1]:
			n=m[1].replace('\d+','').split('del')
			ref=n[1]
			alt=''
			for char in range(len(ref)):
				alt=alt+'O'
		elif 'ins' in m[1]:
			n0=m[1].split('_')
			loc=n0[0]
			n=n0[1].replace('\d+','').split('ins')
			alt=n[1]
			ref=''
			for char in range(len(alt)):
				ref=ref+'O'
		else:
			print(m[1])
		index=('@').join([chr,loc,ref,alt])
		file2.write(('\t').join([chr,loc,ref,alt,gene0+';'+protein0+';'+loc0,index ]) + '\n')



dict={}
dict['gene']=gene
dict['protein']=protein
dict['genomic']=genomic
dict['id']=id
df=pd.DataFrame(dict)
df.drop_duplicates(inplace=True)

df.to_csv('mutations.csv',index=False)

df1=df[df['genomic'] == 'NA']


df1[['gene','protein']].to_csv('mutations_1.csv',sep='\t',index=False)
file2.close()

db=MySQLdb.connect(db="mysql",read_default_file="~/config.cnf")
c=db.cursor()

cmd0='select id,description from cell_lines;'
c.execute(cmd0)
results=c.fetchall()
ids=[]
par=[]
freq=[]
gene=[]
aa=[]
for entry in results:
	m=entry[1]
	n=m.split()
	print(n,'here')
	id=entry[0]
	if len(n)==4:
		ids.append(id)
		par.append(n[0])
		freq.append(n[1])
		if n[2]:
#			print('here',n[2])
			#gene.append(n[2].replace('PI3KÎ±','PIK3CA').replace('PI3KCA','PIK3CA'))
			gene.append(n[2])
#			print(gene.append(n[2]).replace('PI3KÎ±','PIK3CA'))
#			print(n[2],gene.append(n[2]).replace('PI3KÎ±','PIK3CA').replace('PI3KCA','PIK3CA'))
		else:
			gene.append(n[2])
		#aa.append(n[3].replace('185delAG/+','185delAG').replace('V769_D770insASV','V769-D770insASV'))
		aa.append(n[3])

dict_1={}
dict_1['id']=ids
dict_1['par']=par
dict_1['freq']=freq
dict_1['gene']=gene
dict_1['aa']=aa

df1=pd.DataFrame(dict_1)

df1.to_csv('cell_lines.csv')

