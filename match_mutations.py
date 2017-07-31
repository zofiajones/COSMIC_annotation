import matplotlib
matplotlib.use('Agg')
import argparse
import warnings
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib.legend as legend
from matplotlib.backends.backend_pdf import PdfPages
import re
import numpy as np
import pandas as pd
from matplotlib_venn import venn3,venn2
import pylab
import plotly.plotly as py
import plotly
import cufflinks as cf
import plotly.tools as tls
#tls.embed('https://plot.ly/~cufflinks/8')

df0=pd.read_csv('/home/zjones/MultiPlex/mutations.csv',header='infer')
df1=pd.read_csv('multiplex_annotations.txt',header=None,names=['chr','loc','ref','alt','info'],sep='\t')

df11=pd.read_csv('/home/zjones/MultiPlex/cell_lines.csv',header='infer')

#df0=pd.read_csv('test.txt',header='infer')
#df1=pd.read_csv('multiplex_annotations.txt',header=None,names=['chr','loc','ref','alt','info'],sep='\t')

cell=[]

info=list(df1['info'])

gene=[]
cosm=[]
protein=[]
protein1=[]
genomic=[]
cosm_base=[]
i=0
for item in info:
	print(item)
	m=item.split(';')
	gene.append(m[0])
	cosm.append(m[1])
	protein.append(m[2])
	print(m[0],m[2])
	protein1.append(m[2].split('fs')[0])
	#print('here1',m[2].split('fs')[0])
#	print(m[5],'here15')
	if len(m) == 6:
		if 'chr' in m[5]:
			genomic.append(m[5])
			i=i+1
			print('here16',m[5],i)
		else:
			genomic.append('NA')
	else:
		genomic.append('NA')
	if len(m) == 7:
		cosm_base.append(m[6])
	else:
		cosm_base.append('NA')

dict_1={}
dict_1['gene']=gene
dict_1['cosm']=cosm
dict_1['protein']=protein
dict_1['protein1']=protein1
dict_1['genomic']=genomic
dict_1['cosm_base']=cosm_base
df2=pd.DataFrame(dict_1)

print('here13',genomic[:10],i)

cosm=[0]*len(df0['gene'].tolist())
print(len(cosm))
print(cosm[:10])

eng=['NA']*len(df0['gene'].tolist())
cosm_base=[0]*len(df0['gene'].tolist())

corrected=[]
for i in range(len(df0['gene'].tolist())):
	cosm[i]='NS'
	gene=df0['gene'].iloc[i]
	#print(gene)
	protein=df0['protein'].iloc[i]
	new_protein=protein
	print(new_protein,'here16a')
	stra=df11[(df11['gene']==gene) & (protein == df11['aa'])]['id'].tolist()
	strb=''
	for item in stra:
		strb=strb+';'+str(item)
		cell.append(item)
	if not strb:
		strb='NA'
	print('here12',strb,new_protein)
	eng[i]=strb
	loc=df0['genomic'].iloc[i]
	stra=df2[(df2['gene']==gene) & (protein == df2['protein'])]['cosm_base'].tolist()
	strb=''
	c_base=[]
	for item in stra:
		strb=strb+';'+str(item)
		c_base.append(item)
	if not strb:
		strb='NA'
	cosm_base[i]=strb
	#print(gene,protein)
	df3=df2[(gene == df2['gene']) &  (protein == df2['protein'])]
	str0=(';').join(df3['cosm'].tolist())
	if str0:
		cosm[i]=str0;
		print('here12',strb)
		print(new_protein,'here16x')
	else:
		protein1=df0['protein'].iloc[i].split('fs')[0]
		#print(protein1)
		df3=df2[(gene == df2['gene']) &  (protein1 == df2['protein1'])]
		str1=''
		str1=(';').join(df3['cosm'].tolist())
		if 'COSM' in str1:
			cosm[i]=str1;
			new_protein=(';').join(df3['protein'].tolist())
			stra=(';').join(df3['cosm_base'].tolist())
			cosm_base[i]=stra
			print(new_protein,'here16z')
		else:
			print('here11',loc)
			df4=df2[  loc == df2['genomic'] ]
			str2=''
			str2=(';').join(df4['cosm'].tolist())
			new_protein=(';').join(df4['protein'].tolist())
			print(new_protein,'here16y')
			if str2:
				cosm[i]=str2;
				stra=(';').join(df4['cosm_base'].tolist())
				cosm_base[i]=stra
			else:
				print('here1')

	corrected.append(new_protein)
	print(new_protein,'here16b')



df0['corrected']=corrected
df0['eng']=eng
df0['cosm_base']=cosm_base

df0['COSMIC']=cosm
print(cosm[:10])
print(len(df0['gene'].tolist()))
df4=df0[df0['COSMIC'] == 'NS']
df4[['gene' , 'protein']].to_csv('mutations_3.csv',index=False,sep='\t')
df4.to_csv('mutations_4.csv',index=False,sep='\t')

print(len(df4['gene'].tolist()))
df5=df0[df0['COSMIC'] != 'NS']
print(len(df5['gene'].tolist()))

df0.to_csv('mutations_2.csv')

dict_z={}
dict_z['cell']=cell
dfz=pd.DataFrame(dict_z)
dfz.to_csv('mutations_5.csv')

cell0=set(df11['id'].tolist())
cell2=cell0.difference(set(cell))
print(cell2)

for item in cell2:
	#print(item)
	gene=df11[df11['id'] == int(item)]['gene'].tolist()[0].replace('PI3KÎ±','PIK3CA').replace('PI3KCA','PIK3CA')
	aa=df11[df11['id']==item]['aa'].tolist()[0].replace('V769_D770insASV','V769-D770insASV').replace('185delAG','185delAG/+')
	if 'wild' not in aa:
		print(gene,aa,item)
		strb=df0[(df0['gene']==gene) & (df0['protein']==aa)]['eng'].tolist()[0].replace('NA','')
		ind=df0[(df0['gene']==gene) & (df0['protein']==aa)].index.values.tolist()[0]
		print(strb)
		df0.set_value(ind,'eng',strb+';'+str(item))
		print(df0[(df0['gene']==gene) & (df0['protein']==aa)]['eng'].tolist(),strb+';'+str(item))
		#stra=df0['eng'].iloc[strc]
		#df0['eng']
		#if strc:
		#	strb=strc[0]
		#stra=df11[(df11['gene']==gene) & (aa == df11['aa'])]['id'].tolist()
		#for item in stra:
		#	strb=strb+';'+str(item)
		#df0[(df0['gene']==gene) & (df0['protein']==aa)]['eng']=strb
		#print('found',aa,gene,strb,stra)


df0.to_csv('mutations_2.csv')


eng_count=[0]*len(df0['eng'].tolist())
for i in range(len(df0['eng'].tolist())):
	eng_count[i]=len(df0['eng'].tolist()[i].split(';'))

df0['count']=eng_count
df0.to_csv('mutations_2.csv',sep='\t')

pp=PdfPages('venn.pdf')

gene0=list(set(df0['gene'].tolist()))
gene1=[0]*len(gene0)
counts=[0]*len(gene0)

for i in range(len(gene0)):
	print(gene0[i])
	gene1[i]=len(df0[df0['gene'] == gene0[i]]['count'].tolist())
	count=0
	for count0 in df0[df0['gene'] == gene0[i]]['count'].tolist():
		count=count+int(count0)
	counts[i]=count

dict_2={}
dict_2['gene']=gene0
dict_2['count']=gene1
dict_2['cell']=counts
dfx=pd.DataFrame(dict_2)

dfx.to_csv('gene_count.csv')
dfw=pd.DataFrame( dfx.groupby([ 'cell','count'  ]).describe()  )
dfw.to_csv('gene_count_grouped.csv')
